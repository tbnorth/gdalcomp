"""
gdalcomp.py - collect GDAL raster calc. utilities

With thanks to https://pcjericks.github.io/py-gdalogr-cookbook/

Terry Brown, Terry_N_Brown@yahoo.com, Tue Jul 14 10:29:06 2015
"""

import argparse
import os
import random
import sys
import unittest

from collections import namedtuple
from hashlib import sha1

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import numpy as np

BBox = namedtuple("BBox", "l b r t")
Block = namedtuple("Block", "c r w h")
def make_parser():

    parser = argparse.ArgumentParser(
        description="""Raster calcs. with GDAL.
        The first --grid defines the project, extent, cell size, and origin
        for all calculations, all other grids are transformed and resampled
        as needed to match.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--grid", type=str, nargs='*', action='append',
        help="""Supply name, path, and resample method for grid input.
        Resample method is ignored for the first --grid, and defaults
        to 'neared_neighbor'.  Alternatives are 'bilinear' and
        'cubic'"""
    )

    parser.add_argument("--expr", type=str, action='append',
        help="""A statement to execute with grid names as numpy arrays,
        of the form `<var> = <expr>`, e.g. `wl = mwl & (cti >= 10.5)`"""
    )

    parser.add_argument("--output-dir", type=str, default='.',
        help="Base directory for outputs"
    )

    parser.add_argument("--output", type=str, nargs='+', action='append',
        help="Names to output, from --expr (or --grid)", default=[]
    )

    parser.add_argument("--tile-cols", type=int, default=1024,
        help="Tile width."
    )
    parser.add_argument("--tile-rows", type=int, default=None,
        help="Tile height, defaults to tile width."
    )
    parser.add_argument("--tiles-merge", action='store_true', default=False,
        help="Merge output tiles to single .tif."
    )
    parser.add_argument("--tiles-only",
        help="Space separated list of tiles to process, e.g. '1,1 1,2 5,5'"
    )
    parser.add_argument("--overlap-cols", type=int, default=0,
        help="NOT IMPLEMENTED Tile overlap for kernel filters etc."
    )
    parser.add_argument("--overlap-rows", type=int, default=None,
        help="""NOT IMPLEMENTED Tile overlap for kernel filters etc.
        Defaults to `overlap-cols`"""
    )

    parser.add_argument("--import", type=str, action='append', default=[],
        help="NOT IMPLEMENTED `import *` from the named module into the execution environment"
    )


    #parser.add_argument('<|positional(s)|>', type=str, nargs='+',
    #    help="<|help|>"
    #)

    return parser
def grid_annotations(grid):
    """grid_annotations - a dict of useful info for a grid

    :return: dict
    """
    gt = grid.GetGeoTransform()
    return dict(
        rows = grid.RasterYSize,
        cols = grid.RasterXSize,
        left = gt[0],
        top = gt[3],
        sizex = gt[1],
        sizey = -gt[5],
        bottom = gt[3] - -gt[5] * grid.RasterYSize,
        right = gt[0] + gt[1] * grid.RasterXSize,
    )
def annotate_grid(grid):
    """annotate_grid - add useful info to grid as attributes

    :param grid: grid to annotate
    """
    grid.__dict__.update(grid_annotations(grid))
def bbox_for(grid, block):
    """bbox_for - Return BBox(left, bottom, right, top) map unit coords for inputs

    Inverse of cells_for().

    :Parameters:
    - `grid`: grid containing cells, annotated by annotate_grid()
    - `block`: Block of grid cells
    """

    return BBox(
        l = grid.left + block.c * grid.sizex,
        b = grid.top - (block.r+block.h) * abs(grid.sizey),
        r = grid.left + (block.c+block.w) * grid.sizex,
        t = grid.top - block.r * abs(grid.sizey)
    )
def cells_for_naive(grid, bbox):
    """cells_for_naive - identify cells in grid for rectangle in map units
    Return Block(col, row, width, height)

    This "naive" version is clear, but not as efficient as cells_for()

    :Parameters:
    - `grid`: grid containing cells, as annotated by annotate_grid()
    - `bbox`: BBox bounds in map units
    """

    # simple implementation
    #D print (bbox)
    #D print ((bbox.l - grid.left) / grid.sizex)

    lcol = (bbox.l - grid.left) / grid.sizex
    if lcol - int(lcol) > 0.5:
        lcol += 1
    lcol = int(lcol)

    #D print ((bbox.r - grid.left) / grid.sizex)
    rcol = (bbox.r - grid.left) / grid.sizex
    if rcol - int(rcol) < 0.5:
        rcol -= 1
    rcol = int(rcol)

    #D print((grid.top - bbox.t) / grid.sizey)
    trow = (grid.top - bbox.t) / grid.sizey
    if trow - int(trow) > 0.5:
        trow += 1
    trow = int(trow)
    #D print((grid.top - bbox.b) / grid.sizey)
    brow = (grid.top - bbox.b) / grid.sizey
    if brow - int(brow) < 0.5:
        brow -= 1
    brow = int(brow)
    #D print(lcol, rcol, trow, brow)
    # a faster implementation using half grid size steps may be possible

    assert bbox.t >= bbox.b

    if rcol < lcol:
        assert rcol == lcol - 1
        rcol = lcol = int((bbox.l + (bbox.r-bbox.l)/2. - grid.left) / grid.sizex)
    if brow < trow:
        assert brow == trow - 1
        brow = trow = int(grid.top - (bbox.t + (bbox.t-bbox.b)/2.) / grid.sizey)

    assert rcol >= lcol, (rcol, lcol, bbox)
    assert brow >= trow, (brow, trow, bbox)

    return Block(
        c=lcol,
        r=trow,
        w=rcol-lcol+1,
        h=brow-trow+1
    )
def cells_for(grid, bbox):
    """cells_for - identify cells in grid for rectangle in map units
    Return Block(col, row, width, height)

    :Parameters:
    - `grid`: grid containing cells, as annotated by annotate_grid()
    - `bbox`: BBox bounds in map units
    """

    hx = grid.sizex / 2.
    lcol = int((bbox.l - grid.left) / hx + 1) // 2
    rcol = int((bbox.r - grid.left) / hx - 1) // 2
    hy = grid.sizey / 2.
    trow = int((grid.top - bbox.t) / hy + 1) // 2
    brow = int((grid.top - bbox.b) / hy - 1) // 2

    if rcol < lcol:
        rcol = lcol = int((bbox.l + (bbox.r-bbox.l)/2. - grid.left) / grid.sizex)
    if brow < trow:
        brow = trow = int(grid.top - (bbox.t + (bbox.t-bbox.b)/2.) / grid.sizey)

    return Block(
        c=lcol,
        r=trow,
        w=rcol-lcol+1,
        h=brow-trow+1
    )
def cells_from(grid0, grid1, block=None, bbox=None, resample="nearest_neighbor"):
    """cells_from - re-project and re-sample subregion from grid1 to align with grid0

    Must supply either block or bbox.  Grids will be annotate_grid()ed if not already.

    :param grid grid0: reference grid
    :param grid grid1: source grid
    :param Block block: target Block in grid0
    :param BBox bbox: target BBox in grid0
    :param str resample: resample method, "nearest_neighbor", "bilinear", or "cubic"
    :return: grid (annotated)
    """

    if block is None and bbox is None:
        raise Exception("Need either block or bbox for cells_from()")
    elif block is None:
        block = cells_for(grid0, bbox)
    else:
        bbox = bbox_for(grid0, block)

    for grid in grid0, grid1:
        if not hasattr(grid, 'cols'):
            annotate_grid(grid)

    resample = {
        'nearest_neighbor': gdal.GRA_NearestNeighbour,
        'bilinear': gdal.GRA_Bilinear,
        'cubic': gdal.GRA_CubicSpline,
    }[resample]

    driver = gdal.GetDriverByName('MEM')
    outRaster = driver.Create(
        'noname', block.w, block.h,
        grid1.RasterCount, grid1.GetRasterBand(1).DataType)
    outRaster.SetGeoTransform((bbox.l, grid0.sizex, 0, bbox.t, 0, -grid0.sizey))
    outRaster.SetProjection(grid0.GetProjection())
    outRaster.GetRasterBand(1).SetNoDataValue(grid0.GetRasterBand(1).GetNoDataValue())
    outRaster.GetRasterBand(1).Fill(grid0.GetRasterBand(1).GetNoDataValue())
    gdal.ReprojectImage(grid1, outRaster, None, None, resample)
    annotate_grid(outRaster)
    return outRaster
class TileRunner(object):
    """TileRunner - provide tiles for iterating a grid
    """

    def __init__(self, grid, cols=1024, rows=None, overlap_cols=0, overlap_rows=None):
        """
        :param GDAL grid grid: grid to iterate over,
            will be annotate_grid()ed if not already
        :param int xsize: xsize of tiles
        :param int ysize: ysize of tiles, defaults to xsize
        """
        if not rows:
            rows = cols
        if overlap_rows is None:
            overlap_rows = overlap_cols
        self.grid = grid
        self.cols = cols
        self.rows = rows
        self.overlap_cols = overlap_cols
        self.overlap_rows = overlap_rows
        if not hasattr(grid, 'cols'):
            annotate_grid(grid)
        self.n_cols = int(grid.cols // cols)
        if self.n_cols * cols < grid.cols:
            self.n_cols += 1
        self.n_rows = int(grid.rows // rows)
        if self.n_rows * rows < grid.rows:
            self.n_rows += 1

    Tile = namedtuple("Tile", "tile_col tile_row block")

    def __str__(self):
        """__str__ - show properties
        """

        return "cols:{0.cols},rows:{0.rows},n_cols:{0.n_cols},n_rows:{0.n_rows}".format(self)
    def tiles(self):
        """tiles - generator that lists tiles for this grid as Blocks
        """

        for r in range(self.n_rows):
            for c in range(self.n_cols):
                yield self.tile_calc(c, r)
    def tile_calc(self, c, r, use_overlap=True):
        """tile_calc - calculate tile bounds, with or without overlap

        :param int c: *tile* column, not grid column
        :param int r: tile row
        :param bool use_overlap: use or ignore overlap
        :return: tile
        :rtype: Tile
        """

        if use_overlap:
            overlap_cols = self.overlap_cols
            overlap_rows = self.overlap_rows
        else:
            overlap_cols = 0
            overlap_rows = 0

        tc = max(0, c*self.cols-overlap_cols)
        tr = max(0, r*self.rows-overlap_rows)

        return self.Tile(tile_col=c, tile_row=r, block=Block(
            c=tc,
            w=int(min(self.cols+2*overlap_cols, self.grid.cols-tc)),
            r=tr,
            h=int(min(self.rows+2*overlap_rows, self.grid.rows-tr))
        ))

    def trim_datasource(self, datasource, tile):
        """trim_datasource - return a copy of datasource with the overlap edges trimmed off

        :param GDAL datasource: datasource to trim
        :return: GDAL datasource
        """

        overlap = self.tile_calc(tile.tile_col, tile.tile_row, use_overlap=True).block
        trimmed = self.tile_calc(tile.tile_col, tile.tile_row, use_overlap=False).block

        driver = gdal.GetDriverByName('MEM')
        outRaster = driver.Create(
            'noname', trimmed.w, trimmed.h,
            datasource.RasterCount, datasource.GetRasterBand(1).DataType)
        # FIXME: handle bands
        in_col = trimmed.c-overlap.c
        in_row = trimmed.r-overlap.r
        outRaster.SetGeoTransform((
            datasource.left+datasource.sizex*in_col,
            datasource.sizex, 0,
            datasource.top-datasource.sizey*in_row,
            0, -datasource.sizey))
        outRaster.SetProjection(datasource.GetProjection())
        outRaster.GetRasterBand(1).SetNoDataValue(datasource.GetRasterBand(1).GetNoDataValue())
        outRaster.GetRasterBand(1).Fill(datasource.GetRasterBand(1).GetNoDataValue())
        data = datasource.GetRasterBand(1).ReadAsArray(in_col, in_row, trimmed.w, trimmed.h)
        outRaster.GetRasterBand(1).WriteArray(data)
        annotate_grid(outRaster)
        return outRaster
class TestUtils(unittest.TestCase):
    """TestTest - describe class
    """

    tmp_dir = "temp_test_grids"

    @classmethod
    def setUpClass(cls):
        """setUp -
        """

        if not os.path.isdir(cls.tmp_dir):
            os.mkdir(cls.tmp_dir)

        TestGrid = namedtuple('TestGrid', "origx origy pixw pixh cols rows")

        ox = 357000
        oy = 5170000
        cls.grids = []

        driver = gdal.GetDriverByName('GTiff')
        for n, dim in enumerate([
                TestGrid(ox+2127.5, oy+5413.2, 25.5, 22.75, 1000, 1210),
                TestGrid(ox+127, oy+213, 30, 30, 1017, 767),
            ]):

            dat = np.zeros((dim.rows, dim.cols))
            for y in range(dim.rows):
                for x in range(dim.cols):
                    dat[y, x] = (10000*x+y)*(1-2*((x+y)%2))

            newRasterfn = os.path.join(cls.tmp_dir, "%s_%s.tif" % (dim.cols, dim.rows))
            outRaster = driver.Create(newRasterfn, dim.cols, dim.rows, 1, gdal.GDT_Int32)
            outRaster.SetGeoTransform((dim.origx, dim.pixw, 0, dim.origy, 0, -dim.pixh))
            outband = outRaster.GetRasterBand(1)
            outband.WriteArray(dat)
            outRasterSRS = osr.SpatialReference()
            outRasterSRS.ImportFromEPSG(26915)
            outRaster.SetProjection(outRasterSRS.ExportToWkt())
            outband.FlushCache()
            cls.grids.append(outRaster)
            annotate_grid(outRaster)
    @classmethod
    def tearDownClass(cls):
        """setUp -
        """

        print('\n%s can be deleted' % cls.tmp_dir)
    def test_bbox_for(self):
        """test_bbox_for - not much of a test, just a different formulation
        """

        block = Block(r=5, c=7, w=4, h=11)
        for grid in self.grids:
            ans = BBox(
                l=grid.left + grid.sizex * block.c,
                r=grid.left + grid.sizex * block.c + block.w * grid.sizex,
                t=grid.top - grid.sizey * block.r,
                b=grid.top - grid.sizey * block.r - block.h * grid.sizey,
            )

            bbox = bbox_for(grid, block)
            self.assertEqual(bbox, ans)
    def test_grid_match(self):
        """test_grid_match -
        """

        pairs = [
            # 1000_1210,   1017_767
            ( Block(c=202, r=462, w=1, h=1),
              Block(c=238, r=177, w=1, h=1),
              Block(c=202, r=462, w=1, h=1) ),

            ( Block(c=204, r=461, w=1, h=1),
              Block(c=240, r=176, w=1, h=1),
              Block(c=204, r=461, w=1, h=1) ),

            ( Block(c=816, r=1069, w=5, h=5),
              Block(c=760, r=637, w=5, h=4),
              Block(c=816, r=1069, w=6, h=5) ),

            ( Block(c=474, r=737, w=4, h=6),
              Block(c=470, r=386, w=3, h=4),
              Block(c=474, r=738, w=4, h=5) ),

            ( Block(c=474, r=737, w=4, h=7),
              Block(c=470, r=386, w=3, h=5),
              Block(c=474, r=738, w=4, h=6) ),

            ( Block(c=474, r=737, w=5, h=6),
              Block(c=470, r=386, w=4, h=4),
              Block(c=474, r=738, w=5, h=5) ),
        ]

        # dump polys
        outSHPfn = os.path.join(self.tmp_dir, "boxes.shp")
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(outSHPfn):
            shpDriver.DeleteDataSource(outSHPfn)
        outDataSource = shpDriver.CreateDataSource(outSHPfn)
        outLayer = outDataSource.CreateLayer(outSHPfn, geom_type=ogr.wkbPolygon)
        block = 'rcwh'
        box = 'bltr'
        for ch in block:
            fld = ogr.FieldDefn(ch, ogr.OFTInteger)
            outLayer.CreateField(fld)
        for ch in box:
            fld = ogr.FieldDefn(ch+'_', ogr.OFTReal)
            outLayer.CreateField(fld)
        fld = ogr.FieldDefn('n', ogr.OFTInteger)
        outLayer.CreateField(fld)
        fld = ogr.FieldDefn('g', ogr.OFTInteger)
        outLayer.CreateField(fld)
        for np, pair in enumerate(pairs):
            for n, one in enumerate(pair):
                featureDefn = outLayer.GetLayerDefn()
                outFeature = ogr.Feature(featureDefn)
                for ch in block:
                    outFeature.SetField(ch, getattr(one, ch))
                bbox = bbox_for(self.grids[n % 2], one)
                for ch in box:
                    outFeature.SetField(ch+'_', getattr(bbox, ch))
                outFeature.SetField('n', n)
                outFeature.SetField('g', np)
                ring = ogr.Geometry(ogr.wkbLinearRing)
                ring.AddPoint(bbox.l, bbox.b)
                ring.AddPoint(bbox.l, bbox.t)
                ring.AddPoint(bbox.r, bbox.t)
                ring.AddPoint(bbox.r, bbox.b)
                ring.AddPoint(bbox.l, bbox.b)
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ring)
                outFeature.SetGeometry(poly)
                outLayer.CreateFeature(outFeature)
        del outLayer

        for pair in pairs:
            c0 = pair[0]
            bbox0 = bbox_for(self.grids[0], c0)
            if 0:
                print("%s_%s" % (self.grids[0].cols, self.grids[0].rows))
                print(grid_annotations(self.grids[0]))
                print(c0)
                print(bbox0)

            c1 = pair[1]
            bbox1 = bbox_for(self.grids[1], c1)
            if 0:
                print()
                print("%s_%s" % (self.grids[1].cols, self.grids[1].rows))
                print(grid_annotations(self.grids[1]))
                print(c1)
                print(bbox1)
                print()

            c2 = pair[2]

            # these just test the test assumptions
            grid = self.grids[0]
            self.assertEqual(bbox0.l, grid.left + grid.sizex * c0.c)
            self.assertEqual(bbox0.b, grid.top - grid.sizey * c0.r - grid.sizey * c0.h)
            grid = self.grids[1]
            self.assertEqual(bbox1.l, grid.left + grid.sizex * c1.c)
            sep = bbox0.l-bbox1.l
            maxsep = max(self.grids[0].sizex, self.grids[1].sizex)
            self.assertTrue(abs(sep) < maxsep, (sep, maxsep))

            # this is the real test
            self.assertEqual(cells_for(self.grids[0], bbox0), c0)  # ~identity
            self.assertEqual(cells_for(self.grids[1], bbox1), c1)  # ~identity
            self.assertEqual(cells_for(self.grids[1], bbox0), c1)
            self.assertEqual(cells_for(self.grids[0], bbox1), c2)
    def test_grid_match_more(self):
        """test_grid_match - test cells_for() against cells_for_naive()
        """
        minx = min(i.left for i in self.grids)
        maxx = max(i.right for i in self.grids)
        miny = min(i.bottom for i in self.grids)
        maxy = max(i.top for i in self.grids)

        rngx = maxx - minx
        rngy = maxy - miny

        for grid in self.grids:
            for i in range(100000):
                x0 = minx + random.random() * rngx
                x1 = minx + random.random() * rngx
                y0 = miny + random.random() * rngy
                y1 = miny + random.random() * rngy
                if x0 > x1:
                    x0, x1 = x1, x0
                if y0 > y1:
                    y0, y1 = y1, y0
                box = BBox(l=x0, r=x1, b=y0, t=y1)
                self.assertEqual(cells_for(grid, box), cells_for_naive(grid, box))
    def test_TileRunner(self):
        """test_TileRunner -
        """

        for grid in self.grids:
            runner = TileRunner(grid, 200, 198)
            #D print(runner)
            #D print(grid_annotations(grid))
            c = r = n = 0
            for tile in runner.tiles():
                #D print(tile)
                n += 1
                if tile.r == 0:
                    c += tile.w
                if tile.c == 0:
                    r += tile.h
            # big enough
            self.assertTrue(c >= grid.cols, (c, grid.cols))
            self.assertTrue(r >= grid.rows, (r, grid.rows))
            # but no bigger
            self.assertTrue(c-tile.w < grid.cols, (c-tile.w, grid.cols))
            self.assertTrue(r-tile.h < grid.rows, (r-tile.h, grid.rows))
            self.assertEqual(n, runner.n_rows * runner.n_cols)
            self.assertTrue(n*tile.w*tile.h >= grid.cols*grid.rows)
def main():

    opt = make_parser().parse_args()

    grids = {}
    for name_path_resample in opt.grid:
        name, path = name_path_resample[:2]
        if len(name_path_resample) > 2:
            resample = name_path_resample[2]
        else:
            resample = 'nearest_neighbor'
        grid = gdal.Open(path)
        if not grid:
            raise Exception("Failed to open '%s' - '%s'" % (name, path))
        annotate_grid(grid)
        grid.resample = resample
        grids[name] = grid

    ref_grid = grids[opt.grid[0][0]]
    tr = TileRunner(ref_grid, opt.tile_cols, opt.tile_rows,
        overlap_cols=opt.overlap_cols, overlap_rows=opt.overlap_rows)
    driver = gdal.GetDriverByName('GTiff')

    tiles_only = []
    if opt.tiles_only:
        tiles_only = [tuple(map(int, i.split(','))) for i in opt.tiles_only.split()]
        print tiles_only

    for tile in tr.tiles():

        t_c, t_r = tile.tile_col, tile.tile_row
        if tiles_only and (t_c, t_r) not in tiles_only:
            continue

        print tile
        context = {
            'NP': np,
            'NoData': ref_grid.GetRasterBand(1).GetNoDataValue(),
        }
        # FIXME: handle bands
        block = tile.block
        context[opt.grid[0][0]] = ref_grid.GetRasterBand(1).ReadAsArray(
            block.c, block.r, block.w, block.h)

        for name_path_resample in opt.grid[1:]:
            name, path = name_path_resample[:2]
            cells_datasource = cells_from(
                ref_grid, grids[name], block=block, resample=grids[name].resample)
            context[name] = cells_datasource.GetRasterBand(1).ReadAsArray()

        for expr in opt.expr or []:
            name, expr = expr.split('=', 1)
            name = name.strip()
            context[name] = eval(expr, context)

        for output in opt.output:
            name = output[0]
            if len(output) > 1:
                path = os.path.join(output[1], 'tiles')
            else:
                path = os.path.join(opt.output_dir, name, 'tiles')
            cells_datasource.GetRasterBand(1).WriteArray(context[name])

            trimmed = tr.trim_datasource(cells_datasource, tile)

            if not os.path.exists(path):
                os.makedirs(path)
            fn = os.path.join(path, "%s_%04d_%04d.tif" % (
                name, t_c, t_r))
            print fn
            driver.CreateCopy(fn, trimmed)

    # now make VRT
    for output in opt.output:
        name = output[0]
        if len(output) > 1:
            path = os.path.join(output[1])
        else:
            path = os.path.join(opt.output_dir, name)
        os.system("gdalbuildvrt %s/%s.vrt %s/tiles/*.tif" % (
            path, name, path))
        if opt.tiles_merge:
            os.system("gdal_translate %s/%s.vrt %s/%s.tif" % (
                path, name, path, name))
if __name__ == '__main__':
    main()
