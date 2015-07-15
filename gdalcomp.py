"""
gdalcomp.py - collect GDAL raster calc. utilities

With thanks to https://pcjericks.github.io/py-gdalogr-cookbook/

Terry Brown, Terry_N_Brown@yahoo.com, Tue Jul 14 10:29:06 2015
"""

import os
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
def cells_for(grid, bbox):
    """cells_for - identify cells in grid for rectangle in map units
    Return Block(col, row, width, height)

    :Parameters:
    - `grid`: grid containing cells, as annotated by annotate_grid()
    - `bbox`: BBox bounds in map units
    """
    
    # raise Exception("cells_for called")
    
    # mid-line row / col of grid nearest line

    if 0:
        lcol = max(0, min(grid.cols-1, int( (bbox.l - grid.left) / grid.sizex ) ) )
        trow = max(0, min(grid.rows-1, int( (grid.top - bbox.t) / abs(grid.sizey) ) ) )
        rcol = max(0, min(grid.cols-1, int( (bbox.r - grid.left) / grid.sizex ) ) )
        brow = max(0, min(grid.rows-1, int( (grid.top - bbox.b) / abs(grid.sizey) ) ) )

    if 0:
        lcol = max(0, min(grid.cols-1, int( (bbox.l - grid.left) / grid.sizex + 0.5 ) ) )
        trow = max(0, min(grid.rows-1, int( (grid.top - bbox.t) / abs(grid.sizey) + 0.5 ) ) )
        rcol = max(0, min(grid.cols-1, int( (bbox.r - grid.left) / grid.sizex + 0.5 ) ) )
        brow = max(0, min(grid.rows-1, int( (grid.top - bbox.b) / abs(grid.sizey) + 0.5 ) ) )

    # simple implementation
    print (bbox)
    print ((bbox.l - grid.left) / grid.sizex)
    lcol = int((bbox.l - grid.left) / grid.sizex)
    if (bbox.l - (grid.left + lcol * grid.sizex))/2. > 0.5:
        lcol += 1
    print ((bbox.r - grid.left) / grid.sizex)
    rcol = int((bbox.r - grid.left) / grid.sizex)
    if (bbox.r - (grid.left + rcol * grid.sizex))/2. < 0.5:
        rcol -= 1
    print((grid.top - bbox.t) / grid.sizey)
    trow = int((grid.top - bbox.t) / grid.sizey)
    if ((grid.top - trow * grid.sizey) - bbox.t)/2. > 0.5:
        trow += 1
    print((grid.top - bbox.b) / grid.sizey)
    brow = int((grid.top - bbox.b) / grid.sizey)
    if ((grid.top - brow * grid.sizey) - bbox.b)/2. < 0.5:
        brow -= 1
    print(lcol, rcol, trow, brow)
    # a faster implementation using half grid size steps may be possible

    # simple implementation
    print (bbox)
    print ((bbox.l - grid.left) / grid.sizex)
    
    lcol = (bbox.l - grid.left) / grid.sizex
    if lcol - int(lcol) > 0.5:
        lcol += 1
    lcol = int(lcol)
    
    print ((bbox.r - grid.left) / grid.sizex)
    rcol = (bbox.r - grid.left) / grid.sizex
    if rcol - int(rcol) < 0.5:
        rcol -= 1
    rcol = int(rcol)
    
    print((grid.top - bbox.t) / grid.sizey)
    trow = (grid.top - bbox.t) / grid.sizey
    if trow - int(trow) > 0.5:
        trow += 1
    trow = int(trow)
    print((grid.top - bbox.b) / grid.sizey)
    brow = (grid.top - bbox.b) / grid.sizey
    if brow - int(brow) < 0.5:
        brow -= 1
    brow = int(brow)
    print(lcol, rcol, trow, brow)
    # a faster implementation using half grid size steps may be possible
    

    assert bbox.t >= bbox.b
    assert rcol >= lcol, (rcol, lcol)
    assert brow >= trow, (brow, trow)

    return Block(
        c=lcol, 
        r=trow, 
        w=rcol-lcol+1, 
        h=brow-trow+1
    )
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
    def test_bbox_for(self,):
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
            # 1000_1210, 1017_767
            ( Block(c=202, r=462, w=1, h=1), Block(c=238, r=177, w=1, h=1), Block(c=202, r=462, w=1, h=1) ),
            ( Block(c=204, r=461, w=1, h=1), Block(c=240, r=176, w=1, h=1), Block(c=204, r=461, w=1, h=1) ),
            ( Block(c=816, r=1069, w=5, h=5), Block(c=760, r=637, w=5, h=4), Block(c=816, r=1069, w=6, h=5) ),
            ( Block(c=474, r=737, w=4, h=6), Block(c=470, r=386, w=3, h=4), Block(c=474, r=738, w=4, h=5) ),
            ( Block(c=474, r=737, w=4, h=7), Block(c=470, r=386, w=3, h=5), Block(c=474, r=738, w=4, h=6) ),
            ( Block(c=474, r=737, w=5, h=6), Block(c=470, r=386, w=4, h=4), Block(c=474, r=738, w=5, h=5) ),
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
                # outFeature.SetGeometry(point)
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
def main():
    pass

if __name__ == '__main__':
    main()
