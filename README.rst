gdalcomp.py
===========

a command line tool / python library for raster calculation with GDAL
=====================================================================

Compared to gdal_calc.py
------------------------

gdal_calc.py is very useful and inspired gdalcomp.py.

 - gdalcomp does on the fly (OTF) reprojection as needed, raster
   projections and resolutions can differ
 - gdalcomp uses tile based calculation for low memory footprint and
   parallel or distribtuted processing
 - gdalcomp can overlap and then trim tiles to eliminate edge effects
   for methods that consider neighboring cells
 - gdalcomp is easy to integrate with your own code, so you're just
   doing numpy array operations and not worrying about projections,
   grid alignment, etc. etc.

Examples
--------

Cut two tiles out of a dataset for experimentation, development, etc.::

    python gdalcomp.py \
      --grid agrid /path/to/agrid \
      --tiles-only "1,1 1,2" \
      --output-dir /tmp/test \
      --output agrid

creates ``/tmp/test/tiles/agrid_0001_0001.tif``, ``/tmp/test/tiles/agrid_0001_0002.tif`` and ``/tmp/test/agrid.vrt``.

python gdalcomp.py \
      --grid agrid /mnt/usr1/scratch/marsch/marschner \
      --tiles-only "1,1 1,2" \
      --output-dir /tmp/test \
      --output agrid:Float32 \
      --setup "k = GC.hemi_kernel(7)" \
      --do "agrid = GC.apply_kernel(agrid, k, agrid_NoData)" \
      --overlap-cols 5

python gdalcomp.py \
      --grid agrid /mnt/usr1/scratch/marsch/marschner \
      --output-dir /tmp/test \
      --output agrid:Float32 \
      --setup "k = GC.hemi_kernel(5)" \
      --do "agrid = GC.apply_kernel(agrid, k, agrid_NoData)" \
      --overlap-cols 5 \
      --blocks 1 6

python gdalcomp.py \
      --grid agrid /mnt/usr1/scratch/marsch/marschner \
      --output-dir /tmp/test \
      --output agrid \
      --cpus 8
