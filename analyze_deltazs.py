import gispy as gis
from osgeo import gdal
import numpy as np
import time

def testDeltaByPercent(buf, fac, outbase, percents):
    for percent in percents:
        print percent, "running"
        outpath = outbase + str(percent) + ".tif"
        gis.rasterZonesFromVector_delta(buf, fac, outpath, float(percent), minvalue = 125.0)

basepath = "E:/konrad/Projects/usgs/prosper-nhd/data"
facpath = basepath + "/topo/fac/fac_sfboise_5070.tif"
bufhr20path = basepath + "/outputs/sfboise/buf20_hr.shp"
bufmr20path = basepath + "/outputs/sfboise/buf20_mr.shp"
bufmr30path = basepath + "/outputs/sfboise/buf30_mr.shp"
bufmr40path = basepath + "/outputs/sfboise/buf40_mr.shp"

outbasemr = basepath + "/outputs/sfboise/rasterzones/mr40_"

percents = range(100, 1001, 100)
percents.insert(0, 50)
#percents = [100]

########################################################################################################################
# Test delta zonal stats for mr nhd
########################################################################################################################
print "start"
testDeltaByPercent(bufmr40path, facpath, outbasemr, percents)
print "done"