import gisops
from osgeo import gdal, ogr, osr

"""
This section is for merging NHD and PROSPER data within CRB
"""
basepath = "E:/konrad/Projects/usgs/prosper-nhd/data"
#crbpoly = basepath + "/quads/Quads_CRB.shp"
crbquads = basepath + "/quads/Quads_CRB_stats.shp"
ppt = basepath + "/ppt/wateryear/ppt_percentile.tif"
scpdsi = basepath + "/scpdsi/wateryear/scpdsi_checkyear_NA83.tif"
pptcrb = basepath + "/ppt/wateryear/ppt_percentile_crb.tif"
scpdsicrb = basepath + "/scpdsi/wateryear/scpdsi_checkyear_NA83_crb.tif"
gdalwarp = "C:/Program Files/GDAL/gdalwarp"

# gisops.clipRasterWithPolygon(scpdsi, crbquads, scpdsicrb)
# print "scpdsi clipped"
# gisops.clipRasterWithPolygon(ppt, crbquads, pptcrb)
# print "ppt clipped"

from rasterstats import zonal_stats
stats = zonal_stats(crbquads, scpdsicrb, geojson_out=True)
gisops.joinZonalStatsToSHP(crbquads, stats, "id", ["mean"], ["pdsi_mean"], ogr.OFTReal)
print "pdsi zonal stats done"
stats = zonal_stats(crbquads, pptcrb, geojson_out=True)
gisops.joinZonalStatsToSHP(crbquads, stats, "id", ["mean"], ["ppt_mean"], ogr.OFTReal)
print "ppt zonal stats done"


"""
This section is for calculating mean, sd, and percentile for scpdsi and ppt
"""
# #inPath = "E:/konrad/Projects/usgs/prosper-nhd/data/ppt/wateryear/ppt.tif"
# inPath = "E:/konrad/Projects/usgs/prosper-nhd/data/scpdsi/wateryear/scpdsi_wymean.tif"
# yearPath = "E:/konrad/Projects/usgs/prosper-nhd/data/quads/ras/year.tif"
# percentilepath = "E:/konrad/Projects/usgs/prosper-nhd/data/scpdsi/wateryear/scpdsi_percentile_2.tif"
# scorepath = "E:/konrad/Projects/usgs/prosper-nhd/data/scpdsi/wateryear/scpdsi_checkyear_2.tif"
#
# index = gisops.createBandIndex(yearPath, 1896, 2016)
# mask = gisops.createMask(yearPath, 1896, 2016)
# gisops.percentileOfMultibandIndex(inPath, index, percentilepath, scorepath, mask)

