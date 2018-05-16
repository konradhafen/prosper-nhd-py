import gispy as gis
from osgeo import gdal
import numpy as np
import time
import geopandas as gpd
import pandas

def testfunc():
    #get indices where PROSPER has values
    #prob = gis.getRasterBandAsArray(probpath, 1)
    prob = np.array([[0, 0.4, 0, 0, 0, 0],
                     [0, 0.7, 0, 0, 0, 0],
                     [0, 0.05, 0, 0, 0, 0],
                     [0, 0.4, 0, 0, 0, 0],
                     [0, 0.7, 0, 0, 0, 0],
                     [0, 0.05, 0, 0, 0, 0],
                     [0, 0.4, 0, 0, 0, 0],
                     [0, 0.7, 0, 0, 0, 0],
                     [0, 0.05, 0, 0, 0, 0]
                     ])
    pgeot = (0.0, 1.0, 0.0, 0.0, 0.0, -1.0)
    ngeot = (0.0, 2.0, 0.0, 0.0, 0.0, -2.0)
    row, col = np.where(prob > 0)
    print row
    print col
    x, y = gis.coordinatesOfAddress(row, col, pgeot)
    idx = gis.linearIndexOfCoordinates(x, y, ngeot, 5, 3)
    print x, y, idx
    print np.arange(15).reshape((5,3))

def addIDtoNHD(nhdpath):
    gis.createIDField(nhdpath)


###############################################
#delta zonal stats for overall PROSPER average
###############################################
def deltaZonalStatsAverage(bufferpath, facpath, raspath, joinstats, fieldnames):
    zstats = gis.zonalStatisticsDelta(bufferpath, raspath, facpath, deltavalue=200.0, minvalue=125.0)
    print "zonal stats delta done"
    gis.joinZonalStatsToSHP(bufferpath, zstats, "fid", joinstats, fieldnames)
    print "join done"

def deltaZonalStatsByYear(bufferpath, facpath, rasbase, joinstats, fieldnames, fileend = ".tif"):
    for year in range(2004, 2017):
        yearabv = str(year)[2:]
        raspath = rasbase + str(year) + fileend
        zstats = gis.zonalStatisticsDelta(bufferpath, raspath, facpath, deltavalue=200.0, minvalue=125.0)
        #print zstats
        writenames = []
        for i in range(len(fieldnames)):
            writenames.append(yearabv + fieldnames[i])
        gis.joinZonalStatsToSHP(bufferpath, zstats, "fid", joinstats, writenames)
        print year, "done"

def scpdsiDifference(prospath, scppath, scpcheckpath, outpath):
    start = time.time()
    prob = gis.getRasterBandAsArray(prospath)
    prows, pcols, prob_geot = gis.getGeoTransformAndSize(prospath)
    scpyear = gis.getRasterBandAsArray(scppath)
    syrows, sycols, sy_geot = gis.getGeoTransformAndSize(scppath)
    scpcheck = gis.getRasterBandAsArray(scpcheckpath)
    scrows, sccols, sc_geot = gis.getGeoTransformAndSize(scpcheckpath)
    print "everything read"
    row, col = np.where(prob > 0)
    print "rows colums found"
    x, y = gis.coordinatesOfAddress(row, col, prob_geot)
    print "coordinates obtained"
    idx_scpyear = gis.linearIndexOfCoordinates(x, y, sy_geot, syrows, sycols)
    idx_scpcheck = gis.linearIndexOfCoordinates(x, y, sc_geot, scrows, sccols)
    idx_prob = gis.linearIndexOfCoordinates(x, y, prob_geot, prows, pcols)
    print "indices calculated"
    year_vals = np.add(np.take(scpyear, idx_scpyear), 100.0)
    check_vals = np.add(np.take(scpcheck, idx_scpcheck), 100.0)
    write_vals = np.where((year_vals > 0.0) & (check_vals > 0.0), np.subtract(year_vals, check_vals), -9999.0)
    print "values calculated"
    # write_vals = np.take(prob, idx_prob)
    out = np.empty(prob.shape, dtype=np.float32)
    print "output created"
    out.fill(-9999.0)
    print "output filled"
    np.put(out, idx_prob, write_vals)
    print "output values changed"
    del prob
    del scpyear
    del scpcheck
    del year_vals
    del check_vals
    print "time", time.time()-start
    gis.writeArrayAsRaster(outpath, out, prows, pcols, prob_geot, gis.getProjection(prospath),
                           nodata=-9999, nan=-9999, datatype=gdal.GDT_Float32, drivername='GTiff')
    print "done", time.time()-start

def scpdsiDifferenceByYear(prospath, scppath, scpcheckpath, outdir):
    prob = gis.getRasterBandAsArray(prospath)
    prows, pcols, prob_geot = gis.getGeoTransformAndSize(prospath)
    scpcheck = gis.getRasterBandAsArray(scpcheckpath)
    scrows, sccols, sc_geot = gis.getGeoTransformAndSize(scpcheckpath)
    syrows, sycols, sy_geot = gis.getGeoTransformAndSize(scppath)
    print "everything read"
    row, col = np.where(prob > 0)
    out = np.empty(prob.shape, dtype=np.float32)
    del prob
    print "output created"
    x, y = gis.coordinatesOfAddress(row, col, prob_geot)
    print "coordinates calculated"
    idx_prob = gis.linearIndexOfCoordinates(x, y, prob_geot, prows, pcols)
    idx_scpcheck = gis.linearIndexOfCoordinates(x, y, sc_geot, scrows, sccols)
    idx_scpyear = gis.linearIndexOfCoordinates(x, y, sy_geot, syrows, sycols)
    check_vals = np.add(np.take(scpcheck, idx_scpcheck), 100.0)
    print "indices calculated"

    for year in range(2004, 2017):
        print "year", year, "running"
        outpath = outdir + "/diff_" + str(year) + ".tif"
        band = 121 - (2016 - year)
        print "reading array for year"
        scpyear = gis.getRasterBandAsArray(scppath, band)
        year_vals = np.add(np.take(scpyear, idx_scpyear), 100.0)
        write_vals = np.where((year_vals > 0.0) & (check_vals > 0.0), np.subtract(year_vals, check_vals), -9999.0)
        del year_vals
        print "values calculated"

        out.fill(-9999.0)
        print "output filled"
        np.put(out, idx_prob, write_vals)
        print "writing array"

        gis.writeArrayAsRaster(outpath, out, prows, pcols, prob_geot, gis.getProjection(prospath),
                               nodata=-9999, nan=-9999, datatype=gdal.GDT_Float32, drivername='GTiff')


def getSCPDSIforPROSPERYears(scpdsipath, nhdnetwork):
    nhdds = gis.openOGRDataSource(nhdnetwork, 1)
    nhdlyr = nhdds.GetLayer(0)
    newfields = []
    for i in range(2004, 2017):
        newfields.append("scp_"+str(i))
    gis.createFields(nhdlyr, newfields)
    return None

def catSeedStats(catseedpath, statpath, shppath):
    zstats = gis.zonalStatistics_rasterZones(catseedpath, statpath, "GRIDCODE")
    fieldnames = ["max", "min", "mean", "sd", "med", "count"]
    writestats = ["max", "min", "mean", "sd", "median", "count"]
    gis.joinZonalStatsToSHP(shppath, zstats, "GRIDCODE", writestats, fieldnames)

def bufferStats(facpath, statpath, shppath, minindex, maxindex):
    writestats = ["max", "min", "mean", "sd", "median", "count"]
    for bufdist in range(20, 101, 20):
        print "running buffer", bufdist
        bufferpath = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/buf" + str(bufdist) + ".shp"
        fieldnames = ["max", "min", "mean", "sd", "med", "count"]
        for i in range(0, len(writestats)):
            fieldnames[i] = fieldnames[i] + str(bufdist)
        zstats = gis.zonalStatisticsDelta_methodtest(bufferpath, facpath, statpath, deltamin=minindex,
                                                     deltamax=maxindex, minvalue=125.0, idfield='GRIDCODE')
        print "zstats done"
        gis.joinZonalStatsToSHP(shppath, zstats, 'GRIDCODE', writestats, fieldnames)
        print "join done"

# def bufferStatsTest(facpath, statpath, shppath):
#     writestats = ["max", "min", "mean", "sd", "median", "count"]
#     for bufdist in range(40, 50, 20):
#         print "running buffer", bufdist
#         bufferpath = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/buf20.shp"
#         fieldnames = ["max", "min", "mean", "sd", "med", "count"]
#         for i in range(0, len(writestats)):
#             fieldnames[i] = fieldnames[i] + str(bufdist)
#         zstats = gis.zonalStatisticsDelta_methodtest(bufferpath, facpath, statpath, deltamin=-0.95,
#                                                      deltamax=0.9, minvalue=125.0)
#         print "zstats done"
#         gis.joinZonalStatsToSHP(shppath, zstats, 'fid', writestats, fieldnames)
#         print "join done"

basepath = "E:/konrad/Projects/usgs/prosper-nhd/data"
csvpath = basepath + "/outputs/csv/nhd_hr_buf20.shp"
bufferpath = basepath + "/outputs/shp/nhd_hr_buf20.shp"
bufferpath_cat = basepath + "/outputs/shp/nhd_hr_buf20_cat.shp"
bufferpath_prob = basepath + "/outputs/shp/nhd_hr_buf20_prob.shp"
bufferpath_dif = basepath + "/outputs/shp/nhd_hr_buf20_dif.shp"
streamspath = basepath + "/outputs/shp/nhd_stream_network_hr_subset.shp"
catpath = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_MEAN.tif"
sdpath = basepath + "/prosper/RawSPPs/SPP_STD.tif"
facpath = basepath + "/topo/fac/fac_albers83.tif"

prospath = basepath + "/prosper/RawSPPs/SPP_MEAN_wgs84.tif"
prosdir = basepath + "/prosper/RawSPPs"
scpcheck = basepath + "/scpdsi/wateryear/scpdsi_checkyear_crb_wgs84.tif"
scpall = basepath + "/scpdsi/wateryear/scpdsi_wymean_crb_wgs84.tif"
outpath = basepath + "/scpdsi/difference/diff_MEAN.tif"
outdir = basepath + "/scpdsi/difference"

joinstats_cat = ["majority", "mean", "max", "min"]
joinstats_prob = ["sd", "mean", "max", "min", "median"]
#joinstats_prob = ["count", "sd", "mean", "max", "min", "median"]
rasbase_cat = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_"
rasbase_prob = basepath + "/prosper/RawSPPs/SPP_"
rasbase_diff = basepath + "/scpdsi/difference/diff_"
raspath_cat = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_MEAN.tif"
raspath_prob = basepath + "/prosper/RawSPPs/SPP_MEAN.tif"
raspath_diff = basepath + "/scpdsi/difference/diff_MEAN_epsg5070.tif"
fieldnames_prob = ["_sd", "_mean", "_max", "_min", "_med"]
fieldnames_cat = ["_maj", "_mean", "_max", "_min"]
fieldnames_diff = ["_sdd", "_meand", "_maxd", "_mind", "_medd"]

print "running"
#testfunc()
#scpdsiDifference(prospath, scpall, scpcheck, outpath)
#scpdsiDifferenceByYear(prospath, scpall, scpcheck, outdir)
#addIDtoNHD(streamspath)
#deltaZonalStatsAverage(bufferpath_cat, facpath, raspath_diff, joinstats_prob, fieldnames_diff)

########################################################################################################################
#zonal stats on PROSPER categorical CIs
########################################################################################################################
#deltaZonalStatsByYear(bufferpath_cat, facpath, rasbase_cat, joinstats_cat, fieldnames_cat)

########################################################################################################################
#zonal stats on scPDSI difference between PROSPER and NHD
########################################################################################################################
# deltaZonalStatsByYear(bufferpath_dif, facpath, rasbase_diff, joinstats_prob, fieldnames_diff, "_epsg5070.tif")

########################################################################################################################
#zonal stats based on CatSeed grid
########################################################################################################################
catseedpath = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/catseed.tif"
shppath = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/flowline_nd-025-025.shp"
boisefac = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/fac.tif"
boisespp = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/spp.tif"
#gis.maskRasterWithValues(boisefac, 125, 10000000000)
catSeedStats(catseedpath, boisefac, shppath)

########################################################################################################################
#zonal stats based on different buffer distances
########################################################################################################################
#shppath = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/flowline-099-099.shp"
boisefac = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/fac.tif"
boisefaccopy = "E:/konrad/Projects/usgs/prosper-nhd/data/method_dev/nhd/mr/sfboise/facCopy.tif"
bufferStats(boisefac, boisefaccopy, shppath, -0.25, 0.25)