import gispy as gis
from osgeo import gdal
import numpy as np
import time
import geopandas

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

basepath = "E:/konrad/Projects/usgs/prosper-nhd/data"

bufferpath = basepath + "/nhd/MR/buf20_prob.shp"
streamspath = "/outputs/shp/nhd_stream_network_hr.shp"
catpath = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_MEAN.tif"
sdpath = basepath + "/prosper/RawSPPs/SPP_STD.tif"
facpath = basepath + "/topo/fac/fac_albers83.tif"

prospath = basepath + "/prosper/RawSPPs/SPP_MEAN_wgs84.tif"
prosdir = basepath + "/prosper/RawSPPs"
scpcheck = basepath + "/scpdsi/wateryear/scpdsi_checkyear_crb_wgs84.tif"
scpall = basepath + "/scpdsi/wateryear/scpdsi_wymean_crb_wgs84.tif"
outpath = basepath + "/scpdsi/difference/diff_MEAN.tif"
outdir = basepath + "/scpdsi/difference"

print "running"
#testfunc()
#scpdsiDifference(prospath, scpall, scpcheck, outpath)
scpdsiDifferenceByYear(prospath, scpall, scpcheck, outdir)

