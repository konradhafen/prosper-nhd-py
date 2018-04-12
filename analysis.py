import gispy as gis
from osgeo import gdal
import numpy as np
import geopandas

def testfunc():
    #get indices where PROSPER has values
    #prob = gis.getRasterBandAsArray(probpath, 1)
    prob = np.array([[0, 0.4, 0, 0, 0],
                     [0, 0.7, 0, 0, 0],
                     [0, 0.05, 0, 0, 0],
                     [0, 0.4, 0, 0, 0],
                     [0, 0.7, 0, 0, 0],
                     [0, 0.05, 0, 0, 0]
                     ])
    pgeot = (0,1,0,0,0,-1)
    ngeot = (0,5,0,0,0,-5)
    row, col = np.where(prob > 0)
    print row
    print col
    x, y = gis.coordinatesOfAddress(row, col, pgeot)
    idx = gis.linearIndexOfCoordinates(x, y, ngeot, 2, 2)
    print x, y, idx

def scpdsiDifference(prospath, scppath, scpcheckpath, outpath):
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
    write_vals = np.subtract(year_vals, check_vals)
    print "values calculated"
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
    gis.writeArrayAsRaster(outpath, out, prows, pcols, prob_geot, gis.getProjection(prospath),
                           nodata=-9999, nan=-9999, datatype=gdal.GDT_Float32, drivername = 'GTiff')
    print "done"

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
streamspath = "/outputs/shp/nhd_prosper_network.shp"
catpath = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_MEAN.tif"
sdpath = basepath + "/prosper/RawSPPs/SPP_STD.tif"
facpath = basepath + "/topo/fac/fac_albers83.tif"

probpath = basepath + "/prosper/RawSPPs/SPP_MEAN_wgs84.tif"
scpcheck = basepath + "/scpdsi/wateryear/scpdsi_checkyear_crb_wgs84.tif"
scpall = basepath + "/scpdsi/wateryear/scpdsi_wymean_crb_wgs84.tif"
outpath = basepath + "/scpdsi/difference/diff_MEAN.tif"

print "running"

scpdsiDifference(probpath, scpall, scpcheck, outpath)

