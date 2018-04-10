import gispy as gis
import geopandas

basepath = "E:/konrad/Projects/usgs/prosper-nhd/data"

bufferpath = basepath + "/nhd/MR/buf20_prob.shp"
streamspath = "/outputs/shp/nhd_prosper_network.shp"
catpath = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_MEAN.tif"
probpath = basepath + "/prosper/RawSPPs/SPP_MEAN.tif"
sdpath = basepath + "/prosper/RawSPPs/SPP_STD.tif"
facpath = basepath + "/topo/fac/fac_albers83.tif"

def getSCPDSIforPROSPERYears(scpdsipath, nhdnetwork):
    nhdds = gis.openOGRDataSource(nhdnetwork, 1)
    nhdlyr = nhdds.GetLayer(0)
    newfields = []
    for i in range(2004, 2017):
        newfields.append("scp_"+str(i))
    gis.createFields(nhdlyr, newfields)
    return None

