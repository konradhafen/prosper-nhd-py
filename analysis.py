import gispy
import geopandas

basepath = "E:/konrad/Projects/usgs/prosper-nhd/data"

bufferpath = basepath + "/nhd/MR/buf20_prob.shp"
streamspath = "/outputs/shp/nhd_prosper_network.shp"
catpath = basepath + "/prosper/CategoricalSPPs/CategoricalSPP_MEAN.tif"
probpath = basepath + "/prosper/RawSPPs/SPP_MEAN.tif"
sdpath = basepath + "/prosper/RawSPPs/SPP_STD.tif"
facpath = basepath + "/topo/fac/fac_albers83.tif"

