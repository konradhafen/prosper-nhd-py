import os
import numpy as np
from osgeo import gdal
from scipy import stats

def createBandIndex(rasterPath, minValue, maxValue):
    """

    Args:
        rasterPath:
        minValue:
        maxValue:

    Returns:

    """
    band = getRasterBandAsArray(rasterPath)
    print "year", band.shape
    print band[250:255, 600:605]
    if band is not None:
        array = band - minValue
        array = np.where((array < 0) | (array > (maxValue-minValue)), 0, array)
        print "index"
        print array[250:255, 600:605]
        return array
    else:
        return None

def createMask(rasterPath, minValue, maxValue, band=1):
    array = getRasterBandAsArray(rasterPath, band)
    return np.where((array < minValue) | (array > maxValue), 0, 1)

def getRasterAsArray(rasterPath):
    """
    Returns raster as numpy array

    Args:
        rasterPath: Path to input rater

    Returns:
        numpy array if file exists or None if it does not

    """
    if os.path.isfile(rasterPath):
        ds = gdal.Open(rasterPath)
        array = ds.ReadAsArray()
        ds = None
        return array
    else:
        return None

def getRasterBandAsArray(rasterPath, band = 1):
    """
    Returns raster band as 2d numpy array

    Args:
        rasterPath: Path to input raster

    Returns:
        numpy array if band exists or None if band does not exist

    """
    if os.path.isfile(rasterPath):
        ds = gdal.Open(rasterPath)
        bands = ds.RasterCount
        if band > 0 and band <= bands:
            array = ds.GetRasterBand(band).ReadAsArray()
            ds = None
            return array
        else:
            return None
    else:
        return None

def linearTake(values, indices):
    _, nR, nC = values.shape
    idx = nC*nR*indices + nC*np.arange(nR)[:, None] + np.arange(nC)
    return np.take(values, idx), idx

def maskArray(array, mask):
    return array*mask

def percentileMultiband(multi, index):
    """
    Calculate the percentile of a specified band of a multiband raster

    Args:
        multi: Path to multiband raster
        index: Numpy array containing the band index

    Returns:

    """
    mean = np.mean(multi, axis=0)
    print "mean", mean.shape
    print mean[250:255, 600:605]
    sd = np.std(multi, axis=0)
    print "sd", sd.shape
    print sd[250:255, 600:605]
    score, idx = linearTake(multi, index)
    print "index", index.shape
    print index[250:255, 600:605]
    print "idx", idx.shape
    print idx[250:255, 600:605]
    print "score", score.shape
    print score[250:255, 600:605]
    print "values", multi.shape
    print multi[79, 250:255, 600:605]
    result = stats.norm.cdf(score, loc=mean, scale=sd)*100
    print "result", result.shape
    print result[250:255, 600:605]
    return result, mean, sd, score, index

def writeArrayAsGTiff(path, array, rows, cols, geot, srs, nodata=-9999, nan=-9999, type=gdal.GDT_Float32):
    ds = gdal.GetDriverByName("GTiff").Create(path, xsize=cols, ysize=rows, bands=1, eType=type)
    ds.SetProjection(srs)
    ds.SetGeoTransform(geot)
    array = np.where((array==np.nan) | (array==nan), nodata, array)
    ds.GetRasterBand(1).WriteArray(array)
    ds.GetRasterBand(1).SetNoDataValue(nodata)
    ds.GetRasterBand(1).FlushCache()
    ds = None
    return None

def run(datapath, indexpath, outpath,outpath1,outpath2, outpath3, outpath4):
    multi = getRasterAsArray(datapath)
    index = createBandIndex(indexpath, 1896, 2016)
    mask = createMask(indexpath, 1896, 2016)
    print mask.shape
    result, mean, sd, score, index = percentileMultiband(multi, index)

    ds = gdal.Open(datapath)
    writeArrayAsGTiff(outpath, maskArray(result, mask), rows=result.shape[0], cols=result.shape[1],
                      geot=ds.GetGeoTransform(), srs=ds.GetProjection(), nan=0.0)
    ds = None