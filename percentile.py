import os
import numpy as np
from osgeo import gdal
from scipy import stats
import gisops

def createBandIndex(rasterPath, minValue, maxValue):
    """

    Args:
        rasterPath:
        minValue:
        maxValue:

    Returns:

    """
    band = getRasterBandAsArray(rasterPath)
    if band is not None:
        array = band - minValue
        array = np.where((array < 0) | (array > (maxValue-minValue)), 0, array)
        return array
    else:
        return None

def createMask(rasterPath, minValue, maxValue, band=1):
    """
    Create a mask from a raster band

    Args:
        rasterPath: path to raster
        minValue: minimum data value
        maxValue: maximum data value
        band: raster band to use (default: 1)

    Returns:
        Integer array with a value of 0 where the input raster band is < minValue or > maxValue

    """
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
    """
    Get 2d array of band values from a multiband raster of shape (bands, rows, columns)
    Args:
        values: input 3d array containing data values
        indices: 2d array containing indices to bands in the value raster

    Returns:
        2d array where the value in the array corresponds the band value from indices

    """
    _, nR, nC = values.shape
    idx = nC*nR*indices + nC*np.arange(nR)[:, None] + np.arange(nC) #convert 2d indices to linear indices
    return np.take(values, idx), idx

def maskArray(array, mask, nodata=-9999):
    """
    Replace array values with a no data value where a mask is false (0)

    Args:
        array: input array to be masked
        mask: boolean array or integer array with values of 0 (false) and 1 (true)
        nodata: value to write where mask is false (default: -9999)

    Returns:
        array with nodata where mask is false

    """
    return np.where(mask, array, nodata)

def percentileMultiband(multi, index):
    """
    Calculate the percentile of a specified value at a position in a multiband raster

    Args:
        multi: Path to multiband raster
        index: Numpy array containing the band index

    Returns:
        2d arrays result (percentile of value based on all bands), and score (value of requested band)

    """
    mean = np.mean(multi, axis=0) #mean of all bands at each row,col
    sd = np.std(multi, axis=0) #standard deviation of all bands at each row,col
    score, idx = linearTake(multi, index) #value of index band at each row,col
    result = stats.norm.cdf(score, loc=mean, scale=sd)*100 #percentile
    return result, score

def writeArrayAsGTiff(path, array, rows, cols, geot, srs, nodata=-9999, nan=-9999, type=gdal.GDT_Float32):
    """
    Write array to a GeoTiff raster

    Args:
        path: output file for raster
        array: array containing data
        rows: number of rows in array
        cols: number of columns in array
        geot: affine geotransformation for the output raster
        srs: spatial reference for the output raster
        nodata: no data value for the output raster
        nan: value in array that should be written as nodata
        type: gdal data type of output raster (default: GDT_Float32)

    Returns:
        None

    """
    ds = gdal.GetDriverByName("GTiff").Create(path, xsize=cols, ysize=rows, bands=1, eType=type)
    ds.SetProjection(srs)
    ds.SetGeoTransform(geot)
    array = np.where((array==np.nan) | (array==nan), nodata, array)
    ds.GetRasterBand(1).WriteArray(array)
    ds.GetRasterBand(1).SetNoDataValue(nodata)
    ds.GetRasterBand(1).FlushCache()
    ds = None
    return None

def percentileOfMultibandIndex(datapath, index, percentilepath, scorepath=None, mask = None):
    """
    For a multiband raster, calculates the percentile of a band value at each row,col

    Args:
        datapath: input path to multiband raster
        index: 2d array with each row,col containing an index to a band in data path
        percentilepath: output path for 2d raster of percentile values
        scorepath: ouput path for 2d raster of the scores of each band index (optional)
        maskpath: boolean array to mask outputs (optional)

    Returns:
        None

    """
    multi = getRasterAsArray(datapath)
    if index is not None and index.shape == multi.shape[1:]:
        result, score = percentileMultiband(multi, index)
        if mask is None or mask.shape != result.shape:
            mask = np.ones((result.shape))

        ds = gdal.Open(datapath)
        writeArrayAsGTiff(percentilepath, maskArray(result, mask), rows=result.shape[0], cols=result.shape[1],
                          geot=ds.GetGeoTransform(), srs=ds.GetProjection())
        if scorepath is not None:
            writeArrayAsGTiff(scorepath, maskArray(score, mask), rows=result.shape[0], cols=result.shape[1],
                              geot=ds.GetGeoTransform(), srs=ds.GetProjection())
        ds = None
    else:
        print "problem with input index array", index.shape, multi.shape[1:]

    return None
