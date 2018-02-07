import percentile

#inPath = "E:/konrad/Projects/usgs/prosper-nhd/data/ppt/wateryear/ppt.tif"
inPath = "E:/konrad/Projects/usgs/prosper-nhd/data/scpdsi/wateryear/scpdsi_wymean.tif"
yearPath = "E:/konrad/Projects/usgs/prosper-nhd/data/quads/ras/year.tif"
percentilepath = "E:/konrad/Projects/usgs/prosper-nhd/data/scpdsi/wateryear/scpdsi_percentile.tif"
scorepath = "E:/konrad/Projects/usgs/prosper-nhd/data/scpdsi/wateryear/scpdsi_checkyear.tif"

percentile.run(inPath, yearPath, percentilepath, scorepath)

# multi = percentile.getRasterAsArray(pptPath)
# index = percentile.createBandIndex(yearPath, 1896, 2016)
# print index.shape, np.min(index), np.max(index)
# print multi.shape
# result = percentile.percentileMultiband(multi, index)
# print result
# print result.shape
# print np.nanmean(result), np.nanmin(result), np.nanmax(result)

# import numpy as np
#
# array3d = np.array([[[1.0, 2.0, 3.0],
#                      [1.0, 2.0, 3.0],
#                      [1.0, 2.0, 3.0]],
#                     [[1.0, 2.0, 3.0],
#                      [1.0, 2.0, 3.0],
#                      [1.0, 2.0, 3.0]],
#                     [[1.0, 2.0, 3.0],
#                      [1.0, 2.0, 3.0],
#                      [1.0, 2.0, 3.0]]
#                     ])
#
# array3d = np.array([[[1.0, 1.0, 1.0, 1.0],
#                      [1.0, 1.0, 1.0, 1.0],
#                      [1.0, 1.0, 1.0, 1.0]],
#                     [[2.0, 2.0, 2.0, 2.0],
#                      [2.0, 2.0, 2.0, 2.0],
#                      [2.0, 2.0, 2.0, 2.0]],
#                     [[3.0, 3.0, 3.0, 3.0],
#                      [3.0, 3.0, 3.0, 3.0],
#                      [3.0, 3.0, 3.0, 3.0]]
#                     ])
#
#
# array2d = np.array([[1.0, 2.0, 3.0, 4.0],
#                      [1.0, 2.0, 3.0, 4.0],
#                      [1.0, 2.0, 3.0, 4.0],
#                      [1.0, 2.0, 3.0, 4.0]
#                      ])
#
# array2d = np.array([[1.0, 2.0, 3.0],
#                     [1.0, 2.0, 3.0],
#                     [1.0, 2.0, 3.0]
#                     ])
#
# array2d = np.array([[3.0, 0.0, 0.0, 0.0],
#                     [0.0, 3.0, 0.0, 0.0],
#                     [0.0, 0.0, 3.0, 0.0]
#                     ])
#
# final = np.array([[33.33, 33.33, 66.66],
#                   [33.33, 33.33, 66.66],
#                   [33.33, 33.33, 66.66]
#                   ])
#
#
# zindex = percentile.createBandIndexNP(array2d, 1, 3)
# _, nR, nC = array3d.shape
# idx = nC*nR*zindex + nR*np.arange(nR)[:,None] + np.arange(nC)
# print idx
# print np.take(array3d, idx)
# print array3d.ravel()[idx]
# print zindex.choose(array3d)

# index = percentile.createBandIndexNP(array2d, 1, 3)
# print index
# result = percentile.percentileMultibandNP(array3d, index)
# print result