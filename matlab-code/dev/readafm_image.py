import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from numpy import matlib
from scipy import linalg
import seaborn
fname = '/home/arnold/matlab/may12_WORKS.mi'
# fname = '/media/labserver/may12_newgrating.mi'
import ipdb



def split_mi_line(l):
    # all the keys are 14 lines long.
    # split() removes whitespace.
    thekey = l[0:14].split()[0]
    thedata = l[14:]
    return thekey, thedata


def meta2dict(metaDataLst):

    metaData = OrderedDict()
    metaData['global'] = dict()
    frameKey = 'global'
    for l in metaDataLst:
            thekey, thedata = split_mi_line(l)
            if thekey == 'bufferLabel':
                frameKey = thedata.split('\n')[0]
                metaData[frameKey] = dict()
            else:
                metaData[frameKey][thekey] = thedata
    return metaData


def grabBuffer(f, xpix, ypix):
    datmat = np.zeros([ypix, xpix])
    for yrow in range(0, ypix):
        for xcol in range(0, xpix):
            datmat[yrow, xcol] = int(f.readline())
    return datmat


def loadrawafm(fname):
    f = open(fname, 'rU')
    metaDatalst = []
    while True:
        l = f.readline()
        metaDatalst.append(l)
        if l.split()[0] == 'data':
            break

    metaData = meta2dict(metaDatalst)

    xpix = int(metaData['global']['xPixels'])
    ypix = int(metaData['global']['yPixels'])

    thedata = dict()
    for key in metaData.keys():
        if key != 'global':
            thedata[key] = grabBuffer(f, xpix, ypix)
    f.close()
    return metaData, thedata

############################################################
#                                                          #
#                 AFM data plotting helpers                #
#                                                          #
############################################################

def genlabels(arr):
    ticks = []
    for tick in arr:
        ticks.append(str(tick))
    return ticks


def normalizeAxLabels(npix, ax, xy, xy_max, nticks):

    axticks = np.linspace(0, npix, nticks)
    axticks_norm = (axticks/axticks[-1])*xy_max
    axlabs = genlabels(axticks_norm)
    if xy == 'x':
        axticks = ax.set_xticks(axticks)
        ax.set_xticklabels(axlabs)
    elif xy == 'y':
        axticks = ax.set_yticks(axticks)
        ax.set_yticklabels(axlabs, rotation=0)


def unit2num(unit):
    if unit == 'mm':
        num = 1e-3
    elif unit == 'um':
        num = 1e-6
    elif unit == 'nm':
        num = 1e-9
    return num

def flatten_1(datmat, order=1):
    """
    Fit a line to image data along the x-axis and subtract it. This is the
    best fit line to the data in EVERY row.
    """
    xpix, ypix = datmat.shape
    xs = np.linspace(0, xpix, xpix)
    xs = xs.reshape(xs.shape[0], -1)

    # [x1  1][m]  [z1]
    # [x2  1][b] =[z2]
    # [x3  1]     [z3]
    X = matlib.repmat(xs, ypix, 1)
    PHI = np.ones_like(X)
    for n in range(1, order+1):
        PHI = np.hstack((X**n, PHI))
    Z = datmat.reshape(datmat.shape[0]*datmat.shape[1], 1)

    coeffs = linalg.lstsq(PHI, Z)

    coeffs = coeffs[0]
    # best fit polynomial for every row of pixels
    z_fit = np.zeros_like(xs.T)
    for n in range(0, order+1):
        print n, coeffs[n][0]
        z_fit = (xs.T**(order - n)) * coeffs[n][0] + z_fit
    return datmat - z_fit
############################################################
#                                                          #
#                 The main script                          #
#                                                          #
############################################################
plt.ion()

metaData, thedata = loadrawafm(fname)


xlen = float(metaData['global']['xLength'])
ylen = float(metaData['global']['yLength'])
xpix = int(metaData['global']['xPixels'])
ypix = int(metaData['global']['yPixels'])

frame = 'Topography'

datmat = thedata[frame]
bufrange = float(metaData[frame]['bufferRange'])
lenunit = 1e3
datmat = (datmat*bufrange*lenunit)/2**15


datmat2 = flatten_1(datmat)
fig = plt.gcf()
ax = plt.axes()
# plt.set_xticklabels(
# seaborn.heatmap(datmat, cmap=plt.get_cmap('gist_gray'))
pc = ax.pcolor(datmat, cmap=plt.get_cmap('gist_gray'))
normalizeAxLabels(xpix, ax, 'x', 15, 5)
normalizeAxLabels(xpix, ax, 'y', 15, 5)
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
ax.set_title(frame)

cbar = fig.colorbar(pc)
cbar.set_label('nm')

fig.savefig('/home/arnold/matlab/grating.png', format='png', dpi=600)

# ###########################################
# # Deflection
# frame = 'Deflection'
# datmat = thedata[frame]
# bufrange = float(metaData[frame]['bufferRange'])
# lenunit = 1e0
# datmat = (datmat*bufrange*lenunit)/2**15

# fig = plt.figure()
# ax = plt.axes()


# pc = ax.pcolor(datmat, cmap=plt.get_cmap('gist_gray'))
# normalizeAxLabels(xpix, ax, 'x', 15, 5)
# normalizeAxLabels(xpix, ax, 'y', 15, 5)
# ax.set_xlabel('x [$\mu$m]')
# ax.set_ylabel('y [$\mu$m]')
# ax.set_title(frame)

# cbar = fig.colorbar(pc)
# cbar.set_label('V')

# plt.draw()
