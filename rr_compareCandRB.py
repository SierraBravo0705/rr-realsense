"""
Program to compare depth data from camera (C) lying in .csv-files to force data from
respiration belt (RB) lying in .csv-files as well.
Example data can be found under ....

Under "Set parameters" one can choose amongst which datasets to consider, filtering them
as desired.

A plot results of three statistical magnitudes, derived from comparing C to ground truth/RB.
They are PCC, Abs. Error and Rel. Error

created on 2022-07-03 12:33:26.424321
@author: Steffen Brinkmann
"""

## Set up environment
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.patches as mpatches
import numpy as np
import rr_algorithms as rra
import scipy.ndimage
import scipy.signal

## Set parameters
# modify according to what to consider during data evaluation: prob, bpmPacs, distance, method
prob = [1, 2, 4, 5, 6, 7, 8, 9] # probands: 1, 2, (3), 4, 5, 6, 7, 8, 9 exist so far
                                # leave out 3 as there are not all datasets for him
bpmPacs = [10, 15] # paced bpms: 10bpm, 15bpm
distance = [1, 2, 3] # distances 1m, 2m, 3m
freqs = [15] # sampling frequencies: 10fps (only RB), 15fps, 30fps
             # later on automatically newly derived --> freq
             # leave on 15 for example data is only of 15fps
method = ['median'] # methods: mean, median
                    # ONLY consider one at a time
postMedFilt = 14 # filter size for median filter on C data
                 # mean-method: 18 (max. 20)
                 # median-method: 14 (NOT HIGHER)
dec = 3 # leave on 3, decimation filter magnitude in rr_readC was kept 3 constantly

pathC = 'C:/Users/sbrin/Desktop/BA/Data/Processed/'
pathRB = 'C:/Users/sbrin/Desktop/BA/Data/Processed/'

timeScale = 1000 # 1000ms = 1s

# predefine necessary arrays
rAll = np.array([])
errAbsAll = np.array([]) # error in bpm
errRelAll = np.array([]) # error in percent

errAbs1m = np.array([])
errAbs2m = np.array([])
errAbs3m = np.array([])
errRel1m = np.array([])
errRel2m = np.array([])
errRel3m = np.array([])
r1m = np.array([])
r2m = np.array([])
r3m = np.array([])

errAbs10bpm = np.array([])
errAbs15bpm = np.array([])
errRel10bpm = np.array([])
errRel15bpm = np.array([])
r10bpm = np.array([])
r15bpm = np.array([])

probM = [3, 5, 7, 8, 9] # male probands, rest female
errAbsM = np.array([])
errAbsF = np.array([])
errRelM = np.array([])
errRelF = np.array([])
rM = np.array([])
rF = np.array([])

## get errorAbs, errorRel and r for chosen set of signals
for p in prob:
    for bpmPac in bpmPacs:
        for dist in distance:
            for met in method:
                for f in freqs:
                    paramSetC, _ = rra.get_parameterC(p, bpmPac, dist, met, f, dec)
                    paramSetRB, id = rra.get_parameterRB(p, bpmPac, dist)

                    filenameC = pathC+paramSetC+'.csv'
                    filenameRB = pathRB+paramSetRB+'.csv'

                    ## import csv-files with data from respiration belt (RB) and camera (C) respectively
                    tsC, dataC = rra.read_csvC(filenameC)
                    tsRB, dataRB = rra.read_csvRB(filenameRB, id)

                    # get frequency from C data
                    freq = np.int16(np.round(timeScale / tsC[1]))

                    ## align data from RB and C
                    # consider: recording of C and RB data started at different time instants respectively.
                    # Max. deviation is supposed to be <2s.
                    # Truncate shorter signal by 2s at start and end to make sure it lies within the longer signal
                    tsRBI, dataRBI = rra.interpolate(tsRB, dataRB, freq)
                    tsCI, dataCI = rra.interpolate(tsC, dataC, freq)

                    tsCal, dataCal, tsRBal, dataRBal = rra.align(tsCI, dataCI, tsRBI, dataRBI, freq)

                    ## median filter to reduce noise in C signal
                    # for Raw Median Method: size=14
                    # for Raw Mean Method: size=18
                    dataCfilt = scipy.ndimage.median_filter(dataCal, size=postMedFilt)

                    ## get RR from both signals
                    # approach: save ts of peaks, count number of peaks in set window (as big as possible) -->bpm
                    bpmC, bpmRB, errorAbs = rra.get_bpm(dataCfilt, dataRBal, bpmPac, freq)
                    errAbsAll = np.append(errAbsAll, errorAbs)
                    errorRel = errorAbs/bpmRB*100 # in %
                    errRelAll = np.append(errRelAll, errorRel)

                    ## get PCC comparing C to ground truth RB
                    r, _, _, _ = rra.pearsonr_ci(dataCfilt, dataRBal) # caution: using filtered dataCal!
                    rAll = np.append(rAll, r)

                    ## group gained statistical magnitude
                    if dist == 1:
                        errAbs1m = np.append(errAbs1m, errorAbs)
                        errRel1m = np.append(errRel1m, errorRel)
                        r1m = np.append(r1m, r)
                    elif dist == 2:
                        errAbs2m = np.append(errAbs2m, errorAbs)
                        errRel2m = np.append(errRel2m, errorRel)
                        r2m = np.append(r2m, r)
                    elif dist ==3:
                        errAbs3m = np.append(errAbs3m, errorAbs)
                        errRel3m = np.append(errRel3m, errorRel)
                        r3m = np.append(r3m, r)

                    if bpmPac == 10:
                        errAbs10bpm = np.append(errAbs10bpm, errorAbs)
                        errRel10bpm = np.append(errRel10bpm, errorRel)
                        r10bpm = np.append(r10bpm, r)
                    elif bpmPac == 15:
                        errAbs15bpm = np.append(errAbs15bpm, errorAbs)
                        errRel15bpm = np.append(errRel15bpm, errorRel)
                        r15bpm = np.append(r15bpm, r)

                    if p in probM:
                        errAbsM = np.append(errAbsM, errorAbs)
                        errRelM = np.append(errRelM, errorRel)
                        rM = np.append(rM, r)
                    elif p not in probM:
                        errAbsF = np.append(errAbsF, errorAbs)
                        errRelF = np.append(errRelF, errorRel)
                        rF = np.append(rF, r)

## Plot
# fontsizes
sSz = 10
mSz = 14
lSz = 16
plt.rc('font', size=lSz) # default text sizes
plt.rc('axes', titlesize=lSz) # fontsize of axes titles
plt.rc('axes', labelsize=lSz) # fontsize of x, y labels
plt.rc('xtick', labelsize=sSz) # fontsize of xtick labels
plt.rc('ytick', labelsize=sSz) # fontsize of ytick labels
plt.rc('figure', titlesize=lSz) # fontsize of figure title
plt.rc('figure', figsize=(10, 15)) # size of figures to be plotted big enough

# place xtick-labels (and respective data) the following way:
# M   F       1m  2m  3m      10bpm   15bpm       All
#  Sex         Distance            RR
labels = ['M', '\nSex', 'F', '1m', '2m\nDistance', '3m', '10bpm', '\nRR', '15bpm', 'All']
barPos = 1
barPosSex = [barPos, barPos*2]
barPosDist = [barPos*4, barPos*5, barPos*6]
barPosRR = [barPos*8, barPos*9]
barPosAll = [barPos*11]
barWidth = 0.5

labelPos = [barPos, barPos*1.5, barPos*2,
            barPos*4, barPos*5, barPos*6,
            barPos*8, barPos*8.5, barPos*9,
            barPos*11]

colSex = ['darkred', 'red']
colDist = ['darkgreen', 'seagreen','lightgreen']
colRR = ['darkblue', 'lightblue']
colAll = ['yellow']

# 3 subplots, one each for each magnitude, will be stacked one over another
fig = plt.figure()
gs = fig.add_gridspec(3)
(ax1, ax2, ax3) = gs.subplots()


# PCC
bpSex = ax1.boxplot([rM, rF], positions=barPosSex, widths=barWidth, patch_artist=True) # sex
bpDist = ax1.boxplot([r1m, r2m, r3m], positions=barPosDist, widths=barWidth, patch_artist=True) # dist
bpRR = ax1.boxplot([r10bpm, r15bpm], positions=barPosRR, widths=barWidth, patch_artist=True) # RR
bpAll = ax1.boxplot([rAll], positions=barPosAll, widths=barWidth, patch_artist=True) # All

ax1.yaxis.grid(True, color='lightgrey', which='major')
ax1.set(axisbelow=True,
        title='PCC',
        xlabel='Considered proband group',
        ylabel='PCC [-]',
        ylim=(0,1))
ax1.set_xticks(labelPos, labels)
ax1.tick_params(axis='x', bottom=False)
# ax1.set_xticklabels(sublabels, rotation=45)

# set facecols for all boxes, patch_artist needs to be set True!
for patch, color in zip(bpSex['boxes'], colSex):
    patch.set_facecolor(color)
for patch, color in zip(bpDist['boxes'], colDist):
    patch.set_facecolor(color)
for patch, color in zip(bpRR['boxes'], colRR):
    patch.set_facecolor(color)
for patch, color in zip(bpAll['boxes'], colAll):
    patch.set_facecolor(color)

# set median bar to be black, too
for median in bpSex['medians']:
    median.set_color('black')
for median in bpDist['medians']:
    median.set_color('black')
for median in bpRR['medians']:
    median.set_color('black')
for median in bpAll['medians']:
    median.set_color('black')


# Abs Err
bpSex = ax2.boxplot([errAbsM, errAbsF], positions=barPosSex, widths=barWidth, patch_artist=True) # sex
bpDist = ax2.boxplot([errAbs1m, errAbs2m, errAbs3m], positions=barPosDist, widths=barWidth, patch_artist=True) # dist
bpRR = ax2.boxplot([errAbs10bpm, errAbs15bpm], positions=barPosRR, widths=barWidth, patch_artist=True) # RR
bpAll = ax2.boxplot([errAbsAll], positions=barPosAll, widths=barWidth, patch_artist=True) # All

ax2.yaxis.grid(True, color='lightgrey', which='major')
ax2.set(axisbelow=True,
        title='Absolute Error',
        xlabel='Considered proband group',
        ylabel='Abs. error [breaths/min]',)
ax2.set_xticks(labelPos, labels)
ax2.tick_params(axis='x', bottom=False)
ax2.set_ylim([0,8])
# ax2.set_xticklabels(sublabels, rotation=45)

# set facecols for all boxes, patch_artist needs to be set True!
for patch, color in zip(bpSex['boxes'], colSex):
    patch.set_facecolor(color)
for patch, color in zip(bpDist['boxes'], colDist):
    patch.set_facecolor(color)
for patch, color in zip(bpRR['boxes'], colRR):
    patch.set_facecolor(color)
for patch, color in zip(bpAll['boxes'], colAll):
    patch.set_facecolor(color)

# set median bar to be black, too
for median in bpSex['medians']:
    median.set_color('black')
for median in bpDist['medians']:
    median.set_color('black')
for median in bpRR['medians']:
    median.set_color('black')
for median in bpAll['medians']:
    median.set_color('black')


# Rel Err
bpSex = ax3.boxplot([errRelM, errRelF], positions=barPosSex, widths=barWidth, patch_artist=True) # sex
bpDist = ax3.boxplot([errRel1m, errRel2m, errRel3m], positions=barPosDist, widths=barWidth, patch_artist=True) # dist
bpRR = ax3.boxplot([errRel10bpm, errRel15bpm], positions=barPosRR, widths=barWidth, patch_artist=True) # RR
bpAll = ax3.boxplot([errRelAll], positions=barPosAll, widths=barWidth, patch_artist=True) # All

ax3.yaxis.grid(True, color='lightgrey', which='major')
ax3.set(axisbelow=True,
        title='Relative Error',
        xlabel='Considered proband group',
        ylabel='Rel. error [%]',)
ax3.set_xticks(labelPos, labels)
ax3.tick_params(axis='x', bottom=False)
ax3.set_ylim([0,70])
# ax3.set_xticklabels(sublabels, rotation=45)

# set facecols for all boxes, patch_artist needs to be set True!
for patch, color in zip(bpSex['boxes'], colSex):
    patch.set_facecolor(color)
for patch, color in zip(bpDist['boxes'], colDist):
    patch.set_facecolor(color)
for patch, color in zip(bpRR['boxes'], colRR):
    patch.set_facecolor(color)
for patch, color in zip(bpAll['boxes'], colAll):
    patch.set_facecolor(color)

# set median bar to be black, too
for median in bpSex['medians']:
    median.set_color('black')
for median in bpDist['medians']:
    median.set_color('black')
for median in bpRR['medians']:
    median.set_color('black')
for median in bpAll['medians']:
    median.set_color('black')


# Hide x-labels and tick labels for all but bottom plot
for ax in (ax1, ax2, ax3):
    ax.label_outer()

fig.tight_layout()
plt.show()

# calculate medians for my own knowledge
rMedM = np.median(rM)
rMedF = np.median(rF)
rMed1m = np.median(r1m)
rMed2m = np.median(r2m)
rMed3m = np.median(r3m)
rMed10bpm = np.median(r10bpm)
rMed15bpm = np.median(r15bpm)
rMedAll = np.median(rAll)

errAbsMedM = np.median(errAbsM)
errAbsMedF = np.median(errAbsF)
errAbsMed1m = np.median(errAbs1m)
errAbsMed2m = np.median(errAbs2m)
errAbsMed3m = np.median(errAbs3m)
errAbsMed10bpm = np.median(errAbs10bpm)
errAbsMed15bpm = np.median(errAbs15bpm)
errAbsMedAll = np.median(errAbsAll)

errRelMedM = np.median(errRelM)
errRelMedF = np.median(errRelF)
errRelMed1m = np.median(errRel1m)
errRelMed2m = np.median(errRel2m)
errRelMed3m = np.median(errRel3m)
errRelMed10bpm = np.median(errRel10bpm)
errRelMed15bpm = np.median(errRel15bpm)
errRelMedAll = np.median(errRelAll)

print('Operation terminated successfully')