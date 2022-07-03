"""
Algorithms used in the programs
rr_readC.py
rr_compareCandRB

created on 2022-07-03 12:33:26.424321
@author: Steffen Brinkmann
"""
import numpy as np
import scipy.stats
import scipy.signal

def get_parameterC(p, bpmPac, dist, met, freq, dec):
    paramSetC = str(bpmPac)+'bpm_'+str(dist)+'m_'+str(freq)+'fps_'+str(dec)+'dec_'+met+'_'+'prob'+str(p)+'_C'
    paramSetC_bag = str(bpmPac)+'bpm_'+str(dist)+'m_'+str(freq)+'fps_'+'prob'+str(p)
    return paramSetC, paramSetC_bag

def get_parameterRB(p, bpmPac, dist):
    # id = 0: 15bpm_1m_10fps_probX tsRB+dataRB
    # id = 1: 15bpm_2m_10fps_probX tsRB+dataRB
    # id = 2: 15bpm_3m_10fps_probX tsRB+dataRB
    # id = 3: 10bpm_1m_10fps_probX tsRB+dataRB
    # id = 4: 10bpm_2m_10fps_probX tsRB+dataRB
    # id = 5: 10bpm_3m_10fps_probX tsRB+dataRB
    if bpmPac == 15:
        if dist == 1:
            id = 0
        elif dist == 2:
            id = 1
        else:
            id = 2
    elif bpmPac ==10:
        if dist == 1:
            id = 3
        elif dist == 2:
            id = 4
        else:
            id = 5

    paramSetRB = 'prob'+str(p)+'_RBnew'

    return paramSetRB, id

def read_csvC(filenameC):
    with open (filenameC) as fC:
        arrayC = np.loadtxt(filenameC, delimiter=',', skiprows=1)
    tsC = arrayC[:,0] # in ms
    dataC = arrayC[:,1] # depth set in mm
    dataC = (dataC - np.mean(dataC))*-1 # dataC and dataRB act in opposite ways:
                                        # Exhl.: more depth, less force (breast smaller)
                                        # Inhl.: less depth, more force (breast wider)
    return tsC, dataC

def read_csvRB(filenameRB, id):
    ids = [0, 3, 6, 9, 12, 15]
    # col0+col1: 15bpm_1m_10fps_probX tsRB+dataRB
    # col3+col4: 15bpm_2m_10fps_probX tsRB+dataRB
    # col6+col7: 15bpm_3m_10fps_probX tsRB+dataRB
    # col9+col10: 10bpm_1m_10fps_probX tsRB+dataRB
    # col12+col13: 10bpm_2m_10fps_probX tsRB+dataRB
    # col15+col16: 10bpm_3m_10fps_probX tsRB+dataRB
    with open (filenameRB) as fRB:
        arrayRB = np.genfromtxt(filenameRB, delimiter=';', skip_header=2, usecols=(ids[id], ids[id]+1))
    tsRB = arrayRB[:,0]*1000 # s-->ms
    dataRB = arrayRB[:,1] # force set in N
    # removing possible nan-entries
    tsRB = tsRB[~(np.isnan(tsRB))]
    dataRB = dataRB[~(np.isnan(dataRB))]
    dataRB = dataRB - np.mean(dataRB)

    return tsRB, dataRB

def interpolate(ts, data, freq, timeScale=1000):
    tsAligned = np.array(ts) - ts[0]
    timeStep = timeScale/freq
    tsCount = int(tsAligned[-1] / timeStep)
    tsMax = tsCount * timeStep
    tsNew = np.linspace(tsAligned[0], tsMax, tsCount+1)
    dataNew = np.interp(tsNew, tsAligned, data)
    return tsNew, dataNew

def pearsonr_ci(x,y,alpha=0.01):
    '''
    Calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    See: https://zhiyzuo.github.io/Pearson-Correlation-CI-in-Python/
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default, here 0.01
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''

    r, p = scipy.stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = scipy.stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

def align(tsCI, dataCI, tsRBI, dataRBI, freq, timeScale = 1000):
    '''
    Precondition: tsC, tsRB same sample steps, e.g. 0, 666.666... ms for 15fps (cf. interpolation)
    Shift of shorter dataset x-wise only discretely by prementioned steps
    :param tsCI: interpolated
    :param dataCI: interpolated
    :param tsRBI: interpolated
    :param dataRBI: interpolated
    :param freq: higher freq (given before interpolation)
    :return: tsCal, dataCal, tsRBal, dataRBal
    '''
    timeStep = timeScale/freq
    cropL = np.int16(np.ceil(3000/timeStep)) # get crop length ~3s:
    deltaSize = np.abs(tsCI.size - tsRBI.size)

    r_Set = np.zeros(deltaSize + 2*cropL)
    # p_Set = np.zeros(deltaSize + 2*cropL)
    # ci_Set = np.zeros((deltaSize + 2*cropL, 2))

    if tsCI.size <= tsRBI.size:
        dataShort = dataCI[cropL:-cropL]

        for i in range(0, deltaSize + 2*cropL):
            r, _, _, _ = pearsonr_ci(dataRBI[i:dataRBI.size - (deltaSize + 2*cropL) +i], dataShort) # comparison here
            r_Set[i] = r
            # p_Set[i] = p
            # ci_Set[i] = [lo, hi]

        idMax = np.argmax(r_Set) #position of highest correlation coefficient
        tsCI += tsRBI[idMax] - tsCI[cropL] # shift of ts of shorter signal

    else:
        dataShort = dataRBI[cropL:-cropL]

        for i in range(0, deltaSize + 2*cropL):
            r, _, _, _ = pearsonr_ci(dataCI[i:dataCI.size - (deltaSize + 2 * cropL) + i],dataShort)  # comparison here
            r_Set[i] = r
            # p_Set[i] = p
            # ci_Set[i] = [lo, hi]

        idMax = np.argmax(r_Set)
        tsRBI += tsCI[idMax] - tsRBI[cropL]

    # perform alignment and crop both signals to 56s
    tsL = 56000
    if tsCI[0] < tsRBI[0]:
        index = np.argmin(np.abs(tsCI - tsRBI[0]))
        tsCI = np.round(tsCI[index:] - tsRBI[0], decimals=6)
        tsRBI = np.round(tsRBI - tsRBI[0], decimals=6)
        dataCI = dataCI[index:]

    else:
        index = np.argmin(np.abs(tsRBI - tsCI[0]))
        tsRBI = np.round(tsRBI[index:] - tsCI[0], decimals=6)
        tsCI = np.round(tsCI - tsCI[0], decimals=6)
        dataRBI = dataRBI[index:]

    tsCal = tsCI[tsCI <= tsL+timeStep]  # only works with np.arrays
    tsRBal = tsRBI[tsRBI <= tsL+timeStep]
    dataCal = dataCI[:len(tsCal)]
    dataRBal = dataRBI[:len(tsRBal)]

    return tsCal, dataCal, tsRBal, dataRBal

def get_bpm(dataC, dataRB, bpmPac, freq, timeScale=1000):
    '''
    Core consists in memorizing timestamps of RR-peaks inside the given datasets.
    Then time differences between peaks are calculated and averaged --> period time
    --> bpm
    Error is difference between bpmC and bpmRB.

    bpmPac - bpm as specified during paced breathing
    :param dataC:
    :param dataRB:
    :param bpmPac:
    :param freq:
    :param timeScale:
    :return: bpmC, bpmRB, error
    '''
    beatPac = 60*timeScale/bpmPac
    timeStep = timeScale/freq
    spbPac = np.int16(np.round(beatPac/timeStep)) # number of timestamps per paced bpm

    # width = np.int16(np.round(spbPac/4)) # comment out for median method
    distance = np.int16(np.round(spbPac * 0.8))
    peaksRB, _ = scipy.signal.find_peaks(dataRB, distance=distance)
    peaksC, _ = scipy.signal.find_peaks(dataC, distance=distance)

    tsdifRB = np.zeros(np.size(peaksRB)-1)
    tsdifC = np.zeros(np.size(peaksC)-1)
    for i in range(0, peaksRB.size - 1):
        tsdifRB[i] = peaksRB[i + 1] - peaksRB[i]
    for j in range(0, peaksC.size - 1):
        tsdifC[j] = peaksC[j + 1] - peaksC[j]

    spbRB = np.mean(tsdifRB)
    spbC = np.mean(tsdifC)

    beatRB = spbRB * timeStep
    beatC = spbC * timeStep

    bpmRB = (60*timeScale)/beatRB
    bpmC = (60*timeScale)/beatC

    error = np.abs(bpmC-bpmRB)

    return bpmC, bpmRB, error