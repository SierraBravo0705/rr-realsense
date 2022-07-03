"""
Postprocessing of .bag-file containing depth stream into .csv-file containing
mean or median depth of chest area over time

collected with IntelRealSense SDK and D435 depth camera
possible decimation filter implemented to reduce computational time

specifications when recording with the Intel RealSense SDK:
# ONLY depth stream
# set ROI
# resolution: 848x480 (original resolution without decimation)
# framerate: 15fps
# Enable Auto Exposure
# mean intensity setpoint: 1536.000
# ROI as extracted from metadata
# decimation: ON/OFF (if OFF, decimation filter in this script can be used)

parameters that will need to be set by the user:
# if further decimation here desired: decimation parameter 2<=dec<=8
# which datasets to be read and from which storage location

created on 2022-07-03 12:33:26.424321
@author: Steffen Brinkmann
"""

## Set up environment
import numpy as np
import pyrealsense2 as rs
import csv
import matplotlib.pyplot as plt
import rr_algorithms as rra

## Decimation function
def get_decimation(depth_frame, dec):
    decimation = rs.decimation_filter()
    decimation.set_option(rs.option.filter_magnitude, dec)
    dec_depth_frame = decimation.process(depth_frame)
    # consider decimated dec_depth_frame just like depth frame as before
    # permits us to perform a lot more operations with it
    dec_depth_frame = dec_depth_frame.as_depth_frame()
    return dec_depth_frame

## Mean depth function (method 2)
def get_mean_depth(frame, ROI):
    # Get ROI
    ROI_top = ROI[0]
    ROI_bottom = ROI[1]
    ROI_left = ROI[2]
    ROI_right = ROI[3]
    # matrix depth_image containing depth value in mm for every pixel x, y
    depth_image = np.asanyarray(frame.get_data())
    depth_image_ROI = depth_image[ROI_top:ROI_bottom+1, ROI_left:ROI_right+1]   # caution: here indexes row, col; not x, y
    # substitute entries = 0 by nan
    depth_image_ROI = np.where(depth_image_ROI == 0.0, np.nan, depth_image_ROI)
    # mean of all not nan values in array
    mean_depth = np.nanmean(depth_image_ROI) # unit already mm
    return mean_depth

# Median depth function (method 2)
def get_median_depth(frame, ROI):
    # Get ROI
    ROI_top = ROI[0]
    ROI_bottom = ROI[1]
    ROI_left = ROI[2]
    ROI_right = ROI[3]
    # matrix depth_image containing depth value in mm for every pixel x, y
    depth_image = np.asanyarray(frame.get_data())
    depth_image_ROI = depth_image[ROI_top:ROI_bottom + 1,
                      ROI_left:ROI_right + 1]  # caution: here indexes row, col; not x, y
    # substitute entries = 0 by nan
    depth_image_ROI = np.where(depth_image_ROI == 0.0, np.nan, depth_image_ROI)
    # median of all not nan values in array
    median_depth = np.nanmedian(depth_image_ROI)  # unit already mm
    return median_depth

## Set parameters
prob = [1, 2, 4, 5, 6, 7, 8, 9] # probands: 1, 2, (3), 4, 5, 6, 7, 8, 9 so far
                                # leave out 3 as there are not all datasets for him
bpmPacs = [10, 15] # paced bpms: 10bpm, 15bpm
distance = [1, 2, 3] # distances: 1m, 2m, 3m
freq = [15] # sampling frequencies: 10fps (only RB), 15fps, 30fps
method = ['median'] # methods: 'mean', 'median'
                    # ONLY one at a time

pathC_bag = 'C:/Users/sbrin/Desktop/BA/Data/Measurements/' # adding PX later
pathC_csv = 'C:/Users/sbrin/Desktop/BA/Data/Processed/'

# magnitude of decimation, 2<=dec<=8, recommendation: dec=3
dec = 3

for p in prob:
    for bpmPac in bpmPacs:
        for dist in distance:
            for met in method:
                for f in freq:
                    paramSetC_csv, paramSetC_bag = rra.get_parameterC(p, bpmPac, dist, met, f, dec)
                    filenameC_bag = pathC_bag+'P'+str(p)+'/'+paramSetC_bag+'.bag'

                    ## Set up pipeline
                    # Create context object owning handles to all connected realsense devices
                    pipe = rs.pipeline()
                    # Configure streams/Create cfg object
                    cfg = rs.config()
                    # telling cfg that we will use recorded device from file
                    # by pipeline through playback
                    cfg.enable_device_from_file(filenameC_bag, repeat_playback=False)
                    # start streaming from file
                    profile = pipe.start(cfg)
                    # needed so frames are not dropped during processing:
                    playback = profile.get_device().as_playback()
                    playback.set_real_time(False)

                    ## Get ROI and dec_ROI from first frame
                    frame = pipe.wait_for_frames()
                    ROI_top = np.array(frame.get_frame_metadata(rs.frame_metadata_value.exposure_roi_top))
                    ROI_bottom = np.array(frame.get_frame_metadata(rs.frame_metadata_value.exposure_roi_bottom))
                    ROI_left = np.array(frame.get_frame_metadata(rs.frame_metadata_value.exposure_roi_left))
                    ROI_right = np.array(frame.get_frame_metadata(rs.frame_metadata_value.exposure_roi_right))
                    ROI = np.asarray([ROI_top, ROI_bottom, ROI_left, ROI_right], dtype=int)
                    # decimated ROI
                    dec_ROI_top = np.floor(ROI_top/dec)
                    dec_ROI_bottom = np.ceil(ROI_bottom/dec)
                    dec_ROI_left = np.ceil(ROI_left/dec)
                    dec_ROI_right = np.floor(ROI_right/dec)
                    dec_ROI = np.asarray([dec_ROI_top, dec_ROI_bottom, dec_ROI_left, dec_ROI_right], dtype=int)

                    ## Capture frames from depth stream in file and process captured frames
                    timestamp_set = np.array([])
                    depth_set = np.array([])
                    while True:
                        frame_present, frame = pipe.try_wait_for_frames()
                        # playback.pause() # might not be necessary
                        if not frame_present:
                            break
                        depth_frame = frame.get_depth_frame()
                        # getting array "timestamp_set"
                        timestamp = depth_frame.get_timestamp()
                        timestamp_set = np.append(timestamp_set, timestamp)
                        # getting array "depth_set" (with decimated depth frames)
                        dec_depth_frame = get_decimation(depth_frame, dec)
                        depth = get_median_depth(dec_depth_frame, dec_ROI)
                        depth_set = np.append(depth_set, depth)
                    pipe.stop()

                    ## Work timestamp_set and depth_set
                    # get reference for timestamp = 0 to be first frame of captured frameset
                    # each entry is the timestamp in ms
                    timestamp_set = timestamp_set - timestamp_set[0]
                    # set minimal distance in depth_set as reference to 0
                    depth_set = depth_set - np.amin(depth_set)

                    ## Save both arrays, timestamp_set and depth_set, into .csv-file
                    # set name of future csv-file
                    filenameC_csv = pathC_csv+paramSetC_csv+'.csv'
                    # open file in writing mode
                    with open(filenameC_csv, 'w+', newline='') as f:
                        # create new writer-object
                        header = ['timestamp', 'displacement']
                        wrtr = csv.DictWriter(f, fieldnames=header)
                        # write the desired header into csv-file
                        wrtr.writeheader()
                        # write our arrays row per row into csv-file
                        for k in range(0,timestamp_set.size):
                            wrtr.writerow({'timestamp': timestamp_set[k], 'displacement': depth_set[k]})

print('All given files processed')