# rr-realsense
Remote system in python to extract the respiratory rate from depth videos recorded with the D435 Intel RealSense Depth Camera. Evaluation of recorded data against ground truth data included.

The system can work data recorded in the following setting:
  -probands sitting upright, facing depth camera
  -chest as ROI
  -indoors, well lit
  -chair of 0.5m standard height, camera on tripod on height of 1.1m (range 0.9m-1.2m or even more possible) with horizontal view on ROI
  -10bpm or 15bpm paced
  -three distances: 1m (min.), 2m, 3m (max.)
  -two paced respiratory rates (RRs): 10bpm, 15bpm
    ->without pacing: suppose generic RR (12bpm-20bpm for healthy adults) and recalculate RR iteratively, adjusting supposed RR for an optimal RR recognition.
  -8 probands: 1, 2, 4, 5, 6, 7, 8, 9 - male: 3, 5, 7, 8, 9
  -Go Direct respiration belt from Vernier around chest as ground truth
The sytem can easily be modified to other settings.

Exemplary experimental data and evaluation plots are attached.
Mean approach:
![Mean](https://user-images.githubusercontent.com/108615772/177041469-6acff0a4-31dc-4ff4-936a-655522ef91b5.png)
Median approach:
![Median](https://user-images.githubusercontent.com/108615772/177041473-11b20144-7074-4b60-b536-6a504ed4bc9f.png)


Files that build up the remote respiratory rate (RR) system included in this repository:
  rr_readC
  rr_compareCandRB
  rr_algorithms
  
rr_readC filters out necessary data from .bag-files created with IntelRealSense Viewer (settings cf. code). It can be chosen to either employ the median or the mean approach. Both average depths of pixels lying inside the ROI specified during recording for a certain time instant.
The mean approach has proven to deliver more precise results in  between 1m and 2m distance between proband and camera.
A .csv-file is created with the depth data and their respective timestamp.

rr_compareCandRB can read the camera (C) .csv-files and extract the RR. Same happens with .csv-files of the respiration belt (RB).
The Pearson Correlation Coefficient (PCC), Absolute Error (Abs Err) and Relative Error (Rel Err) comparing the C-RR and RB-RR as ground truth are calculated and plotted.

rr_algorithms contains all functions called during the exertion of the two previous python files.


Credits to Intel RealSense Team and the pyrealsense2 - library.
