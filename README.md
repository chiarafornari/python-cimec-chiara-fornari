# A simple resting state EEG pre-processing pipeline

This repository contains the final project for the course **Python for (open) Neuroscience** taught by Luigi Petrucco at the Doctoral school in Cognitive and Brain Sciences at CIMeC, University of Trento.

***
The project consists of a simple pre-processing pipeline for EEG data collected in humans. In the spirit of the course title, the script uses the BIDS format to organize the data.  
The assigment is built in Spyder (5.4.3) using Python (3.11.5).

It includes:
- a file with a custom Python module that contains the whole script (`python_assigment_final_BIDS_Fornari.py`).
- an example of resting state (rs) EEG data of one older adults subject (ID: 045).

## Dependencies

| Package   | Version used | 
|-----------|--------------|
| [MNE-Python](https://mne.tools/stable/index.html)     | 1.7.1        | 
| [MNE-BIDS](https://mne.tools/mne-bids/stable/index.html)  | 0.15.0       |
| [Matplotlib](https://matplotlib.org/)| 3.7.2        |
| [NumPy](https://numpy.org/)    | 1.23.5       | 
| [SciPy](https://scipy.org/)    | 1.11.1       |


## Installation

For the analysis of EEG data is necessary to install MNE package [[1](https://doi.org/10.3389/fnins.2013.00267)] and MNE-BIDS [[2](https://doi.org/10.21105/joss.01896),[3](https://doi.org/10.1038/s41597-019-0104-8)] to save and convert the data in the BIDS format.

Run the followings in a command line terminal: `pip install mne` and `pip install mne-bids`

***

## EEG Data info and processing 

#### Data collection:
The raw EEG data are collected using BrainVision Recorder using a 64 channel active EasyCap with a 1000 Hz sampling rate frequency. The impedence of the eletrocdes was keept below 20 kΩ using a saline gel.
The online reference electrode was placed on TP9 and the ground electrode on Fpz.
The signal is recorded for 3 minutes in a resting state condition asking to the participant to keep the eyes closed.

#### Data filtering:
Raw data were downsampled to 250 Hz and both highpass (0.1 Hz) and lowpass (40 Hz) filters were applied, as well as notch filter (48-52 Hz) to reduce the line noise.

#### Data quality:
A visual inspection of the data is provided in order to check the quality of the signal and the presence of bad channels that in case are interpolated. HEOG, VEOG and ECG channel are created to control/clean for eye movement and heart rate.

#### Data ri-referencing:
In order to reduce the noise, the average electrodes offline reference was applied.

#### Independent Component Analysis (ICA):
The ICA is applied and the bad components removed.

#### Power spectral density (PSD):
A general PSD plot is provided to check the occipital alpha (8-12 Hz) and frontal theta (4-8 Hz).

***

## Initializate the environment

The raw EEG data (with the three extentions: `.eeg` `.vhdr` `.vmrk`) provided as example should be saved in your local directory (variable name in the `.py` file: bids_root). 
Rememer to change it accordlying before running the script! When you have to change information a comment `# CHANGE` appears.
Once defined the directory, the script automatically reshape the data and results in the BIDS format. 

To keep track of all the processess, `.fif` files are saved in sub-folder (Filtereddata) that is automatically generated. I prefer to save also this data, so it is easier to check the quality and so on.
To summarize, a report file is also saved for each participant in a sub-folder (Report) that is automatically generated.

***

## References
[1] Gramfort, et al. (2013). MEG and EEG data analysis with MNE-Python. _Frontiers in Neuroscience_, 7(267), 1–13. https://doi.org/10.3389/fnins.2013.00267.

[2] Appelhoff, et al. (2019). MNE-BIDS: Organizing electrophysiological data into the BIDS format and facilitating their analysis. _Journal of Open Source Software_ 4: (1896). https://doi.org/10.21105/joss.01896

[3] Pernet, et al. (2019). EEG-BIDS, an extension to the brain imaging data structure for electroencephalography. _Scientific Data_, 6, 103. https://doi.org/10.1038/s41597-019-0104-8
