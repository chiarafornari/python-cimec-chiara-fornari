# A simple resting state EEG pre-processing pipeline

This repository contains the final project for the course **Python for (open) Neuroscience** taught by Luigi Petrucco at the Doctoral school in Cognitive and Brain Sciences at CIMeC, University of Trento.

***
The project consists of a simple pre-processing pipeline for EEG data collected in humans. In the spirit of the course title, the script uses the BIDS format to organize the data.  
The assigment is built in Spyder (5.4.3) using Python (3.11.5).

It includes:
- a file with a custom Python module that contains the whole script (`python_assigment_final_BIDS_Fornari.py`).
- an example of resting state (rs) EEG data for one older subject (ID: 045).

## Dependencies

| Package   | Version used | 
|-----------|--------------|
| [MNE-Python](https://mne.tools/stable/index.html)     | 1.7.1        | 
| [MNE-BIDS](https://mne.tools/mne-bids/stable/index.html)  | 0.15.0       |
| [Matplotlib](https://matplotlib.org/)| 3.7.2        |
| [NumPy](https://numpy.org/)    | 1.23.5       | 
| [SciPy](https://scipy.org/)    | 1.11.1       |


## Installation

For the analysis of EEG data is necessary to install MNE package [[1].(https://doi.org/10.3389/fnins.2013.00267)] and MNE-BIDS [[2].(https://doi.org/10.21105/joss.01896),[3].(https://doi.org/10.1038/s41597-019-0104-8)] to save and convert the data in the BIDS format.

Run the followings in a command line terminal: `pip install mne` and `pip install mne-bids`

***

## EEG Data info and processing 

#### Data collection:
The raw EEG data are collected using BrainVision Recorder using a 64 channel active EasyCap with a 1000 Hz sampling frequency. The impedence of the eletrocdes was keept below 20 kΩ using a saline gel.
The online reference electrode was TP9 and the ground electrode was Fpz.
The signal is recorded for 3 minutes in a resting state condition asking to the participant to keep the eyes closed.

#### Data filtering:
Raw data were downsampled to 250 Hz and both highpass (0.1 Hz) and lowpass (40 Hz) were applied, as well as notch filter (48-52 Hz) to reduce line noise.

#### Data quality:
A visual inspection of the data is provided in order to check the quality of the signal and the presence of bad channels that in case are interpolated. HEOG, VEOG ans ECG channel are created to control for eye movement and heart rate.

#### Data ri-referencing:
In order to reduce the noise, the average electrodes references was applied.

#### Independent Component Analysis (ICA):
The ICA is applied and the bad components removed.

#### Power spectral density (PSD):
A general PSD plot is provided to check the occipital alpha (8-12 Hz) and frontal theta (4-8 Hz)

***

## Initializate the environment

The raw EEG data (with the three extentions: `.eeg` `.vhdr` `.vmrk`) provided as example should be seved in your local directory (variable name: bids_root). 
Rememer to change it accordlying before running the script! When you have to change information a comment `# CHANGE` appears.
Once defined the directory, the script automatically reshape the data and results in the BIDS format. 

To keep trak of all the processess, `.fif` files are saved in sub-folder (Filtereddata) that are automatically generated. I prefer to save also this data, so it is easier to check the quality and so on.
To summarize, a report file is also saved for each participant in a sub-folder (Report) that are automaticallu generated.

***

## References
[1] Alexandre Gramfort, Martin Luessi, Eric Larson, Denis A. Engemann, Daniel Strohmeier, Christian Brodbeck, Roman Goj, Mainak Jas, Teon Brooks, Lauri Parkkonen, and Matti S. Hämäläinen. MEG and EEG data analysis with MNE-Python. Frontiers in Neuroscience, 7(267):1–13, 2013. https://doi.org/10.3389/fnins.2013.00267.

[2] Appelhoff, S., Sanderson, M., Brooks, T., Vliet, M., Quentin, R., Holdgraf, C., Chaumon, M., Mikulan, E., Tavabi, K., Höchenberger, R., Welke, D., Brunner, C., Rockhill, A., Larson, E., Gramfort, A. and Jas, M. (2019). MNE-BIDS: Organizing electrophysiological data into the BIDS format and facilitating their analysis. Journal of Open Source Software 4: (1896). https://doi.org/10.21105/joss.01896

[3] Pernet, C. R., Appelhoff, S., Gorgolewski, K. J., Flandin, G., Phillips, C., Delorme, A., Oostenveld, R. (2019). EEG-BIDS, an extension to the brain imaging data structure for electroencephalography. Scientific Data, 6, 103. https://doi.org/10.1038/s41597-019-0104-8



