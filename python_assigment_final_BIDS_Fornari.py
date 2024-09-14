# -*- coding: utf-8 -*-
"""
Created on Mon July 08 2024

@author: chiara.fornari
Source: https://mne.tools/stable/auto_tutorials/preprocessing/index.html

Title: a quick pre-processing pipline resting state EEG (eyes-closed) 

Dependencies (versions):
    python (3.11.5)
    mne (1.7.1)
    spyder (5.4.3)
    matplotlib (3.7.2)
    numpy (1.23.5)

"""

#%%############################################################################
""" IMPORT THE NECESSATY PYTHON MODULES """
###############################################################################

import numpy as np
#print(f"NumPy Version: {np.__version__}")
import os
import matplotlib
matplotlib.use('Qt5Agg')

# Checking the backends
#backends = ['Qt5Agg', 'TkAgg', 'Agg']
#for backend in backends:
#    try:
#        matplotlib.use(backend)
#        import matplotlib.pyplot as plt
#        # Add a simple plot to test
#        plt.plot([1, 2, 3], [1, 4, 9])
#        plt.title(f'Testing {backend}')
#        plt.show()
#        print(f'{backend} works!')
#    except Exception as e:
#        print(f'{backend} does not work: {e}')
        
import matplotlib.pyplot as plt
from pathlib import Path

import mne
#print(f"Version MNE: {mne.__version__}")
from mne_bids import write_raw_bids, BIDSPath
from scipy import signal


#%%############################################################################
""" Convert to BIDS format and DEFINE THE PATHS """ 
###############################################################################

# Variables
ID = '045' ## CHANGE
sID = 'sub-' + ID  
session = 'ses-01'  # if this is simpy "ses-{run}", consider define it from the run variable to keep them consistent and aligned
task = 'restEC'
run = '01'
suffix = 'eeg'

# Path
bids_root = Path() ## CHANGE based on your directory
eeg_file = bids_root / f'{ID}_RS_C_1.vhdr' # based on the name in which you saved the raw data

# Load EEG data
raw = mne.io.read_raw_brainvision(str(eeg_file))

# meta-data
bids_path = BIDSPath(subject=ID, session=run, task=task, run=run, root=bids_root, suffix=suffix)

# BIDS format
write_raw_bids(raw, bids_path, overwrite=True)

# sub-folders for a better organizazione of the pre-processing
source_data_dir = bids_root / sID / session / suffix
interp_data_dir = bids_root / sID / session / suffix / "Filtereddata" 
# Check if the folder exists, if not create it
if not interp_data_dir.exists():
    interp_data_dir.mkdir(parents=True)  # Create the folder and any missing parent directories
reports = bids_root / sID / session / suffix / "Report" 
if not reports.exists():  # LP: you can use the exist_ok parameter of the mkdir to avoid the if statement
    reports.mkdir(parents=True)  

# LP : don't leave around commented code
#code_dir = bids_root / "Code"
#if not code_dir.exists():
#    code_dir.mkdir(parents=True)  
  
source_data_dir.mkdir(parents=True, exist_ok=True)
interp_data_dir.mkdir(parents=True, exist_ok=True)
reports.mkdir(parents=True, exist_ok=True)

#%%############################################################################
""" Update README """ 
###############################################################################

# Define the path to your README file  # LP: define all custom variables at the beginning of the script
readme_file = bids_root / 'README.txt'  # Update with your actual path
# LP: I am not really following the logic of this readme update
# Information to add to the README file
additional_info = """

### Usage

This pipeline can be useful for general pre-processing resting-state eeg data. 
Please read the [GitHub](https://github.com/chiarafornari/python-cimec-chiara-fornari) README.md for additional information.

### Contact

For inquiries or issues, contact [Chiara](mailto:chiara.fornari@unitn.it).
"""

# Check if README file exists
# LP: My suggestion would be to just use pathlib classes and methods to check if the file exists, not mix with os functions
if os.path.exists(readme_file):
    # Open the README file in append mode and add additional information
    with open(readme_file, 'a') as f:
        f.write(additional_info)
    print(f"Additional information added to '{readme_file}'")
else:
    # Create the README file and add the information
    with open(readme_file, 'w') as f:
        f.write("# README\n\n")
        f.write("Initial content of README file.\n")
        # LP: I guess this is leftover code
        f.write(additional_info)
    print(f"README file '{readme_file}' created and additional information added.")
    

#%%############################################################################
""" ADD info and check the vhdr DATA """ 
###############################################################################

# add other info
raw.info['experimenter'] = 'Chiara Fornari'
# LP: I would suggest keeping every single variable a user might want to change at the beginning of the script

# Check the data structure
print(raw)
print(raw.info)

print("RAW DATA UPLOADED")     

""" NOTE
For each plot I create a new "named variable" so then I can include the figures in the final report. 
For each big pre-processing step I create a new dataset so it is possible to make a step back easly
"""

#%%############################################################################
""" CHANNEL LOCATION """ 
###############################################################################

# load the channel information for EasyCap - set montage
montage = mne.channels.make_standard_montage('easycap-M1')
raw.set_montage(montage, verbose=False)

# check the montage
montage.plot(kind='3d', show=True) #plot 3D representation of the electrode montage
raw.plot_sensors(show=True) #plot 2D montage

print("CHANNEL LOCATION UPLOADED")  

#%%############################################################################
""" PLOT vhdr DATA """ 
###############################################################################

# set parameters for plots and CHANGE if necessary
start = 0
duration=10
n_channels=64
scalings=50e-6

# plot raw data at 100 microV, here CHANGE 'start' depending on the data
rw = raw.plot(start=start, duration=duration, n_channels=n_channels, scalings=scalings, block=False) 
# plot the PSD (0-80 Hz)
rw_psd = raw.compute_psd(fmin=0, fmax=80).plot()
 
# create and plot the ccovariance matrix
cov = mne.compute_raw_covariance(raw, tmin=0, tmax=20) 

# Get the covariance values and channel names
cov_data = cov.data
ch_names = raw.ch_names

# Plot covariance matrix with channel names
cov_matrix = plt.figure(figsize=(10, 8))  # Adjust figure size as needed
plt.imshow(cov_data, origin='lower', cmap='RdBu_r', vmin=np.min(cov_data), vmax=np.max(cov_data))
plt.colorbar(label='Covariance')
plt.xticks(np.arange(len(ch_names)), ch_names, rotation=90)
plt.yticks(np.arange(len(ch_names)), ch_names)
plt.xlabel('Channels')
plt.ylabel('Channels')
plt.title('Covariance Matrix with Channel Names')

plt.tight_layout()
plt.show()

#%%############################################################################
""" PREPROCESSING """ 
###############################################################################

# downsampling: 250 Hz
dwnsamp_raw = raw.copy().resample(sfreq=250) 

# plot the psd of the dwnsampled db and the raw one. 
# Dwnsamp produces a low-pass filter at 1/2 of the dwnsamp value, in this case at 125 Hz, because aliasing.
for data, title in zip([raw, dwnsamp_raw], ["Origianl", "Downsampled"]):
    fig = data.compute_psd().plot(average=True, picks="data", exclude="bads")
    fig.subplots_adjust(top=0.9)
    fig.suptitle(title)
    plt.setp(fig.axes, xlim=(0,250))
dwn = dwnsamp_raw.plot(start=start, duration=start, n_channels=n_channels, bgcolor='w', color='y', scalings=scalings)
print("DOWNSAMPLING DONE")    

# filter the data: high-pass at 0.1 Hz and low-pass at 40 Hz
l_freq = 0.1 
h_freq = 40
filt_raw = dwnsamp_raw.copy()
filt_raw.load_data().filter(l_freq, h_freq) 
# plot 
filt = filt_raw.plot(start=start, duration=duration, n_channels=n_channels, bgcolor='w', color='b', scalings=scalings) 
filt_psd = filt_raw.compute_psd(fmin=0, fmax=80).plot()  
print("FILTERS DONE")

# Notch Filter: 48-52 Hz
notch_raw = filt_raw.copy()
notch_raw.notch_filter(freqs = [48,52])
# LP: For all plots, I would try to find a way to save them without showing them, so one check them out only if needed
notch = notch_raw.plot(start=start, duration=duration, n_channels=n_channels, bgcolor='w', color='g', scalings=scalings) 
print("NOTCH DONE")

#Create the final db after dwnsamp and filters
EEG_preproc = notch_raw.copy()
print(EEG_preproc)
print(EEG_preproc.info)
print("PRE-PROCESSED DB CREATED")

# save the db
EEG_preproc.save(str(interp_data_dir / f'{sID}_EEG_preprocessed.fif'), overwrite=True)
print("SAVED")

#%%############################################################################
"""REMOVE AND INTERPOLATE BAD CHANNELS"""
###############################################################################
# before marking any channels, the list is empty:
bc = EEG_preproc.copy()
print(bc.info['bads'])

# Visualize
fig = bc.plot(butterfly=True, color='#00000022', bad_color='r', title = 'bads_marked', scalings=scalings)

# Interpolate the bad channels
interp_bc = bc.copy()
interp_bc.load_data().interpolate_bads(reset_bads=False)

#%%############################################################################
"""CHANNELS REMOVED""" 
###############################################################################

# Removed channels if very very very bad.
ch_removed = 0 # CHANGE!!!!#

#%%############################################################################
"""CREATING HEOG, ECG, VEOG """
###############################################################################

heog_ch_name = 'HEOG'
# LP: some hardcoding here as well
mne.set_bipolar_reference(inst = interp_bc, anode = ['FT9'], cathode = ['FT10'], ch_name= heog_ch_name, drop_refs= False, copy = False)
veog_evoked = mne.channels.combine_channels(interp_bc, groups= dict(VEOG = [0,1]), method= 'mean', keep_stim= True)
ecg_evoked = mne.channels.combine_channels(interp_bc, groups= dict(ECG = [35,35]), method= 'mean', keep_stim= True)
print(interp_bc.info['ch_names'][35])
veog_evoked.plot()
interp_bc.add_channels([veog_evoked], force_update_info = True)
interp_bc.add_channels([ecg_evoked], force_update_info = True)
interp_bc.set_channel_types({'ECG':'ecg', 'VEOG' :'eog', heog_ch_name: 'eog'})
print(interp_bc.info['ch_names'])

#%%#############################################################################
""" SAVE THE DATA """ 
###############################################################################
# original raw data file which is filtered and interpolated
interp_bc.save(str(interp_data_dir / f'{sID}_EEG_filtered_interpolated_.fif'), overwrite=True)
interp_bc.plot(start=start, duration=duration, n_channels=n_channels, bgcolor='w', color='b', scalings=scalings)

#%%############################################################################
""" SAVE REPORT """ 
###############################################################################

# recall for saving in the report but not showing it
rw = raw.plot(start=start, duration=duration, n_channels=n_channels, scalings=scalings) 
plt.close()  # LP: this is what I meant - maybe have a variable at the beginning of the script that you can set to True or False to show the plots
dwn = dwnsamp_raw.plot(start=start, duration=duration, n_channels=n_channels, bgcolor='w', color='y', scalings=scalings)
plt.close() 
filt = filt_raw.plot(start=start, duration=duration, n_channels=n_channels, bgcolor='w', color='b', scalings=scalings) 
plt.close() 
notch = notch_raw.plot(start=start, duration=duration, n_channels=n_channels, bgcolor='w', color='g', scalings=scalings) 
plt.close() 

report0 = mne.Report(title= 'Report - Preprocessing')
report0.add_raw(raw = interp_bc, title = sID)
report0.add_figure(title ='Raw', fig= rw)
report0.add_figure(title ='Raw PSD', fig= rw_psd)
report0.add_figure(title ='Raw channel covariance matrix', fig= cov_matrix)
report0.add_figure(title ='Downsampling', fig= dwn)
report0.add_figure(title ='H-pass and L-pass filters', fig= filt)
report0.add_figure(title ='Psd after filtering', fig= filt_psd)
report0.add_figure(title ='Notch', fig= notch)

#save 
report0.save(reports / f'{sID}_EEG_preprocessing.html', overwrite=True, open_browser=False)
print("REPORT SAVED")

print("THE PREPROCESSING IS FINISHED :)")

#%%############################################################################
""" RE-REFERENCING """ 
###############################################################################

# remove useless channels 
re_ref = interp_bc.copy()
re_ref = re_ref.drop_channels(ch_names=['HEOG', 'VEOG', 'ECG'], on_missing='raise')

print("HEOG, VEOG, ECG REMOVED")  

# set the online ref channel
re_ref = mne.add_reference_channels(re_ref, ref_channels=["TP10"])
reref = re_ref.plot(start=start, duration=duration, n_channels=n_channels, scalings=scalings, title = 'Add ref channel') 

# re-referencing using the average
ave_ref = re_ref.copy()
ave_ref.set_eeg_reference(ref_channels='average').apply_proj() 
averaged = ave_ref.plot(start=start, duration=duration, n_channels=68, scalings=scalings, title ='rereferenced average')

print("RE-REFERECING DONE")

#Add referencing info in the report
report0.add_figure(title ='Referencing', fig= averaged)
report0.save(reports / f'{sID}_EEG_preprocessing.html', overwrite=True, open_browser=False)
print("REPORT UPDATED")

# save
ave_ref.save(str(interp_data_dir / f'{sID}_EEG_referencing.fif'), overwrite=True)
print("SAVED")

#%%############################################################################
""" RUN ICA """ 
###############################################################################

# processing ICA
ica_database = ave_ref.copy()
# LP: these (n_components, max_iter...) is exactly the kind of parameters that you want to keep up center in the beghinning of the script,
# maybe even uppercase! 
ica = mne.preprocessing.ICA(n_components= 62, max_iter="auto", random_state=97) # Extended ICA: Infomax; Defaults ICA: fastica
ica.fit(ica_database)

# Plot ICA components
components = ica.plot_components(inst = ica_database)

#%%############################################################################
""" Removing components after ICA """ 
###############################################################################

# Specify components to remove (adjust as needed based on inspection)
ica.exclude = [0, 1, 2]   ### CHANGE !!!
excluded = ica.plot_properties(ica_database, picks=ica.exclude)

# Apply artifact rejection
cleaned = ica_database.copy()
ica.apply(cleaned)

# inspect the signal
cleaned_plot = cleaned.plot(start=start, duration=duration, n_channels=n_channels, scalings=scalings, color='purple', title = 'Cleaned database after ICA')
print("ICA DONE") 

# save
cleaned.save(str(interp_data_dir / f'{sID}_ICA_cleaned.fif'), overwrite=True)
print("SAVED")

#Add ICA info in the report
report0.add_ica(ica = ica, title = 'ICA', inst = cleaned)
report0.add_figure(title ='ICA components removed', fig=excluded)
report0.add_figure(title ='Cleaned database', fig= cleaned_plot)

report0.save(reports / f'{sID}_EEG_preprocessing.html', overwrite=True, open_browser=False)

print("REPORT UPDATED")

#%%############################################################################
""" Power Spectral Density (PSD) """
###############################################################################

# Example of posterior electrode selection (adjust based on your specific montage)
posterior_electrodes = ['O1', 'O2', 'P3', 'P4', 'Pz', 'POz', 'PO7', 'PO8', 'PO3', 'PO4']  
frontal_electrodes = ['Fz', 'F3', 'F4', 'F1', 'F2', 'F7', 'F8', 'F6', 'F5']

# Select EEG channels and posterior/frontal electrodes
picks = mne.pick_types(cleaned.info, meg=False, eeg=True, eog=False, stim=False, exclude='bads')
posterior_picks = mne.pick_channels(cleaned.info['ch_names'], include=posterior_electrodes)
frontal_picks = mne.pick_channels(cleaned.info['ch_names'], include=frontal_electrodes)

# Calculate power spectral density using Welch's method for posterior electrodes
fs = cleaned.info['sfreq']  # Sampling frequency
freqs, psd_posterior = signal.welch(cleaned.get_data(picks=posterior_picks), fs=fs, nperseg=fs*2, scaling='density')

# Calculate power spectral density using Welch's method for frontal electrodes
freqs, psd_frontal = signal.welch(cleaned.get_data(picks=frontal_picks), fs=fs, nperseg=fs*2, scaling='density')

# Calculate average PSD across electrodes
avg_psd_posterior = np.mean(psd_posterior, axis=0)
avg_psd_frontal = np.mean(psd_frontal, axis=0)

# Plot power spectrum and topography
psd_final = plt.figure(figsize=(10, 6))

# Plot Power Spectrum
plt.subplot(1, 2, 1)
plt.plot(freqs, 10 * np.log10(avg_psd_posterior), color='purple', label='Posterior (Alpha)')
plt.plot(freqs, 10 * np.log10(avg_psd_frontal), color='green', label='Frontal (Theta)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density (dB/Hz)')
plt.title('Power Spectrum - Alpha and Theta')
plt.xlim([1, 20])  # Adjust xlim based on your frequency range of interest
plt.legend()
plt.grid(True)

# Plot Topography of Alpha and Theta Band Power
avg_alpha_power_topo = np.zeros(len(cleaned.info['ch_names']))
avg_alpha_power_topo[posterior_picks] = np.mean(psd_posterior[:, (freqs >= 8) & (freqs <= 12)], axis=1)

avg_theta_power_topo = np.zeros(len(cleaned.info['ch_names']))
avg_theta_power_topo[frontal_picks] = np.mean(psd_frontal[:, (freqs >= 4) & (freqs <= 8)], axis=1)

avg_power_topo = avg_alpha_power_topo + avg_theta_power_topo

mne.viz.plot_topomap(avg_power_topo, cleaned.info, names=cleaned.info['ch_names'])

plt.tight_layout()
plt.show()

#Add PSD info in the report
report0.add_figure(title ='Power Spectral Density (PSD)', fig=psd_final)
report0.save(reports / f'{sID}_EEG_preprocessing.html', overwrite=True, open_browser=False)

print("REPORT UPDATED")