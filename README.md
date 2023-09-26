# +plot #
Package for basic MATLAB data plots.

## Usage ##
See function descriptions, or `help plot` in the MATLAB editor for usage details.

### Git ###
To add this file as a submodule package to a MATLAB project:
1. Open a git bash in the desired project folder.
2. Use the following syntax to add the submodule as a package folder:  
```(git)
git submodule add git@github.com:Neuro-Mechatronics-Interfaces/MATLAB_Data_Plots.git +plot
```

## Contents ##
### General Plotting Utilities ###  
 + [correlation](correlation.m) - Plot correlation (optionally specifying an axes handle).  
 + [parameters](parameters.m) - Return parameters struct, which sets default values for things like epoch durations etc.  
 + [raincloud](raincloud.m) - Plots a combination of half-violin, boxplot, and raw.  
 + [raster](raster.m) - Plot raster pulse trains from cell array of sample times.  
 + [time_amp_axes](time_amp_axes.m) - Return figure and axes handle for time-vs-amplitude axes.  
 + [time_trial_axes](time_trial_axes.m) - Return figure and axes handle for time-vs-trials axes.  
 + [tt_spline](tt_spline.m) - Plot smoothing spline fitted model from time table.  

### EMG-Specific Plots ###  
 + [emg_array](emg_array.m) - Plot emg array given samples as they come from TMSi SDK.  
 + [emg_averages](emg_averages.m) - Handle plotting gridded averages (TMSi/HD-EMG array).  
 + [emg_averages__bipolar](emg_averages__bipolar.m) - EMG processing for bipolar TMSi channels.  
 + [emg_averages__imu](emg_averages__imu.m) - EMG processing for TMSi IMU channels.  
 + [emg_averages__unipolar_array](emg_averages__unipolar_array.m) - EMG processing for HD-EMG TMSi Array.  
 + [emg_rms](emg_rms.m) - Plot RMS heatmap with filled contour lines.  
 + [emg_stack](emg_stack.m) - Create stack of EMG individual trial traces.  
 + [emg_waterfall](emg_waterfall.m) - Create waterfall plots to show individual stimulus trial responses.  
 + [impedances_hdemg](impedances_hdemg.m) - Generate impedance plot in grid orientation for HD-EMG.  
 + [raw_tmsi_channel](raw_tmsi_channel.m) - Handle plotting single-channel raw TMSi data.  
 + [synergy](synergy.m) - Plot muscle synergies from nnmf results.  

### Abstract-Concept Plots ###  
 + [direct_form_I](direct_form_I.m) - Plot Direct Form-I implementation block diagram and pole/zero plot. 
 + [direct_form_II](direct_form_II.m) - Plot Direct Form-II implementation block diagram and pole/zero plot.  

### Experiment-Specific Plots ###  
 + [diff_simulated_fields](diff_simulated_fields.m) - Plots difference between simulated fields from Block_A and Block_B.  
 + [pattern_summary](pattern_summary.m) - Plot summary of a given set of patterns, superimposed on approximation of sulcus.  
 + [potentiometer_event_alignments](potentiometer_event_alignments.m) - Handle plotting potentiometers from wrist task against parsed events for a given exported alignment.  
 + [response_curves](response_curves.m) - Creates set of response curves for parsed "response" metric data.  
 + [sba](sba.m) - Plot Scalable Brain Atlas for macaque.  
 + [target_classification_tracking](target_classification_tracking.m) - Plot tracking through time for F1-score.  
 + [trial_durations](trial_durations.m) - Plot bar graphs showing distribution of trial durations.  
 + [trigger_deltas](trigger_deltas.m) - Handle plotting stem display of times between each "trigger" (sync) sample.  
