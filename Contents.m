% +PLOT  Package with code for generating different types of custom plots.
%
% Labeling Utilities
%   add_multi_ylabel - Adds multiple y-axis labels to an axis with customized alignment, color, spacing, and font options.
%   add_scale_bar - Adds a scale bar to an existing axis, with customizable units, color, and position.
%   add_titles - Adds a title and subtitle to an axis with customizable font and color options.
%
% General Plotting Utilities
%   correlation                    - Plot correlation (optionally specifying an axes handle)
%   get_parent_figure              - Returns the parent figure even if it's not directly 'Parent' in the tree of the provided graphics object.
%   parameters                     - Return parameters struct, which sets default values for things like epoch durations etc.
%   raincloud                      - Plots a combination of half-violin, boxplot, and raw
%   raster                         - Plot raster pulse trains from cell array of sample times
%   time_amp_axes                  - Return figure and axes handle for time-vs-amplitude axes
%   time_trial_axes                - Return figure and axes handle for time-vs-trials axes
%   tt_spline                      - Plot smoothing spline fitted model from time table.
%
% EMG-Specific Plots
%   emg_array                      - Plot emg array given samples as they come from TMSi SDK.
%   emg_averages                   - Handle plotting gridded averages (TMSi/HD-EMG array)
%   emg_averages__bipolar          - EMG processing for bipolar TMSi channels
%   emg_averages__imu              - EMG processing for TMSi IMU channels
%   emg_averages__unipolar_array   - EMG processing for HD-EMG TMSi Array
%   emg_rms                        - Plot RMS heatmap with filled contour lines.
%   emg_stack                      - Create stack of EMG individual trial traces.
%   emg_waterfall                  - Create waterfall plots to show individual stimulus trial responses.
%   impedances_hdemg               - Generate impedance plot in grid orientation for HD-EMG
%   raw_tmsi_channel               - Handle plotting single-channel raw TMSi data.
%   synergy                        - Plot muscle synergies from nnmf results
%
% Abstract-Concept Plots
%   direct_form_I                  - Plot Direct Form-I implementation block diagram and pole/zero plot.
%   direct_form_II                 - Plot Direct Form-II implementation block diagram and pole/zero plot.
%
% Experiment-Specific Plots
%   diff_simulated_fields          - Plots difference between simulated fields from Block_A and Block_B
%   pattern_summary                - Plot summary of a given set of patterns, superimposed on approximation of sulcus.
%   potentiometer_event_alignments - Handle plotting potentiometers from wrist task against parsed events for a given exported alignment.
%   response_curves                - Creates set of response curves for parsed "response" metric data.
%   sba                            - Plot Scalable Brain Atlas for macaque
%   target_classification_tracking - Plot tracking through time for F1-score
%   trial_durations                - Plot bar graphs showing distribution of trial durations.
%   trigger_deltas                 - Handle plotting stem display of times between each "trigger" (sync) sample.
