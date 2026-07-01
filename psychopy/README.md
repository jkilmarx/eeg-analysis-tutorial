# PsychoPy and LSL experiment notes

The analysis assumes a PsychoPy task with 20 eyes-open trials and 20 eyes-closed trials, each lasting 15 seconds. At each trial onset, the task should push a short, consistent LSL marker such as `eyes_open` or `eyes_closed`.

Before recording participants:

1. Confirm the EEG and marker streams are visible in LabRecorder.
2. Start LabRecorder before the first task marker and stop it after the final marker.
3. Test marker labels, trial counts, and timing with a short pilot recording.
4. Save the experiment version and marker vocabulary with the study documentation.

The notebooks intentionally allow several label variants because real experiments often use `open`, `EO`, `closed`, or `EC`. Consistent labels make analysis safer.
