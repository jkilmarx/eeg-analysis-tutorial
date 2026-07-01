# EEG Analysis Tutorial: Eyes Open vs. Eyes Closed

[Python](https://www.python.org/) · [MNE-Python](https://mne.tools/) · [pyxdf](https://github.com/xdf-modules/pyxdf) · [Jupyter](https://jupyter.org/)

This beginner-friendly tutorial follows an EEG recording from an XDF file through preprocessing and condition-level spectral analysis. The example experiment uses PsychoPy to present eyes-open and eyes-closed trials, Lab Streaming Layer (LSL) to send event markers, and LabRecorder to store synchronized EEG and marker streams in XDF.

## Scientific question

Can we detect stronger alpha-band (8–12 Hz) activity during eyes-closed than eyes-open EEG, especially over posterior scalp channels?

The expected group-level pattern is increased posterior alpha power during eyes closed. An individual recording may differ because of data quality, sensor placement, participant state, referencing, or preprocessing choices.

## Learning objectives

By completing the notebooks, you will learn to:

- inspect an XDF recording without assuming its stream layout;
- convert an EEG stream into an MNE `Raw` object;
- identify bad channels, filter, notch, interpolate, and re-reference EEG;
- align LSL task markers with EEG and create 15-second epochs; and
- compare spectra and alpha topographies between conditions.

## Expected experiment structure

| Item | Expected value |
|---|---:|
| Conditions | Eyes open and eyes closed |
| Trials per condition | 20 |
| Trial duration | 15 seconds |
| Event timing | Markers sent through LSL |
| Recording container | XDF produced by LabRecorder |

Treat these values as an experimental expectation, not a guarantee. Always inspect the streams and marker counts before analysis.

## Repository structure

```text
eeg-analysis-tutorial/
├── README.md
├── environment.yml
├── .gitignore
├── data/                 # Place local XDF recordings here
├── outputs/              # Generated FIF files and figures
├── psychopy/             # Notes for connecting the experiment
├── notebooks/            # Three sequential tutorial notebooks
├── src/                  # Reusable, readable helper functions
└── images/               # Optional documentation images
```

## Setup

1. Install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or another Conda-compatible distribution.
2. From this repository's root, create the environment:

   ```bash
   conda env create -f environment.yml
   conda activate eeg-analysis-tutorial
   python -m ipykernel install --user --name eeg-analysis-tutorial --display-name "Python (EEG tutorial)"
   ```

3. Put a de-identified recording at `data/sample_recording.xdf`, or edit `XDF_PATH` in the notebooks.
4. Start Jupyter from the repository root:

   ```bash
   jupyter lab
   ```

5. Open `notebooks/01_data_loading_and_raw_eeg.ipynb` and select the `Python (EEG tutorial)` kernel.

Run the notebooks in order. They use relative paths and are designed to be run with the notebook working directory set to `notebooks/`. The analysis cells are guided skeletons: read each instruction, replace the `...` placeholders, and complete every `TODO` before moving to the next step.

## Notebook guide

1. **Data Loading and Raw EEG** — implement XDF inspection, stream selection, MNE conversion, marker alignment, and saving.
2. **EEG Preprocessing** — decide and implement bad-channel handling, interpolation, filtering, referencing, and quality-control comparisons.
3. **Task Comparison: PSDs and Topomaps** — implement marker normalization and 15-second epoching, then use MNE's PSD plots and a simple mean-alpha summary to create topomaps.

## Data are not included

Sample EEG data are intentionally not included by default. XDF and common EEG data formats are ignored by Git so that a local recording is not accidentally committed.

Before sharing any human-participant data, follow the consent language, institutional review requirements, and data-use agreements for the study. Remove direct and indirect identifiers, inspect stream metadata for names or device identifiers, and share only through an approved access-controlled or public repository. De-identification reduces risk but does not automatically make a dataset safe to publish.

## Reproducibility notes

- Record preprocessing decisions, rejected trials, bad channels, software versions, and any deviations from the expected protocol.
- Keep raw data unchanged; write derived files to `outputs/`.
- Inspect PSDs, topographies, marker counts, and retained trial counts before interpreting a condition difference.
