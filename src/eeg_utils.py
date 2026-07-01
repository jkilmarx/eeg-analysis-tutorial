"""Readable helper functions for the eyes-open/eyes-closed EEG tutorial.

The functions favor clear checks and useful error messages over clever abstractions.
XDF metadata can vary between acquisition systems, so students should inspect the
stream summary before selecting streams or interpreting channel information.
"""

from __future__ import annotations

from collections.abc import Mapping

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix


def _first(value, default=None):
    """Return the first item from XDF's commonly used one-item lists."""
    if isinstance(value, (list, tuple, np.ndarray)):
        return value[0] if len(value) else default
    return default if value is None else value


def _stream_info_value(stream, key, default=None):
    """Read one metadata value while tolerating missing or unusual XDF fields."""
    info = stream.get("info", {})
    return _first(info.get(key), default)


def summarize_xdf_streams(streams):
    """Return one summary row per XDF stream as a pandas DataFrame.

    Parameters
    ----------
    streams : sequence of dict
        The stream list returned by ``pyxdf.load_xdf``.
    """
    rows = []
    for index, stream in enumerate(streams):
        samples = np.asarray(stream.get("time_series", []))
        observed_channels = samples.shape[1] if samples.ndim == 2 else 1
        rows.append(
            {
                "index": index,
                "name": _stream_info_value(stream, "name", "unknown"),
                "type": _stream_info_value(stream, "type", "unknown"),
                "channel_count": int(
                    float(_stream_info_value(stream, "channel_count", observed_channels))
                ),
                "nominal_srate_hz": float(
                    _stream_info_value(stream, "nominal_srate", 0.0)
                ),
                "sample_count": int(len(stream.get("time_stamps", []))),
            }
        )
    return pd.DataFrame(rows)


def find_stream_by_type_or_name(streams, keyword):
    """Return streams whose name or type contains ``keyword`` (case-insensitive).

    A list is returned because an XDF file can contain multiple matching streams.
    Use :func:`summarize_xdf_streams` to decide which match is appropriate.
    """
    keyword = str(keyword).casefold()
    matches = []
    for stream in streams:
        name = str(_stream_info_value(stream, "name", "")).casefold()
        stream_type = str(_stream_info_value(stream, "type", "")).casefold()
        if keyword in name or keyword in stream_type:
            matches.append(stream)
    return matches


def get_stream_channel_names(stream):
    """Extract channel labels from common XDF metadata layouts.

    Returns an empty list if labels are absent. Students may then provide labels
    manually or use generated names such as ``EEG001``.
    """
    try:
        channels = stream["info"]["desc"][0]["channels"][0]["channel"]
    except (KeyError, IndexError, TypeError):
        return []

    names = []
    for channel in channels:
        label = _first(channel.get("label")) if isinstance(channel, Mapping) else None
        names.append(str(label) if label else f"EEG{len(names) + 1:03d}")
    return names


def create_mne_raw_from_xdf_stream(eeg_stream, channel_names=None, ch_types="eeg"):
    """Convert a regularly sampled XDF EEG stream to ``mne.io.RawArray``.

    XDF stores samples as time-by-channel; MNE expects channel-by-time. The
    function converts common microvolt recordings to volts only when the XDF
    channel metadata explicitly report ``microvolts`` or ``uV``.
    """
    data = np.asarray(eeg_stream.get("time_series", []), dtype=float)
    if data.ndim == 1:
        data = data[:, np.newaxis]
    if data.ndim != 2 or data.shape[0] == 0:
        raise ValueError("The selected EEG stream has no two-dimensional sample data.")

    sampling_rate = float(_stream_info_value(eeg_stream, "nominal_srate", 0.0))
    if sampling_rate <= 0:
        timestamps = np.asarray(eeg_stream.get("time_stamps", []), dtype=float)
        if timestamps.size < 2:
            raise ValueError("Sampling rate is missing and cannot be inferred from timestamps.")
        sampling_rate = 1.0 / np.median(np.diff(timestamps))

    n_channels = data.shape[1]
    if channel_names is None:
        channel_names = get_stream_channel_names(eeg_stream)
    if not channel_names:
        channel_names = [f"EEG{number:03d}" for number in range(1, n_channels + 1)]
    if len(channel_names) != n_channels:
        raise ValueError(
            f"Found {n_channels} data columns but {len(channel_names)} channel names."
        )

    if isinstance(ch_types, str):
        ch_types = [ch_types] * n_channels
    if len(ch_types) != n_channels:
        raise ValueError("ch_types must be one value or one value per channel.")

    # XDF does not enforce units. Convert only when the metadata are explicit.
    try:
        channel_metadata = eeg_stream["info"]["desc"][0]["channels"][0]["channel"]
        units = [str(_first(item.get("unit"), "")).casefold() for item in channel_metadata]
    except (KeyError, IndexError, TypeError):
        units = []
    if units and len(units) == n_channels and all(
        unit in {"uv", "µv", "microvolt", "microvolts"} for unit in units
    ):
        data = data * 1e-6

    info = mne.create_info(channel_names, sfreq=sampling_rate, ch_types=ch_types)
    return mne.io.RawArray(data.T, info, verbose=False)


def extract_marker_annotations(marker_stream, raw_start_time=None):
    """Convert an XDF marker stream to MNE annotations.

    Parameters
    ----------
    marker_stream : dict
        An XDF stream containing marker samples and timestamps.
    raw_start_time : float | None
        XDF timestamp corresponding to time zero in the EEG ``Raw`` object. If
        omitted, the first marker becomes time zero; using the first EEG timestamp
        is usually the scientifically correct choice.
    """
    timestamps = np.asarray(marker_stream.get("time_stamps", []), dtype=float)
    samples = marker_stream.get("time_series", [])
    if len(timestamps) != len(samples):
        raise ValueError("Marker samples and marker timestamps have different lengths.")
    if len(timestamps) == 0:
        return mne.Annotations(onset=[], duration=[], description=[])

    if raw_start_time is None:
        raw_start_time = float(timestamps[0])

    descriptions = []
    for sample in samples:
        value = _first(sample, sample)
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        descriptions.append(str(value))

    return mne.Annotations(
        onset=timestamps - float(raw_start_time),
        duration=np.zeros(len(timestamps)),
        description=descriptions,
    )


def compute_bandpower_from_epochs(epochs, fmin, fmax, picks="eeg"):
    """Compute mean power in a frequency band for every epoch and channel.

    Returns
    -------
    bandpower : ndarray, shape (n_epochs, n_channels)
        Welch power spectral density averaged from ``fmin`` through ``fmax``.
    """
    if fmin >= fmax:
        raise ValueError("fmin must be smaller than fmax.")
    spectrum = epochs.compute_psd(
        method="welch", fmin=fmin, fmax=fmax, picks=picks, verbose=False
    )
    frequencies = spectrum.freqs
    psd = spectrum.get_data()
    # Integrating accounts for frequency-bin spacing and yields band power.
    return trapezoid(psd, frequencies, axis=-1)


def extract_bandpower_features(epochs, bands, picks="eeg"):
    """Create a labeled DataFrame of per-channel band-power features.

    ``bands`` should map a readable name to ``(fmin, fmax)``, for example
    ``{"theta": (4, 8), "alpha": (8, 12), "beta": (13, 30)}``.
    """
    if not isinstance(bands, Mapping) or not bands:
        raise ValueError("bands must be a non-empty mapping of name to (fmin, fmax).")

    picked_indices = mne.pick_types(epochs.info, eeg=picks == "eeg") if picks == "eeg" else None
    picked_names = (
        [epochs.ch_names[index] for index in picked_indices]
        if picked_indices is not None
        else epochs.copy().pick(picks).ch_names
    )

    columns = {}
    for band_name, limits in bands.items():
        if len(limits) != 2:
            raise ValueError(f"Band {band_name!r} must contain (fmin, fmax).")
        values = compute_bandpower_from_epochs(epochs, *limits, picks=picks)
        for channel_index, channel_name in enumerate(picked_names):
            columns[f"{band_name}_{channel_name}"] = values[:, channel_index]
    return pd.DataFrame(columns)


def plot_confusion_matrix_from_predictions(y_true, y_pred, labels=None):
    """Plot and return a confusion matrix from predicted class labels."""
    if labels is None:
        labels = list(dict.fromkeys(np.concatenate([np.asarray(y_true), np.asarray(y_pred)])))
    matrix = confusion_matrix(y_true, y_pred, labels=labels)
    display = ConfusionMatrixDisplay(confusion_matrix=matrix, display_labels=labels)
    display.plot(cmap="Blues", colorbar=False)
    display.ax_.set_title("Confusion matrix")
    plt.tight_layout()
    return display
