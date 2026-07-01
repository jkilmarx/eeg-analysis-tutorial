"""Small helper functions for the EEG analysis tutorial."""

from .eeg_utils import (
    compute_bandpower_from_epochs,
    create_mne_raw_from_xdf_stream,
    extract_bandpower_features,
    extract_marker_annotations,
    find_stream_by_type_or_name,
    get_stream_channel_names,
    plot_confusion_matrix_from_predictions,
    summarize_xdf_streams,
)

__all__ = [
    "summarize_xdf_streams",
    "find_stream_by_type_or_name",
    "get_stream_channel_names",
    "create_mne_raw_from_xdf_stream",
    "extract_marker_annotations",
    "compute_bandpower_from_epochs",
    "extract_bandpower_features",
    "plot_confusion_matrix_from_predictions",
]
