from .seq_annotation_objects import annotation_seq
from .seq_io import seq_io
from .seq_objects import run_sequence_demo
from .seq_plots import seq_plots


def run_sequence_examples() -> None:
    run_sequence_demo()
    annotation_seq()
    # seq_plots()
    seq_io()


__all__ = ["run_sequence_demo", "annotation_seq", "seq_plots", "seq_io", "run_sequence_examples"]
