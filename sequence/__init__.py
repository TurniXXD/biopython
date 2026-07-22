from .seq_annotation_objects import annotation_seq
from .seq_objects import run_sequence_demo


def run_sequence_examples() -> None:
    run_sequence_demo()
    annotation_seq()


__all__ = ["run_sequence_demo", "annotation_seq", "run_sequence_examples"]
