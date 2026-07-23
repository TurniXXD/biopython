from collections import defaultdict
from gzip import GzipFile
from io import TextIOWrapper
from pathlib import Path
from urllib.request import urlopen

import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

from sequence.utils import trim_ends

STATIC_DIR = Path(__file__).resolve().parent.parent / "static"

yersinia_pestis_biovar_fasta_url_fna = (
    "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna"
)

yersinia_pestis_biovar_fasta_url_gb = (
    "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb"
)

url_fastq = (
    "https://raw.githubusercontent.com/biopython/biopython/master/Tests/Quality/example.fastq.gz"
)


def seq_plots() -> None:
    # Fasta files are seq only
    with urlopen(yersinia_pestis_biovar_fasta_url_fna) as response:
        with TextIOWrapper(response, encoding="utf-8") as handle:
            record = SeqIO.read(handle, "fasta")

    seq = record.seq

    if seq is None:
        raise ValueError("Record has no sequence")

    window_size = 100

    positions = []
    g_fractions = []

    for start in range(0, len(seq), window_size):
        window = seq[start : start + window_size]

        positions.append(start + 1)
        g_fractions.append(window.count("G") / len(window))

    fig, ax = plt.subplots(figsize=(16, 9))
    ax.plot(positions, g_fractions)

    ax.set_xlabel("Sequence position")
    ax.set_ylabel("G fraction per 100 bases")

    # plt.show()

    url = "https://raw.githubusercontent.com/biopython/biopython/master/Tests/Quality/example.fastq.gz"
    qual_pos = defaultdict(list)

    with urlopen(url) as response:
        with GzipFile(fileobj=response) as gz:
            with TextIOWrapper(gz, encoding="utf-8") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    qualities = record.letter_annotations["phred_quality"]

                    trimmed = trim_ends(record, min_quality=20)

                    if len(trimmed) >= 20:
                        print(trimmed.id, len(record), "->", len(trimmed))

                    for index, quality in enumerate(qualities, start=1):
                        # if index <= 25 or quality == 40:
                        #     continue

                        qual_pos[index].append(quality)

    if not qual_pos:
        raise ValueError("No quality values matched the filter")

    positions = sorted(qual_pos)
    values_per_position = [qual_pos[position] for position in positions]

    fig, ax = plt.subplots(figsize=(16, 9))
    sns.boxplot(data=values_per_position, ax=ax)

    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels([str(position) for position in positions])

    ax.set_xlabel("Read position")
    ax.set_ylabel("Phred quality")

    # plt.show()
