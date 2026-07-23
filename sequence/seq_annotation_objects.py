from io import TextIOWrapper
from pathlib import Path
from urllib.request import urlopen

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

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


def annotation_seq() -> None:
    simple_seq = Seq("GATC")
    simple_seq_record = SeqRecord(simple_seq, id="AC1")
    print(simple_seq_record.seq)
    simple_seq_record.annotations["evidence"] = "You better believe I have some."
    # Phred quality score estimates how confident a DNA sequencer is that it identified a nucleotide correctly.
    simple_seq_record.letter_annotations["phred_quality"] = [40, 40, 38, 30]
    print(simple_seq_record.annotations)

    # Fasta files are seq only
    with urlopen(yersinia_pestis_biovar_fasta_url_fna) as response:
        with TextIOWrapper(response, encoding="utf-8") as handle:
            record = SeqIO.read(handle, "fasta")

    seq = record.seq
    if seq is None:
        raise ValueError("The FASTA record has no sequence")

    print(record)
    print(len(seq))
    print(seq[:50])

    # GenBank files are seq + SeqFeature annotations
    with urlopen(yersinia_pestis_biovar_fasta_url_gb) as response:
        with TextIOWrapper(response, encoding="utf-8") as handle:
            record_gb = SeqIO.read(handle, "genbank")

    for feature in record_gb.features:
        print(feature)
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", ["Unknown"])[0]
            translation = feature.qualifiers.get("translation", ["Missing"])[0]
            print(f"Product: {product}")
            print(f"Protein: {translation}")

    seq = record.seq

    if seq is None:
        raise ValueError("Record has no sequence")

    # Gene is read from the opposite strand
    feature = SeqFeature(SimpleLocation(5, 18, strand=-1), type="gene")

    location = feature.location

    if location is None:
        raise ValueError("Feature has no location")

    feature_seq = seq[location.start : location.end].reverse_complement()
    print(feature_seq)
    extracted_feature_seq = feature.extract(seq)
    print(extracted_feature_seq)

    unformatted_record = SeqRecord(
        Seq(
            "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
            "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
            "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
            "SSAC"
        ),
        id="gi|14150838|gb|AAK54648.1|AF376133_1",
        description="chalcone synthase [Cucumis sativus]",
    )
    print(unformatted_record.format("fasta"))
