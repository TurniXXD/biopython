from io import TextIOWrapper
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez, SeqIO
from Bio.SeqUtils.CheckSum import seguid


def seq_io() -> None:
    for seq_record in SeqIO.parse("./static/ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    identifiers = [r.id for r in SeqIO.parse("./static/ls_orchid.gbk", "genbank")]
    print(identifiers)

    record_iterator = SeqIO.parse("./static/ls_orchid.fasta", "fasta")
    first_record = next(record_iterator)
    print(first_record.id)
    print(first_record.description)

    records = list(SeqIO.parse("./static/ls_orchid.gbk", "genbank"))
    print("Found %i records" % len(records))
    print(repr(records[-1].seq))
    print(records[-1].annotations["organism"])
    print(records[0].annotations["organism"])

    with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="6273291") as handle:
        seq_record = SeqIO.read(handle, "fasta")

    print("%s with %i features" % (seq_record.id, len(seq_record.features)))

    with Entrez.efetch(db="nucleotide", rettype="db", retmode="text", id="6273291") as handle:
        seq_record = SeqIO.read(handle, "gb")

    print("%s with %i features" % (seq_record.id, len(seq_record.features)))

    with Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text", id="6273291,6273290,6273289"
    ) as handle:
        for seq_record in SeqIO.parse(handle, "gb"):
            print("%s %s..." % (seq_record.id, seq_record.description[:50]))
            print(
                "Sequence length %i, %i features, from: %s"
                % (
                    len(seq_record),
                    len(seq_record.features),
                    seq_record.annotations["source"],
                )
            )

    url_uniprot = "https://rest.uniprot.org/uniprotkb/O23729.txt"

    try:
        with urlopen(url_uniprot, timeout=15) as response:
            with TextIOWrapper(response, encoding="UTF-8") as handle:
                seq_record = SeqIO.read(handle, "swiss")

        print(seq_record.id)
        print(seq_record.annotations.get("keywords", []))

    except HTTPError as error:
        print(f"UniProt returned HTTP {error.code}: {error.reason}")

    except URLError as error:
        print(f"Could not connect to UniProt: {error.reason}")

    for record in SeqIO.parse("./static/ls_orchid.gbk", "genbank"):
        # Sequence Globally Unique Identifier
        # SEGUID is a checksum used to create a compact identifier from a biological sequence.
        # Two records have different ID but might have same seq, with SEGUID we can detect duplicates.
        print(record.id, seguid(record.seq))

    records_by_checksum = SeqIO.to_dict(
        SeqIO.parse("./static/ls_orchid.gbk", "genbank"),
        key_function=lambda record: seguid(record.seq),
    )

    for checksum, record in records_by_checksum.items():
        print(checksum, record.id)
