from pathlib import Path
from typing import Any, cast

from Bio import Entrez, SeqIO
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq, Seq, back_transcribe, reverse_complement, transcribe, translate
from Bio.SeqUtils import gc_fraction

STATIC_DIR = Path(__file__).resolve().parent.parent / "static"


def run_sequence_demo() -> None:
    seq = Seq("ACGT")
    print(f"Hello from biopython! {seq}")
    print(f"Hello from biopython! {seq.complement()}")
    print(f"Hello from biopython! {seq.complement_rna()}")
    print(f"Hello from biopython! {seq.reverse_complement()}")

    print("\n\n")

    for seq_record in SeqIO.parse(STATIC_DIR / "ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    print("\n\n")

    for seq_record in SeqIO.parse(STATIC_DIR / "ls_orchid.gbk", "genbank"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    handle = Entrez.esearch(
        db="nucleotide",
        term="TP53[Gene] AND Homo sapiens[Organism]",
        retmax=1,
    )
    record = cast(dict[str, Any], Entrez.read(handle))
    handle.close()

    print(record["IdList"])

    # DNA with high GC content is more stable.
    print(f"My seq gc fraction: {gc_fraction(seq)}")
    print("Seq: %s" % seq)

    coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    messenger_rna = coding_dna.transcribe()
    protein_seq = messenger_rna.translate()

    print(f"Protein seq: {protein_seq}")
    # Translating nucleotides up to the first in-frame stop codon and then stop.
    print(
        "Stop codon translation "
        f"{coding_dna.translate(table='Vertebrate Mitochondrial', to_stop=True)}"
    )

    gene = Seq(
        "GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA"
        "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT"
        "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT"
        "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT"
        "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA"
    )

    print(gene.translate(table="Bacterial"))

    # In the bacterial genetic code GTG is a valid start codon, and while it does
    # normally encode Valine, if used as a start codon it should be translated as
    # methionine. This happens if you tell Biopython your sequence is a complete CDS (coding seq).
    print(gene.translate(table="Bacterial", cds=True))

    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

    print(
        f"Translation tables: \n\n{standard_table}, {mito_table.start_codons}, {mito_table.stop_codons}"
    )

    # Partially defined sequence that is 159345973 bytes long with starting position at 117512683
    partial_seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length=159345973)

    # Finding sub seq
    search_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
    # subseq can be a string, bytes, bytearray, Seq or MutableSeq
    search_seq.index("ATGGGCCGC")
    search_seq.index(b"ATGGGCCGC")
    search_seq.index(bytearray(b"ATGGGCCGC"))
    search_seq.index(Seq("ATGGGCCGC"))
    search_seq.index(MutableSeq("ATGGGCCGC"))
    search_seq.find("ACTG")
    # Search starting from the right
    search_seq.rfind("ACTG")

    for index, sub in search_seq.search(["CC", "GGG", "CC"]):
        print(index, sub)

    # Working with strings
    str_seq = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
    reverse_complement(str_seq)
    transcribe(str_seq)
    back_transcribe(str_seq)
    translate(str_seq)
