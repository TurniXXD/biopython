import os
from typing import Any, cast

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.getenv("ENTREZ_EMAIL", "")


def main():
    if not Entrez.email:
        raise RuntimeError("Set ENTREZ_EMAIL in your .env file before running this script.")

    seq = Seq("ACGT")
    print(f"Hello from biopython! {seq}")
    print(f"Hello from biopython! {seq.complement()}")
    print(f"Hello from biopython! {seq.complement_rna()}")
    print(f"Hello from biopython! {seq.reverse_complement()}")

    print("\n\n")

    for seq_record in SeqIO.parse("./static/ls_orchid.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    print("\n\n")

    for seq_record in SeqIO.parse("./static/ls_orchid.gbk", "genbank"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    handle = Entrez.esearch(db="nucleotide", term="TP53[Gene] AND Homo sapiens[Organism]", retmax=1)

    record = cast(dict[str, Any], Entrez.read(handle))
    handle.close()

    print(record["IdList"])

    # DNA with high GC content is more stable
    print(f"My seq gc fraction: {gc_fraction(seq)}")
    print("Seq: %s" % seq)

    coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    template_dna = coding_dna.reverse_complement()
    messenger_rna = coding_dna.transcribe()
    protein_seq = messenger_rna.translate()

    print(f"Protein seq: {protein_seq}")
    # Translating nucleotides up to the first in frame stop codon and then stop
    print(
        f"Stop codon translation {
            coding_dna.translate(table='Vertebrate Mitochondrial', to_stop=True)
        }"
    )


if __name__ == "__main__":
    main()
