import os

from Bio import Entrez
from dotenv import load_dotenv

from sequence import run_sequence_examples

_ = load_dotenv()

Entrez.email = os.getenv("ENTREZ_EMAIL", "")


def main() -> None:
    if not Entrez.email:
        raise RuntimeError("Set ENTREZ_EMAIL in your .env file before running this script.")

    run_sequence_examples()


if __name__ == "__main__":
    main()
