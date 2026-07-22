from Bio.SeqRecord import SeqRecord


def trim_ends(record: SeqRecord, min_quality: int = 20) -> SeqRecord:
    qualities = record.letter_annotations["phred_quality"]

    start = 0
    end = len(qualities)

    while start < end and qualities[start] < min_quality:
        start += 1

    while end > start and qualities[end - 1] < min_quality:
        end -= 1

    return record[start:end]
