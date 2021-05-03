from pathlib import Path
import sys

sys.path.insert(1, "./")

from scripts.LONG_Genome_assembly_as_shortest_superstring import (
    read_sequences,
    find_overlaps,
)

sample_dictionary = {
    "Rosalind_56": "ATTAGACCTG",
    "Rosalind_57": "CCTGCCGGAA",
    "Rosalind_58": "AGACCTGCCG",
    "Rosalind_59": "GCCGGAATAC",
}


def test_read_sequences():
    # Make sure the sample dataset works as intended
    sequence_dictionary = read_sequences(
        "data/Example_genome_assembly_as_shortest_superstring.fasta"
    )

    assert len(sequence_dictionary) == 4
    assert sequence_dictionary == sample_dictionary


def test_find_overlaps():
    # Make sure the sample dataset works as intended
    assert find_overlaps(sample_dictionary) == "ATTAGACCTGCCGGAATAC"


def test_rosalind_exercise():
    # Make sure that the script produces the right answer as judged by Rosalind
    test_file = "data/rosalind_long5.txt"
    answer_file = "results/long5.txt"

    with open(answer_file, "r") as read_file:
        answer = read_file.readline().strip()

    sequence_dictionary = read_sequences(test_file)
    assert find_overlaps(sequence_dictionary) == answer
