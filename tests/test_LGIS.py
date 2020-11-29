import sys

sys.path.insert(1, "./")

from scripts.LGIS_Longest_increasing_subsequence import (
    read_length_and_sequence,
    longest_subsequence,
)

# function import method thanks to addicted on stackoverflow:
# https://stackoverflow.com/questions/4383571/importing-files-from-different-folder#comment104909466_40612922


def test_read_length_and_sequence():
    with open("data/Example_longest_increasing_subsequence.txt", "r") as read_file:
        length = int(read_file.readline())
        sequence = list(read_file.readline().split())

    assert length == 5
    assert sequence == ["5", "1", "4", "2", "3"]


def test_longest_subsequence():
    test_sequence = ["5", "1", "4", "2", "3"]

    assert longest_subsequence(
        seq=test_sequence, mode="strictly", order="increasing"
    ) == ["1", "2", "3"]
    assert longest_subsequence(
        seq=test_sequence, mode="strictly", order="decreasing"
    ) == ["5", "4", "3"] or longest_subsequence(
        seq=test_sequence, mode="strictly", order="decreasing"
    ) == [
        "5",
        "4",
        "2",
    ]
    # In this example, there are two possibilities for the longest decreasing subsequence.

