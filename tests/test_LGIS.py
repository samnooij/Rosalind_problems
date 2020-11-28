import sys

sys.path.insert(1, "./")

from scripts.LGIS_Longest_increasing_subsequence import read_length_and_sequence

# function import method thanks to addicted on stackoverflow:
# https://stackoverflow.com/questions/4383571/importing-files-from-different-folder#comment104909466_40612922


def test_read_length_and_sequence():
    with open("data/Example_longest_increasing_subsequence.txt", "r") as read_file:
        length = int(read_file.readline())
        sequence = list(read_file.readline().split())

    assert length == 5
    assert sequence == ["5", "1", "4", "2", "3"]
