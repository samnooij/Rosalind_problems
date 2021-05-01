# Rosalind problem scripts

These are the Python scripts that I used to solve the 
exercises on [rosalind.info](https://www.rosalind.info).

_Note: before 'LGIS' (#24) I used to work with Jupyter notebooks._  
_After that, I started developing separate Python scripts that_
_are unit-tested with pytest._

## 24. LGIS: Longest Increasing Subsequence

I started looking for this term online and found that it seems to
be a well-known computer science problem 
([LIS on wikipedia](https://en.wikipedia.org/wiki/Longest_increasing_subsequence)).
There is even a code example on Wikipedia, which is not in Python,
so I cannot copy-paste it.

When looking a bit further, I found a wonderful Python implementation by
[arekolek](https://stackoverflow.com/users/1916449/arekolek) on Stackoverflow: 
https://stackoverflow.com/a/38337443.  
That function is so nice that I could not have done it better myself and
I copied it into my own script. It also solves both increasing and decreasing
subsequences!
And it has a nice docstring from which I can learn a thing or two.

## 25. LONG: Genome Assembly as Shortest Superstring

This is also a common problem for which different code solutions
are available online. See for example:
https://leetcode.com/problems/find-the-shortest-superstring/solution/

I decided to try and make my own implementation this time.
Not because I think I can make a faster, more optimised solution,
but because I want to learn how to do this and understand
exactly how the code works. That way, I might also be able to easily add
alternative functionality.

_Suggested bonus functionality:_
 1. handle reverse complement sequences
 2. set a minimum length for overlaps
 3. handle datasets with not exclusively overlapping sequences
    (return multiple sequences)