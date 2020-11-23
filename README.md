[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/samnooij/Rosalind_problems/master)

# Rosalind_problems
Scripts and notebooks used to solve exercises from rosalind.info

---

This is a collection of scripts I have written to solve the bioinformatics exercises at [rosalind.info](http://rosalind.info).
They are written in Python, using Jupyter notebooks to mix my personal notes and comments with code.

I anticipate this will mostly be a record for myself, but feel free to copy any part that is of use! 
This project is licenced under the [BSD-2-clause](LICENSE) licence.

---

_Update: 2020-11-23_

From here on, I would like to try a different workflow.  
Instead of developing the scripts as Jupyter notebook, I am going to use 
VSCodium to write the code as separate scripts/applications.

These scripts will then use the [Black](https://github.com/psf/black) code
format, and functions will be tested using `pytest`, which can also be
implemented as automated test on GitHub (with Actions).

Later on, I might add documentation (of the functions) using 
[Sphinx](https://www.sphinx-doc.org/en/master/) and 
[readthedocs](https://readthedocs.org/).  
And - but that might as well be in a separate repository -
put all useful functions together in one or a few scripts
as my own sort of bioinformatics toolbox.