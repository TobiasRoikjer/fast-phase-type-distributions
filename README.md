# Fast, scalable graph-based phase-type distributions for populations genetics
This library implements phase-type distributions in population genetics used for my masters thesis in Bioinformatics. [Link to pdf](thesis.pdf). It implements algorithms I have developed for a graph-based definitions of phase-type distributions such that properties and transformations can be computed efficiently. From this, we can scale the applications to a larger sample size and more population genetic models.

Note that this repository is a) there is no focus on freeing memory or failing gracefully (yet) and b) the command line part is not mature at all

## Building
This project uses CMake to build and compile. It requires the GNU scintific library to be installed (GSL).

## Platforms
This proejct supports Linux, OSX, and Windows Subsystem for Linux, and likely Windows as well, depending on the compiler used.

## Contact
tobiasgitemail@gmail.com
