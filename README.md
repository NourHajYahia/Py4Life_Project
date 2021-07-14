# Py4Life_Project

The main context of my project is about absorption
spectroscopy of DNA nucleotides. In our lab we try to Construct an efficient in
vitro expansion methods of DNA repeats, is paramount for understanding the
structures and functions of such sequences. 
In addition, expansions and deletion generally begins with realignment
of the primer/template complex leading to slipped-strand mispairing during DNA
replication. Studies that used different polymerases in vitro in order to
understand the expansion mechanism. Long
products of a polyG:polyC and polyA:polyT duplex has been synthesized using
Klenow fragment of DNA polymerase 1 by our lab. For this field of research we rely a lot on
spectroscopy absorption measurement of DNA nucleotide for concentration
evaluating at wavelength 260 nm in particular, as result of the huge reports
and exact calculations of molar absorption coefficient factor at this wavelength.

Considering this, it is also important to have efficient
methods for processing; parsing and analyzing this absorption spectra data,
with this python program can help a lot for faster and better building of the
needed primers.

The main goal of the proposed python script is to give a
right prediction of the type of given DNA nucleotide absorption spectra by
relying on fed reference absorptions spectra. 


Data:

Two types of data for the analysis introduced to the program
using python:

1. References data: which represent our predictive data, that we try to fit our sample to it and compare the similarity to our samples data, in other words to try to find the nearest reference to our sample.

2. Samples data: which represent our true data, the data we are trying to predict and calculating its concentration. In addition, finding the possible pairings of double strand in our data, which can construct polyApolyT or polyGpolyC primer. 
