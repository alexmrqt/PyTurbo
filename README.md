# PyTurbo

A work-in-progress set of Python class (actually, Python wrappers for C++ class)
implementing several basic turbo-algorithms (turbo-decoding, turbo-equalization, etc.).

## Currelently Implements
* The Viterbi Algorithm
* The Log BCJR Algorithm (sometimes referred as log-MAP or log-forward/backward algorithm).
 
# Installation
## Dependencies
You will need Cython, as well as a C++ compiler.

## Command line
```
python3 setup.py install
```

# Based on
* Viterbi algorithm implementation is taken from the gr-lazyviterbi GNURadio OOT module (https://github.com/alexmrqt/gr-lazyviterbi).
* Trellis description is taken for the gr-trellis module of GNURadio (https://github.com/gnuradio/gnuradio).
