README - Network-based Frequency Analysis (NFA)

Derrible S, Ahmad N (2015) Network-Based and Binless Frequency Analyses. PLoS ONE 10(11): e0142108. doi:10.1371/journal.pone.0142108 Available at: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142108

version 1.10.1 - for pandas 0.18.0 and later.

See version 1.00 if you have a pandas version earlier than 0.18.0.


What’s new

- minor bug fix in the NFA_functions.py code with new pandas version. Added “.iloc” after “self” on lines 27, 28, 37.

- minor bug fix in NFA_main.py. The headers do not have to be numbers now to be able to plot the final figure. Changed “int” to “eval” at line 170. 

- added “# -*- coding: utf-8 -*-“ on top of every code to prevent some bugs when using some characters.
