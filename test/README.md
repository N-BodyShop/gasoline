##Description##
The full suite in doall.py takes roughly 3 hours to complete on 1 core of a
reasonably new workstation.  Running multiple threads should decrease the run time,
especially for the Sedov blast

##Test Scripts##
* _shocktube.py_ A simple Sod Shock tube.  One of the simplest hydro tests.
* _sedovblast.py_ A Sedov blast. A harsher hydro test, also useful for testing
multistepping.
* _onestar.py_  A single star in a virialized ball of gas.  Checks that the
energy output from feedback is correct.
* _doall.py_ The full suite of tests, runs all of the gasoline test scripts.
##Utility Scripts
These scripts don't test gasoline's performance or accuracy, but are more useful
for general development tasks.
* _checkdefines.sh_ Checks that all of the defines throughout the sourcecode are
included in the Makefile.
