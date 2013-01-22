#!/usr/bin/python
import os
import time
import shutil
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
from testutil import run_gasoline, build_gasoline

#Test scripts
import shocktube
import sedovblast
import onestar

if __name__ == '__main__':
	exe = "./gasoline"
	start_time = time.time()
#Build gasoline
	build_gasoline('data/Makefile.gas', exe)
#Run the shock tube
	testdir='shocktube'
	files = ["data/shocktube.std", "data/shocktube.param"]
	run_gasoline(testdir, files, 'shocktube.param', exe, args="-sz 1")
	#shocktube.make_plots(testdir)
#Run the Sedov blast
	testdir='sedov'
	files = ["data/sedov.std", "data/sedov.param"]
	run_gasoline(testdir, files, 'sedov.param', exe, args="-sz 1")
	sedovblast.make_plots(testdir)
#Onestar needs the starformation stuff, so let's build it
	build_gasoline('data/Makefile.stars', exe)
#Run the onestar
	testdir='onestar'
	files = ["data/onestar.tbin", "data/onestar.param"]
	run_gasoline(testdir, files, 'onestar.param', exe, args="-sz 1")
	onestar.make_plots(testdir)
	end_time = time.time()
	print "\n\n\nFinished all Tests."
	print "Total runtime was %d seconds" % int(end_time-start_time)
