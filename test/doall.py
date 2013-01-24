#!/usr/bin/python
import os
import time
import shutil
import argparse
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
from testutil import run_gasoline, build_gasoline

#Test scripts
import shocktube
import sedovblast
import onestar
def parse_param(paramfile):
	params = open(paramfile).readlines()
	makefile = ''
	args = ''
	paramdict = {}
	for i in params:
		if i.split()[0] == 'MAKEFILE':
			makefile = i.split()[-1]
		elif i.split()[0] == 'ARGS':
			rags = i.split()[-1]
		else:
			paramdict[i.split()[0]] = i
			
	if makefile == '':
		print "No makefile specified by paramfile!"
		exit(1)
	return (makefile, args, paramdict)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run the gasoline test suite.')
	parser.add_argument("paramfile", nargs='?', default='data/default.test',
			help="The makefile you want to build gasoline with")
	args = parser.parse_args()
	exe = "./gasoline"
	makefile, args, paramdict = parse_param(args.paramfile)
	start_time = time.time()
#Build gasoline
	build_gasoline(makefile, exe)
#Run the shock tube
	testdir='shocktube'
	files = ["data/shocktube.std", "data/shocktube.param"]
	run_gasoline(testdir, files, 'shocktube.param', exe, args=args, params=paramdict)
	#shocktube.make_plots(testdir)
#Run the Sedov blast
	testdir='sedov'
	files = ["data/sedov.std", "data/sedov.param"]
	run_gasoline(testdir, files, 'sedov.param', exe, args=args, params=paramdict)
	sedovblast.make_plots(testdir)
#Run the onestar
	testdir='onestar'
	files = ["data/onestar.tbin", "data/onestar.param"]
	run_gasoline(testdir, files, 'onestar.param', exe, args=args, params=paramdict)
	onestar.make_plots(testdir)
	end_time = time.time()
	print "\n\n\nFinished all Tests."
	print "Total runtime was %d seconds" % int(end_time-start_time)
