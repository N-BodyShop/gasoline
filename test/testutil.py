#!/usr/bin/python
import os
import subprocess
import glob
import shutil
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
import pynbody.plot.sph as p_sph


def run_gasoline(testdir, files, param, exe="gasoline", args="", paramdict={}):
	try:
		os.mkdir(testdir)
	except OSError:
		print "Output directory already exists, not running gasoline"
		return
	for fname in files:
		shutil.copy(fname, testdir)
		shutil.copy(fname, testdir)
	shutil.copy(exe, testdir+'/gasoline')
	print "Starting gasoline..."
	os.chdir(testdir)
	if len(paramdict) > 0:
		shutil.copy(param, "tmp.param")
		os.remove(param)
		paramlines = open("tmp.param").readlines()
		for line in paramlines():
			if line.split()[0] in paramdict:
				line = paramdict[line.split()[0]]
				del paramdict[line.split()[0]]
		[paramlines.append(line) for line in paramdict.itervalues()]
		outfile = open(param, 'w')
		[outfile.write("%s\n" for line in paramlines]
		outfile.close()
		os.remove("tmp.param")
	exitcode = subprocess.call('./gasoline '+args+' '+param , shell=True)
	os.chdir('..')
	if exitcode != 0:
		print "gasoline failed with exit code: " + str(exitcode)

def build_gasoline(makefile, obj="gasoline", mdl="pthread"):
	os.mkdir("buildtmp")
	shutil.copy(makefile, "buildtmp/Makefile")
	shutil.copytree("../mdl/", "buildtmp/mdl", ignore=shutil.ignore_patterns("*.o"))
	for i in glob.iglob("../*c"):
		shutil.copy(i, "buildtmp/")
	for i in glob.iglob("../*h"):
		shutil.copy(i, "buildtmp/")
	os.chdir("buildtmp")
	exitcode = subprocess.call(['make', mdl])
	if exitcode != 0:
		print "make failed with exit code: " + str(exitcode)
	shutil.copy(obj, "../gasoline")
	os.chdir("..")
	shutil.rmtree("buildtmp")
