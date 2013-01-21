#!/usr/bin/python
import os
import shutil
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
from testutil import run_gasoline

#Test scripts
import shocktube
import sedovblast
import onestar

if __name__ == '__main__':
#Run the shock tube
	testdir='shocktube'
	files = ["data/shocktube.std", "data/shocktube.param"]
	exe = "../gasoline.ce92ce93"
	run_gasoline(testdir, files, 'shocktube.param', exe)
	shocktube.make_plots(testdir)
#Run the Sedov blast
	testdir='sedov'
	files = ["data/sedov.std", "data/sedov.param"]
	exe = "../gasoline.ce92ce93"
	run_gasoline(testdir, files, 'sedov.param', exe)
	sedovblast.make_plots(testdir)
#Run the onestar
	testdir='onestar'
	files = ["data/onestar.tbin", "data/onestar.param"]
	exe = "../gasoline.ce92ce93"
	run_gasoline(testdir, files, 'onestar.param', exe)
	onestar.make_plots(testdir)
