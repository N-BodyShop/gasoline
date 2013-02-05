#!/usr/bin/python
import os
import shutil
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
from testutil import run_gasoline

def plot_density(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['rho'], 'r.')
	plt.xlabel('position')
	plt.ylabel(r'$\rho$')

def plot_temperature(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['temp'], 'r.')
	plt.xlabel('position')
	plt.ylabel('temperature')

def plot_velocity(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['vx'], 'r.')
	plt.xlabel('position')
	plt.ylabel('velocity')
	
def plot_pressure(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['rho']*sim.g['temp'], 'r.')
	plt.xlabel('position')
	plt.ylabel(r'$\rho$')

def make_plots(testdir):
	for infile, step in [("shocktube.std", "000"), ("shocktube.00100", "100"), 
			("shocktube.00200", "200"), ("shocktube.00300", "300"),
			("shocktube.00400", "400")]:
		print "Plotting Step " + step
		plt.subplot(221)
		plot_density(testdir, infile, step)
		plt.subplot(222)
		plot_temperature(testdir, infile, step)
		plt.subplot(223)
		plot_velocity(testdir, infile, step)
		plt.subplot(224)
		plot_pressure(testdir, infile, step)
		plt.tight_layout(pad=0.1)
		plt.savefig(testdir+"/shocktube_step%s.png" %step, dpi=150)
		plt.clf()

if __name__ == '__main__':
	testdir='shocktube'
	files = ["data/shocktube.std", "data/shocktube.param"]
	exe = "../gasoline.ce92ce93"
	run_gasoline(testdir, files, 'shocktube.param', exe)
	make_plots(testdir)
