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
	plt.savefig(testdir+"/onestar_density_step"+step+".png")
	plt.clf() 

def plot_entropy(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], np.power(sim.g['temp'], 1.5)/sim.g['rho'], 'r.')
	plt.xlabel('position')
	plt.ylabel(r'S')
	plt.savefig(testdir+"/onestar_entropy_step"+step+".png")
	plt.clf() 

def plot_temperature(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['temp'], 'r.')
	plt.xlabel('position')
	plt.ylabel('temperature')
	plt.savefig(testdir+"/onestar_temperature_step"+step+".png")
	plt.clf() 

def plot_velocity(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['vx'], 'r.')
	plt.xlabel('position')
	plt.ylabel('velocity')
	plt.savefig(testdir+"/onestar_velocity_step"+step+".png")
	plt.clf() 
	
def plot_pressure(testdir, infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['rho']*sim.g['temp'], 'r.')
	plt.xlabel('position')
	plt.ylabel(r'$\rho$')
	plt.savefig(testdir+"/onestar_pressure_step"+step+".png")
	plt.clf() 

def make_plots(testdir):
	for infile, step in [("onestar.std", "0"), ("onestar.00100", "100"), 
			("onestar.00200", "200"), ("onestar.00300", "300"),
			("onestar.00400", "400")]:
		print "Plotting Step " + step
		plot_density(testdir, infile, step)
		plot_entropy(testdir, infile, step)
		plot_temperature(testdir, infile, step)
		plot_velocity(testdir, infile, step)
		plot_pressure(testdir, infile, step)

if __name__ == '__main__':
	testdir='onestar'
	files = ["data/onestar.tbin", "data/onestar.param"]
	exe = "./gasoline"
	run_gasoline(testdir, files, 'onestar.param', exe)
	make_plots(testdir)
