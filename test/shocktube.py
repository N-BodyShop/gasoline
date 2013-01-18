#!/usr/bin/python
import os
import shutil
import pynbody as pyn
import matplotlib.pyplot as plt
testdir='shocktube'

def plot_density(infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['rho'], 'r.')
	plt.xlabel('position')
	plt.ylabel(r'$\rho$')
	plt.savefig(testdir+"/shocktube_density_step"+step+".png")
	plt.clf() 

def plot_temperature(infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['temp'], 'r.')
	plt.xlabel('position')
	plt.ylabel('temperature')
	plt.savefig(testdir+"/shocktube_temperature_step"+step+".png")
	plt.clf() 

def plot_velocity(infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['vx'], 'r.')
	plt.xlabel('position')
	plt.ylabel('velocity')
	plt.savefig(testdir+"/shocktube_velocity_step"+step+".png")
	plt.clf() 
	
def plot_pressure(infile, step):
	sim = pyn.load(testdir+"/"+infile)
	plt.plot(sim.g['x'], sim.g['rho']*sim.g['temp'], 'r.')
	plt.xlabel('position')
	plt.ylabel(r'$\rho$')
	plt.savefig(testdir+"/shocktube_pressure_step"+step+".png")
	plt.clf() 

def make_plots():
	for infile, step in [("shocktube.std", "0"), ("shocktube.00100", "100"), 
			("shocktube.00200", "200"), ("shocktube.00300", "300"),
			("shocktube.00400", "400")]:
		print "Plotting Step " + step
		plot_density(infile, step)
		plot_temperature(infile, step)
		plot_velocity(infile, step)
		plot_pressure(infile, step)

def run_gasoline(exe, args):
	try:
		os.mkdir(testdir)
	except OSError:
		print "Output directory already exists, not running gasoline"
		return
	shutil.copy('data/shocktube.std', testdir)
	shutil.copy('data/shocktube.param', testdir)
	shutil.copy(exe, testdir+'/gasoline')
	print "Starting gasoline..."
	os.chdir(testdir)
	exitcode = os.spawnl(os.P_WAIT, './gasoline', args, 'shocktube.param')
	os.chdir('..')
	if exitcode != 0:
		print "gasoline failed with exit code: " + str(exitcode)

def readconfig():
	output = {'EXE':'gasoline', 'ARGS':''}
	configfile = open('test.cfg')
	for line in configfile:
		if line[0] == '#':
			continue
		if line.split()[0] == 'EXE:':
			output['EXE'] = line.split()[1]
		if line.split()[0] == 'ARGS:':
			output['ARGS'] = line.split()[1]
	return output

if __name__ == '__main__':
	config = readconfig()
	run_gasoline(config['EXE'], config['ARGS'])
	make_plots()
