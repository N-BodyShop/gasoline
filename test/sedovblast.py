#!/usr/bin/python
import os
import shutil
import pynbody as pyn
import matplotlib.pyplot as plt
import pynbody.plot.sph as p_sph
testdir='sedov'

def plot_density():
	rshock=[0.8995, 1.1864,1.3958]
	time_10 = pyn.load(testdir+'/sedov.00010')
	time_20 = pyn.load(testdir+'/sedov.00020')
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.plot(time_10['r'], time_10['rho'], ',r')
	plt.plot(time_20['r'], time_20['rho'], ',g')
	plt.plot(time_30['r'], time_30['rho'], ',b')
	for i in rshock:
		plt.axvline(i, color='black')
	plt.xlabel("Radius $(kpc)$")
	plt.ylabel("Density $(M_\odot/kpc^3)$")
	plt.xlim((0,2))
	plt.title(title)

def plot_density_image():
	rshock=[0.8995, 1.1864,1.3958]
	rshock10x = rshock[0]*np.cos(np.arange(0,2*np.pi,0.0001))
	rshock10y = rshock[0]*np.sin(np.arange(0,2*np.pi,0.0001))
	rshock20x = rshock[1]*np.cos(np.arange(0,2*np.pi,0.0001))
	rshock20y = rshock[1]*np.sin(np.arange(0,2*np.pi,0.0001))
	rshock30x = rshock[2]*np.cos(np.arange(0,2*np.pi,0.0001))
	rshock30y = rshock[2]*np.sin(np.arange(0,2*np.pi,0.0001))
	plt.set_cmap('gray')
	time_10 = pyn.load(testdir+'/sedov.00010')
	time_20 = pyn.load(testdir+'/sedov.00020')
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.subplot(311)
	p_sph.image(time_10)
	plt.plot(rshock10x, rshock10y, 'r.')
	plt.subplot(312)
	p_sph.image(time_20)
	plt.plot(rshock20x, rshock20y, 'r.')
	plt.subplot(313)
	p_sph.image(time_30)
	plt.plot(rshock30x, rshock30y, 'r.')
	plt.savefig('sedov_density_image.png')

def plot_velocity():
	time_10 = pyn.load(testdir+'/sedov.00010')
	time_20 = pyn.load(testdir+'/sedov.00020')
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.plot(time_10['r'], time_10['vr'], ',r')
	plt.plot(time_20['r'], time_20['vr'], ',g')
	plt.plot(time_30['r'], time_30['vr'], ',b')
	plt.xlabel("Radius $(kpc)$")
	plt.ylabel("Velocity $(km/s)$")
	plt.xlim((0,2))
	plt.ylim((-10,40))
	plt.title(title)

def plot_temperature():
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.semilogy(time_30['r'], time_30['temp'], ',k')
	plt.xlabel("Radius $(kpc)$")
	plt.ylabel("Temperature $(K)$")
	plt.xlim((0,2))
	plt.ylim((1,1e8))
	plt.title(title)

def plot_entropy():
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.semilogy(time_30['r'], np.power(time_30['temp'], 1.5)/time_30['rho'], '.k')
	plt.semilogy([0,0.1], [5.13e6, 5.13e6], 'r-')
	plt.xlabel("Radius $(kpc)$")
	plt.ylabel("$A(S) K^{1.5} cc/g$")
	plt.xlim((0,2))
	plt.ylim((6e3,6e7))
	plt.title(title)


def make_plots():
	plt.subplot(221)
	plot_density()
	plt.subplot(222)
	plot_temperature()
	plt.subplot(223)
	plot_velocity()
	plt.subplot(224)
	plot_entropy()
	plt.savefig(testdir+'/sedov.png')
	plot_density_image()

def run_gasoline(exe, args):
	try:
		os.mkdir(testdir)
	except OSError:
		print "Output directory already exists, not running gasoline"
		return
	shutil.copy('data/sedov.std', testdir)
	shutil.copy('data/sedov.param', testdir)
	shutil.copy(exe, testdir+'/gasoline')
	print "Starting gasoline..."
	os.chdir(testdir)
	exitcode = os.spawnl(os.P_WAIT, './gasoline', args, 'sedov.param')
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
