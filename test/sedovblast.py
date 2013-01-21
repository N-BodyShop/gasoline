#!/usr/bin/python
import os
import shutil
from testutil import run_gasoline
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
import pynbody.plot.sph as p_sph

def plot_density(testdir):
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

def plot_density_image(testdir):
	plt.ion()
	rshock=[0.8995, 1.1864,1.3958]
	rshock10x = rshock[0]*np.cos(np.arange(0,2*np.pi,0.0001))
	rshock10y = rshock[0]*np.sin(np.arange(0,2*np.pi,0.0001))
	rshock20x = rshock[1]*np.cos(np.arange(0,2*np.pi,0.0001))
	rshock20y = rshock[1]*np.sin(np.arange(0,2*np.pi,0.0001))
	rshock30x = rshock[2]*np.cos(np.arange(0,2*np.pi,0.0001))
	rshock30y = rshock[2]*np.sin(np.arange(0,2*np.pi,0.0001))
	time_10 = pyn.load(testdir+'/sedov.00010')
	time_20 = pyn.load(testdir+'/sedov.00020')
	time_30 = pyn.load(testdir+'/sedov.00030')
	f = plt.figure(figsize=(9,3))
	p = f.add_subplot(131)
	p_sph.image(time_10, width=4, show_cbar=False, subplot=p, clear=False, cmap='gray')
	plt.plot(rshock10x, rshock10y, 'r,')
	p = f.add_subplot(132)
	p_sph.image(time_20, width=4, show_cbar=False, subplot=p, clear=False, cmap='gray')
	plt.plot(rshock20x, rshock20y, 'r,')
	p = f.add_subplot(133)
	p_sph.image(time_30, width=4, show_cbar=False, subplot=p, clear=False, cmap='gray')
	plt.plot(rshock30x, rshock30y, 'r,')
	plt.tight_layout(pad=0.1)
	plt.savefig(testdir+'/sedov_density_image.png', dpi=100)

def plot_velocity(testdir):
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

def plot_temperature(testdir):
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.semilogy(time_30['r'], time_30['temp'], '.k')
	plt.xlabel("Radius $(kpc)$")
	plt.ylabel("Temperature $(K)$")
	plt.xlim((0,2))
	plt.ylim((1,1e8))

def plot_entropy(testdir):
	time_30 = pyn.load(testdir+'/sedov.00030')
	plt.semilogy(time_30['r'], np.power(time_30['temp'], 1.5)/time_30['rho'], '.k')
	plt.semilogy([0,0.1], [5.13e6, 5.13e6], 'r-')
	plt.xlabel("Radius $(kpc)$")
	plt.ylabel("$A(S) K^{1.5} cc/g$")
	plt.xlim((0,2))
	plt.ylim((6e3,6e7))


def make_plots(testdir):
	plt.subplot(221)
	plot_density(testdir)
	plt.subplot(222)
	plot_temperature(testdir)
	plt.subplot(223)
	plot_velocity(testdir)
	plt.subplot(224)
	plot_entropy(testdir)
	plt.tight_layout(pad=0.1)
	plt.savefig(testdir+'/sedov.png', dpi=150)
	plot_density_image(testdir)

if __name__ == '__main__':
	testdir='sedov'
	files = ["data/sedov.std", "data/sedov.param"]
	exe = "../gasoline.ce92ce93"
	run_gasoline(testdir, files, 'sedov.param', exe)
	make_plots(testdir)
