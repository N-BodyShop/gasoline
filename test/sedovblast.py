#!/usr/bin/python
import os
import shutil
from testutil import run_gasoline, test_result
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
import pynbody.plot.sph as p_sph

#Parameters for the test
testdir='sedov'
files = ["data/sedov.std", "data/sedov.param"]
exe = "./gasoline"
testname = "Sedov Blast"

def test_radius():
    description = """ 
    This test ensures that the blastwave is growing at the
    correct rate (that the radius at times 10, 20, 30 are correct).  Each peak
    density must be within 10% of the theoretical result on average to pass."""
    success = False
    rshock=[0.8995, 1.1864,1.3958]
    threshold = 0.1
    time_10 = pyn.load(testdir+'/sedov.00010')
    time_20 = pyn.load(testdir+'/sedov.00020')
    time_30 = pyn.load(testdir+'/sedov.00030')
    fig = plt.figure()
    plt.plot(time_10['r'], time_10['rho'], ',r', figure=fig)
    plt.plot(time_20['r'], time_20['rho'], ',g', figure=fig)
    plt.plot(time_30['r'], time_30['rho'], ',b', figure=fig)
    peak_10 = time_10['r'][np.where(time_10['rho'] == np.max(time_10['rho']))]
    peak_20 = time_20['r'][np.where(time_20['rho'] == np.max(time_20['rho']))]
    peak_30 = time_30['r'][np.where(time_30['rho'] == np.max(time_30['rho']))]
    value = 1.0-np.mean((peak_10/rshock[0], peak_20/rshock[1], peak_30/rshock[2]))
    if np.abs(value) < threshold:
        success = True
        plt.title('PASS', color='green')
    else:
        plt.title('FAIL', color='red')
    for i in rshock:
        plt.axvline(i, color='black', linewidth=2)
    plt.xlabel("Radius $(kpc)$")
    plt.ylabel("Density $(M_\odot/kpc^3)$")
    plt.xlim((0,2))
    return (fig, success, value, threshold, description)


def test_entropy():
    description = """ 
    This test ensures that average entropy of the inner 8 particles is within
    10% of the theoretical value at the final step."""
    success = False
    entropy = 5.13e6
    threshold = 0.1
    time_30 = pyn.load(testdir+'/sedov.00030')
    fig = plt.figure()
    plt.semilogy(time_30['r'], np.power(time_30['temp'], 1.5)/time_30['rho'], '.k')
    plt.semilogy([0,0.1], [5.13e6, 5.13e6], 'r-')
    center_idx = np.argsort(time_30['r'])[:8]
    value = 1-entropy/np.mean(np.power(time_30['temp'][center_idx], 1.5)/time_30['rho'][center_idx])
    if np.abs(value) < threshold:
        success = True
        plt.title('PASS', color='green')
    else:
        plt.title('FAIL', color='red')
    plt.xlabel("Radius $(kpc)$")
    plt.ylabel("$A(S) K^{1.5} cc/g$")
    plt.xlim((0,1))
    plt.ylim((6e3,6e7))
    return (fig, success, value, threshold, description)


def run():
    run = run_gasoline(testdir, files, 'sedov.param', exe)
    results = test_result(run, testname)
    if run == 0:
        results.add_result(*test_radius())
        results.add_result(*test_entropy())
    return results

if __name__ == '__main__':
    results = run()
    if results.run:
        for i in range(len(results.success)):
            print "TEST:"
            print results.description[i]
            print "Pass\tValue\tThreshold\t"
            print "%r\t%e\t%e\t" % (results.success[i], results.value[i], results.threshold[i])
            results.plots[i].savefig(testname+str(i)+".png")
