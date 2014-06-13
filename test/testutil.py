#!/usr/bin/python
import os
import subprocess
import glob
import shutil
import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
import pynbody.plot.sph as p_sph


def run_gasoline(testdir, files, param, exe="gasoline", args="", params={}):
    try:
        os.mkdir(testdir)
    except OSError:
        print "Output directory already exists, not running gasoline"
        return 0
    for fname in files:
        shutil.copy(fname, testdir)
        shutil.copy(fname, testdir)
    shutil.copy(exe, testdir+'/gasoline')
    print "Starting gasoline..."
    os.chdir(testdir)
    if len(params) > 0:
        shutil.copy(param, "tmp.param")
        os.remove(param)
        paramlines = open("tmp.param").readlines()
        for line in paramlines():
            if line.split()[0] in params:
                line = params[line.split()[0]]
                del params[line.split()[0]]
        [paramlines.append(line) for line in params.itervalues()]
        outfile = open(param, 'w')
        [outfile.write("%s\n" % line) for line in paramlines]
        outfile.close()
        os.remove("tmp.param")
    exitcode = subprocess.call('./gasoline '+args+' '+param , shell=True)
    os.chdir('..')
    if exitcode != 0:
        print "gasoline failed with exit code: " + str(exitcode)
    return exitcode

def build_gasoline(makefile, obj="gasoline", mdl="pthread"):
    os.mkdir("buildtmp")
    shutil.copy(makefile, "buildtmp/Makefile")
    try:
        shutil.copytree("../mdl/", "buildtmp/mdl", ignore=shutil.ignore_patterns("*.o"))
    except OSError:
        shutil.copytree("../../mdl/", "buildtmp/mdl", ignore=shutil.ignore_patterns("*.o"))
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
    return exitcode

class test_result():
    def __init__(self, exitcode, name=""):
        if exitcode == 0:
            self.run = True
        else:
            self.run = False
        self.plots = []
        self.success = []
        self.value = []
        self.threshold = []
        self.description = []
        self.name = name
    def add_result(self, plot, success, value, threshold, description):
        self.plots.append(plot)
        self.success.append(success)
        self.value.append(value)
        self.threshold.append(threshold)
        self.description.append(description)
