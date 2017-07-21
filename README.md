```
 	  ▄████  ▄▄▄        ██████  ▒█████   ██▓     ██▓ ███▄    █ ▓█████
	 ██▒ ▀█▒▒████▄    ▒██    ▒ ▒██▒  ██▒▓██▒    ▓██▒ ██ ▀█   █ ▓█   ▀
	▒██░▄▄▄░▒██  ▀█▄  ░ ▓██▄   ▒██░  ██▒▒██░    ▒██▒▓██  ▀█ ██▒▒███   
	░▓█  ██▓░██▄▄▄▄██   ▒   ██▒▒██   ██░▒██░    ░██░▓██▒  ▐▌██▒▒▓█  ▄
	░▒▓███▀▒ ▓█   ▓██▒▒██████▒▒░ ████▓▒░░██████▒░██░▒██░   ▓██░░▒████▒
	 ░▒   ▒  ▒▒   ▓▒█░▒ ▒▓▒ ▒ ░░ ▒░▒░▒░ ░ ▒░▓  ░░▓  ░ ▒░   ▒ ▒ ░░ ▒░ ░
	  ░   ░   ▒   ▒▒ ░░ ░▒  ░ ░  ░ ▒ ▒░ ░ ░ ▒  ░ ▒ ░░ ░░   ░ ▒░ ░ ░  ░
	░ ░   ░   ░   ▒   ░  ░  ░  ░ ░ ░ ▒    ░ ░    ▒ ░   ░   ░ ░    ░   
	      ░       ░  ░      ░      ░ ░      ░  ░ ░           ░    ░  ░
```                                                                  

#  Gasoline:  Particle Hydrodynamics Have Never Been Smoother

[![Build Status](http://imp.mcmaster.ca:8080/buildStatus/icon?job=gasoline-n-bodyshop/master)](http://imp.mcmaster.ca:8080/job/gasoline-n-bodyshop/job/master/)

Gasoline is a modern TREESPH code for solving the equations of gravity and
hydrodynamics in astrophysical problems.  To use the code, please check out the
[home page](http://gasoline-code.com) or the documentation on the
[wiki](https://github.com/N-BodyShop/gasoline/wiki).  If you run into a problem,
please open an [issue](https://github.com/N-BodyShop/gasoline/issues).  

# Quickstart

* Pull gasoline and the [MDL](https://github.com/N-BodyShop/mdl) helper library:
```
	git clone git@github.com:N-BodyShop/mdl.git
	git clone git@github.com:N-BodyShop/gasoline.git
```
* Copy the `Makefile.default` to `Makefile` and edit it to enable the feature
  you want:
```
	cp Makefile.default Makefile
```
* Compile the code with `make ARCH`, where `ARCH` is the architecture you want
  to use (for example, basic MPI):
```
	make mpi
```
* Run the code on for your simulation, using the parameter file as the argument
  (again, with MPI as the example):
```
	mpirun -np 1024 /path/to/gasoline_repo/gasoline simulation.param
```


# License

Gasoline is released under a GNU Public License v2.  Please take a look at
LICENSE.md for the terms.
