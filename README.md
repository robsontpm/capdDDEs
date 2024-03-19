RIGOROUS INTEGRATOR FOR DELAY DIFFERENTIAL EQUATIONS
====================================================

A library to do rigorous numerics for Delay Differential Equations
based on [^phd-thesis], [^rignum-ddes], [^high-order-ddes].

TL;DR
-----

If you have ```Debian``` based system (tested: Ubuntu 20.04),
then to see the project in action, you can run:

```bash
sudo apt-get install libgmp-dev libmpfr-dev libboost-all-dev git cmake autoconf libtool
mkdir capdDDEs
cd capdDDEs
git clone git@github.com:robsontpm/capdDDEs.git .
chmod a+x tldr.sh
./tldr.sh
```

It might take some time to compile everything. 
The next sections describes what is done in ```tldr.sh``` script. 

Requirements
------------

The requirements are the same as for [CAPD library](http://capd.ii.uj.edu.pl), which are:

 - gcc, g++ or clang++ (must support at least C++17 standard).
 - make, cmake (version 3.13.0 or newer)
 - pkg-config
 - git
 - gmp, mpfr (for multiprecision, some tools in capdDDEs uses those)
 - Boost library (header files and binaries)

In Debian based systems you can run:

```bash
sudo apt-get install libgmp-dev libmpfr-dev libboost-all-dev git cmake autoconf libtool
```

The documentation and the compilation system uses for now the ```g++``` compiler.
You can change it in **.mk** files in the root. 

_TODO: get the ```$CXX``` from environment vars (?)_

Compilation
-----------

Assuming you have all the needed requirements, 
to start using this project you can do the following steps 
(we assume you have -nix type of system, Windows is unsupported):

First clone the repo:

```bash
mkdir capdDDEs
cd capdDDEs
git clone git@github.com:robsontpm/capdDDEs.git .
```

Next, go to external, add permissions and build CAPD 
library (it is a long process):

```bash
cd external
chmod a+x capd/configure
chmod a+x capd-build.sh
./capd-build.sh
```

Then, after capd is build without errors, go back to root 
directory and compile some exemplary programs using make:

```bash
cd ..
make examples/mackey_glass_stable
```

Compiled programs will go to **./bin** directory:

```bash
ls bin/examples/
ls bin/examples/mackey_glass_stable

# sample output of the above ls commands are:
#
#    mackey_glass_stable
#
#    find_and_prove_n6
#
```

Finally, run proofs or demos, e.g. :

```bash
cd bin/examples/mackey_glass_stable
./find_and_prove_n6
```

You can see list of available programs and demos by running
in the main directory of the repository the following command:

```bash
make list
```

If you are interested in more documentation, see ./docs

Contact
-------

In case of any problems, please contact me at
robert.szczelina@uj.edu.pl - I will provide as much support as
possible, especially I will contact people from CAPD project if there
would be any problems with building it. Please attach as much
information as possible (e.g. output of build.sh / configure / make
programs),version of your OS, maybe some screenshots. 

Aditional programs (proofs) are usually attached to the recent 
manuscripts, see the References section below.


DESCRIPTION OF THIS PROJECT
---------------------------

The repository contains source codes for the rigorous 
integration library for Delay Differential Equations (DDEs) 
of the form[^1]:

```
x'(t) = f(x(t), x(t-d1), ..., x(t-dm))
```

Where d1,..., dm are m delays. It is build upon ideas
from CAPD library for ODEs (capd.ii.uj.edu.pl).
By rigorous methods we understand routines that compute the 
enclosure (interval bounds) on the **true** solutions to the 
problems. Those in turn might be applied in Computer Assisted Proofs,
that is, in mathematical proofs, that verify some assumptions
by appropriate computer programs.
Please reffer to [CAPD](capd.ii.uj.edu.pl) documentation for 
general papers on the topic. See **References** for the 
papers specyfic to DDEs [^1][^2][^3][^4].

This project requires the [CAPD library](capd.ii.uj.edu.pl),
but the repository contains compatible CAPD version
in './external/' directory. You can use your own version
of the CAPD, but in that case you need to modify ```CAPDBINDIR``` and ```CAPDLIBDIR``` variables in the 
template Makefile scripts ```build-conf.mk``` and ```build-conf-mp.mk```. 
As of 2022 the integrator uses CAPD DynSys version 5.1.2. 

_TODO: change it to have CAPD as a submodule (?) or downloaded from git separately._

The work is open-source, under GPLv3 license. Please consider
citing some of the papers listed in **References**, when using this 
work in any scientific project.

PROJECT STRUCTURE (UPDATED 2023)
--------------------------------

A short description of available folders / files. A detailed README.md
file can be found in some of them, when necessary. The project is in 
constant evolution, so changes can be made. When the core part of 
library would become part of CAPD, then the structure would be more
or less constant.

 - ```include/```
	As the library is mostly based on C++ templates, most
	of the code is here. 
 
 - ```include/ddes```
 	the core of the rigorous integrator for DDEs
 
 - ```include/ddeshelper```
	some extra routines to help with programs. 
	Might not be the part of the CAPD library, when 
	ddes are included. 
 
 - ```src/```
	As not all functions are templated, their implementation 
	need to go to .cpp files in here. 
 
 - ```programs/```
 	this is a directory to put any programs to be compiled with
 	the build-in compilation system. See **CUSTOM PROGRAMS** later
 	to learn how you can use this structure. 
 
 	_TODO: describe in more detail._
 
 - ```Makefile```
 	this file may be used to build programs inside the library.
 
 	Just call:
 
	```bash
	make _program-name_
	```
 
 	where _program-name_ is equal to the name of any directory under
 	```programs/```. The resulting executable file will go to the
 	bin/_program-name_ directory. 
 
	```bash
	make list
	```
 
 	will print out the list of available programs. 
 
 - ```external/``` 
 	You should build CAPD from this folder. 
 	The instructions are somwhere in this README.md file.  
 
 - ```scripts/```
 	some usefull scripts in languages other than C++ (bash / Python / tex)
 	to be used in some proofs to generate human-readable output
 	(e.g. the formulation of existence theorems with nice estimates)
 	This folder may not be present.
 
 - ```bin/```
 	this directory is initially empty (or non-existent). 
 	If you will compile any program which accompany the library 
 	(programs folder), then it will be most probably placed here. 
 	Also compiled version of bundled CAPD will be placed here. 
 
 - ```build-common.mk```
 	```build-conf.mk```
 	```build-conf-mp.mk```
 	```build-conf-svn.mk```
 
 	those files serve to configure compilation 
 	(make, Makefile files in programs, etc.), see above
 	for the instruction on compilation. The ```-mp``` version is used
 	for multiprecision compilation (see for example converter programs), ```-svn``` version 
 	is for using your own (external) compilation of CAPD.
 
 - ```README.md```
 	this file

The most important part is the rigorous integrator for DDEs. 
It is contained in ```include/capd/ddes```
As of 2022, I have decided to reorganize it to be more
complaint with the CAPD structure. I have followed example of ```external/capd/capdDynSys4/include/capd/pdes``` 
folder structure in this manner. I hope capdDDEs can become a part of CAPD library 
at some point, when the library is mature enough.

YOUR OWN CAPD LIBRARY
---------------------

If you want to use your own version of CAPD, e.g. from ```svn``` repo, 
then you need to modify file ```build-conf.mk``` 
(and other ```.mk``` probably) first. Follow the instructions in that 
file to make the ```Makefile``` work. 

In case of any errors in compiling CAPD library, please check requirements
on http://capd.ii.uj.edu.pl, or contact me at robert.szczelina@uj.edu.pl.


### CAPD v. 5 AND SHARED LIBRARIES (UPDATED 2022)

In a CAPD v5 there is a problem with compiled program, when the library
was installed localy. The system does not see the dynamically linked 
libraries. 

From 'make install' of CAPD you get info similar to that:

```
> Libraries have been installed in:
>  ~/workspace/capdDDEs5.1.2/bin/capd/lib
> 
> If you ever happen to want to link against installed libraries
> in a given directory, LIBDIR, you must either use libtool, and
> specify the full pathname of the library, or use the '-LLIBDIR'
> flag during linking and do at least one of the following:
>   - add LIBDIR to the 'LD_LIBRARY_PATH' environment variable
>     during execution
>   - add LIBDIR to the 'LD_RUN_PATH' environment variable
>     during linking
>   - use the '-Wl,-rpath -Wl,LIBDIR' linker flag
>   - have your system administrator add LIBDIR to '/etc/ld.so.conf'
>
> See any operating system documentation about shared libraries for
> more information, such as the ld(1) and ld.so(8) manual pages.
```

For now, I do the following:

 - create file '/etc/ld.so.conf.d/CAPDv5.conf' with the LIBDIR in it:
 
 	```sudo cat '/home/robson/workspace/capdDDEs5.1.2/bin/capd/lib' > /etc/ld.so.conf.d/CAPDv5.conf```
 
 - reload ldconfig:
 	
 	```sudo ldconfig```

The above solution requires admin rights. I think the following option: 

 - use the ```'-Wl,-rpath -Wl,LIBDIR'``` linker flag

	could work better, it can be implemented with ```Makefile```. 

	_TODO: Try to implement this method instead of ldconfig._

REFERENCES
----------

[^high-order-ddes]: Szczelina, R., Zgliczyński, P. 
  "High-Order Lohner-Type Algorithm for Rigorous Computation of Poincaré Maps in Systems of Delay Differential Equations with Several Delays", 
  Found Comput Math, Accepted, published online, (2023). 
  [doi:s10208-023-09614-x](https://doi.org/10.1007/s10208-023-09614-x).

[^rignum-ddes]: Szczelina, R., Zgliczyński, P. 
  "Algorithm for Rigorous Integration  of Delay Differential Equations and the Computer-Assisted Proof of Periodic Orbits in the Mackey–Glass Equation.", 
  Found Comput Math 18, 1299–1332 (2018). 
  [doi:10.1007/s10208-017-9369-5](https://doi.org/10.1007/s10208-017-9369-5)

[^logistic-dde]:	Szczelina, R. 
  "A computer assisted proof of multiple periodic orbits in some first order non-linear delay differential equation", 
  Electronic Journal of Qualitative Theory of Differential Equations, No. 83, 1-19, (2016) 
  [doi:10.14232/ejqtde.2016.1.83](https://doi.org/10.14232/ejqtde.2016.1.83)

[^phd-thesis]: Szczelina, R. 
  "Rigorous Integration of Delay Differential Equations", 
  PhD Thesis, (2015)