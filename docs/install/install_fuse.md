# FUSE Installation

We have successfully installed FUSE on a number of Unix-like (\*nix) operating
systems, including Linux and Darwin (Mac OS X). 

To compile FUSE, change into the `build/` directory inside your FUSE installation and run `make`:
```
cd /path/to/fuse/build
make
```

Test by:
`/path/to/fuse/bin/fuse.exe`

# Dependencies

To compile FUSE, you will need:
 
### A Fortran compiler

  We have successfully used the intel Fortran compiler (`ifort`, version 17.x) and the GNU Fortran compiler (`gfortran`, version 6 or higher), the latter of which is freely available. Since we do not use any compiler-specific extensions, you should be able to compile FUSE with other Fortran compilers as well. If you do not have a Fortran compiler, you can install `gfortran` for free. The easiest way is to use a package manager (e.g., Homebrew). Note that `gfortran` is installed as part of the `gcc` compiler suite (for Homebrew, `brew install gcc`).

### The NetCDF libraries

  [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) or the Network Common Data Format is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. For Homebrew, you can install NetCDF library as `brew install netcdf-fortran`. Most \*nix package managers include a NetCDF port. Note that you need to ensure that:

  - You have NetCDF version 4.x;
  - The NetCDF libraries are compiled with the same compiler as you plan to use for compiling FUSE (if you installed NetCDF via Homebrew, and you compile FUSE using the Homebrew gfortran, you’re almost always consistent); and
  - You have the NetCDF Fortran library installed (`libnetcdff.*`) and not just the C-version.
 
### A copy of the FUSE source code from [this repo](https://github.com/CH-EARTH/fuse)

  You have a number of options:

  - If you just want to use the latest stable release of FUSE, then simply look for the [latest release](https://github.com/CH-EARTH/fuse/releases);
  - If you want the latest and greatest (and potentially erroneous), download a copy of the [development branch](https://github.com/CH-EARTH/fuse/tree/develop) (or clone it);
  - If you may want to do FUSE development, then fork the repo on github and start editing your own copy.

### pkg-config

  `pkg-config` is a command-line tool that helps software builds find the right compiler and linker flags for installed libraries (like HDF5, netCDF, etc.). After it’s installed, you can use `pkg-config` in build systems (Makefiles, CMake, configure scripts) to automatically discover the correct -I include paths and -L/-l library flags, instead of you having to set those paths manually. In FUSE `pkg-config` is used in the Makefile.
