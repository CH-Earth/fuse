#========================================================================
# Makefile to compile FUSE
#========================================================================
#
#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides
F_MASTER = your/path/to/fuse/repository/

# Core directory that contains FUSE source code
F_KORE_DIR = $(F_MASTER)build/FUSE_SRC/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin/

#========================================================================
# PART 1: Define the libraries, driver programs, and executables
#========================================================================

# Define the fortran compiler. You have a few options: i) set the FC
# variable in your environment, ii) set it when you compile this
# Make file (e.g. make FC=ifort), iii) or if don't define it, the compiler
# specified below is used
ifndef FC
	#FC = ifort
	FC = gfortran
endif

# Define the NetCDF and HDF5 libraries. Use the libraries associated with
# the compiler you selected above. Note that these paths are machine-dependent
# - consider asking the person in charge of your machine where the libraries are
# installed. The folders should contain a lib directory - see LIBRARIES below.
ifeq "$(FC)" "ifort"
	NCDF_PATH = ${ISCA_NCDF_PATH_ifort}
	HDF_PATH  = ${ISCA_HDF_PATH_ifort}
endif

NCDFF_LIB_PATH = /your/path/to/netcdf-fortran
NCDF_LIB_PATH = /your/path/to/netcdf
HDF_LIB_PATH = /your/path/to/hdf5

LIBRARIES = -L$(NCDFF_LIB_PATH)/lib -lnetcdff -L$(NCDF_LIB_PATH)/lib -lnetcdf -L$(HDF_LIB_PATH)/lib -lhdf5_hl -lhdf5
INCLUDE = -I$(NCDFF_LIB_PATH)/include -I$(NCDF_LIB_PATH)/include -I$(HDF_LIB_PATH)/include

# Define the driver program and associated subroutines for the fidelity test
FUSE_DRIVER = \
  sobol.f90 \
  fuse_metric.f90 \
  functn.f90 \
  fuse_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(FUSE_DRIVER))

# Define the executables
DRIVER_EX = fuse.exe

#========================================================================
# PART 2: Assemble all of the FUSE sub-routines
#========================================================================

# Define directories
NUMREC_DIR = $(F_KORE_DIR)FUSE_NR
HOOKUP_DIR = $(F_KORE_DIR)FUSE_HOOK
DRIVER_DIR = $(F_KORE_DIR)FUSE_DMSL
NETCDF_DIR = $(F_KORE_DIR)FUSE_NETCDF
ENGINE_DIR = $(F_KORE_DIR)FUSE_ENGINE
SCE_DIR    = $(F_KORE_DIR)FUSE_SCE
TIME_DIR   = $(F_KORE_DIR)FUSE_TIME

# Utility modules
FUSE_UTILMS= \
	kinds_dmsl_kit_FUSE.f90 \
	utilities_dmsl_kit_FUSE.f90 \
	fuse_fileManager.f90
UTILMS = $(patsubst %, $(HOOKUP_DIR)/%, $(FUSE_UTILMS))

# Numerical Recipes utilities
FUSE_NRUTIL= \
	nrtype.f90 \
	nr.f90 nrutil.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NRUTIL))

# Data modules
FUSE_DATAMS=    \
  model_defn.f90 \
	model_defnames.f90 \
	multiconst.f90 \
	multiforce.f90 \
	multibands.f90 \
	multiparam.f90 \
	multistate.f90 \
	multi_flux.f90 \
	multiroute.f90 \
	multistats.f90 \
	model_numerix.f90
DATAMS = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_DATAMS))

# Time I/O modules
FUSE_TIMEMS=    \
  time_io.f90
TIMUTILS = $(patsubst %, $(TIME_DIR)/%, $(FUSE_TIMEMS))

# Information modules
FUSE_INFOMS= \
	metaoutput.f90 \
	metaparams.f90 \
	meta_stats.f90 \
	selectmodl.f90 \
	putpar_str.f90 \
	getpar_str.f90 \
	par_insert.f90 \
	parextract.f90 \
	varextract.f90 \
	sumextract.f90 \
	str_2_xtry.f90 \
	xtry_2_str.f90
INFOMS = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_INFOMS))

# Numerical Recipes
FUSE_NR_SUB= \
	ludcmp.f90 lubksb.f90 svbksb.f90 svdcmp.f90 pythag.f90 \
	gammln.f90 gammp.f90 gcf.f90 gser.f90
NR_SUB = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NR_SUB))

# Model guts
FUSE_MODGUT=\
	mod_derivs.f90 \
	update_swe.f90 \
	qrainerror.f90 \
	qsatexcess.f90 \
	evap_upper.f90 \
	evap_lower.f90 \
	qinterflow.f90 \
	qpercolate.f90 \
	q_baseflow.f90 \
	q_misscell.f90 \
	logismooth.f90 \
	mstate_eqn.f90 \
	fix_states.f90 \
	meanfluxes.f90 \
	wgt_fluxes.f90 \
	updatstate.f90 \
	q_overland.f90
MODGUT = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_MODGUT))

# Solver
FUSE_SOLVER= \
	interfaceb.f90 \
	limit_xtry.f90 \
	viol_state.f90 \
	fuse_deriv.f90 \
	fmin.f90 fdjac_ode.f90 flux_deriv.f90 disaggflux.f90 \
	fuse_sieul.f90 \
	newtoniter.f90 lnsrch.f90
SOLVER = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_SOLVER))

# Define routines for FUSE preliminaries
FUSE_PRELIM= \
	ascii_util.f90 \
	uniquemodl.f90 \
	getnumerix.f90 \
	getparmeta.f90 \
	assign_stt.f90 \
	assign_flx.f90 \
	assign_par.f90 \
	adjust_stt.f90 \
	par_derive.f90 \
	bucketsize.f90 \
	mean_tipow.f90 \
	qbsaturatn.f90 \
	qtimedelay.f90 \
	init_stats.f90 \
	init_state.f90
PRELIM = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_PRELIM))

FUSE_MODRUN= \
        metrics.f90 \
	conv_funcs.f90 \
	force_info.f90 \
	clrsky_rad.f90 \
	getPETgrid.f90 \
	get_mbands.f90 \
	get_time_indices.f90\
	initfluxes.f90 \
	set_all.f90 \
	ode_int.f90 \
	fuse_solve.f90 \
	comp_stats.f90 \
	mean_stats.f90
MODRUN = $(patsubst %, $(ENGINE_DIR)/%, $(FUSE_MODRUN))

# Define NetCDF routines
FUSE_NETCDF = \
    handle_err.f90 \
    extractor.f90 juldayss.f90 caldatss.f90 \
    get_gforce.f90 \
    get_smodel.f90 \
    get_fparam.f90 \
    def_params.f90 \
    def_output.f90 \
    def_sstats.f90 \
    put_params.f90 \
    put_output.f90 \
    put_sstats.f90
NETCDF = $(patsubst %, $(NETCDF_DIR)/%, $(FUSE_NETCDF))

SCE = \
	sce_16plus.o

# ... and stitch it all together...
FUSE_ALL = $(UTILMS) $(NRUTIL) $(DATAMS) $(TIMUTILS) $(INFOMS) \
	$(NR_SUB) $(MODGUT) $(SOLVER) $(PRELIM) $(MODRUN) \
	$(NETCDF) $(SCE)

#========================================================================
# PART 3: Compile the puppy
#========================================================================

# Define flags - these are compiler-dependent
ifeq "$(FC)" "ifort"
 FLAGS_NORMA = -O3 -FR -auto -fltconsistency -fpe0 -fpp
 FLAGS_DEBUG = -O0 -p -g -debug -warn all -check all -FR -auto -WB -traceback -fltconsistency -fpe0 -fpp
 FLAGS_DEBUG = -O0 -p -g -debug -check all -FR -auto -WB -traceback -fltconsistency -fpe0 -fpp # without -warn all
 FLAGS_FIXED = -O2 -c -fixed
endif
ifeq "$(FC)" "gfortran"
 FLAGS_NORMA = -O3 -ffree-line-length-none -fmax-errors=0 -cpp -fallow-argument-mismatch
 FLAGS_DEBUG = -p -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds -cpp
 FLAGS_DEBUG = -p -Og -ffree-line-length-none -fmax-errors=0 -fcheck=bounds -cpp # without -warn all
 FLAGS_FIXED = -O2 -c -ffixed-form
endif

# select a set of flags
FLAGS = $(FLAGS_NORMA)

# MPI: FUSE with MPI has been compiled successfully with mpif90 and mpiifort.
# Note: override must be specifed to enable FC passed as argument to be overridden
# https://www.gnu.org/software/make/manual/html_node/Override-Directive.html#Override-Directive

ifeq "$(MODE)" "distributed"

	ifeq "$(FC)" "ifort"
	 override FC = mpiifort
	endif

	ifeq "$(FC)" "gfortran"
		override FC = mpif90
	endif

	MPI_FLAGS = -D__MPI__
	DRIVER_EX = fuse_mpi.exe

endif

# Compile SCE code written in Fortran 77
sce_16plus.o: $(SCE_DIR)/sce_16plus.f
	$(FC) -O2 -c -ffixed-form -o $@ $<

# Compile
all: compile install clean

# compile FUSE into a static library
compile: sce_16plus.o
	$(FC) -O3 -ffree-line-length-none -fmax-errors=0 -cpp -fallow-argument-mismatch \
	$(FUSE_ALL) $(DRIVER) \
	$(LIBRARIES) $(INCLUDE) $(MPI_FLAGS) -o $(DRIVER_EX)

# Remove object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *__genmod.f90

# Copy the executable to the bin directory
install:
	mkdir -p $(EXE_PATH)
	mv $(DRIVER_EX) $(EXE_PATH)
