#========================================================================
# Makefile to compile FUSE
#========================================================================
#
#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides
F_MASTER = ${HOME}/Documents/analysis/diffModel/FUSE/source/fuse/

# Core directory that contains FUSE source code
F_KORE_DIR = $(F_MASTER)build/FUSE_SRC/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)build/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin/

#========================================================================
# PART 1: Define the libraries, driver programs, and executables
#========================================================================

.DEFAULT_GOAL := all

# Define the fortran compiler. 
#FC = ifort
FC = gfortran

# Define the NetCDF and HDF5 libraries. Use the libraries associated with
# the compiler you selected above. Note that these paths are machine-dependent
# - consider asking the person in charge of your machine where the libraries are
# installed. The folders should contain a lib directory - see LIBRARIES below.
ifeq "$(FC)" "ifort"
	NCDF_PATH = ${ISCA_NCDF_PATH_ifort}
	HDF_PATH  = ${ISCA_HDF_PATH_ifort}
endif

ifeq "$(FC)" "gfortran"
	INC_NETCDF := $(shell nf-config --fflags)
	LIB_NETCDF := $(shell nf-config --flibs) $(shell nc-config --libs)
endif

LIBRARIES = $(LIB_NETCDF)
INCLUDE = $(INC_NETCDF)

GENINC := $(F_MASTER)build/generated
VERSIONFILE := $(GENINC)/fuseversion.inc

# Use /bin/sh-compatible commands
VERSION    := $(shell git describe --tags --abbrev=0 2>/dev/null || echo "no-tag")
BUILDTIME  := $(shell date -u +"%Y-%m-%dT%H:%M:%SZ")
GITBRANCH  := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "detached")
GITHASH    := $(shell git rev-parse HEAD 2>/dev/null || echo "unknown")

GENINC := $(F_MASTER)build/generated
VERSIONFILE := $(GENINC)/fuseversion.inc

$(GENINC):
	mkdir -p $@

$(VERSIONFILE): | $(GENINC)
	@{ \
	  echo "! Auto-generated: do not edit"; \
	  echo "integer, parameter :: FUSE_VERSION_LEN = 64"; \
	  echo "integer, parameter :: FUSE_BUILDTIME_LEN = 32"; \
	  echo "integer, parameter :: FUSE_GITBRANCH_LEN = 64"; \
	  echo "integer, parameter :: FUSE_GITHASH_LEN = 64"; \
	  printf "character(len=FUSE_VERSION_LEN),   parameter :: FUSE_VERSION    = '%s'\n" "$(VERSION)"; \
	  printf "character(len=FUSE_BUILDTIME_LEN), parameter :: FUSE_BUILDTIME  = '%s'\n" "$(BUILDTIME)"; \
	  printf "character(len=FUSE_GITBRANCH_LEN), parameter :: FUSE_GITBRANCH  = '%s'\n" "$(GITBRANCH)"; \
	  printf "character(len=FUSE_GITHASH_LEN),   parameter :: FUSE_GITHASH    = '%s'\n" "$(GITHASH)"; \
	} > $@

#========================================================================
# PART 2: Assemble all of the FUSE sub-routines
#========================================================================

# Define directories
NUMREC_DIR  = $(F_KORE_DIR)numrec
HOOKUP_DIR  = $(F_KORE_DIR)hookup
DRIVER_DIR  = $(F_KORE_DIR)driver
NETCDF_DIR  = $(F_KORE_DIR)netcdf
DSHARE_DIR  = $(F_KORE_DIR)dshare
PRELIM_DIR  = $(F_KORE_DIR)prelim
RUNTIME_DIR = $(F_KORE_DIR)runtime
PHYSICS_DIR = $(F_KORE_DIR)physics
OLDPHYS_DIR = $(F_KORE_DIR)physics_orig
SOLVER_DIR  = $(F_KORE_DIR)solver_orig
UTILMS_DIR  = $(F_KORE_DIR)util
SCE_DIR     = $(F_KORE_DIR)sce
TIME_DIR    = $(F_KORE_DIR)netcdf

# Define the executables
DRIVER_EX = fuse.exe

# Define the driver program and associated subroutines for the fidelity test
FUSE_DRIVER = \
	get_fuse_prelim.f90 \
	fuse_rmse.f90 functn.f90 \
	sce_driver.f90 \
  fuse_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(FUSE_DRIVER))

# Manager modules
FUSE_HOOKUP= \
	kinds_dmsl_kit_FUSE.f90 \
	utilities_dmsl_kit_FUSE.f90 \
	fuse_fileManager.f90
HOOKUP = $(patsubst %, $(HOOKUP_DIR)/%, $(FUSE_HOOKUP))

# Numerical Recipes utilities
FUSE_NRUTIL= \
	nrtype.f90 \
	nr.f90 nrutil.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NRUTIL))

# Data modules
FUSE_DATAMS=    \
  model_defn.f90 \
	data_types.f90 \
	model_defnames.f90 \
	globaldata.f90 \
	multiconst.f90 \
	multiforce.f90 \
	multibands.f90 \
	multiparam.f90 \
	multistate.f90 \
	multi_flux.f90 \
	multiroute.f90 \
	multistats.f90 \
	model_numerix.f90
DATAMS = $(patsubst %, $(DSHARE_DIR)/%, $(FUSE_DATAMS))

# Time I/O modules
FUSE_TIMEMS=    \
  time_io.f90
TIMUTILS = $(patsubst %, $(TIME_DIR)/%, $(FUSE_TIMEMS))

# Utility modules
FUSE_UTILMS= \
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
UTILMS = $(patsubst %, $(UTILMS_DIR)/%, $(FUSE_UTILMS))

# Numerical Recipes
FUSE_NR_SUB= \
	ludcmp.f90 lubksb.f90 svbksb.f90 svdcmp.f90 pythag.f90 \
	gammln.f90 gammp.f90 gcf.f90 gser.f90
NR_SUB = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NR_SUB))

# FUSE physics (differentiable model)
FUSE_PHYSICS= \
	smoothers.f90 \
	get_parent.f90 \
	update_swe_diff.f90 \
	qsatexcess_diff.f90 \
	evap_upper_diff.f90 \
	evap_lower_diff.f90 \
	qinterflow_diff.f90 \
	qpercolate_diff.f90 \
	q_baseflow_diff.f90 \
	q_misscell_diff.f90 \
	mstate_rhs_diff.f90 \
	mod_derivs_diff.f90 \
	conserve_clamp.f90 \
	fix_ovshoot.f90 \
	implicit_solve.f90
PHYSICS = $(patsubst %, $(PHYSICS_DIR)/%, $(FUSE_PHYSICS))

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
	mstate_eqn.f90 \
	fix_states.f90 \
	meanfluxes.f90 \
	wgt_fluxes.f90 \
	updatstate.f90 \
	q_overland.f90
MODGUT = $(patsubst %, $(OLDPHYS_DIR)/%, $(FUSE_MODGUT))

# Solver
FUSE_SOLVER= \
	interfaceb.f90 \
	limit_xtry.f90 \
	viol_state.f90 \
	fuse_deriv.f90 \
	fmin.f90 fdjac_ode.f90 flux_deriv.f90 disaggflux.f90 \
	fuse_sieul.f90 \
	newtoniter.f90 lnsrch.f90
SOLVER = $(patsubst %, $(SOLVER_DIR)/%, $(FUSE_SOLVER))

# Define routines for FUSE preliminaries
FUSE_PRELIM= \
	parse_command_args.f90 \
	ascii_util.f90 \
	uniquemodl.f90 \
	getnumerix.f90 \
	force_info.f90 \
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
PRELIM = $(patsubst %, $(PRELIM_DIR)/%, $(FUSE_PRELIM))

FUSE_RUNTIME= \
	conv_funcs.f90 \
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
RUNTIME = $(patsubst %, $(RUNTIME_DIR)/%, $(FUSE_RUNTIME))

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
FUSE_ALL = $(HOOKUP) $(NRUTIL) $(DATAMS) $(TIMUTILS) $(UTILMS) \
	$(NR_SUB) $(PHYSICS) $(MODGUT) $(SOLVER) $(PRELIM) $(RUNTIME) \
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

 FLAGS_NORMA = -O3 -ffree-line-length-none -fmax-errors=0 -cpp
 #FLAGS_DEBUG = -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds -cpp
 FLAGS_DEBUG = -O0 -ffree-line-length-none -fmax-errors=0 -fcheck=bounds -cpp # without -warn all
 FLAGS_FIXED = -O2 -c -ffixed-form

endif

# select a set of flags
#FLAGS = $(FLAGS_NORMA) -I$(GENINC)
FLAGS = $(FLAGS_DEBUG) -I$(GENINC)

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

# Compile
all: compile install clean

# compile FUSE into a static library
compile: sce_16plus.o | $(GENINC)
	@rm -f $(VERSIONFILE)
	@$(MAKE) $(VERSIONFILE)
	$(FC) $(FLAGS) $(FUSE_ALL) $(DRIVER) \
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

# describe how to compile SCE code written in Fortran 77
sce_16plus.o: $(SCE_DIR)/sce_16plus.f
	$(FC) $(FLAGS_FIXED) $(SCE_DIR)/sce_16plus.f
