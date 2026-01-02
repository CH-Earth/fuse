#========================================================================
# Makefile to compile FUSE
#========================================================================

.DEFAULT_GOAL := all

#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides
FUSE_ROOT = $(abspath $(CURDIR)/..)

# Core directory that contains FUSE source code
FUSE_SOURCE_DIR = $(FUSE_ROOT)/build/FUSE_SRC/

# Location of the compiled modules
MOD_PATH = $(FUSE_ROOT)/build/

# Define the directory for the executables
EXE_PATH = $(FUSE_ROOT)/bin/

#========================================================================
# PART 1: Define the libraries, driver programs, and executables
#========================================================================

# default Fortran compiler is set to `gfortran`
# other options: ifort
FC = gfortran

# find HDF5
# Check if pkg-config is available for HDF5
# Check if all required environment variables are set
ifneq ($(HDF5_LIB_DIR)$(HDF5_INCLUDE_DIR),)
  # If all variables are set, proceed
  :
else
  # check if the `pkg-config' is available
  HAS_PKG_CONFIG := $(shell pkg-config --exists hdf5 && echo yes)

  ifeq ($(HAS_PKG_CONFIG),yes)
    HDF5_INCLUDE_DIR := $(shell pkg-config --variable=includedir hdf5)
    HDF5_LIB_DIR := $(shell pkg-config --variable=libdir hdf5)
  else
    $(error "pkg-config for HDF5 is not available. Set HDF5_INCLUDE_DIR and HDF5_LIB_DIR environment variables manually")
  endif # pkg-config ... hdf5
endif # HDF5

# netcdf-c
# Use the nc-config to set the proper flags
# Check if nc-config is available
ifneq ($(NETCDF_C_PREFIX),)
    NETCDF_C_INCLUDE := $(NETCDF_C_PREFIX)/include
    NETCDF_C_LIB := $(NETCDF_C_PREFIX)/lib
else
  ifneq ($(shell which nc-config 2>/dev/null),)
    NETCDF_C_INCLUDE := $(shell nc-config --includedir)
    NETCDF_C_LIB := $(shell nc-config --libdir)
  else
    $(error nc-config not found. Please install netcdf-c library, or set NETCDF_C_PREFIX manually)
  endif
endif

# netcdf-fortran
# Use the nf-config to set the proper flags
# Check if nf-config is available
ifneq ($(NETCDF_F_PREFIX),)
    NETCDF_F_INCLUDE := $(NETCDF_F_PREFIX)/include
    NETCDF_F_LIB := $(NETCDF_F_PREFIX)/lib
else
  ifneq ($(shell which nf-config 2>/dev/null),)
    # Seems problematic for arm64 macs
    # NETCDF_F_INCLUDE := $(shell nf-config --includedir)
    # NETCDF_F_LIB := $(shell nf-config --prefix)/lib
    # Using pkgconf (pkg-config) for the time being
    NETCDF_F_INCLUDE := $(shell pkgconf --variable=includedir netcdf-fortran)
    NETCDF_F_LIB := $(shell pkgconf --variable=libdir netcdf-fortran)
  else
    $(error nf-config not found. Please install netcdf-fortran library, or set NETCDF_F_PREFIX manually)
  endif
endif

# define appropriate flags
INCLUDES += -I$(HDF5_INCLUDE_DIR) -I$(NETCDF_C_INCLUDE) -I$(NETCDF_F_INCLUDE)
LIBS += -L$(HDF5_LIB_DIR) -lhdf5 -lhdf5_hl -L$(NETCDF_F_LIB) -lnetcdff -L$(NETCDF_C_LIB) -lnetcdf

$(info INCLUDES are $(INCLUDES))
$(info LIBS are $(LIBS))

#========================================================================
# PART 2: Assemble all of the FUSE sub-routines
#========================================================================

# Define directories
NUMREC_DIR  = $(FUSE_SOURCE_DIR)numrec
HOOKUP_DIR  = $(FUSE_SOURCE_DIR)hookup
DRIVER_DIR  = $(FUSE_SOURCE_DIR)driver
NETCDF_DIR  = $(FUSE_SOURCE_DIR)netcdf
DSHARE_DIR  = $(FUSE_SOURCE_DIR)dshare
PRELIM_DIR  = $(FUSE_SOURCE_DIR)prelim
RUNTIME_DIR = $(FUSE_SOURCE_DIR)runtime
PHYSICS_DIR = $(FUSE_SOURCE_DIR)physics
OLDPHYS_DIR = $(FUSE_SOURCE_DIR)physics_orig
SOLVER_DIR  = $(FUSE_SOURCE_DIR)solver_orig
UTILMS_DIR  = $(FUSE_SOURCE_DIR)util
SCE_DIR     = $(FUSE_SOURCE_DIR)sce
TIME_DIR    = $(FUSE_SOURCE_DIR)netcdf

# Define the executables
DRIVER_EX = fuse.exe

# Define the driver program and associated subroutines
FUSE_DRIVER =
#FUSE_DRIVER += setup_domain.f90
#FUSE_DRIVER += setup_model_definition.f90
FUSE_DRIVER += fuse_metric.f90 functn.f90
#FUSE_DRIVER += sce_driver.f90
FUSE_DRIVER += fuse_driver.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(FUSE_DRIVER))

# Manager modules
FUSE_HOOKUP =
FUSE_HOOKUP += kinds_dmsl_kit_FUSE.f90
FUSE_HOOKUP += utilities_dmsl_kit_FUSE.f90
HOOKUP = $(patsubst %, $(HOOKUP_DIR)/%, $(FUSE_HOOKUP))

# Numerical Recipes utilities
FUSE_NRUTIL =
FUSE_NRUTIL += nrtype.f90
FUSE_NRUTIL += nr.f90 nrutil.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NRUTIL))

# Data modules
FUSE_DATAMS =
FUSE_DATAMS += model_defn.f90
#FUSE_DATAMS += data_types.f90
FUSE_DATAMS += model_defnames.f90
#FUSE_DATAMS += globaldata.f90
FUSE_DATAMS += multiconst.f90
FUSE_DATAMS += multiforce.f90
FUSE_DATAMS += multibands.f90
FUSE_DATAMS += multiparam.f90
FUSE_DATAMS += multistate.f90
FUSE_DATAMS += multi_flux.f90
FUSE_DATAMS += multiroute.f90
FUSE_DATAMS += multistats.f90
FUSE_DATAMS += model_numerix.f90
DATAMS = $(patsubst %, $(DSHARE_DIR)/%, $(FUSE_DATAMS))

# Time I/O modules
FUSE_TIMEMS =
FUSE_TIMEMS += time_io.f90
TIMUTILS = $(patsubst %, $(TIME_DIR)/%, $(FUSE_TIMEMS))

# Utility modules
FUSE_UTILMS =
FUSE_UTILMS += fuse_fileManager.f90
#FUSE_UTILMS += alloc_domain.f90
#FUSE_UTILMS += alloc_scratch.f90
FUSE_UTILMS += metaoutput.f90
FUSE_UTILMS += metaparams.f90
FUSE_UTILMS += meta_stats.f90
FUSE_UTILMS += selectmodl.f90
FUSE_UTILMS += putpar_str.f90
FUSE_UTILMS += getpar_str.f90
FUSE_UTILMS += par_insert.f90
FUSE_UTILMS += parextract.f90
FUSE_UTILMS += varextract.f90
FUSE_UTILMS += sumextract.f90
FUSE_UTILMS += str_2_xtry.f90
FUSE_UTILMS += xtry_2_str.f90
UTILMS = $(patsubst %, $(UTILMS_DIR)/%, $(FUSE_UTILMS))

# Numerical Recipes
FUSE_NR_SUB =
FUSE_NR_SUB += ludcmp.f90 lubksb.f90 svbksb.f90 svdcmp.f90 pythag.f90
FUSE_NR_SUB += gammln.f90 gammp.f90 gcf.f90 gser.f90
NR_SUB = $(patsubst %, $(NUMREC_DIR)/%, $(FUSE_NR_SUB))

# FUSE physics (differentiable model)
FUSE_PHYSICS =
FUSE_PHYSICS += smoothers.f90
FUSE_PHYSICS += get_parent.f90
FUSE_PHYSICS += update_swe_diff.f90
FUSE_PHYSICS += qsatexcess_diff.f90
FUSE_PHYSICS += evap_upper_diff.f90
FUSE_PHYSICS += evap_lower_diff.f90
FUSE_PHYSICS += qinterflow_diff.f90
FUSE_PHYSICS += qpercolate_diff.f90
FUSE_PHYSICS += q_baseflow_diff.f90
FUSE_PHYSICS += q_misscell_diff.f90
FUSE_PHYSICS += mstate_rhs_diff.f90
FUSE_PHYSICS += mod_derivs_diff.f90
FUSE_PHYSICS += conserve_clamp.f90
FUSE_PHYSICS += fix_ovshoot.f90
FUSE_PHYSICS += implicit_solve.f90
PHYSICS = $(patsubst %, $(PHYSICS_DIR)/%, $(FUSE_PHYSICS))

# Old physics (model "guts")
FUSE_MODGUT =
FUSE_MODGUT += mod_derivs.f90
FUSE_MODGUT += update_swe.f90
FUSE_MODGUT += qrainerror.f90
FUSE_MODGUT += qsatexcess.f90
FUSE_MODGUT += evap_upper.f90
FUSE_MODGUT += evap_lower.f90
FUSE_MODGUT += qinterflow.f90
FUSE_MODGUT += qpercolate.f90
FUSE_MODGUT += q_baseflow.f90
FUSE_MODGUT += q_misscell.f90
FUSE_MODGUT += mstate_eqn.f90
FUSE_MODGUT += fix_states.f90
FUSE_MODGUT += meanfluxes.f90
FUSE_MODGUT += wgt_fluxes.f90
FUSE_MODGUT += updatstate.f90
FUSE_MODGUT += q_overland.f90
MODGUT = $(patsubst %, $(OLDPHYS_DIR)/%, $(FUSE_MODGUT))

# Old solver
FUSE_SOLVER =
FUSE_SOLVER += interfaceb.f90
FUSE_SOLVER += limit_xtry.f90
FUSE_SOLVER += logismooth.f90
FUSE_SOLVER += viol_state.f90
FUSE_SOLVER += fuse_deriv.f90
FUSE_SOLVER += fmin.f90 fdjac_ode.f90 flux_deriv.f90 disaggflux.f90
FUSE_SOLVER += fuse_sieul.f90
FUSE_SOLVER += newtoniter.f90 lnsrch.f90
SOLVER = $(patsubst %, $(SOLVER_DIR)/%, $(FUSE_SOLVER))

# FUSE preliminaries
FUSE_PRELIM =
#FUSE_PRELIM += parse_command_args.f90
FUSE_PRELIM += ascii_util.f90
FUSE_PRELIM += uniquemodl.f90
FUSE_PRELIM += getnumerix.f90
FUSE_PRELIM += force_info.f90
FUSE_PRELIM += getparmeta.f90
FUSE_PRELIM += assign_stt.f90
FUSE_PRELIM += assign_flx.f90
FUSE_PRELIM += assign_par.f90
FUSE_PRELIM += adjust_stt.f90
FUSE_PRELIM += par_derive.f90
FUSE_PRELIM += bucketsize.f90
FUSE_PRELIM += mean_tipow.f90
FUSE_PRELIM += qbsaturatn.f90
FUSE_PRELIM += qtimedelay.f90
FUSE_PRELIM += init_stats.f90
FUSE_PRELIM += init_state.f90
PRELIM = $(patsubst %, $(PRELIM_DIR)/%, $(FUSE_PRELIM))

# FUSE runtime
FUSE_RUNTIME =
FUSE_RUNTIME += metrics.f90
FUSE_RUNTIME += conv_funcs.f90
FUSE_RUNTIME += clrsky_rad.f90
FUSE_RUNTIME += getPETgrid.f90
#FUSE_RUNTIME += get_time_windows.f90
FUSE_RUNTIME += get_time_indices.f90
FUSE_RUNTIME += initfluxes.f90
FUSE_RUNTIME += set_all.f90
FUSE_RUNTIME += ode_int.f90
FUSE_RUNTIME += fuse_solve.f90
FUSE_RUNTIME += comp_stats.f90
FUSE_RUNTIME += mean_stats.f90
RUNTIME = $(patsubst %, $(RUNTIME_DIR)/%, $(FUSE_RUNTIME))

# NetCDF routines
FUSE_NETCDF =
FUSE_NETCDF += handle_err.f90
FUSE_NETCDF += extractor.f90 juldayss.f90 caldatss.f90
#FUSE_NETCDF += domain_decomp.f90
FUSE_NETCDF += get_gforce.f90
FUSE_NETCDF += get_mbands.f90
FUSE_NETCDF += get_smodel.f90
FUSE_NETCDF += get_fparam.f90
FUSE_NETCDF += def_params.f90
FUSE_NETCDF += def_output.f90
FUSE_NETCDF += def_sstats.f90
FUSE_NETCDF += put_params.f90
FUSE_NETCDF += put_output.f90
FUSE_NETCDF += put_sstats.f90
NETCDF = $(patsubst %, $(NETCDF_DIR)/%, $(FUSE_NETCDF))

SCE = sce_16plus.o

# ... and stitch it all together...
FUSE_ALL = 
FUSE_ALL += $(HOOKUP)
FUSE_ALL += $(NRUTIL)
FUSE_ALL += $(DATAMS)
FUSE_ALL += $(UTILMS)
FUSE_ALL += $(TIMUTILS)
FUSE_ALL += $(NR_SUB)
#FUSE_ALL += $(PHYSICS)
FUSE_ALL += $(MODGUT)
FUSE_ALL += $(SOLVER)
FUSE_ALL += $(PRELIM)
FUSE_ALL += $(RUNTIME)
FUSE_ALL += $(NETCDF) 
FUSE_ALL += $(SCE)

#=====================
# PART 3: Compile fuse
#=====================

# Define flags based on specified compiler
ifeq ($(FC),ifort)
  FFLAGS_NORMA = -O3 -FR -auto -fltconsistency -fpe0 -fpp
  FFLAGS_DEBUG = -O0 -g -debug -warn all -check all -FR -auto -WB -traceback -fltconsistency -fpe0 -fpp
  FFLAGS_FIXED = -O2 -c -fixed
endif

ifeq ($(FC),gfortran)
  FFLAGS_NORMA = -O3 -ffree-line-length-none -fmax-errors=0 -cpp
  FFLAGS_DEBUG = -g -Wall -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=bounds -cpp
  FFLAGS_FIXED = -O2 -c -ffixed-form
endif

# Default flags
FFLAGS = $(FFLAGS_NORMA)

# Target-specific flags for 'debug' target
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: compile

# Special provision for gcc>13
ifeq ($(FC),gfortran)

	GFORTRAN_VERSION := $(shell gfortran -dumpversion | cut -d. -f1)
  $(info compiler version is $(GFORTRAN_VERSION))
  GFORTRAN_GT_13 := $(shell expr $(GFORTRAN_VERSION) \>= 13)

  ifeq ($(GFORTRAN_GT_13),1)
    FFLAGS += -fallow-argument-mismatch
  endif
endif

# MPI: FUSE with MPI has been compiled successfully with mpif90 and mpiifort.
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

# compile target
compile: sce_16plus.o
	$(FC) $(FUSE_ALL) $(DRIVER) \
	$(FFLAGS) $(LIBS) $(INCLUDES) -o $(DRIVER_EX)

# Remove object files
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *__genmod.f90

# Copy the executable to the bin directory
install:
	mkdir -p $(EXE_PATH)
	mv $(DRIVER_EX) $(EXE_PATH)

.PHONY: all compile install clean debug
