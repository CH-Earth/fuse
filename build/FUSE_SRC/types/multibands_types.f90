MODULE multibands_types

 ! Created by Brian Henn to allow multi-band snow modeling, 6/2013
 ! Based on module MULTIFORCE by Martyn Clark

 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026

 USE nrtype

 implicit none
 private

 public :: BANDS, BANDS_INFO, BANDS_VAR


 TYPE BANDS ! for catchment scale modeling
  INTEGER(I4B)                         :: NUM             ! band number (-)
  REAL(SP)                             :: Z_MID           ! band mid-point elevation (m)
  REAL(SP)                             :: AF              ! fraction of basin area in band (-)
  REAL(SP)                             :: SWE             ! band snowpack water equivalent (mm)
  REAL(SP)                             :: SNOWACCMLTN     ! new snow accumulation in band (mm day-1)
  REAL(SP)                             :: SNOWMELT        ! snowmelt in band (mm day-1)
  REAL(SP)                             :: DSWE_DT         ! rate of change of band SWE (mm day-1)
 ENDTYPE BANDS

 ! for distributed modeling MBANDS is split between time-independent and time-dependent charactertistics

 TYPE BANDS_INFO ! invariant characteristics
  REAL(SP)                             :: Z_MID           ! band mid-point elevation (m)
  REAL(SP)                             :: AF              ! fraction of basin area in band (-)
 ENDTYPE BANDS_INFO

 TYPE BANDS_VAR ! time-dependent characteristics
  REAL(SP)                             :: SWE             ! band snowpack water equivalent (mm)
  REAL(SP)                             :: SNOWACCMLTN     ! new snow accumulation in band (mm day-1)
  REAL(SP)                             :: SNOWMELT        ! snowmelt in band (mm day-1)
  REAL(SP)                             :: DSWE_DT         ! rate of change of band SWE (mm day-1)
 ENDTYPE BANDS_VAR

END MODULE multibands_types
