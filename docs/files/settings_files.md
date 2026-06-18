The `settings` directory must contain the following files (provided for the catchment case study):

1. A file called `input_info.txt` providing metadata for the NetCDF input file. It defines the name and units of the variables in the input file.

1. The file `M_DECISIONS` (called `fuse_zDecisions_902.txt` in the case studies) describes the different options available in the FUSE modeling framework. These modeling decisions are described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). See below for a preliminary overview that maps code onto equations in the paper.

2. The file `CONSTRAINTS` (called `fuse_zConstraints_snow.txt` in the case studies) defines in particular the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).

3. The file `MOD_NUMERIX` (called `fuse_zNumerix.txt` in the case studies) defines decisions regarding the numerical solution technique. Examples of the impact of these decisions are described by [Clark and Kavetski (WRR 2010)](http://dx.doi.org/10.1029/2009WR008894) and [Kavetski and Clark (WRR 2010)](http://dx.doi.org/10.1029/2009WR008896).

```
FUSE has a number of modeling decisions, with the following options:

- RFERR: rainfall error
	- additive_e: [Eq. n/a*]       additive rainfall error
	- multiplc_e: [Eq. n/a*]       multiplicative rainfall error
 
- ARCH1: upper-zone architecture
	- onestate_1: [Eq. 1a, 2b]     upper layer defined by a single state variable, enables lower zone evap
	- tension1_1: [Eq. 1b, 2b]     upper layer broken up into tension and free storage, enables lower zone evap
	- tension2_1: [Eq. 1c, 2a]     tension storage sub-divided into recharge and excess, disables lower zone evap

- ARCH2: lower-zone architecture
	- unlimfrc_2: [Eq. 2a/2b, 6a]  baseflow reservoir of unlimited size (0-HUGE), fractional rate
	- tens2pll_2: [Eq. 2c,    6b]  tension reservoir plus two parallel tanks; baseflow from two tanks only
	- fixedsiz_2: [Eq. 2a/2b, 6c]  baseflow reservoir of fixed size, power recession
	- topmdexp_2: [Eq. n/a**]      baseflow reservoir of fixed size, exponential recession
	- unlimpow_2: [Eq. 2a/2b, 6d]  baseflow reservoir of unlimited size (0-HUGE), power recession

- QSURF: surface runoff
	- prms_varnt: [Eq. 9a]         PRMS variant (fraction of upper tension storage)
	- arno_x_vic: [Eq. 9b]         ARNO/Xzang/VIC parameterization (upper zone control)
	- tmdl_varnt: [Eq. 9c]         TOPMODEL parameterization (only valid for TOPMODEL qb)

- QPERC: percolation
	- perc_f2sat: [Eq. 4a]         water from (field capacity to saturation) available for percolation
	- perc_w2sat: [Eq. 4b]         water from (wilting point to saturation) available for percolation
	- perc_lower: [Eq. 4c, 4d]     perc defined by moisture content in lower layer (SAC)

- ESOIL: evaporation
	- sequential: [Eq. 3a, 3b]     sequential evaporation model
	- rootweight: [Eq. 3c, 3d]     root weighting

- QINTF: interflow
	- intflwnone: [Eq. 5a]         no interflow
	- intflwsome: [Eq. 5b]         interflow

- Q_TDH: time delay
	- no_routing: [Eq. n/a*]       no routing
	- rout_gamma: [Eq. 13a, 13b]   use a Gamma distribution with shape parameter = 2.5

- SNOWM: snow model
	- no_snowmod: [Eq. n/a***]     no snow model
	- temp_index: [Eq. n/a***]     temperature-index snow model

Other flows:
- bucket overflow: [Eq. 12a-g]      overflow from buckets with a maximum value, not a decision


Notes:
*    Not explicitly described in Clark et al. (2008). Look at QRAINERROR.F90 and UPDATE_SWE.F90 for rainfall error options. No routing is self-explanatory.
**   Not explicitly described in Clark et al. (2008). Possibly added later?
***  Not described in Clark et al. (2008) but introduced later in Henn et al. (2015).
```