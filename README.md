# RRTMG_LW: Longwave Radiative Transfer Model for GCMs
---
**Contents**

1. [Introduction](#intro)
2. [Releases](#releases)
3. [Column Version](#column)
4. [GCM Version](#gcm)
5. [Contact](#contact)
6. [References](#ref)

## Introduction <a name="intro"></a>

This package contains the source code and sample makefiles necessary to run the latest version of RRTMG\_LW, a correlated *k*-distribution longwave radiative transfer model developed at AER for application to GCMs. This version of RRTMG\_LW has been modified from the standard RRTM\_LW distributed by AER to enhance its performance for use within general circulation models. This code has also been modified to utilize updated FORTRAN coding features. Two modes of operation are possible: 

1) RRTMG_LW can be run as a [column model](https://github.com/AER-RC/RRTMG_LW#rrtmg_lw--column-version) using the [input files](https://github.com/AER-RC/RRTMG_LW#input-data) and [source modules](https://github.com/AER-RC/RRTMG_LW#source-code) described in the column version section, or 
2) it can be implemented as a [subroutine into an atmospheric general circulation model or single column model](https://github.com/AER-RC/RRTMG_LW#rrtmg_lw--gcm-version). 

The version of RRTMG\_LW provided here has been modified from the standard RRTM\_LW to enhance performance with little effect on the accuracy. The total number of *g*-points used has been reduced from 256 to 140. Fluxes are accurate to within 0.5 W m<sup>-2</sup> and cooling rate within 0.1 K day<sup>-1</sup> relative to the standard RRTM\_LW, which is itself accurate to within 1 W m<sup>-2</sup> of the data-validated line-by-line radiative transfer model, LBLRTM. Required absorption coefficient input data can be read in either from data stored within the code or from an external netCDF file as selected in the makefile. 

This model can also utilize McICA, the Monte-Carlo Independent Column Approximation, to represent sub-grid scale cloud variability such as cloud fraction and cloud overlap. If the McICA option is selected to model a cloudy profile in column mode, then the model will run stochastically, and the output fluxes and heating rates will be an average over 200 samples. In GCM mode, the code will calcualte a single column per profile, and the statistical basis is provided by the spatial and temporal dimensions of the 3-D calculations. Several cloud overlap methods are available for partial cloudiness including maximum-random, exponential, and exponential-random. 

The model includes an optional feature to provide simultaneously with a normal forward calculation the change in upward flux with respect to surface temperature for each model level. This option is controlled by the input flag, `idrv`. Setting this flag to 1 will output dF/dT for total sky and clear sky in GCM mode in new output arrays `duflx_dt` and `duflxc_dt`. These can be utilized to approximate the change in upward flux for a change in surface temperature only at time intervals between full radiation calls. In single column mode, setting idrv to 1 requires the extra input of a dT change in surface temperature relative to the input surface temperature, and the provided dT will be applied to the flux derivative to output a modified upward flux profile for that dT change in surface temperature. The default `idrv` setting of 0 provides the original forward radiative transfer calculation.

For more information on the model, see the [Wiki Description Page](https://github.com/AER-RC/RRTMG_LW/wiki/Description)

## Releases <a name="releases"></a>

[Version 5.0 is the latest version of the model](https://github.com/AER-RC/RRTMG_LW/releases/tag/v5.0)

Releases before Version 5.0 are not publicly available.

## RRTMG_LW : Column Version <a name="column"></a>

### DOCUMENTATION
The following text files (in the `doc` directory), along with this `README` provide information on release updates and on using and running RRTMG\_LW:

| Filename | Description |
| :--- | :--- |
| `release_notes.txt` | Code archive update information |
| `rrtmg_lw_instructions.txt` | Input instructions for files INPUT_RRTM, IN_CLD_RRTM and IN_AER_RRTM |

### SOURCE CODE
The following source files (in the `src` directory) must be used to run RRTMG\_LW in stand-alone mode as a column model (the utility files are stored separately in the `aer_rt_utils` directory):

| Filename | Description |
| :--- | :--- |
| `rrtmg_lw.1col.f90` | RRTMG_LW main module |
| `rrtmg_lw_cldprop.f90` | Calculation of cloud optical properties |
| `rrtmg_lw_cldprmc.f90` | Calculation of cloud optical properties (McICA) |
| `rrtmg_lw_init.f90` | RRTMG_LW initialization routine; reduces *g*-intervals from 256 to 140 |
| `rrtmg_lw_k_g.f90` | Absorption coefficient data file |
| `rrtmg_lw_read_nc.f90` | Optional absorption coefficient data netCDF input |
| `rrtmg_lw_rtrn.f90` | Calculation of clear and cloudy radiative transfer using random cloud overlap |
| `rrtmg_lw_rtrnmr.f90` | Calculation of clear and cloudy radiative transfer using maximum-random cloud overlap |
| `rrtmg_lw_rtrnmc.f90` | Calculation of clear and cloudy radiative transfer using McICA (with user-selectable overlap method) |
| `rrtmg_lw_setcoef.f90` | Set up routine |
| `rrtmg_lw_taumol.f90` | Calculation of optical depths and Planck fractions for each spectral band |
| `mcica_random_numbers.f90` | Random number generator for McICA |
| `mcica_subcol_gen_lw.1col.f90` | Sub-column generator for McICA |
| `rrtatm.f` | Process user-defined input data files |
| `extra.f` | Process input data files |
| `util_**.f` | Utilities (available for multiple platforms) |

The following module files (in the `modules` directory) must be used to run RRTMG\_LW in stand-alone mode as a column model (these must be compiled before the source code files):

| Filename | Description |
| :--- | :--- |
| `parkind.f90` | real and integer kind type parameters |
| `parrrtm.f90` | main configuration parameters |
| `rrlw_cld.f90` | cloud property coefficients |
| `rrlw_con.f90` | constants |
| `rrlw_kg**.f90` | absorption coefficient arrays for 16 spectral bands |
| `rrlw_ncpar.f90` | parameters for netCDF input data option |
| `rrlw_ref.f90` | reference atmosphere data arrays |
| `rrlw_tbl.f90` | exponential look up table arrays |
| `rrlw_vsn.f90` | version number information |
| `rrlw_wvn.f90` | spectral band and *g*-interval array information |

### INPUT DATA
The following file (in the `data` directory) is the optional netCDF input file containing absorption coefficient and other input data for the model. The file is used if keyword `KGSRC` is set for netCDF input in the makefile. 

| Filename | Description |
| :--- | :--- |
| `rrtmg_lw.nc` | Optional netCDF input data file |

### MAKEFILES
The following files (in the `build/makefiles` directory) can be used to compile RRTMG\_LW in stand-alone mode as a column model on various platforms.  Link one of these into the `build` directory to compile. 

| Filename | Description |
| :--- | :--- |
| `make_rrtmg_lw_sgi` | Sample makefile for SGI |
| `make_rrtmg_lw_sun` | Sample makefile for SUN |
| `make_rrtmg_lw_linux_pgi` | Sample makefile for LINUX (PGI compiler) |
| `make_rrtmg_lw_aix_xlf90` | Sample makefile for AIX (XLF90 compiler) |
| `make_rrtmg_lw_OS_X_g95` | Sample makefile for OS_X (G95 compiler) |
| `make_rrtmg_lw_OS_X_ibm_xl` | Sample makefile for OS_X (IBM XL compiler) |

### SAMPLE INPUT/OUTPUT
Several sample input and output files are included in the `run_examples_std_atm` directory. Note that user-defined profiles may be used for as many as 200 layers.

| Filename | Description |
| :--- | :--- |
| `INPUT_RRTM` | Required input file for (clear sky) atmospheric specification |
| `IN_CLD_RRTM` | Required input file for cloud specification if clouds are present |
| `IN_AER_RRTM` | Required input file for aerosol specification if aerosols are present |
| `OUTPUT_RRTM` | Main output file for atmospheric fluxes and heating rates |
| `input_rrtm_ICRCCM_sonde` | Sample radiosonde-style input profile for clear sky |
| `input_rrtm.MLS-clr` | Sample 51 layer mid-latitude summer standard atmosphere for clear sky |
| `input_rrtm.MLS-cld-imca0-icld2` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and maximum-random cloud overlap selected (without McICA) |
| `input_rrtm.MLS-cld-imca1-icld2` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and maximum-random cloud overlap selected (with McICA) |
| `input_rrtm.MLS-cld-imca1-icld4-idcor0` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and exponential cloud overlap and constant decorrelation length selected (with McICA) |
| `input_rrtm.MLS-cld-imca1-icld5-idcor0` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and exponential-random cloud overlap and constant decorrelation length selected (with McICA) |
| `input_rrtm.MLS-cld-imca1-icld5-idcor1` | Sample 51 layer mid-latitude summer standard atmosphere with cloud flag turned on and exponential-random cloud overlap and varying decorrelation length selected (with McICA) |
| `input_rrtm.MLS-clr-aer12` | Sample 51 layer mid-latitude summer standard atmosphere with aersol flag set |
| `input_rrtm.MLS-clr-xsec` | Sample 51 layer mid-latitude summer standard atmosphere with cross-section input (CFCs, etc.) |
| `input_rrtm.MLS-clr-idrv1` | Sample 51 layer mid-latitude summer standard atmosphere with derivative option set to provide modified upward fluxes for the provided change in surface temperature |
| `input_rrtm.MLW-clr` | Sample 51 layer mid-latitude winter standard atmosphere |
| `input_rrtm.SAW-clr` | Sample 51 layer sub-arctic winter standard atmosphere |
| `input_rrtm.TROP-clr` | Sample 51 layer tropical standard atmosphere |
| `in_cld_rrtm-cld5` | Sample cloud input file |
| `in_cld_rrtm-cld7` | Sample cloud input file |
| `in_aer_rrtm-aer12` | Sample aerosol input file |
| `script.run_std_atm | UNIX script for running the full suite of example cases, which will put the output into similarly named files prefixed with `output_rrtm*` |

### INSTRUCTIONS FOR COMPILING AND RUNNING THE COLUMN MODEL
1) In the `build` directory, link one of the makefiles from the `makefile` sub-directory into `build/make.build`. To use the optional netCDF input file, switch the keyword `KGSRC` in the makefile from `dat` to `nc`. Compile the model with `make -f make.build`
2) Link the executable to, for example, `rrtmg_lw` in the `run_examples_std_atm` directory
3) If the optional netCDF input file was selected when compiling, link the file `data/rrtmg_lw.nc` into the `run_examples_std_atm` directory.  
4) In the `run_examples_std_atm` directory, run the UNIX script `./script.run_std_atm` to run the full suite of example cases. To run a single case, modify `INPUT_RRTM` following the instructions in `doc/rrtmg_lw_instructions.txt`, or copy one of the example `input_rrtm*` files into `INPUT_RRTM`. If clouds are selected (`ICLD` > 0), then modify `IN_CLD_RRTM` or copy one of the `in_cld_rrtm*` files into `IN_CLD_RRTM`. If aerosols are selected (`IAER` > 0), then modify `IN_AER_RRTM` or set it to the sample file `in_aer_rrtm-aer12`.
5) In column mode, if McICA is selected (`IMCA`=1) with partial cloudiness defined, then RRTMG\_LW will run the case 200 times to derive adequate statistics, and the average of the 200 samples will be written to the output file, `OUTPUT_RRTM`. 

## RRTMG_LW : GCM version <a name="gcm"></a>

### SOURCE CODE:
The following source files (in the `src` directory) must be used to run RRTMG\_LW as a callable subroutine:

| Filename | Description |
| :--- | :--- |
| `rrtmg_lw_rad.f90` | RRTMG_LW main module (with McICA) |
| `rrtmg_lw_rad.nomcica.f90` | Optional RRTMG_LW main module (without McICA only) |
| `rrtmg_lw_cldprop.f90` | Calculation of cloud optical properties |
| `rrtmg_lw_cldprmc.f90` | Calculation of cloud optical properties (McICA) |
| `rrtmg_lw_init.f90` | RRTMG_LW initialization routine; reduces *g*-intervals from 256 to 140; (This has to run only once and should  be installed in the GCM initialization section) |
| `rrtmg_lw_k_g.f90` | Absorption coefficient data file |
| `rrtmg_lw_read_nc.f90` | Optional absorption coefficient data netCDF input |
| `rrtmg_lw_rtrn.f90` | Calculation of clear and cloudy radiative transfer using random cloud overlap |
| `rrtmg_lw_rtrnmr.f90` | Calculation of clear and cloudy radiative transfer using maximum-random cloud overlap |
| `rrtmg_lw_rtrnmc.f90` | Calculation of clear and cloudy radiative transfer using McICA (with selectable overlap method) |
| `rrtmg_lw_setcoef.f90` | Set up routine |
| `rrtmg_lw_taumol.f90` | Calculation of optical depths and Planck fractions for each spectral band |
| `mcica_random_numbers.f90` | Random number generator for McICA |
| `mcica_subcol_gen_lw.f90` | Sub-column generator for McICA (must be called in GCM just before call to RRTMG) |

**NOTE**: Only one of `rrtmg_lw_k_g.f90` or `rrtmg_lw_read_nc.f90` is required. 

The following module files (in the `modules` directory) must be used to run RRTMG\_LW as a callable subroutine (these must be compiled before the source code)

| Filename | Description |
| :--- | :--- |
| `parkind.f90` | real and integer kind type parameters |
| `parrrtm.f90` | main configuration parameters |
| `rrlw_cld.f90` | cloud property coefficients |
| `rrlw_con.f90` | constants |
| `rrlw_kg**.f90` | absorption coefficient arrays for 16 spectral bands |
| `rrlw_ncpar.f90` | parameters for netCDF input data option |
| `rrlw_ref.f90` | reference atmosphere data arrays |
| `rrlw_tbl.f90` | look up table arrays |
| `rrlw_vsn.f90` | version number information |
| `rrlw_wvn.f90` | spectral band and *g*-interval array information |

### INPUT DATA:
The following file (in the `data` directory) is the optional netCDF file containing absorption coefficient and other input data for the model. The file is used if source file `rrtmg_lw_read_nc.f90` is used in place of `rrtmg_lw_k_g.f90` (only one or the other is required). 

| Filename | Description |
| :--- | :--- |
| `rrtmg_lw.nc` | Optional netCDF input data file |

### NOTES ON RUNNING THE GCM (SUBROUTINE) VERSION OF THE CODE:

1) The module `rrtmg_lw_init.f90` is the initialization routine that has to be called only once.  The call to this subroutine should be moved to the initialization section of the host model if RRTMG\_LW is called by a GCM or SCM. 
2) The number of model layers and the number of columns to be looped over should be passed into RRTMG\_LW through the subroutine call along with the other model profile arrays.  
3) To utilize McICA, the sub-column generator (`mcica_subcol_gen_lw.f90`) must be implemented in the GCM so that it is called just before RRTMG\_LW. The cloud overlap method is selected using the input flag, icld. If either exponential (`ICLD`=4) or exponential-random (`ICLD`=5) cloud overlap is selected, then the subroutine `get_alpha` must be called prior to calling `mcica_subcol_lw` to define the vertical correlation parameter, `alpha`, needed for those overlap methods. Also for those methods, use the input flag `idcor` to select the use of either a constant or latitude-varying decorrelation length. If McICA is utilized, this will run only a single statistical sample per model grid box. There are two options for the random number generator used with McICA, which is selected with the variable irnd in `mcica_subcol_gen_lw.f90`. When using McICA, then the main module is `rrtmg_lw_rad.f90`. If McICA is not used, then the main module is `rrtmg_lw_rad.nomcica.f90` and the cloud overlap method is selected by setting flag `icld`. 

## Maintenance and Contact Info <a name="contact"></a>
Atmospheric and Environmental Research, 
131 Hartwell Avenue, Lexington, MA 02421

Original version:   E. J. Mlawer, et al. (AER)
Revision for GCMs:  Michael J. Iacono (AER)

Contact:   Michael J. Iacono   (E-mail: miacono@aer.com)

## References <a name="ref"></a>

* [AER Radiative Transfer Models Documentation](https://www.rtweb.aer.com)
* [Github Wiki](https://github.com/AER-RC/RRTMG_LW/wiki)
* **RRTMG_LW, RRTM_LW**
    * Iacono, M.J., J.S. Delamere, E.J. Mlawer, M.W. Shephard, S.A. Clough, and W.D. Collins, Radiative forcing by long-lived greenhouse gases: Calculations with the AER radiative transfer models, *J. Geophys. Res.*, 113, D13103, doi:10.1029/2008JD009944, 2008.
    * Clough, S.A., M.W. Shephard, E.J. Mlawer, J.S. Delamere, M.J. Iacono, K. Cady-Pereira, S. Boukabara, and P.D. Brown, Atmospheric radiative transfer modeling: a summary of the AER codes, *J. Quant., Spectrosc. Radiat. Transfer*, 91, 233-244, 2005.
    * Iacono, M.J., J.S. Delamere, E.J. Mlawer, and S.A. Clough, Evaluation of upper tropospheric water vapor in the NCAR Community Climate Model (CCM3) using modeled and observed HIRS radiances. *J. Geophys. Res.*, 108(D2), 4037, doi:10.1029/2002JD002539, 2003.
    * Iacono, M.J., E.J. Mlawer, S.A. Clough, and J.-J. Morcrette, Impact of an improved longwave radiation model, RRTM, on the energy budget and thermodynamic properties of the NCAR Community Climate Model, CCM3, *J. Geophys. Res.*, 105, 14873-14890, 2000.
    * Mlawer, E.J., S.J. Taubman, P.D. Brown, M.J. Iacono, and S.A. Clough:  Radiative transfer for inhomogeneous atmospheres: RRTM, a validated correlated-k model for the longwave.  *J. Geophys. Res.*, 102, 16663-16682, 1997.
* **McICA**
    * Pincus, R., H. W. Barker, and J.-J. Morcrette, A fast, flexible, approximation technique for computing radiative transfer in inhomogeneous cloud fields, *J. Geophys. Res.*, 108(D13), 4376, doi:10.1029/2002JD003322, 2003.
*  **Latitude-Varying Decorrelation Length**
    *  Oreopoulos, L., D. Lee, Y.C. Sud, and M.J. Suarez, Radiative impacts of cloud heterogeneity and overlap in an atmospheric General Circulation Model, *Atmos. Chem. Phys.*, 12, 9097-9111, doi:10.5194/acp-12-9097-2012, 2012.
* [Full list of references](https://github.com/AER-RC/RRTMG_LW/wiki/References)
