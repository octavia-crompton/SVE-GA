Things to fix:
ncol and nrow
j and k

# Things to talk about

- presenting fluxes in a non-confusing way
- what to do with GW.
- library of tests.
- making the animations functionality usable.
  - see simulations. 


### List of notebook files:

- `Check_fluxes.ipynb` :
  - Check mass balance in the SVE model


- `GW_GA_validation.ipynb`
  - Checks SVE-GA model simulations against GW solutions
  - links to: `sve_GA_validate_sims/GW_validation_sims`
 

- `image_read.ipynb`
	- code to (i) read an image file, and (ii) regrid if necessary
    - consolidate with other image handling code
 
- `SVE_sims-GA_infl_test.ipynb`
	- checks that the SVE model correctly implements Green-Ampt
	- links to: `GA_infl_test`
  - TODO: check larger range of parameters


# TODO:  
## rename variables:

"dt_p" : "dt_print"
"rainD" : "rain_depth_cm"
"fluxD" : "flux_depth_cm"
"zinflc" : "infl_cm"
"influx" : ""
"flux1" ... : "flux1_m3" ()

## add new inputs:
"min_vol" : min volume at which to trunctate the simulation?


##  add new summary variables:
"fluxD" -->  "flux1_depth_cm", ...  ?
