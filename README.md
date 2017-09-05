# postprocessing_tools

## Description

This directory contains Python code for postprocessing results from the [GENE](http://genecode.org) gyrokinetic code.
The code provides the following tools to read its raw output into dictionary:
* `ParIO.py` for reading `parameters` files
* `fieldlib.py` for reading `field` files which contain phi, apar, bpar
* `momlib.py` for reading `mom` files with dens, temp, upar, etc
* `read_write_geometry.py` for reading `tracer_efit` files with geometric quantities

Using these code to convert raw data into dictionaries, I write more code for specific research interest. Some of the often used code are written as methods in `__Helper.py` code:
* `parIOHelper.py' has two methods: 
   * `init_read_parameters` makes use of `ParIO.py`
   * `otherRef` calculates derived reference value
* `momHelper.py` has seven methods:
   * `momen_step_time` plots step time for entire simulation
   * `global_moments` returns time, deln, tperp for given time index and space indices
   * `momen_xz` returns time, dens_xz, tperp_xz for given time which adds all the ky in the simulation at zeta = 0
   * `momen_tx` returns time, deln_tx, tperp_tx for given ky and z indices
   * `momen_rms` returns the RMS given a time series signal
   * `momen_tky` returns time, deln_tky, tperp_tky for given x and z indices
   * `radiometer` returns the integral in given frequency range given signal in the frequency domain
## Usage

TODO: reconstruct zgrid and find the value of z corresponding to the decay factor

TODO: time averaging, kx, kz averaging or selecting
