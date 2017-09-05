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
   * `init\_read\_parameters` makes use of `ParIO.py`
   * `otherRef` calculates derived reference value


## Usage

TODO: reconstruct zgrid and find the value of z corresponding to the decay factor

TODO: time averaging, kx, kz averaging or selecting
