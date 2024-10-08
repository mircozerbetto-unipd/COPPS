#!/bin/bash
../fitFB2potential ./angles_trj.dat 100 100 3 3 -180.0 -180.0 0 298.15
# 100 100 : x-bins and y-bins.
# 3 3     : xNmax, yNmax expansion parameters.
# -180.0 -180.0 : initial values of x and y. Final values are assumed
#                 to be +360 deg with respect to the initial values.
