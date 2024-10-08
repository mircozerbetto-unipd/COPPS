#!/bin/bash
../fitFB1potential ./angle.trj 100 10 -180.0
# 100    : n-bins
# 10     : nMax expansion parameter
# -180.0 : initial value of x. Final value is assumed to be +360 deg with respect to the initial value.
