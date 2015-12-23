#!/usr/bin/python
from bs_func import *

def main():
    """
    This script runs a SINGLE stellar model + Likelihood calculation for testing purposes. It does not use multinest.
    """

    ################################
    # TO TEST THE PRIOR FUNCTIONS:
    ################################
    cube = [ 0.5, 0.2, 0.005 ]
    ndim = len(cube)
    nparams = len(cube)
    # calc_allpriors transforms normalized cube to real units inputs with prior distribution:
    calc_allpriors( cube, ndim, nparams )

    ########################################
    # TO TEST THE MODELLING AND LIKELIHOOD:
    ########################################
    # You can define here the input parameters of the model (real units, not normalized):
    # (or comment the next line and to just use the cube calculated by calc_allpriors)
    cube = [ 1.003e+00, 2.407e-02, 3.330e+01 ]
    # calculates Likelihood:
    print calc_Likelihood4stmd( cube, ndim, nparams )

if __name__ == "__main__":
    main()
