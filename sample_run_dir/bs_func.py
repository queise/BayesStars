from bs_class import *

def calc_Likelihood4stmd(cube, ndim, nparams):
    mass  = cube[0]
    metal = cube[1]
    age   = cube[2]

    # Initializes stellar model with given inputs:
    stmod = cl_stmod(mass, metal, age)
    print stmod.modelname

    # Prepares and runs DMP2k code for stellar evolution:
    stmod.run_st_evolution()
    if not stmod.OK:
        stmod.print_results(Errflag='DMP-Err')
        return -np.inf

    # Prepares and runs ADIPLS code for stellar oscillations:
    stmod.run_oscillations()
    if not stmod.OK:
        stmod.print_results(Errflag='OSC-Err')
        return -np.inf

    # Prepares and runs SepCalc small code that calculates seismic parameters and Chi2:
    stmod.calc_seismicparameters()

    # Calculates likelihood:
    stmod.calc_Likelihood()

    # Prints results in summary file:
    stmod.print_results()

    return stmod.Like


def calc_allpriors(cube, ndim, nparams):
    """
    Transforms normalized cube (0-1):
		cube[0] : mass
		cube[1] : metallicity
		cube[2] : final age
    to real units, taking into account a prior distribution defined here.

    When defining the prior functions, use:
        cube[n] = min + cube[n]*(max-min) for linear priors
        cube[n] = min + m.pow(cube[n]/((max-min)^(-1/2.35)),2.35) for a Salpeter initial mass function
    """
    cube[0] = 0.7 + m.pow(cube[0]/1.04585450769,2.35) 	# Salpeter,  range is approx. [0.6, 1.6] for m010.zams
    cube[1] = 0.0005 + cube[1]*(0.0400-0.0005)	     	# Linear, range allowed by opa tables is [0.0005, 0.0404]
    cube[2] = 0.0 + cube[2]*13800.			# Linear, range [0, 13800.]

    # intervals February 2015:
    # cube[0]       [0.7, 1.6]
    # cube[1]       [0.0005, 0.0400]
    # cube[2]       [0.0, 13800.]

    # intervals January 2015:
    # cube[0]	[1.0, 1.6]
    # cube[1]	[0.0095, 0.0334]
    # cube[2]	[0.0, 8000.]
