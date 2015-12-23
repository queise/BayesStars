#!/usr/bin/python
from bs_class import *
import pymultinest

def main():
    """
    *************
     Bayes Stars
    *************

    Calculates posterior probabilities (and best fit) on a set of parameters comparing the results
    of stellar modelling (dmp2k + adipls) with stellar observations.

    The observational parameters to be compared with the model must be introduced in bs_obsdata.py, and some functions
    in bs_class.py must be modified accordingly (i.e. load_model_chars, calc_Chi_and_Factor, etc.).

    Requirements:
        - dmp2k version >= Nov14
            * recommended: use version >= Feb15 to use specific SRs for cluster runs (in src/for_cluster_runs/)
            * compile with silent=.TRUE. in src/dmpmain/mod_dmpmain.f
        - Adipls 2012
        - pymultinest
        - OpenMPI
        - (optional: SepCalc v2.2.6, not needed because now the seism parameters are calculated with python functions)

    Structure:
        - run_dir/ : to be created by user, is a copy of the sample_run_dir/ directory provided (name can be changed), with:
		       bs_main.py, bs_class.py, bs_func.py, bs_func2.py, dmp2k.out, sun.don
		       and, if running in cluster with cue system: condor.sub  openmpiscript
        - inputs/  : initially empty, put it in /local if running in cluster
        - outputs/ : initially empty, put it in /local if running in cluster
        - default_input_files/ m010.zams, adipls.c.in, agsm.l5bi.d.15.p2, redistrb.in

    Usage:
        - First, set:
            - A name for your run in the var run_name in bs_main.py.
            - The options for the stellar model run and the paths (NOM_OUTDIR and NOM_CHEMIN) in sun.don.
            - The options for the oscillations in adipls.c.in.
            - The priors for the sampled parameters in the function allpriors in bs_func.py.
            - The observational data in bs_obsdata.py
            - The paths for defaults/inputs/outputs in the variables self.path* in bs_class.py
            - If using condor cue system, set the paths in condor.sub and in openmpiscript
        - Run:
            - using only OpenMPI (no cueing system):
              nohup mpirun -np 32 ./bs_main.py > /dev/null &
            - same, but all output written to output.txt (the file usually gets extremely big, use only for debugging):
              nohup mpirun -np 32 ./bs_main.py > output.txt &
            - using Condor cueing system:
              condor_submit condor.sub

    """

    # Please give a name to your run, it will be used as a root name for multinest files:
    run_name = 'NoDM'

    # The "cube" array is that of the parameters to be sampled:
    cube = [ 0.9, 0.5, 0.1 ] # initial values are not used, only the number of params (dimension of cube) is relevant
    nparams = len(cube)

    # Multinest run:
    pymultinest.run(calc_Likelihood4stmd, calc_allpriors, nparams, importance_nested_sampling = False, resume = False,
                    verbose = True, n_live_points=32, outputfiles_basename=run_name+"_", sampling_efficiency=0.02,
                    const_efficiency_mode=True, init_MPI=False)


if __name__ == "__main__":
    main()
