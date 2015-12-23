
def load_observed_data():
    """
    Please fill here the observational data of your star of choice.
    Note that the variables that are used for the Likelihood calculation are defined in the function calc_Chi_and_Factor
    in bs_class.py. The number of variables defined and returned here has no effect on how the likelihood is calculated.
    """

    # log(g):
    obs_logg = 4.3			    # cgs
    obs_std_logg = 0.22		    # cgs

    # Teff:
    obs_teff = 6200			    # K
    obs_std_teff = 200		    # K

    # Z/Xs:
    obs_ZXs = 0.019
    obs_std_ZXs = 0.006

    # Min and max observed frequencies (to calculate the separations of the model within this range)
    f_min = 1.250               # mHz
    f_max = 2.275               # mHz

    # Large separation: (see Tassoul ApJSS 43 (1980), Lopes & Turck-Chieze A&A 290 (1994))
    obs_large_sep = 0.0878421	# mHz
    obs_std_large_sep = 0.0006	# mHz

    # Small separation:
    obs_small_sep = 0.0065780	# mHz
    obs_std_small_sep = 0.0006	# mHz

    # Ratios:
    #   the mean and slope are calculated for a specific frequency range:
    f_min_sl = 1.6			    # mHz
    f_max_sl = 1.9			    # mHz
    #   mean r02 and r010: (see Roxburgh & Vorontsov A&A 411 (2003), Cunha & Metcalfe ApJ 666 (2007), arXiv: 0705.3009)
    obs_mean_r02 = 0.0762410
    obs_std_mean_r02 = 0.003
    obs_mean_r010 = 0.0438890
    obs_std_mean_r010 = 0.0025
    #   slope of r02, r010 and Lsr010 (for Lsr010 see Brandao, Cunha & Christensen-Dalsgaard (2014), arXiv: 1311.7600)
    obs_sl_r02 = -0.0336458
    obs_std_sl_r02 = 0.005
    obs_sl_r010 = -0.0366931
    obs_std_sl_r010 = 0.0068747
    obs_sl_Lsr010 = -0.0031958
    obs_std_sl_Lsr010 = 0.0006037

    return obs_logg, obs_std_logg, obs_teff, obs_std_teff, obs_ZXs, obs_std_ZXs, \
           f_min, f_max, obs_large_sep, obs_std_large_sep, obs_small_sep, obs_std_small_sep, f_min_sl, f_max_sl
