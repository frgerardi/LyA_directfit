# LyA_directfit

This code is mainly based on the use of [Vega](https://github.com/andreicuceu/vega), previously developed by A.Cuceu, and [Cobaya](https://github.com/CobayaSampler/cobaya). 

Given a data vector, the whole workflow runs an ad-hoc likelihood function (likelihood.py), which is called by Cobaya, via the info file, where all the parameters to sample are also specified. Inside this likelihood, for a given step of the sampling process, cosmological quantities of interest are computated, via CAMB, and the corresponding correlation functions are calculated using a slightly different version of Vega, where all amendments are included within the [vega/vega_interface_mod.py](https://github.com/frgerardi/LyA_directfit/blob/main/vega/vega_interface_mod.py) module.

To run the whole code, after specifying the correlations you are interested in via the [config_files/main_auto.ini](https://github.com/frgerardi/LyA_directfit/blob/main/config_files/main_auto.ini), it is simply necessary to run the following command line in your conda environment

<cobaya-run info.yaml>

