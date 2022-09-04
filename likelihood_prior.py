from cobaya.likelihood import Likelihood
import numpy as np
import pandas as pd
import os
import fnmatch
import matplotlib.pyplot as plt
from vega.vega_interface_mod import VegaInterface
import scipy
from scipy.interpolate import interp1d,InterpolatedUnivariateSpline
from typing import Optional



class Likelihood(Likelihood):
    
    effective_redshift: Optional[float]
    correlation_type: Optional[str]
    derived_beta_tracer2: Optional[bool] 

    def initialize(self, **params_values):
        
        print('Running prior likelihood')


    def get_requirements(self):
        
        reqs = {'H0': None, 'ombh2': None , 'omch2': None, 'omnuh2': None, 'omk': None,'As': None, 'ns': None,\
                'Hubble':{'z': [self.effective_redshift]}, 'omegal': None,'fsigma8':{'z': [self.effective_redshift]}, 
                'rdrag': None,  'omegam': None,'sigma8_z':{'z': [self.effective_redshift]}
                }
       
        reqs['log_biasLYA']=None
        #reqs['bias_LYA']=None  
        reqs['log_biasQSO']=None
        #reqs['bias_QSO']=None
        
        reqs['beta_LYA']=None

        reqs['sigma_velo_disp_lorentz_QSO']=None
            
        if self.correlation_type=='AUTO+CROSS':
            #reqs['bias_QSO']=None
            if self.derived_beta_tracer2 == False:
                reqs['beta_QSO']=None
        
        return reqs

    def logp(self, **params_values):

        growth_rate = (self.provider.get_fsigma8(self.effective_redshift))/(self.provider.get_sigma8_z(self.effective_redshift))

        params_values["growth_rate"] = growth_rate[0]

        params_values["_derived"]["growth_rate"] = growth_rate[0]
        params_values["_derived"]["f_sigma8"] = self.provider.get_fsigma8(self.effective_redshift)[0]

        return 0
