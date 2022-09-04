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
    
    debugging: Optional[bool]
    correlation_type: Optional[str]
    effective_redshift: Optional[float]
    D_M_fid: Optional[float]
    D_H_fid: Optional[float]
    derived_beta_tracer2: Optional[bool]
    
    def initialize(self, **params_values):
        
        if self.debugging:
            print('effective redshift = ', self.effective_redshift)
        
        if self.correlation_type=='AUTO':
            self.vega = VegaInterface('config_files/main_auto.ini')
        elif self.correlation_type=='AUTO+CROSS':
            if self.derived_beta_tracer2 == True:
                self.vega = VegaInterface('config_files/main_auto+cross_bifQSO.ini')
            else:
                self.vega = VegaInterface('config_files/main_auto+cross.ini')
        else:
            print('----------------------------CHANGE INI SETTINGS-------------------------------------------')
        self.k_grid = np.logspace(-4,2,700)
        self.vega.fiducial['z_fiducial'] = self.effective_redshift
        
        self.k_max = 1


    def get_requirements(self):
        
        reqs = {'H0': None, 'ombh2': None , 'omch2': None, 'omnuh2': None, 'omk': None,'As': None, 'ns': None,\
                'Pk_grid': {'z': [self.effective_redshift], 'k_max':self.k_max, 'nonlinear':[False]},
                'angular_diameter_distance':{'z': [self.effective_redshift]}, 'Hubble':{'z': [self.effective_redshift]}, 'omegal': None,
                'sigma8_z':{'z': [self.effective_redshift]}, 'fsigma8':{'z': [self.effective_redshift]}, 'rdrag': None,  'omegam': None, 
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
        
        scale_factor = 1/(1+self.effective_redshift)
        h = params_values['H0']/100
        D_M = (self.provider.get_angular_diameter_distance(self.effective_redshift)/scale_factor)*h
        D_H = ((scipy.constants.c/1000)/self.provider.get_Hubble(self.effective_redshift))*h
        params_values['ap_full'] = D_H[0]/self.D_H_fid
        params_values['at_full'] = D_M[0] /self.D_M_fid
       
        params_values['bias_LYA'] = -10**(params_values['log_biasLYA'])
        params_values['bias_QSO'] = 10**(params_values['log_biasQSO'])
        
              
        growth_rate = (self.provider.get_fsigma8(self.effective_redshift))/(self.provider.get_sigma8_z(self.effective_redshift))
        
        params_values["growth_rate"] = growth_rate[0]

        params_values["_derived"]["growth_rate"] = growth_rate[0]
        params_values["_derived"]["f_sigma8"] = self.provider.get_fsigma8(self.effective_redshift)[0]

             
        params_values["debugging"] = self.debugging
        
        if self.debugging:
            tracers_dict = {'bias_LYA': params_values['bias_LYA']}
            tracers_dict['beta_LYA'] = params_values['beta_LYA']
            if self.correlation_type == 'AUTO+CROSS':
                tracers_dict['bias_QSO'] = params_values['bias_QSO']
                if self.derived_beta_tracer2 == False:
                    tracers_dict['beta_QSO'] = params_values['beta_QSO']
                else:
                    tracers_dict['beta_QSO'] =  (growth_rate/params_values["bias_QSO"])[0]
            print('########### in likelihood properties are ', tracers_dict)

        

        
        k_Mpc, z, pk_Mpc = self.provider.get_Pk_grid(nonlinear = False)
        k = k_Mpc/h
        pk = pk_Mpc[0]*(h**3)
#        np.save('test_pk.npy',np.column_stack((k,pk)))        
        
        cs = interp1d(np.log(k), np.log(pk), kind='cubic')
        fit = np.polyfit(np.log(k[-10:]), np.log(pk[-10:]),deg=1)
        cs_linear = np.poly1d(fit)
        idx = np.where((self.k_grid<=self.k_max))
        pk_hMpc = np.zeros(self.k_grid.shape)
        pk_hMpc[idx] = np.exp(cs(np.log(self.k_grid[idx])))
        pk_hMpc[idx[0][-1]:] = np.exp(cs_linear(np.log(self.k_grid[idx[0][-1]:])))
               
        
#        plot_isotropic_pk = 'isotropic_pk.npy'
#        if plot_isotropic_pk is not None:
#            print('Saving CAMB power spectra')
#            self.vega.save_sample(sample = {'H0':[params_values['H0']],\
#                                            'omch2':[round(params_values['omch2'], 4)],\
#                                            'ombh2':[round(params_values['ombh2'], 4)],\
#                                            'k_hMpc': [self.k_grid.tolist()], 'pk_hMpc': [pk_hMpc.tolist()]} , 
#                                  filename=plot_isotropic_pk)
        
        if 'pk_full' not in self.vega.fiducial.keys():
            self.vega.fiducial['k'] = self.k_grid
            self.vega.fiducial['pk_full'] = pk_hMpc
        self.pk_full = pk_hMpc
            
        self.vega.fiducial['Omega_m'] = self.provider.get_param("omegam")
        self.vega.fiducial['Omega_de'] = self.provider.get_param("omegal")
        
        chi2 = self.vega.chi2(params_values, direct_pk = self.pk_full)
        return -chi2 / 2
