# Code to invert surface observables for best-fit viscosity profile
import sys
import lmfit
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from os.path import join
from math import erf

def read_observed_fields():
  real_DT=np.loadtxt('REFERENCE/DynTop')
  real_geoid=np.loadtxt('REFERENCE/Ge-oid')
  #real_grav=np.loadtxt('REFERENCE/')
  real_CMB=np.loadtxt('REFERENCE/CMBTOP')
  return real_DT, real_geoid, real_CMB

def read_model_fields():
#    model_DT=np.loadtxt('DT_model.llz')
#    model_geoid=np.loadtxt('geoid_model.llz')
#    model_grav=np.loadtxt('grav_model.llz')
#    model_CMB=np.loadtxt('CMB_model.llz')
#    model_vel=np.loadtxt('vel_model.llz')
    model_DT=np.loadtxt('OUTPUT/DynTop')
    model_geoid=np.loadtxt('OUTPUT/Ge-oid')
#    model_grav=np.loadtxt('grav_model.llz')
    model_CMB=np.loadtxt('OUTPUT/CMBTOP')
#    model_vel=np.loadtxt('OUTPUT/VelTOP')
    return model_DT, model_geoid, model_CMB

#def callbackf(pin):
  #global Nfeval
  #print Nfeval
  #Nfeval += 1
   
def objective_func(pin):
  global Nfeval
  
  v = pin.valuesdict()
  moth_ipath = './INPUT/THOt.00_00_edit'
  moth_opath = './OUTPUT/'
  moth_rpath = './REFERENCE/'
  vis_all= './INPUT/visc_edit_'+repr(Nfeval)+'_iter'

  # top 5 slices (~100 km) = litho
  visc=np.zeros(nv)
  visc[0:33]=v['a']
  visc[33:65]=v['um']
  visc[65:130]=v['lm']

  np.savetxt(vis_all,visc)

  proc=sp.Popen('/space2/fdr22/Desktop/Viscosity_Inversion/KERNELS '+moth_ipath+' '+moth_opath+\
    ' '+vis_all+' '+moth_rpath,   shell=True, stdout=sp.PIPE)
  while True:
    line = proc.stdout.readline()
    if line != '':
      #the real code does filtering here
      misfit=line.rstrip()
    else:
      break
  Nfeval += 1
  print Nfeval-1
  return misfit

def final_misfit(oDT,oG,oCMB):
 
  moth_ipath = './INPUT/THOt.00_00_edit'
  moth_opath = './OUTPUT/'
  vis_all= './INPUT/visc_edit_'+repr(Nfeval-1)+'_iter'
  # initiate viscosity
  # top 5 slices (~100 km) = litho
  visc=np.zeros(nv)
  visc[0:33]=a_inv
  visc[33:65]=um_inv
  visc[65:130]=lm_inv

  sp.call('/space2/fdr22/Desktop/Viscosity_Inversion/KERNELS_OUT '+moth_ipath+' '+moth_opath+\
    ' '+vis_all,   shell=True)

  mDT,mg,mCMB=read_model_fields()
  
  misfit_DT=np.sqrt((1./nx)*np.sum(np.power(oDT[:,2]-mDT[:,2],2)))
  misfit_geoid=np.sqrt((1./nx)*np.sum(np.power(og[:,2]-mg[:,2],2)))
  misfit_grav=0.
  misfit_CMB=np.sqrt((1./nx)*np.sum(np.power(oCMB[:,2]-mCMB[:,2],2)))
  misfit_vel=0.
  umisfit=misfit_DT+misfit_geoid+misfit_grav+misfit_CMB+misfit_vel
  return umisfit
  
#number of viscosity values
nv=130
#depths slices ~22.5 km
Nfeval = 1

#number of latlon values
oDT,og,oCMB=read_observed_fields()
nx=len(oDT)

#weightings
w1=1.
w2=1.
w3=1.
w4=1.
w5=1.

#initial guesses
# lithosphere, 2 upper mantle layers and 3 lower mantle layers
a=1000.
um=1000.
lm=1000.

p1 = lmfit.Parameters()
p1.add('a',value=um,min=1.,max=1000000.)
p1.add('um',value=um,min=1.,max=1000000.)
p1.add('lm',value=lm,min=1.,max=1000000.)

mi_pwl=lmfit.minimize(objective_func, p1, args=(oDT,og,oCMB,), method='powell')

mean_res_pwl=objective_func(mi_pwl.params,oDT,og,oCMB)

a_inv=mi_pwl.params.valuesdict().values()[0]
um_inv=mi_pwl.params.valuesdict().values()[1]
lm_inv=mi_pwl.params.valuesdict().values()[2]

# read in files
oDT,og,oCMB=read_observed_fields()
um=final_misfit(oDT,og,oCMB)

f=open('visc_struct_3l_'+repr(w1)+'_'+repr(w2)+'_'+repr(w3)+'_'+repr(w4)+'_'+repr(w5)+'.txt',"w")
f.write("%s %s %s %s %s\n" % ("a", "um","lm", "misfit", "unweighted misfit"))
f.write("%.4g %.4g %.4g %.4g %.4g\n" % (a_inv, um_inv, lm_inv, mean_res_pwl,um))


