import netCDF4
import numpy

variables1 = ['ssu_m','ssv_m','sst_m','sss_m','ssh_m','frq_m','rnf_b','rnf_hc_b']
variables2 = ['rnf_sc_b','utau_b','vtau_b','qns_b','emp_b','sfx_b','en','avt']
variables3 = ['avm','avmu','avmv','dissl','sbc_hc_b','sbc_sc_b']
variables4 = ['qsr_hc_b','fraqsr_1lev','gcx','gcxb','ub','vb','tb','sb']
variables5 = ['rotb','hdivb','sshb','un','vn','tn','sn','rotn','hdivn','sshn','rhop']

variables = variables1 + variables2 + variables3 + variables4 + variables5

rest_pcg = "RESTART_PCG/O2L3_LONG_00001000_restart.nc"
rest_pet = "RESTART_PETSC/O2L3_LONG_00001000_restart.nc"

f1 = netCDF4.Dataset(rest_pcg,'r')
f2 = netCDF4.Dataset(rest_pet,'r')

for i in variables:
	print i
	var_pcg = f1.variables[i][0,:]
	var_pet = f2.variables[i][0,:]
	print i+' shape',var_pcg.shape
	rmse = numpy.sqrt(((var_pet - var_pcg) ** 2).mean())
	print 'rmse '+i+':', rmse

f1.close()
f2.close()

