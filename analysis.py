import sys
import numpy as np
import pathlib
import h5py
import matplotlib.pyplot as plt
import scipy.optimize
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
import csv
import os
from mpi4py import MPI

basedir = os.path.dirname(os.path.realpath(__file__))

dfile = basedir+'/Runfiles/fields_s1_lambda_rho=0.4_new.h5'
#dfile1 = basedir+'/mergedpot_lambda_rho=0.4.h5'
#dfile2 = basedir+'/fields_s1_1_lambda_rho=1.h5'
data = h5py.File(str(dfile), "r")
#data1 = h5py.File(str(dfile1), "r")
#data2 = h5py.File(str(dfile2), "r")

rho = data['tasks/rho'][0, 14, 12, :40]
#rho1 = data1['tasks/rho'][5, 14, 12, :]
#rho2 = data2['tasks/rho'][5, 14, 12, :]

z = data['scales/z/1.0'][:40]
print(z)
print(rho)
plt.plot(z, rho, 'r.', label = "lambda_rho = 0.4")
#plt.semilogy(z, rho1, 'g.', label = "lambda_rho = 0.4")
#plt.semilogy(z, rho2, 'b.', label = "lambda_rho = 1")

plt.title("initial density distribution")
plt.legend(loc='lower left')
plt.xlabel("z")
plt.ylabel("rho")
plt.savefig('analysis1.png')
