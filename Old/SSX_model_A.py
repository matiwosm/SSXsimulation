"""SSX_model_A.py
    
    This is the *simplest* model we will consider for modelling spheromaks evolving in the SSX wind tunnel.
    
    Major simplificiations fall in two categories
    
    Geometry
    --------
    We consider a square duct using parity bases (sin/cos) in all directions.
    
    Equations
    ---------
    The equations themselves are those from Schaffner et al (2014), with the following simplifications
    
    * hall term off
    * constant eta instead of Spitzer
    * no wall recycling term
    * no mass diffusion
    
    For this first model, rather than kinematic viscosity nu and thermal
    diffusivitiy chi varying with density rho as they should, we are here
    holding them *constant*. This dramatically simplifies the form of the
    equations in Dedalus.
    
    We use the vector potential, and enforce the Coulomb Gauge, div(A) = 0.
    Usage:
    SSX_model_A.py  <config_file>
    
    """

import os
import sys
import time
import pathlib
import numpy as np


import dedalus.public as de
from dedalus.extras import flow_tools

from configparser import ConfigParser

from matplotlib import pyplot
import matplotlib as mpl

from spheromak import spheromak_A

import logging
logger = logging.getLogger(__name__)

runconfig = ConfigParser()
config_file = pathlib.Path(sys.argv[-1])
runconfig.read(str(config_file))
logger.info("Using config file {}".format(config_file))

# for optimal efficiency: nx should be divisible by mesh[0], ny by mesh[1], and
# nx should be close to ny. Bridges nodes have 28 cores, so mesh[0]*mesh[1]
# should be a multiple of 28.
# Domain Setup
domain_setup = runconfig['domain']
nx = domain_setup.getint('nx')
ny = domain_setup.getint('ny')
nz = domain_setup.getint('nz')
r = domain_setup.getint('r')
length = domain_setup.getint('length')

# for 3D runs, you can divide the work up over two dimensions (x and y).
# The product of the two elements of mesh *must* equal the number
# of cores used.
# mesh = None
mesh = [14,12]
params = runconfig['params']
# kappa is heat conductivity
kappa = params.getfloat('kappa')
# mu is viscosity
mu = params.getfloat('mu')
# eta is resistivity
eta = params.getfloat('eta')

rho0 = params.getfloat('rho0')
# gamma is the adiabatic constant
gamma = params.getfloat('gamma')

x = de.SinCos('x', nx, interval=(-r, r))
y = de.SinCos('y', ny, interval=(-r, r))
z = de.SinCos('z', nz, interval=(0,length))

domain = de.Domain([x,y,z],grid_dtype='float', mesh=mesh)


SSX = de.IVP(domain, variables=['lnrho','T', 'vx', 'vy', 'vz', 'Ax', 'Ay', 'Az', 'phi'])
###########################################################################################################################################
#-----------------------------------------------Meta_Parameters---------------------------------------------------------------------------#
###########################################################################################################################################
SSX.meta['T','lnrho']['x', 'y', 'z']['parity'] = 1
SSX.meta['phi']['x', 'y', 'z']['parity'] = -1

SSX.meta['vx']['y', 'z']['parity'] =  1
SSX.meta['vx']['x']['parity'] = -1
SSX.meta['vy']['x', 'z']['parity'] = 1
SSX.meta['vy']['y']['parity'] = -1
SSX.meta['vz']['x', 'y']['parity'] = 1
SSX.meta['vz']['z']['parity'] = -1

SSX.meta['Ax']['y', 'z']['parity'] =  -1
SSX.meta['Ax']['x']['parity'] = 1
SSX.meta['Ay']['x', 'z']['parity'] = -1
SSX.meta['Ay']['y']['parity'] = 1
SSX.meta['Az']['x', 'y']['parity'] = -1
SSX.meta['Az']['z']['parity'] = 1
#########################################################################################################################################
#-----------------------------------------------Parameters------------------------------------------------------------------------------#
#########################################################################################################################################
SSX.parameters['mu'] = mu
# chi is normalized thermal diffusivity
SSX.parameters['chi'] = kappa/rho0
# nu normalized viscosity
SSX.parameters['nu'] = mu/rho0
# eta is resistivity
SSX.parameters['eta'] = eta
# gamma is the adiabatic constant
SSX.parameters['gamma'] = gamma
###########################################################################################################################################
#-----------------------------------------------Substitutions-----------------------------------------------------------------------------#
###########################################################################################################################################
SSX.substitutions['divv'] = "dx(vx) + dy(vy) + dz(vz)"
SSX.substitutions['vdotgrad(A)'] = "vx*dx(A) + vy*dy(A) + vz*dz(A)"
SSX.substitutions['Bdotgrad(A)'] = "Bx*dx(A) + By*dy(A) + Bz*dz(A)"
SSX.substitutions['Lap(A)'] = "dx(dx(A)) + dy(dy(A)) + dz(dz(A))"
SSX.substitutions['Bx'] = "dy(Az) - dz(Ay)"
SSX.substitutions['By'] = "dz(Ax) - dx(Az)"
SSX.substitutions['Bz'] = "dx(Ay) - dy(Ax)"

# Coulomb Gauge implies J = -Laplacian(A)
SSX.substitutions['jx'] = "-Lap(Ax)"
SSX.substitutions['jy'] = "-Lap(Ay)"
SSX.substitutions['jz'] = "-Lap(Az)"
SSX.substitutions['J2'] = "jx**2 + jy**2 + jz**2"
SSX.substitutions['rho'] = "exp(lnrho)"
# CFL substitutions
SSX.substitutions['Va_x'] = "Bx/sqrt(rho)"
SSX.substitutions['Va_y'] = "By/sqrt(rho)"
SSX.substitutions['Va_z'] = "Bz/sqrt(rho)"
SSX.substitutions['Cs'] = "sqrt(gamma*T)"
########################################################################################################################################
#-----------------------------------------------Equations------------------------------------------------------------------------------#
########################################################################################################################################
# Continuity
SSX.add_equation("dt(lnrho) + divv = -vdotgrad(lnrho)")

# Momentum
SSX.add_equation("dt(vx) + dx(T) - nu*Lap(vx) = T*dx(rho) - vdotgrad(vx) + (jy*Bz - jz*By)/rho")
SSX.add_equation("dt(vy) + dy(T) - nu*Lap(vy) = T*dy(rho) - vdotgrad(vy) + (jz*Bx - jx*Bz)/rho")
SSX.add_equation("dt(vz) + dz(T) - nu*Lap(vz) = T*dz(rho) - vdotgrad(vz) + (jx*By - jy*Bx)/rho")

# MHD equations: A
# Resistive Ohm's Law: A-form
SSX.add_equation("dt(Ax) + eta*jx + dx(phi) = vy*Bz - vz*By")
SSX.add_equation("dt(Ay) + eta*jy + dy(phi) = vz*Bx - vx*Bz")
SSX.add_equation("dt(Az) + eta*jz + dz(phi) = vx*By - vy*Bx")
# Coulomb Gauge for A
# Away from walls
SSX.add_equation("dx(Ax) + dy(Ay) + dz(Az) = 0", condition="(nx != 0) or (ny != 0) or (nz != 0)")
# Density condition at the orgin
# Electric Potential at the walls
SSX.add_equation("phi = 0", condition="(nx == 0) and (ny == 0) and (nz == 0)")


# Energy
SSX.add_equation("dt(T) - (gamma - 1) * chi*Lap(T) = - (gamma - 1) * T * divv  - vdotgrad(T) + (gamma - 1)*eta*J2")
##########################################################################################################################################
#-----------------------------------------------Solver_Setup-----------------------------------------------------------------------------#
##########################################################################################################################################
solver = SSX.build_solver(de.timesteppers.RK443)
logger.info('Solver built')

run = runconfig['run']
# Initial timestep
dt = run.getfloat('dt')

# Integration parameters
if run.getfloat('stop_sim_time'):
    solver.stop_sim_time = run.getfloat('stop_sim_time')
else:
    solver.stop_sim_time = np.inf
if run.getfloat('stop_wall_time'):
    solver.stop_wall_time = run.getfloat('stop_wall_time')
else:
    solver.stop_wall_time = np.inf
if run.getfloat('stop_iteration'):
    solver.stop_iteration = run.getfloat('stop_iteration')
else:
    solver.stop_iteration = np.inf


# Initial conditions
Ax = solver.state['Ax']
Ay = solver.state['Ay']
Az = solver.state['Az']
lnrho = solver.state['lnrho']
T = solver.state['T']
vz = solver.state['vz']

x = domain.grid(0)
y = domain.grid(1)
z = domain.grid(2)
fullGrid = x*y*z

# Initial condition parameters
initial = runconfig['initial']
R = r
L = R
#lambda_rho = 0.4 # half-width of transition region for initial conditions
lambda_rho = initial.getfloat('lambda_rho')
rho_min = initial.getfloat('rho_min')
T0 = initial.getfloat('T0')
##########################################################################################################################################
#-----------------------------------------------Spheromak Initial Conditions-------------------------------------------------------------#
##########################################################################################################################################
## Spheromak initial condition
aa_x, aa_y, aa_z = spheromak_A(domain, center=(0,0, 0), R=R, L=L)
# The vector potential is subject to some perturbation.
Ax['g'] = aa_x*(1 + 0.5*x*np.exp(-z**2))
Ay['g'] = aa_y*(1 + 0.5*x*np.exp(-z**2))
Az['g'] = aa_z*(1 + 0.5*x*np.exp(-z**2))

limits = runconfig['limits']
max_vel = limits.getfloat('max_vel')

vz['g'] = -np.tanh(8*z - 15)*max_vel/2 + max_vel/2


# Plasma number density initialization.
for i in range(x.shape[0]):
    xVal = x[i,0,0]
    for j in range(y.shape[1]):
        yVal = y[0,j,0]
        for k in range(z.shape[2]):
            zVal = z[0,0,k]
            if(zVal <= 1- lambda_rho):
                fullGrid[i][j][k] = 1
            elif(zVal<=1+lambda_rho and zVal >= 1-lambda_rho and (np.sqrt(xVal**2 + yVal**2)<R)):
                fullGrid[i][j][k] = (1 + rho_min)/2 + (1 - rho_min)/2*np.sin((1-zVal) * np.pi/(2*lambda_rho)) #rho_min + rho_min*np.cos(zVal*np.pi/(2*lambda_rho))
            else:
                fullGrid[i][j][k] = rho_min


rho0 = domain.new_field()
rho0['g'] = fullGrid

lnrho['g'] = np.log(rho0['g'])
T['g'] = T0 * rho0['g']**(gamma - 1)


# analysis output
#data_dir = './'+sys.argv[0].split('.py')[0]
file_handling = runconfig['file_handling']
wall_dt_checkpoints = file_handling.getfloat('wall_dt_checkpoints')  
output_cadence = file_handling.getfloat('output_cadence')  # This is in simulation time units
#output_cadence = .05 # This is in simulation time units

'''checkpoint = solver.evaluator.add_file_handler('checkpoints2', max_writes=1, wall_dt=wall_dt_checkpoints, mode='overwrite')
    checkpoint.add_system(solver.state, layout='c')'''
field_max_writes = file_handling.getint('field_max_writes')
field_writes = solver.evaluator.add_file_handler('fields', max_writes=field_max_writes, sim_dt = output_cadence, mode='overwrite')
field_writes.add_task('vx')
field_writes.add_task('vy')
field_writes.add_task('vz')
field_writes.add_task('Bx')
field_writes.add_task('By')
field_writes.add_task('Bz')
field_writes.add_task("exp(lnrho)", name='rho')
field_writes.add_task('T')

parameter_max_writes = file_handling.getint('parameter_max_writes')
parameter_writes = solver.evaluator.add_file_handler('parameters', max_writes=parameter_max_writes, sim_dt=output_cadence, mode='overwrite')
parameter_writes.add_task('mu')
parameter_writes.add_task('eta')
parameter_writes.add_task('nu')
parameter_writes.add_task('chi')
parameter_writes.add_task('gamma')

load_max_writes = file_handling.getint('load_max_writes')
load_writes = solver.evaluator.add_file_handler('load_data', max_writes=load_max_writes, sim_dt = output_cadence, mode='overwrite')
load_writes.add_task('vx')
load_writes.add_task('vy')
load_writes.add_task('vz')
load_writes.add_task('Ax')
load_writes.add_task('Ay')
load_writes.add_task('Az')
load_writes.add_task('lnrho')
load_writes.add_task('T')
load_writes.add_task('phi')



# Flow properties
cfl_params = runconfig['cfl_params']

flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("sqrt(vx*vx + vy*vy + vz*vz) / nu", name='Re_k')
flow.add_property("sqrt(vx*vx + vy*vy + vz*vz) / eta", name='Re_m')
flow.add_property("sqrt(vx*vx + vy*vy + vz*vz)", name='flow_speed')
flow.add_property("sqrt(vx*vx + vy*vy + vz*vz) / sqrt(T)", name='Ma')

char_time = cfl_params.getfloat('char_time')
CFL_safety = cfl_params.getfloat('cfl_safety') # this should be set to a characteristic time in the problem (the alfven crossing time of the tube, for example)

CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=1, safety=CFL_safety,
                     max_change=1.5, min_change=0.5, max_dt=output_cadence, threshold=0.05)
CFL.add_velocities(('vx', 'vy', 'vz'))
CFL.add_velocities(('Va_x', 'Va_y', 'Va_z'))
CFL.add_velocities(( 'Cs', 'Cs', 'Cs'))


good_solution = True
# Main loop
try:
    logger.info('Starting loop')
    logger_string = 'kappa: {:.3g}, mu: {:.3g}, eta: {:.3g}, dt: {:.3g}'.format(kappa, mu, eta, dt)
    logger.info(logger_string)
    start_time = time.time()
    while solver.ok and good_solution:
        dt = CFL.compute_dt()
        solver.step(dt)
        
        if (solver.iteration-1) % 1 == 0:
            logger_string = 'iter: {:d}, t/tb: {:.2e}, dt/tb: {:.2e}, sim_time: {:.4e}, dt: {:.2e}'.format(solver.iteration, solver.sim_time/char_time, dt/char_time, solver.sim_time, dt)
            Re_k_avg = flow.grid_average('Re_k')
            Re_m_avg = flow.grid_average('Re_m')
            v_avg = flow.grid_average('flow_speed')
            logger_string += ' Max Re_k = {:.2g}, Avg Re_k = {:.2g}, Max Re_m = {:.2g}, Avg Re_m = {:.2g}, Max vel = {:.2g}, Avg vel = {:.2g}, Max Ma = {:.1g}'.format(flow.max('Re_k'), Re_k_avg, flow.max('Re_m'),Re_m_avg, flow.max('flow_speed'), v_avg, flow.max('Ma'))
            logger.info(logger_string)
            if not np.isfinite(Re_k_avg):
                good_solution = False
                logger.info("Terminating run.  Trapped on Reynolds = {}".format(Re_k_avg))
            if not np.isfinite(Re_m_avg):
                good_solution = False
                logger.info("Terminating run. Trapped on magnetic Reynolds = {}".format(Re_m_avg))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))
    logger.info('Iter/sec: {:g}'.format(solver.iteration/(end_time-start_time)))


