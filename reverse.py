"""
Reverses calculations from test to find coefficient of friction
"""
import math
import numpy as np
import sympy as sp
import utils
from plotter import Plotter
from limit_surface import LimitSurface
from kinematic_car import KinematicCar

# Parameters
params = utils.load_params('config.yaml')

# Inputs
max_steering_angle = 0.10
min_turning_radius = params['wheelbase']/np.tan(max_steering_angle)
max_push = params['speed']/min_turning_radius
min_push = -1 * max_push
maxp_twist = np.zeros((6,))
minp_twist = np.zeros((6,))
maxp_twist[4:] = [params['speed'], max_push]
minp_twist[4:] = [params['speed'], min_push]

# Motion of kinematic_car
kcar = KinematicCar(params)
print(maxp_twist)
twist_converted = kcar.c2b_twists(maxp_twist, params['dist'])

# Reverse limit surface
params['radius'] = np.sqrt(2*(params['dist']**2)) # block radius in meters (assuming it is a circle)
r, theta = sp.symbols('r theta') # can use these variables in pressure distribution
params['pressure_distribution']= params['normal_force']/(math.pi * params['radius']**2) # circular pressure distribution
ls = LimitSurface(params)  
wrench = ls.find_wrench(twist_converted)
cof = abs(wrench[3] / wrench[4])

utils.print_params(params)

# Motion of kinematic_car
max_car = kcar.get_car_limits()

# Convert limits into car terms
print('Estimated COF: {}'.format(cof))
print('Pushing min turning radius: {}'.format(round(min_turning_radius,3)))
print('Pushing max steering angle: {}'.format(round(max_steering_angle,3)))

# Plot limits of car
plotter = Plotter()
plotter.set_ls_params(ls.ft_max, ls.m_max)
plotter.plot_vec(max_car, 'vs', normalize=True, color='black')
max_car[2] = max_car[2]*-1 # it's even
plotter.plot_vec(max_car, 'vs', normalize=True, color='black')

# Plot limits of pushing
plotter.plot_vec(minp_twist[3:], 'vs', normalize=True, color='purple')
plotter.plot_vec(maxp_twist[3:], 'vs', normalize=True, color='purple')

# Plot wrench
plotter.plot_vec(wrench[3:], 'ls', normalize=True, color='purple')

plotter.plot_all(ls)
