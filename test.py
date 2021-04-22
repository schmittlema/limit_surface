import math
import numpy as np
import sympy as sp
import utils
from plotter import Plotter
from limit_surface import LimitSurface
from stable import STABLE
from kinematic_car import Kinematic_Car

# Parameters
params = utils.load_params('config.yaml')

# Block definitions
params['radius'] = np.sqrt(2*(params['dist']**2)) # block radius in meters (assuming it is a circle)
r, theta = sp.symbols('r theta') # can use these variables in pressure distribution
# TODO square pressure distribution
params['pressure_distribution']= params['normal_force']/(math.pi * params['radius']**2) # circular pressure distribution
r1 = np.array([-params['dist'], params['dist']])
r2 = np.array([params['dist'], params['dist']])
r3 = np.array([-params['dist'], -params['dist']])
r4 = np.array([params['dist'], -params['dist']])

utils.print_params(params)

# Pushing constraints
ls = LimitSurface(params) 
cone = ls.create_friction_cone(params['coefficient'], r3, r4)
twist_cone = ls.find_valid_twists(cone)

# Motion of kinematic_car
kcar = Kinematic_Car(params)
dists = np.ones(twist_cone.shape[0]) * params['dist']
twist_cone_converted = np.array(list(map(kcar.b2c_twists,twist_cone, dists)))
max_car, min_ls, max_ls = kcar.find_velocity_limits(twist_cone_converted, params['dist'])

# Convert limits into car terms
min_turning_radius = params['speed']/min(max_ls[-1], max_car[-1])
max_steering_angle = np.arctan(params['wheelbase']/min_turning_radius)
print('Pushing min turning radius: {}'.format(round(min_turning_radius,3)))
print('Pushing max steering angle: {}'.format(round(max_steering_angle,3)))

# STABLE -- alternative to friction cone through limit surface
stable = STABLE()
noslip_leftc, noslip_rightc = stable.calculate_noslip_ics(cone, np.array([r1, r2, r3, r4]))
min_stable, max_stable = kcar.stable_limits(params['dist'], noslip_leftc, noslip_rightc)
stable_min_turning_radius = params['speed']/min(max_stable[-1], max_car[-1])
stable_max_steering_angle = np.arctan(params['wheelbase']/stable_min_turning_radius)
print('STABLE Pushing min turning radius: {}'.format(round(stable_min_turning_radius,3)))
print('STABLE Pushing max steering angle: {}'.format(round(stable_max_steering_angle,3)))

# Plot limits of car
plotter = Plotter()
plotter.set_ls_params(ls.ft_max, ls.m_max)
plotter.plot_vec(max_car, 'vs', normalize=True, color='black')
max_car[2] = max_car[2]*-1 # it's even
plotter.plot_vec(max_car, 'vs', normalize=True, color='black')

# Plot limits of pushing
plotter.plot_vec(min_ls, 'vs', normalize=True, color='purple')
plotter.plot_vec(max_ls, 'vs', normalize=True, color='purple')

# Plot limits of pushing with STABLE
plotter.plot_vec(min_stable, 'vs', normalize=True, color='pink')
plotter.plot_vec(max_stable, 'vs', normalize=True, color='pink')

plotter.plot_all(ls, cone, twist_cone_converted)
