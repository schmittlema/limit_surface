import math
import numpy as np
import sympy as sp
from limit_surface import LimitSurface
from kinematic_car import Kinematic_Car

# Pushing params
normal_force = 5.0 # N block weight
dist = 0.05 # COM to edge of block
radius = np.sqrt(2*dist**2) # block radius in meters (assuming it is a circle)
ls_coefficient = 0.5 # coefficient of friction between block and floor
coefficient = 0.5 # coefficient of friction between block and car

print()
print('BLOCK PARAMETERS')
print()
print('Normal Force: {}'.format(normal_force))
print('COM to Edge Distance: {}'.format(dist))
print('Radius of Contact Circle Approx.: {}'.format(radius))
print('Block & Floor Friction Coefficient: {}'.format(ls_coefficient))
print('Block & Pusher Friction Coefficient: {}'.format(coefficient))
print()

r, theta = sp.symbols('r theta') # can use these variables in pressure distribution

# TODO square pressure distribution
pressure_distribution = normal_force/(math.pi * radius**2) # circular press. distribution

# Car params
wheelbase = 0.225 # meters
bumper_to_front_axle = 0.05 # meters
max_steering_angle = 0.35 # radians
speed = 0.5 # m/s constant speed to push at

print()
print('CAR PARAMETERS')
print()
print('Wheelbase: {}'.format(wheelbase))
print('Bumper to front axle: {}'.format(bumper_to_front_axle))
print('Maximum Steering Angle: {}'.format(max_steering_angle))
print('Speed: {}'.format(speed))
print()

# Pushing constraints
ls = LimitSurface(ls_coefficient, pressure_distribution, radius, normal_force)
cone = ls.create_friction_cone(coefficient, np.array([-dist, -dist, 0.0]), np.array([dist, -dist, 0.0]))
twist_cone = ls.find_valid_twists(cone)

# Motion of kinematic_car
kcar = Kinematic_Car(wheelbase, bumper_to_front_axle, max_steering_angle)
dists = np.ones(twist_cone.shape[0]) * dist
twist_cone_converted = np.array(list(map(kcar.b2c_twists,twist_cone, dists)))
max_car, min_ls, max_ls = kcar.find_velocity_limits(twist_cone_converted, speed, dist)

# Convert limits into car terms
min_turning_radius = speed/min(max_ls[-1], max_car[-1])
max_steering_angle = np.arctan(wheelbase/min_turning_radius)
print('Pushing min turning radius: {}'.format(min_turning_radius))
print('Pushing max steering angle: {}'.format(max_steering_angle))

# Plot limits of car
ls.plot_vec(max_car, 'vs', normalize=True, color='black')
max_car[2] = max_car[2]*-1
ls.plot_vec(max_car, 'vs', normalize=True, color='black')

# Plot limits of pushing
ls.plot_vec(min_ls, 'vs', normalize=True, color='purple')
ls.plot_vec(max_ls, 'vs', normalize=True, color='purple')

ls.plot_all(cone, twist_cone_converted)
