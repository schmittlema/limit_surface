import math
import numpy as np
import sympy as sp
import limit_surface as ls
from primitives import Primitives

# Pushing params
normal_force = 5.0 # N block weight
dist = 0.05 # COM to edge
radius = np.sqrt(2*dist**2) # block radius in meters (assuming it is a circle)
ls_coefficient = 0.8 # coefficient of friction between block and floor
coefficient = 0.8 # coefficient of friction between block and car

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
car_length = 0.33 # meters
max_steering_angle = 0.34 # radians
speed = 0.5 # m/s constant speed to push at

print()
print('CAR PARAMETERS')
print()
print('Car Length: {}'.format(car_length))
print('Maximum Steering Angle: {}'.format(max_steering_angle))
print('Speed: {}'.format(speed))
print()

# Pushing constraints
ls = ls.LimitSurface(ls_coefficient, pressure_distribution, radius, normal_force)
cone = ls.create_friction_cone(coefficient, np.array([-dist, -dist, 0.0]), np.array([dist, -dist, 0.0]))
twist_cone = ls.find_valid_twists(cone)

# Motion Primitives
mprims = Primitives(car_length, max_steering_angle)
dists = np.ones(twist_cone.shape[0]) * dist
twist_cone_converted = np.array(list(map(mprims.b2c_twists,twist_cone, dists)))
max_car, min_ls, max_ls = mprims.find_velocity_limits(twist_cone_converted, speed, dist)

# Plot limits of car
ls.plot_vec(max_car, 'vs', normalize=True, color='blue')
max_car[2] = max_car[2]*-1
ls.plot_vec(max_car, 'vs', normalize=True, color='blue')

# Plot limits of pushing
ls.plot_vec(min_ls, 'vs', normalize=True, color='purple')
ls.plot_vec(max_ls, 'vs', normalize=True, color='purple')

ls.plot_all(cone, twist_cone_converted)

