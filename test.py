import math
import numpy as np
import sympy as sp
from limit_surface import LimitSurface
from kinematic_car import Kinematic_Car

# Pushing params
normal_force = 5.0 # N block weight
dist = 0.5 # COM to edge 
radius = np.sqrt(2*dist**2) # block radius in meters (assuming it is a circle)
ls_coefficient = 0.5 # coefficient of friction between block and floor
coefficient = 0.4 # coefficient of friction between block and car

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
print('Bumer to front axle: {}'.format(bumper_to_front_axle))
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

# Plot limits of car
ls.plot_vec(max_car, 'vs', normalize=True, color='black')
max_car[2] = max_car[2]*-1
ls.plot_vec(max_car, 'vs', normalize=True, color='black')

# Plot limits of pushing
ls.plot_vec(min_ls, 'vs', normalize=True, color='purple')
ls.plot_vec(max_ls, 'vs', normalize=True, color='purple')

ls.plot_all(cone, twist_cone_converted)

# SBPL motion primitive calculation
"""
#WIP
timestep = 1.0 # seconds
resolution = 0.1 # meters. map resolution
num_angles = 16 # number of discretized angles 
angular_velocity = min_ls[-1]
angular_velocity = max_ls[-1]
#angular_velocity = 0 
start_angle = np.radians(22.5)
angle_resolution = 2*np.pi/num_angles

if angular_velocity != 0:
    turning_radius = speed/angular_velocity
    delta_theta = angular_velocity * timestep
    delta_x = turning_radius*(np.sin(delta_theta+start_angle) - np.sin(start_angle))
    delta_y = turning_radius*(-np.cos(delta_theta+start_angle) + np.cos(start_angle))
else:
    turning_radius = 1e5
    delta_theta = 0
    delta_x = speed * np.cos(start_angle) 
    delta_y = speed * np.sin(start_angle) 

print("SBPL PRIMITIVE INFO")
print("timestep: {}".format(timestep))
print("speed: {}".format(speed))
print("start angle: {}".format(start_angle))
print("angular_velocity: {}".format(angular_velocity))
print("turning radius: {}".format(turning_radius))
print()
print("delta_x: {}".format(delta_x))
print("delta_y: {}".format(delta_y))
print("delta_theta: {}".format(delta_theta))
print("res_delta_x: {}".format(delta_x/resolution))
print("res_delta_y: {}".format(delta_y/resolution))
print("res_delta_theta: {}".format(delta_theta/angle_resolution))
print("length: {}".format(np.sqrt(delta_x**2 + delta_y**2)))
"""
