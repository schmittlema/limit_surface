# The Pushing Car Motion Primitives
# Primitives respect limit surface and kinematic constraints
# Schmittle

import numpy as np


class Primitives:

    def __init__(self, car_length, max_steering_angle):
        """
        Constructor
        """
        self.length = car_length
        self.max_steering_angle = max_steering_angle
        self.w = self.length/np.tan(self.max_steering_angle) # H + W/2 in kinematic car model

    def find_velocity_limits(self, twist_cone, speed, radius):
        """
        Find the angular velocity limits pusher given twist cone * max steering angle 

        Params: 
            twist_cone ([4x6x1] ndarray): Four vectors making the composite twist cone
            speed (double): fixed speed of car
            radius (double): object radius
        Returns:
            ???
        """

        # car's limits
        min_turning_radius = np.sqrt(self.w**2 + (self.length+radius)**2)
        max_angular_velocity = speed/min_turning_radius

        # pushing limits

        # find maximum angular velocity vecs
        max_omega_idx = np.argmax(twist_cone[:,-1])
        t_cone = np.zeros((twist_cone.shape))
        t_cone[:] = twist_cone
        t_cone[max_omega_idx,-1] = -1e5
        second_max_omega_idx = np.argmax(t_cone[:,-1])

        # define plane between them
        norm_vec = np.cross(twist_cone[max_omega_idx,3:], twist_cone[second_max_omega_idx,3:])

        # use plane equation to find max angular vel for a given vx = 0, vy
        z = twist_cone[max_omega_idx,5] + ((norm_vec[0] * twist_cone[max_omega_idx,3] - norm_vec[1] * (speed - twist_cone[max_omega_idx, 4]))/norm_vec[2])
        print(z, max_angular_velocity)



        

    
