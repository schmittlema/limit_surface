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

    def find_velocity_limits(self, twist_cone, speed, dist):
        """
        Find the angular velocity limits pusher given twist cone * max steering angle 

        Params: 
            twist_cone ([4x6x1] ndarray): Four vectors making the composite twist cone
            speed (double): fixed speed of car
            dist (double): COM to pusher 
        Returns:
           max_car ([3x1] ndarray): un-normalized vector for maximum steering angle
           max_ls ([3x1] ndarray): un-normalized vector for maximum angular velocity
        """

        # car's limits
        min_turning_radius = self.length/np.tan(self.max_steering_angle)
        max_angular_velocity = speed/min_turning_radius
        max_car = np.array([0, speed, max_angular_velocity])

        # pushing limits
        min_ls, max_ls = self.pushing_limits(twist_cone, speed)

        return max_car, min_ls, max_ls

    def b2c_twists(self, block_twist, dist):
        """
        Convert block twist to car twist
        v_car = v_block + w x r

        Params: 
            block_twist ([6x1] ndarray): Vector twist for block
            dist (double): COM to pusher 
        Returns:
            car_twist ([6x1] ndarray): Vector twist for car (rear axle)
        """

        r = np.array([0, -(self.length + dist), 0]) # vector from COM of block to center back axle
        omega = np.array([0,0,block_twist[-1]]) # angular velocity
        diff_twist = np.cross(omega, r)
        twist_car_short = block_twist[3:] + diff_twist # add velocities

        # create full vec
        car_twist = np.zeros((6,))
        car_twist[3:] = twist_car_short

        return car_twist

    def pushing_limits(self, twist_cone, speed):
        """
        Find the velocity limits of pushing 

        Params: 
            twist_cone ([4x6x1] ndarray): Four vectors making the composite twist cone
            speed (double): fixed speed of car
        Returns:
           min_ls ([3x1] ndarray): un-normalized vector for minimum angular velocity
           max_ls ([3x1] ndarray): un-normalized vector for maximum angular velocity
        """
        # TODO fix this so that limit vectors can be on the side edges too

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

        max_ls = np.array([0, speed, z])

        # find minimum angular velocity vecs
        min_omega_idx = np.argmin(twist_cone[:,-1])
        t_cone = np.zeros((twist_cone.shape))
        t_cone[:] = twist_cone
        t_cone[min_omega_idx,-1] = 1e5
        second_min_omega_idx = np.argmin(t_cone[:,-1])

        # define plane between them
        norm_vec = np.cross(twist_cone[min_omega_idx,3:], twist_cone[second_min_omega_idx,3:])

        # use plane equation to find max angular vel for a given vx = 0, vy
        z = twist_cone[min_omega_idx,5] + ((norm_vec[0] * twist_cone[min_omega_idx,3] - norm_vec[1] * (speed - twist_cone[min_omega_idx, 4]))/norm_vec[2])

        min_ls = np.array([0, speed, z])

        return min_ls, max_ls
