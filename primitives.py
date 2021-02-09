# The Pushing Car Motion Primitives
# Primitives respect limit surface and kinematic constraints
# Schmittle

import numpy as np


class Primitives:

    def __init__(self, car_length, bumper2frontaxle, max_steering_angle):
        """
        Constructor
        """
        self.length = car_length
        self.bump2_front_axle = bumper2frontaxle 
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

        r = np.array([0, -(self.length + dist + self.bump2_front_axle), 0]) # vector from COM of block to center back axle
        omega = np.array([0, 0, block_twist[-1]]) # angular velocity
        diff_twist = np.cross(omega, r)
        twist_car_short = block_twist[3:] + diff_twist # add velocities

        # create full vector
        car_twist = np.zeros((6,))
        car_twist[3:] = twist_car_short

        """
        print("block twist")
        print(block_twist)
        print("omega")
        print(omega)
        print("r")
        print(r)
        print("diff")
        print(diff_twist)
        print("car twist")
        print(twist_car_short)
        """

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
        # probably could have written this better...
        
        # find vectors the vx=0 plane would fall between
        vec_01 = twist_cone[1,:] - twist_cone[0,:] + np.tile(twist_cone[0,3:],2)
        vec_02 = twist_cone[2,:] - twist_cone[0,:] + np.tile(twist_cone[0,3:],2)
        vec_23 = twist_cone[3,:] - twist_cone[2,:] + np.tile(twist_cone[2,3:],2)
        vec_31 = twist_cone[1,:] - twist_cone[3,:] + np.tile(twist_cone[3,3:],2)

        diff_01 = np.sign(vec_01[0]*vec_01[3])
        diff_02 = np.sign(vec_02[0]*vec_02[3])
        diff_23 = np.sign(vec_23[0]*vec_23[3])
        diff_31 = np.sign(vec_31[0]*vec_31[3])
        diffs = [diff_01, diff_02, diff_23, diff_31]

        defining_vecs = []
        for idx, diff in enumerate(diffs):
            if diff < 0.0:
                # one of the diffs we want
                if idx == 0:
                   defining_vecs.append(twist_cone[0]) 
                   defining_vecs.append(twist_cone[1]) 
                if idx == 1:
                   defining_vecs.append(twist_cone[0]) 
                   defining_vecs.append(twist_cone[2]) 
                if idx == 2:
                   defining_vecs.append(twist_cone[2]) 
                   defining_vecs.append(twist_cone[3]) 
                if idx == 3:
                   defining_vecs.append(twist_cone[3]) 
                   defining_vecs.append(twist_cone[1]) 


        defining_vecs = np.array(defining_vecs)

        norm_vec = np.cross(defining_vecs[0,3:], defining_vecs[1,3:])
        z = defining_vecs[0,5] + ((norm_vec[0] * defining_vecs[0,3] - norm_vec[1] * (speed - defining_vecs[0, 4]))/norm_vec[2])

        max_ls = np.array([0, speed, z])

        norm_vec = np.cross(defining_vecs[2,3:], defining_vecs[3,3:])
        z = defining_vecs[2,5] + ((norm_vec[0] * defining_vecs[2,3] - norm_vec[1] * (speed - defining_vecs[2, 4]))/norm_vec[2])

        min_ls = np.array([0, speed, z])

        return min_ls, max_ls
