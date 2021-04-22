# The Pushing Car
# Controls that respect limit surface and kinematic constraints
# Assuming rear-axle Kinematic car model
# Schmittle

import numpy as np

class Kinematic_Car:

    def __init__(self, params): 
        """
        Constructor
        """
        self.speed = params['speed']
        self.wheelbase = params['wheelbase']
        self.bump2_front_axle = params['bumper2front_axle']
        self.max_steering_angle = params['max_steering_angle']
        self.w = self.wheelbase/np.tan(self.max_steering_angle) # H + W/2 in kinematic car model

    def find_velocity_limits(self, twist_cone, dist):
        """
        Find the angular velocity limits pusher given twist cone * max steering angle 

        Params: 
            twist_cone ([4x6x1] ndarray): Four vectors making the composite twist cone
            dist (double): COM to pusher 
        Returns:
           max_car ([3x1] ndarray): un-normalized vector for maximum angular velocity steering angle
           max_ls ([3x1] ndarray): un-normalized vector for pushing limits (max angular velocity) given twist cone 
           min_ls ([3x1] ndarray): un-normalized vector for pushing limits (min angular velocity) given twist cone
        """

        # car's limits
        min_turning_radius = self.wheelbase/np.tan(self.max_steering_angle)
        max_angular_velocity = self.speed/min_turning_radius
        max_car = np.array([0, self.speed, max_angular_velocity])

        # pushing limits
        min_ls, max_ls = self.pushing_limits(twist_cone)

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

        r = np.array([0.0, -(self.wheelbase + dist + self.bump2_front_axle), 0]) # vector from COM of block to center back axle
        omega = np.array([0.0, 0.0, block_twist[-1]]) # angular velocity
        diff_twist = np.cross(omega, r)
        twist_car_short = block_twist[3:] + diff_twist # add velocities

        # create full vector
        car_twist = np.zeros((6,))
        car_twist[3:] = twist_car_short

        return car_twist

    def stable_limits(self, dist, left_cone, right_cone):
        """
        Find angular velocity limits for stable pushing from output of STABLE
        Params: 
            dist (float): distance from block COM to edge
            left_cone ([2x4] ndarray): int_x, int_y, vec_x, vec_y. Cone that spans no-slip ICs
            right_cone ([2x4] ndarray): int_x, int_y, vec_x, vec_y 
        Returns:
           min_stable ([3x1] ndarray): un-normalized vector for minimum twist 
           max_stable ([3x1] ndarray): un-normalized vector for maximum twist 
        """ 

        # find where car COR line (rear axle) intersects bottom of cone
        # top of cone not possible with a car
        car_ic_line = -1*(self.wheelbase + self.bump2_front_axle + dist)
        c_l = (car_ic_line - left_cone[0,1])/left_cone[0,3]
        c_r = (car_ic_line - right_cone[0,1])/right_cone[0,3]
        left_int = c_l*left_cone[0,2:] + left_cone[0,:2]
        right_int = c_r*right_cone[0,2:] + right_cone[0,:2]

        # convert x values (turning radius) to max angular velocities
        left_w = self.speed/left_int[0]
        right_w = self.speed/right_int[0]

        min_stable = np.array([0, self.speed, left_w])
        max_stable = np.array([0, self.speed, right_w])

        return min_stable, max_stable

    def pushing_limits(self, twist_cone):
        """
        Find the velocity limits of pushing using block limits and car model

        Params: 
            twist_cone ([4x6x1] ndarray): Four vectors making the composite twist cone
        Returns:
           min_ls ([3x1] ndarray): un-normalized vector for minimum twist 
           max_ls ([3x1] ndarray): un-normalized vector for maximum twist 
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
        z = defining_vecs[0,5] + ((norm_vec[0] * defining_vecs[0,3] - norm_vec[1] * (self.speed - defining_vecs[0, 4]))/norm_vec[2])

        max_ls = np.array([0, self.speed, z])

        norm_vec = np.cross(defining_vecs[2,3:], defining_vecs[3,3:])
        z = defining_vecs[2,5] + ((norm_vec[0] * defining_vecs[2,3] - norm_vec[1] * (self.speed - defining_vecs[2, 4]))/norm_vec[2])

        min_ls = np.array([0, self.speed, z])

        return min_ls, max_ls
