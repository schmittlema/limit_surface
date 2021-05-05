# The Limit Surface Class
# Does not contain car specific params
# Schmittle
import math
import sympy as sp
import numpy as np

import utils

class LimitSurface:

    def __init__(self, params): 
        """ 
        Constructor, creates initial limit surface
        """
        self.fc = params['coefficient']
        self.pd = params['pressure_distribution']
        self.radius = params['radius']
        self.nf = params['normal_force']

        self.create_surface()
    
    def create_surface(self):
        """
        Create limit surface ellipsoid approx

        Params: 
            None
        Returns:
            None
        """
        # Tangential max
        self.ft_max = self.fc * self.nf

        # Moment max
        r, theta = sp.symbols('r theta')
        integral_1 = sp.integrate(self.fc * self.pd * r**2, (r, 0, self.radius))
        self.m_max = sp.integrate(integral_1, (theta, 0, 2*math.pi))

    def find_normal(self, wrench):
        """
        Returns normal vector to limit surface from the point on the surface. 
        The point on the surface is defined by a wrench from the origin to the surface.

        Params: 
            wrench ([6x1] ndarray): Force wrench from 0,0,0 to limit surface
        Returns:
            twist ([6x1] ndarray): Velocity twist normal to limit surface
        """
        # take derivative of ellipsoid eq. at given wrench and use the coefficeints as vector
        return 2*wrench/np.array([1,1,1, self.ft_max**2,self.ft_max**2,self.m_max**2]) + np.array([wrench[3], wrench[4], wrench[5], 0, 0 ,0]) 

    def find_wrench(self, twist):
        """
        Given twist vector, find limit surface point (wrench) that the twist is normal to
        Params: 
            twist ([6x1] ndarray): twist we wish to find point on limit surface to 
        Returns:
            wrench ([6x1] ndarray): Force wrench from 0,0,0 to limit surface
        """
        # rearange find_normal to solve for wrench from twist
        wrench = twist / ((2/np.array([1,1,1,self.ft_max**2,self.ft_max**2,self.m_max**2])) + 1)
        return wrench

    def find_valid_twists(self, cone):
        """
        Find valid velocities given a frictioncone and limit surface

        Params: 
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
        Returns:
            twist_cone ([4x3x1] ndarray): Four vectors making the composite twist cone
        """
        ft_maxes = [self.ft_max]*cone.shape[0]
        m_maxes = [self.m_max]*cone.shape[0]
        cone = np.array(list(map(utils.fit_vec_to_surface, ft_maxes, m_maxes, cone)))
        twist_cone = map(self.find_normal, cone)
        return np.array(list(twist_cone))

    def create_friction_cone(self, cof, r_1, r_2):
        """
        Creates four vectors defining the composite wrench cone.
        Assumes two contact points

        Params: 
            cof (double): coefficient of friction

            r_1 ([2x1] ndarray): vector from COM to contact point 1

            r_2 ([2x1] ndarray): vector from COM to contact point 2

        Returns:
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
        """
        # define friction cones. f_3 & f_4 are currently equivalent to f_1, f_2
        f_1 = np.array([0, 0, 0, -cof * self.ft_max, self.ft_max, 0])
        f_2 = np.array([0, 0, 0, cof * self.ft_max, self.ft_max, 0])
        
        # moments
        m_1 = np.cross(r_1, f_1[3:])[2]
        m_2 = np.cross(r_1, f_2[3:])[2]
        m_3 = np.cross(r_2, f_1[3:])[2]
        m_4 = np.cross(r_2, f_2[3:])[2]

        # friction cone vectors
        F_1 = f_1 + np.array([0,0,0,0,0,m_1])
        F_2 = f_2 + np.array([0,0,0,0,0,m_2])
        F_3 = f_1 + np.array([0,0,0,0,0,m_3])
        F_4 = f_2 + np.array([0,0,0,0,0,m_4])

        return np.array([F_1, F_2, F_3, F_4])

