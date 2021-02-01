# The Limit Surface Class
# Schmittle

import math
import sympy as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class LimitSurface:

    def __init__(self, coefficient, pressure, radius, normal_force):
        """ 
        Constructor, creates initial limit surface
        """
        self.fc = coefficient
        self.pd = pressure # Assuming circular pressure distribution for now
        self.radius = radius
        self.nf = normal_force
        
        self.create_surface()
    
    def create_surface(self):
        """
        Create limit surface ellipse approx

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
            wrench ([3x1] ndarray): Force wrench from 0,0,0 to limit surface
        Returns:
            twist ([3x1] ndarray): Velocity twist normal to limit surface
        """
        # take derivative of ellipsoid eq. at given wrench and use the coefficeints as vector
        return 2*wrench/np.array([1,1,1, self.ft_max**2,self.ft_max**2,self.m_max**2]) + np.array([wrench[3], wrench[4], wrench[5], 0, 0 ,0]) 

    def find_valid_twists(self, cone):
        """
        Find valid velocities given a frictioncone and limit surface

        Params: 
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
        Returns:
            twist_cone ([4x3x1] ndarray): Four vectors making the composite twist cone
        """
        cone = np.array(list(map(self.fit_vec_to_surface, cone)))
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
        m_1 = np.cross(f_1[3:], r_1)[2]
        m_2 = np.cross(f_2[3:], r_1)[2]
        m_3 = np.cross(f_1[3:], r_2)[2]
        m_4 = np.cross(f_2[3:], r_2)[2]

        # friction cone vectors
        F_1 = f_1 + np.array([0,0,0,0,0,m_1])
        F_2 = f_2 + np.array([0,0,0,0,0,m_2])
        F_3 = f_1 + np.array([0,0,0,0,0,m_3])
        F_4 = f_2 + np.array([0,0,0,0,0,m_4])
        return np.array([F_1, F_2, F_3, F_4])

    def fit_vec_to_surface(self, vec):
        """
        Adjusts vector to lie on the surface for easy viz

        Params: 
            vec ([3x1] ndarray): Vector
        Returns:
            vec ([3x1] ndarray): Vector fit to limit surface 
        """

        v_hat = vec / np.linalg.norm(vec)

        # find constant multiple and fit
        # solving for c*F_1_hat that fits ellipsoid equation
        fit = lambda v_hat: v_hat * math.sqrt((1/(v_hat[3]**2/self.ft_max**2 + v_hat[4]**2/self.ft_max**2 + v_hat[5]**2/self.m_max**2)))

        return fit(v_hat)

    def plot_ls(self, ax=None, cone=None):
        """
        Plots limit surface

        Params: 
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
        Returns:
            None
        """

        if ax is None:
            fig = plt.figure(figsize=(10,8))  # Square figure
            ax = fig.add_subplot(111, projection='3d')

        coefs = np.array([1/(self.ft_max**2), 1/(self.ft_max**2), 1/self.m_max**2]
                , dtype=np.float)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1 

        # Radii corresponding to the coefficients:
        rx, ry, rz = 1/np.sqrt(coefs)

        # Set of all spherical angles:
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        # Cartesian coordinates that correspond to the spherical angles:
        # (this is the equation of an ellipsoid):
        x = rx * np.outer(np.cos(u), np.sin(v))
        y = ry * np.outer(np.sin(u), np.sin(v))
        z = rz * np.outer(np.ones_like(u), np.cos(v))

        # Vectors (6D because not all vectors are from origin)

        # example vector
        v = np.array([[0, 0, 0, x[20][20], y[20][20], z[20][20]]])
        X, Y, Z, U, V, W = zip(*v)
        ax.quiver(X, Y, Z, U, V, W, color = (92/255,38/255,134/255,1.0), arrow_length_ratio=0.2, linewidths=4.0)

        # example normal
        normal_vector = self.find_normal(v[0])
        X, Y, Z, U, V, W = zip(*[normal_vector])
        ax.quiver(X, Y, Z, U, V, W, color = (54/255,205/255,196/255,1.0), arrow_length_ratio=0.3, linewidths=3.0)

        # plot friction cone
        if cone is not None:
            self.plot_cone(cone, ax, rx, ry, rz)
            self.plot_twists(cone, ax)

        # Plot:

        ax.plot_surface(x, y, z,  rstride=5, cstride=5, color= (92/255, 90/255, 90/255, 0.09))
        ax.set_xlabel(r'$f_x$')
        ax.set_ylabel(r'$f_y$')
        ax.set_zlabel(r'$m_z$')
        ax.zaxis.set_rotate_label(False) 
        ax.yaxis.set_rotate_label(False) 
        ax.xaxis.set_rotate_label(False) 
        max_radius = max(rx, ry, rz) + 1.0
        ax.set_xlim(-max_radius, max_radius)
        ax.set_ylim(max_radius, -max_radius)
        ax.set_zlim(-max_radius, max_radius)

    def plot_cone(self, cone, ax, rx, ry, rz):
        """
        Plots friction cone on limit surface

        Params: 
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone

            ax (matplotlib axis): matplotlib graph axis
            
            rx, ry, rz (doubles): coefficients of limit surface
        Returns:
            None
        """
        cone = np.array(list(map(self.fit_vec_to_surface, cone)))

        X, Y, Z, U, V, W = zip(*cone)
        ax.quiver(X, Y, Z, U, V, W, color = (244/255,214/255,118/255,1.0), arrow_length_ratio=0.2, linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[0][3:], cone[1][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= (244/255, 214/255, 118/255, 1.0), linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[0][3:], cone[2][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= (244/255, 214/255, 118/255, 1.0), linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[2][3:], cone[3][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= (244/255, 214/255, 118/255, 1.0), linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[3][3:], cone[1][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= (244/255, 214/255, 118/255, 1.0), linewidths=4.0)

    def calculate_connector(self, p1, p2, rx, ry, rz):
        """
        Calculates line on limit surface connecting two points

        #WARNING this is not a 100% accurate drawing becuase of the way I am drawing the horizontal component.

        Params: 
            p1, p2 ([3x1] ndarray): Vectors defining the points to connect
            
            rx, ry, rz (doubles): coefficients of limit surface
        Returns:
            x_l, y_l, z_l ([3x100] ndarrays): Positions of points along connection 
        """

        # vector between the points
        connect_vec = p2 - p1
        t = np.linspace(0, 1.0, 100)
        

        # indexs along connect_vec
        vecs = np.outer(t,connect_vec) + np.tile(p1,(t.shape[0],1))

        # format
        full_vecs = np.zeros((vecs.shape[0], 6))
        full_vecs[:,3:] = vecs

        # fit to surface
        vecs_fit = np.array(list(map(self.fit_vec_to_surface, full_vecs)))

        # separate components
        x_l = np.array([vecs_fit[:,3]])
        y_l = np.array([vecs_fit[:,4]])
        z_l = np.array([vecs_fit[:,5]])

        return x_l, y_l, z_l

    def plot_twists(self, cone, ax):
        """
        Plots friction cone on limit surface

        Params: 
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone

            ax (matplotlib axis): matplotlib graph axis
            
        Returns:
            None
        """
        cone = np.array(list(map(self.fit_vec_to_surface, cone)))

        twist_cone = self.find_valid_twists(cone)

        X, Y, Z, U, V, W = zip(*twist_cone)
        ax.quiver(X, Y, Z, U, V, W, color = (54/255,205/255,196/255,1.0), arrow_length_ratio=0.2, linewidths=2.0)


    def plot_all(self, cone=None, vecs=None):
        """
        Plots limit surface & twist sphere with cones and/or additional vectors

        Params: 
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
            vecs ([-1x3x1] ndarray): Vectors to plot 
        Returns:
            None
        """
        fig = plt.figure(figsize=(20,8))  # Square figure
        ax_ls = fig.add_subplot(111, projection='3d')
        ax_vs = fig.add_subplot(122, projection='3d')

        self.plot_ls(ax_ls, cone)

        plt.show()


