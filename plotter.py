# Plotting - Helps separate logic from ploting stuff
# Schmittle

import numpy as np
import matplotlib
import warnings
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D

import utils

class Plotter:

    def __init__(self):
        self.ls_plot_vectors = []
        self.vs_plot_vectors = []
        self.colors = {'purple': (92/255,38/255,134/255,1.0), \
                        'blue': (54/255,205/255,196/255,1.0), \
                        'yellow': (244/255,214/255,118/255,1.0), \
                        'gray': (92/255, 90/255, 90/255, 0.09), \
                        'black':(20/255,20/255,20/255,1.0),\
                        'pink':(255/255,22/255,114/255,1.0)}

    def plot_ls(self, ls, ax=None, cone=None):
        """
        Plots limit surface

        Params: 
            ls (LimitSurface): limit surface instance
            ax (matplotlib axis): axis to plot on
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
        Returns:
            None
        """

        # If ax is none, create a figure
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
        for v in self.ls_plot_vectors:
            X, Y, Z, U, V, W, color = zip(*[v])
            ax.quiver(float(X[0]), float(Y[0]), float(Z[0]), float(U[0]), float(V[0]), float(W[0]), color = self.colors[color[0]], arrow_length_ratio=0.2, linewidths=4.0)

            normal_vector = ls.find_normal(np.array(v[:-1], dtype=np.float))
            X, Y, Z, U, V, W = zip(*[normal_vector])
            ax.quiver(float(X[0]), float(Y[0]), float(Z[0]), float(U[0]), float(V[0]), float(W[0]), color = self.colors['blue'], arrow_length_ratio=0.3, linewidths=3.0)

        # plot friction cone
        if cone is not None:
            self.plot_cone(cone, ax, rx, ry, rz)
            twist_cone = np.array(ls.find_valid_twists(cone), dtype=np.double)
            self.plot_twists_on_ls(twist_cone, ax)

        # plot limit surface
        ax.plot_surface(x, y, z,  rstride=5, cstride=5, color= self.colors['gray'])

        # Plot design
        self.style_axis(ax, 'Limit Surface', '$f_x$', '$f_y$', '$m_z$', max(rx,ry,rz))

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
        ft_maxes = [self.ft_max]*cone.shape[0]
        m_maxes = [self.m_max]*cone.shape[0]
        cone = np.array(list(map(utils.fit_vec_to_surface, ft_maxes, m_maxes, cone)))

        X, Y, Z, U, V, W = zip(*cone)
        ax.quiver(X, Y, Z, U, V, W, color = self.colors['black'], arrow_length_ratio=0.1, linewidths=2.0)

        self.connect_cone(ax, cone, rx, ry, rz)

    def connect_cone(self, ax, cone, rx, ry, rz, color=None):
        """
        Calculates line on limit surface connecting two points

        Params: 
            ax (matplotlib axis): axis to plot on

            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone

            rx, ry, rz (doubles): coefficients of ellipsoid

            color (str): one of the color strings from self.color 
        Returns:
            None
        """
        # deprecation warning from numpy with matplotlib functions
        warnings.filterwarnings("ignore")

        if color is None:
            color = 'black'

        x_l, y_l, z_l = self.calculate_connector(cone[0][3:], cone[1][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= self.colors[color], linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[0][3:], cone[2][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= self.colors[color], linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[2][3:], cone[3][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= self.colors[color], linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(cone[3][3:], cone[1][3:], rx, ry, rz)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= self.colors[color], linewidths=4.0)

    def calculate_connector(self, p1, p2, rx, ry, rz):
        """
        Calculates line on limit surface connecting two points

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

        rx = np.ones(full_vecs.shape[0])*rx
        ry = np.ones(full_vecs.shape[0])*ry
        rz = np.ones(full_vecs.shape[0])*rz
        # fit to surface
        ft_maxes = [self.ft_max]*full_vecs.shape[0]
        m_maxes = [self.m_max]*full_vecs.shape[0]
        vecs_fit = np.array(list(map(utils.fit_vec_to_surface, ft_maxes, m_maxes, full_vecs, rx, ry, rz)))

        # separate components
        x_l = np.array([vecs_fit[:,3]])
        y_l = np.array([vecs_fit[:,4]])
        z_l = np.array([vecs_fit[:,5]])

        return x_l, y_l, z_l

    def plot_twists_on_ls(self, twist_cone, ax):
        """
        Plots friction cone normal twists on limit surface

        Params: 
            twist_cone ([4x3x1] ndarray): Four vectors making the composite twist cone

            ax (matplotlib axis): matplotlib graph axis
            
        Returns:
            None
        """
        ft_maxes = [self.ft_max]*twist_cone.shape[0]
        m_maxes = [self.m_max]*twist_cone.shape[0]
        cone = np.array(list(map(utils.fit_vec_to_surface, ft_maxes, m_maxes, twist_cone)))

        normalize = lambda vec: [vec[0], vec[1], vec[2], vec[3]/np.linalg.norm(vec[3:]), vec[4]/np.linalg.norm(vec[3:]), vec[5]/np.linalg.norm(vec[3:])]
        twist_cone = list(map(normalize, twist_cone))

        X, Y, Z, U, V, W = zip(*twist_cone)
        ax.quiver(X, Y, Z, U, V, W, color = self.colors['blue'], arrow_length_ratio=0.2, linewidths=2.0)

    def plot_vs(self, ax, twist_cone=None, vecs=None):
        """
        Plots twist sphere with cones and/or additional vectors

        Params: 
            twist_cone ([4x3x1] ndarray): Four vectors making the composite twist cone
            vecs ([-1x3x1] ndarray): Vectors to plot 
        Returns:
            None
        """

        # Twist cone
        if twist_cone is not None:
            abc = np.ones(twist_cone.shape[0])
            ft_maxes = [self.ft_max]*twist_cone.shape[0]
            m_maxes = [self.m_max]*twist_cone.shape[0]
            twist_fit = np.array(list(map(utils.fit_vec_to_surface, ft_maxes, m_maxes, twist_cone, abc, abc, abc)))

            X, Y, Z, U, V, W = zip(*twist_fit)
            ax.quiver(X, Y, Z, U, V, W, color = self.colors['blue'], arrow_length_ratio=0.1, linewidths=4.0)

            self.connect_cone(ax, twist_cone, 1.0, 1.0, 1.0, color='blue')

        # Vectors
        for v in self.vs_plot_vectors:
            X, Y, Z, U, V, W, color = zip(*[v])
            ax.quiver(float(X[0]), float(Y[0]), float(Z[0]), float(U[0]), float(V[0]), float(W[0]), color = self.colors[color[0]], arrow_length_ratio=0.2, linewidths=4.0)

        x_l, y_l, z_l = self.calculate_connector(np.array(self.vs_plot_vectors[-3][3:6], dtype=np.double), np.array(self.vs_plot_vectors[-4][3:6], dtype=np.double), 1.0,1.0,1.0)
        ax.plot_wireframe(x_l, y_l, z_l,  rstride=1, cstride=1, color= self.colors['purple'], linewidths=4.0)

        # Unit velocity sphere
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones_like(u), np.cos(v))

        ax.plot_surface(x, y, z,  rstride=5, cstride=5, color= self.colors['gray'])

        # Plot design
        self.style_axis(ax, 'Unit Velocity Sphere', '$v_x$', '$v_y$', '$\omega$', 1.0)

    def plot_all(self, ls, cone=None, twist_cone=None, vecs=None):
        """
        Plots limit surface & twist sphere with cones and/or additional vectors

        Params: 
            ls (LimitSurface): limit surface instance
            cone ([4x3x1] ndarray): Four vectors making the composite wrench cone
            vecs ([-1x3x1] ndarray): Vectors to plot 
        Returns:
            None
        """
        fig = plt.figure(figsize=(20,8))  # Square figure
        ax_ls = fig.add_subplot(121, projection='3d')
        ax_vs = fig.add_subplot(122, projection='3d')

        self.plot_ls(ls, ax_ls, cone)
        self.plot_vs(ax_vs, twist_cone)

        plt.show()

    def set_ls_params(self, ft_max, m_max):
        " Sets ft_max and m_max which is needed throughout for fit_vec_to_surface"
        self.ft_max = ft_max
        self.m_max = m_max

    def style_axis(self, ax, title, x_title, y_title, z_title, limit):
        """
        Styles a given subplot with axis labels limits etc.

        Params: 
            ax (matplotlib axis): matplotlib graph axis

            title (str): title of plot

            x_title (str): x-axis title

            y_title (str): y-axis title

            z_title (str): z-axis title

            limit (double): limit of axes

        Returns:
            None
        """
        if title == 'Unit Velocity Sphere':
            pop_a = mpatches.Patch(color=self.colors['blue'], label='Stable Velocity Twists')
            pop_b = mpatches.Patch(color=self.colors['black'], label="Car's limits")
            pop_c = mpatches.Patch(color=self.colors['purple'], label="Pushing limits")
            pop_c = mpatches.Patch(color=self.colors['pink'], label="STABLE limits")
            ax.legend(handles=[pop_a,pop_b, pop_c], loc='lower left', bbox_to_anchor=(0.0, -0.1))
        else:
            pop_a = mpatches.Patch(color=self.colors['blue'], label='Velocity Twists')
            pop_b = mpatches.Patch(color=self.colors['black'], label='Composite Wrench Cone')
            ax.legend(handles=[pop_a,pop_b], loc='lower left', bbox_to_anchor=(0.0, -0.1))

        ax.set_title(title, fontdict={'fontsize': 30}, pad=50)
        ax.set_xlabel(r'{}'.format(x_title))
        ax.set_ylabel(r'{}'.format(y_title))
        ax.set_zlabel(r'{}'.format(z_title))
        ax.zaxis.set_rotate_label(False) 
        ax.yaxis.set_rotate_label(False) 
        ax.xaxis.set_rotate_label(False) 
        ax.set_xlim(-limit, limit)
        ax.set_ylim(limit, -limit)
        ax.set_zlim(-limit, limit)

    def plot_vec(self, vec, plot_name, normalize=False, color='purple'):
        """
        Plots vector on given plot

        Params: 
            vec ([3x1] ndarray): vector to plot

            plot_name (str): name of plot (ls or vs)

            normalize (bool): whether or not to normalize vector

        Returns:
            None
        """

        if plot_name == 'vs':
            if normalize:
                vec = utils.fit_vec_to_surface(self.ft_max, self.m_max, np.array([0,0,0,vec[0],vec[1],vec[2]], dtype=np.double), 1.0, 1.0, 1.0)
                vec = np.append(vec, color)
                self.vs_plot_vectors.append(vec)
            else:
                self.vs_plot_vectors.append([0,0,0,vec[0], vec[1], vec[2], color])

        if plot_name == 'ls':
            if normalize:
                vec = utils.fit_vec_to_surface(self.ft_max, self.m_max, np.array([0,0,0,vec[0],vec[1],vec[2]], dtype=np.double))
                vec = np.append(vec, color)
                self.ls_plot_vectors.append(vec)
            else:
                self.ls_plot_vectors.append([0,0,0,vec[0], vec[1], vec[2], color])
