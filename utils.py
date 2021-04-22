# Utility functions
# Schmittle

import os, yaml, math
import numpy as np

def load_params(filename):
    '''
    Loads params from filename 
    **Parameters**
      - **filename (string)** filename
    **Returns**
      - **dict** paramters
    '''
    src_dir = os.path.dirname(os.path.realpath(__file__))
    with open(src_dir+'/'+filename, 'r') as stream:
        config = yaml.safe_load(stream)
    return config

def print_params(params):
    ''' Prints parameters '''
    print()
    print('BLOCK PARAMETERS')
    print()
    print('Normal Force: {}'.format(params['normal_force']))
    print('COM to Edge Distance: {}'.format(params['dist']))
    print('Radius of Contact Circle Approx.: {}'.format(params['radius']))
    print('Block & Floor Friction Coefficient: {}'.format(params['ls_coefficient']))
    print('Block & Pusher Friction Coefficient: {}'.format(params['coefficient']))
    print()

    print()
    print('CAR PARAMETERS')
    print()
    print('Wheelbase: {}'.format(params['wheelbase']))
    print('Bumper to front axle: {}'.format(params['bumper2front_axle']))
    print('Maximum Steering Angle: {}'.format(params['max_steering_angle']))
    print('Speed: {}'.format(params['speed']))
    print()

def fit_vec_to_surface(ft_max, m_max, vec, a=None, b=None, c=None):
    """
    Adjusts vector to lie on the surface for easy viz

    Params: 
        vec ([6x1] ndarray): Vector
        a,b,c (double): coeff of ellipsoid
    Returns:
        vec ([6x1] ndarray): Vector fit to limit surface 
    """
    v_hat = vec / np.linalg.norm(vec)

    if a is None:
        a = ft_max
    if b is None:
        b = ft_max
    if c is None:
        c = m_max

    # find constant multiple and fit
    # solving for c*v_hat_hat that fits ellipsoid equation
    fit = lambda v_hat: v_hat * math.sqrt((1/(v_hat[3]**2/a**2 + v_hat[4]**2/b**2 + v_hat[5]**2/c**2)))

    return fit(v_hat)



