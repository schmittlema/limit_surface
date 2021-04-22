# Utility functions
# Schmittle

import os, yaml 

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
