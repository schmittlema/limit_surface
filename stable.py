# Implementation of STABLE
# Schmittle

import numpy as np

class STABLE:
    def __init__(self):
        """ 
        Constructor
        """
        pass

    def compute_valid_ics():
        """
        Computes stable instaneous centers of rotation for a given object
        Params: 
            None
        Returns:
            None
        """
        #TODO
        pass

    def calculate_noslip_ics(self, comp_fcone, box):
        """
        Calculates no slip ICs using vertical limit
        Params: 
            comp_fcone ([4x3x1] ndarray): Four vectors making the composite wrench cone
            box ([4x2] ndarray): Four points of the box top_left, top_right, bottom_left, bottom_right
        Returns:
            None
        """

        # right edge ICs
        right_vecs = self.calculate_vertical_limit(comp_fcone[3], box[1:3,:])

        # left edge ICs
        left_vecs = self.calculate_vertical_limit(comp_fcone[0], np.array([box[0], box[3]]))

        # left intersection point
        c2 = (left_vecs[0,1] - ((left_vecs[0,0]-right_vecs[1,0])*right_vecs[1,4])/right_vecs[1,3]) \
                /(((left_vecs[0,3]*right_vecs[1,4])/right_vecs[1,3])-left_vecs[0,4])
        c1 = (c2*left_vecs[0,3] + left_vecs[0,0] - right_vecs[1,0])/ right_vecs[1,3]
        print(c1, c2)
        c1 = (c2*left_vecs[0,4] + left_vecs[0,1] - right_vecs[1,1])/ right_vecs[1,4]
        print(c1, c2)
        # TODO fix my algrebra
        #print(right_vecs[1,3:], left_vecs[0,3:])
        #print(c1*right_vecs[1,3:], c2*left_vecs[0,3:])
        #print(right_vecs[1,:3], left_vecs[0,:3])
        #print(c1*right_vecs[1,3:]+right_vecs[1,:3], c2*left_vecs[0,3:]+left_vecs[0,:3])

    def calculate_vertical_limit(self, vector, edge_pts):
        """
        Calculates vertical strip limit for a given friction cone and object.
        Params: 
            vec ([6x1] ndarray): Vector
            edge_pts ([2x2] ndarray): x,y points on edge of object 
        Returns:
            edge_vecs ([2x6] ndarray): vectors defining edges of vectical limit 
        """
        # get perpendicular vector in x,y plane
        perp_vec = np.zeros(vector.shape)
        perp_vec[:] = vector
        perp_vec[-2] = -1*vector[-3]
        perp_vec[-3] = vector[-2]

        # position vecs on edges of box
        edge_vec_1 = perp_vec + np.array([edge_pts[0,0], edge_pts[0,1], 0, 0, 0, 0])
        edge_vec_2 = perp_vec + np.array([edge_pts[1,0], edge_pts[1,1], 0, 0, 0, 0])
        return np.array([edge_vec_1, edge_vec_2])

    def calculate_bisector(self):
        """
        Calculates bisector for a contact point and center of friction 
        Params: 
            None
        Returns:
            None
        """
        #TODO
        pass

    def calculate_tipline(self):
        """
        Calculates tipline from circumscribed disk and contact point
        Params: 
            None
        Returns:
            None
        """
        #TODO
        pass
