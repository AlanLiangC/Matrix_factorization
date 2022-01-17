from factorization.matrixtools import Mattools
from factorization.qr_factorization import householder
import numpy as np
import copy


class URV(Mattools):
    def __init__(self, mat):
        super().__init__(mat)

    def URV_factorization(self):
        '''
        For the URV factorization of any matrix

        Example
        ========
        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[4,8,12,-8],[2,3,2,1],[-3,-1,1,-4]])
        >>>test_object = Matrix(a)
        >>>test_object.URV_factorization()

        array([[ 0.97372899,  0.20462466,  0.02942979, -0.09546984],
        [-0.18257419,  0.9636314 , -0.03927159,  0.19115159],
        [ 0.06085806, -0.06428777, -0.97949468,  0.18097945],
        [ 0.12171612, -0.15940626,  0.19540167,  0.95999636]]),
         array([[ 1.35918623e+01, -2.70100613e-16,  9.30409510e-16,
                  7.99440812e-16],
                [ 1.67486275e+01,  7.45697467e+00, -1.70048074e-17,
                 -1.09712297e-16],
                [ 6.36668910e+00,  4.29004562e+00,  2.53348854e+00,
                 -7.21644966e-16],
                [-3.61921181e-01, -6.10522177e-01, -1.50162168e-02,
                 -5.50147470e-01]]),
         array([[ 0.30223373,  0.50417063,  0.70148075, -0.4029783 ],
                [-0.67882758, -0.26245638,  0.09059456, -0.67978113],
                [ 0.3899655 , -0.82256059,  0.41357157, -0.01671969],
                [ 0.54385181, -0.01796399, -0.57330266, -0.61255739]])

        See alse
        ========

        URV(Mattools)
        '''
        result_mat = copy.deepcopy(self.mat)
        P, R1 = householder(result_mat)
        r = self.rank()
        B = R1[:r, :]
        Q, R2 = householder(B.T)
        T = R2[:r, :r]
        R = np.zeros(result_mat.shape)
        R[:r, :r] = T.T
        return P.T, R, Q


if __name__ == "__main__":
    pass
