import copy
import numpy as np
from factorization.matrixtools import rank, rref
from factorization.lu_factorization import PLU
from factorization.qr_factorization import QR
from factorization.urv_factorization import URV


def det(matrix):
    '''
    Get the det of a square matrix using the method of PLU factorization

    Example
    ========

    >>>from factorization.matrix import det
    >>>import numpy as np
    >>>a = np.array([[1,0,0],[3,0,-4],[4,11,-2]])
    >>>det(a)

    44.000001311302185
    '''
    mat = copy.deepcopy(matrix)
    det_object = PLU(mat)
    _, _, U = det_object.PLU_factorization()
    if det_object.ncol != det_object.nrow:
        error = "The Matrix is not a square matrix , There no det"
        return error
    if det_object.rank() < det_object.ncol:
        mat_det = 0
        return mat_det
    else:
        mat_det = 1
        for i in range(det_object.ncol):
            mat_det *= U[i, i]
        return (-1) ** det_object.transtime * mat_det


def linear_equation(matrix, b):
    '''

    For any liner system's result.
    - When its a nonsingular matrix ,then use PLU factorization
    - Else, use Gauss-Jordan to get its rref style
        - When b = 0,return the null space of the matrix
        - Else, find the general result and augmented matrix last column.

    Example
    ========

    >>>from factorization.matrix import linear_equation
    >>>import numpy as np
    >>>a = np.array([[1,1,2,2,1],[2,2,4,4,3],[2,2,4,4,2],[3,5,8,6,5]])
    >>>b = np.array([1,1,2,3])
    >>>linear_equation(a,b)

    [array([ 1.,  1.,  0.,  0., -1.]), array([[-1., -2.],
        [-1., -0.],
        [ 1.,  0.],
        [ 0.,  1.],
        [ 0.,  0.]])]
    '''
    mat = copy.deepcopy(matrix)
    linear_object = PLU(mat)
    if np.all(b == 0):
        return linear_object.nullspace()
    else:
        augmented_mat = np.c_[mat, b]
        augmented_mat_rank = rank(augmented_mat)
        mat_rank = linear_object.rank()
        if augmented_mat_rank == mat_rank:
            if mat_rank == linear_object.ncol == linear_object.nrow:
                return linear_object.linear_equation(b)
            else:
                if linear_object.nrow >= linear_object.ncol:
                    const = rref(augmented_mat)[:, -1]
                    general = linear_object.nullspace()
                    return const, general
                else:
                    general = linear_object.nullspace()
                    pivots = linear_object.stables
                    const = rref(augmented_mat)[:, -1]
                    anser = np.zeros(linear_object.ncol)
                    for i in range(len(pivots)):
                        anser[pivots[i]] = const[i]
                    return [anser, general]
        else:
            print("无解")
            return 0


class Matrix(PLU, QR, URV):
    def __init__(self, mat):
        super().__init__(mat)

    def det(self):
        '''
        Get the det of a square matrix using the method of PLU factorization

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[1,0,0],[3,0,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.det()

        44.000001311302185

        See also
        ========

        def det(matrix)
        '''
        mat = copy.deepcopy(self.mat)
        return det(mat)

    def linear_equation(self, b):
        '''

        For any liner system's result.
        - When its a nonsingular matrix ,then use PLU factorization
        - Else, use Gauss-Jordan to get its rref style
            - When b = 0,return the null space of the matrix
            - Else, find the general result and augmented matrix last column.

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[1,1,2,2,1],[2,2,4,4,3],[2,2,4,4,2],[3,5,8,6,5]])
        >>>b = np.array([1,1,2,3])
        >>>test_object = Matrix(a)
        >>>test_object.linear_equation(b)

        [array([ 1.,  1.,  0.,  0., -1.]), array([[-1., -2.],
            [-1., -0.],
            [ 1.,  0.],
            [ 0.,  1.],
            [ 0.,  0.]])]

        See also
        ========

        def linear_equation(matrix, b)
        '''

        return linear_equation(self.mat, b)


if __name__ == "__main__":
    pass
