from factorization.matrixtools import Mattools
import numpy as np
import copy


def gram_schmidt(matrix):
    '''
    QR for Gram_Schmidt_orthogonalization

    Examples
    ========

    >>>from factorization.qr_factorization import gram_schmidt
    >>>import numpy as np
    >>>a = np.array([[0,-20,-14],[3,27,-4],[4,11,-2]])
    >>>gram_schmidt(a)

    array([[ 0.  , -0.8 , -0.6 ],
        [ 0.6 ,  0.48, -0.64],
        [ 0.8 , -0.36,  0.48]]), array([[ 5., 25., -4.],
        [ 0., 25., 10.],
        [ 0.,  0., 10.]])

    '''
    mat = copy.deepcopy(matrix)
    nrow, ncol = mat.shape
    R = np.zeros([ncol, ncol])
    Q = np.zeros(mat.shape)
    try:
        for i in range(ncol):
            a = mat[:, i]
            if i == 0:
                v = np.sqrt(np.sum(a ** 2))
                R[i, i] = v
                q = a / v
                Q[:, i] = q
            else:
                fourier = np.zeros([nrow])
                for j in range(i):
                    R[j, i] = np.dot(Q[:, j].T, a)
                    fourier += R[j, i] * Q[:, j]
                q = a - fourier
                v = np.sqrt(np.sum(q ** 2))
                R[i, i] = v
                q = q / v
                Q[:, i] = q
        return Q, R
    except:
        print("Some errors happend,Please check the matrix to ensure its column independent")


def get_reflector(submatix):
    '''
    For reflector of a matrix (line)
    '''
    mat = copy.deepcopy(submatix)
    col = mat.size
    e = np.zeros(col)
    e[0] = 1
    u = submatix - np.sqrt(np.sum(submatix ** 2)) * e
    if np.all(u == 0):
        return np.eye(col)
    else:
        u = u.reshape(col, 1)
        R = np.eye(col) - 2 * (u @ u.T) / (u.T @ u)
        return R


def householder(matrix):
    '''
    For Householder_reduction using reflectors

    :return P: orthogonal matrix
            R: upper triangular matrix

    Example
    ========

    >>>from factorization.qr_factorization import householder
    >>>import numpy as np
    >>>a = np.array([[0,-20,-14],[3,27,-4],[4,11,-2]])
    >>>householder(a)

    array([[ 0.  ,  0.6 ,  0.8 ],
        [-0.8 ,  0.48, -0.36],
        [-0.6 , -0.64,  0.48]]),
     array([[ 5.00000000e+00,  2.50000000e+01, -4.00000000e+00],
            [ 0.00000000e+00,  2.50000000e+01,  1.00000000e+01],
            [ 0.00000000e+00, -1.33226763e-15,  1.00000000e+01]])

    '''
    mat = copy.deepcopy(matrix)
    nrow, ncol = mat.shape
    P = np.eye(nrow)
    submat = copy.deepcopy(mat)
    if nrow > ncol:
        for i in range(min(nrow, ncol)):
            submat_this = submat[i:, i:]
            R_this = np.eye(nrow)
            if np.all(submat_this[:, 0] == 0):
                continue
            else:
                sub_R_this = get_reflector(submat_this[:, 0])
            R_this[i:, i:] = sub_R_this
            submat = R_this @ submat
            P = R_this @ P
        return P, submat
    else:
        for i in range(min(nrow, ncol) - 1):
            submat_this = submat[i:, i:]
            R_this = np.eye(nrow)
            if np.all(submat_this[:, 0] == 0):
                continue
            else:
                sub_R_this = get_reflector(submat_this[:, 0])
            R_this[i:, i:] = sub_R_this
            submat = R_this @ submat
            P = R_this @ P
        return P, submat


def get_rotation(submatrix):
    '''
    For rotation of a matrix (line)
    '''
    mat = copy.deepcopy(submatrix)
    ncol = mat.size
    P = np.eye(ncol)
    for i in range(1, ncol):
        P_this = np.eye(ncol)
        noml = np.sqrt(mat[0] ** 2 + mat[i] ** 2)
        if noml != 0:
            c = mat[0] / noml
            s = mat[i] / noml
            P_this[0, 0] = c
            P_this[0, i] = s
            P_this[i, 0] = -s
            P_this[i, i] = c
            mat = P_this @ mat
            P = P_this @ P
        else:
            continue
    return P


def givens(matrix):
    '''
    For Givens_reduction using reflectors

    :return P: orthogonal matrix
            R: upper triangular matrix

    Example
    ========

    >>>from factorization.qr_factorization import householder
    >>>import numpy as np
    >>>a = np.array([[0,-20,-14],[3,27,-4],[4,11,-2]])
    >>>householder(a)

    array([[ 0.  ,  0.6 ,  0.8 ],
        [-0.8 ,  0.48, -0.36],
        [-0.6 , -0.64,  0.48]]),
     array([[ 5.00000000e+00,  2.50000000e+01, -4.00000000e+00],
            [ 2.66453526e-16,  2.50000000e+01,  1.00000000e+01],
            [-3.55271368e-16, -3.10862447e-16,  1.00000000e+01]])
    '''
    mat = copy.deepcopy(matrix)
    nrow, ncol = mat.shape
    P = np.eye(nrow)
    submat = copy.deepcopy(mat)
    for i in range(min(nrow, ncol)):
        submat_this = submat[i:, i:]
        P_this = np.eye(nrow)
        sub_P_this = get_rotation(submat_this[:, 0])
        P_this[i:, i:] = sub_P_this
        submat = P_this @ submat
        P = P_this @ P
    return P, submat


class QR(Mattools):
    def __init__(self, mat):
        super().__init__(mat)

    def Gram_Schmidt_orthogonalization(self):
        '''

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,-20,-14],[3,27,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.Gram_Schmidt_orthogonalization()

        array([[ 0.  , -0.8 , -0.6 ],
        [ 0.6 ,  0.48, -0.64],
        [ 0.8 , -0.36,  0.48]]), array([[ 5., 25., -4.],
        [ 0., 25., 10.],
        [ 0.,  0., 10.]])

        See alse
        ========

        def gram_schmidt(matrix)
        '''
        return gram_schmidt(self.mat)

    def Householder_reduction(self):
        '''

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,-20,-14],[3,27,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.Householder_reduction()

        array([[ 0.  ,  0.6 ,  0.8 ],
        [-0.8 ,  0.48, -0.36],
        [-0.6 , -0.64,  0.48]]),
         array([[ 5.00000000e+00,  2.50000000e+01, -4.00000000e+00],
                [ 0.00000000e+00,  2.50000000e+01,  1.00000000e+01],
                [ 0.00000000e+00, -1.33226763e-15,  1.00000000e+01]])

        See alse
        ========

        def gram_schmidt(matrix)
        '''

        return householder(self.mat)

    def Givens_reduction(self):
        '''

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,-20,-14],[3,27,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.Givens_reduction()

        array([[ 0.  ,  0.6 ,  0.8 ],
        [-0.8 ,  0.48, -0.36],
        [-0.6 , -0.64,  0.48]]),
         array([[ 5.00000000e+00,  2.50000000e+01, -4.00000000e+00],
                [ 0.00000000e+00,  2.50000000e+01,  1.00000000e+01],
                [ 0.00000000e+00, -1.33226763e-15,  1.00000000e+01]])

        See alse
        ========

        def gram_schmidt(matrix)
        '''

        return givens(self.mat)


if __name__ == "__main__":
    pass