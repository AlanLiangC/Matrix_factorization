import copy
import numpy as np


def rowmul(matrix, index, k):
    """
    Multiplies index row with k
    """
    matrix[index, :] = k * matrix[index, :]
    return matrix


def rowadd(matrix, index1, index2, k):
    """
    Adds the index1 row with index2 row which in turn is multiplied by k
    """
    matrix[index1, :] += k * matrix[index2, :]
    return matrix


def update_mat(mat, row):
    '''
    Put the pivots to the first line
    '''
    mat_this = abs(mat)
    max_index = np.argmax(mat_this[row:, row]) + row
    if max_index == row:
        return mat
    else:
        change = copy.deepcopy(mat[max_index, :])
        mat[max_index, :] = mat[row, :]
        mat[row, :] = change
    return mat


def rref(matrix):
    """
    Returns the reduced row echelon form of a Matrix.

    Examples
    ========

    >>>from factorization.matrixtools import rref
    >>>import numpy as np
    >>> a = np.array([[3,7,4],[2,4,5],[6,2,3]])
    >>> rref(a)

    array([[1., 0., 0.],
       [0., 1., 0.],
       [0., 0., 1.]])

    """
    result_matrix = copy.deepcopy(matrix)
    result_matrix = result_matrix.astype(np.float)
    nrow, ncol = result_matrix.shape
    for i in range(min(nrow, ncol)):
        result_matrix = update_mat(result_matrix, i)
        if np.all(np.round(result_matrix[i, :], decimals=4) == 0):
            continue
        frist_nzero = np.where(np.round(result_matrix[i, :], decimals=4) != 0)[0][0]
        if result_matrix[i, frist_nzero] != 1:
            rowmul(result_matrix, i, 1 / result_matrix[i, frist_nzero])
        rows = [k for k in range(nrow)]
        rows.remove(i)
        for j in rows:
            if (np.round(result_matrix[j, frist_nzero], decimals=4) != 0):
                rowadd(result_matrix, j, i, -result_matrix[j, frist_nzero])
    return result_matrix


def rank(matrix):
    '''
    Return the rank of the matrix

    Examples
    ========

    >>>from factorization.matrixtools import rank
    >>>import numpy as np
    >>> a = np.array([[3,7,4],[2,4,5],[6,2,3]])
    >>> rank(a)

    3
    '''
    result_matrix = copy.deepcopy(matrix)
    rref_state = rref(result_matrix)
    rank_num = 0
    for i in range(rref_state.shape[0]):
        if np.any(np.round(rref_state[i, :], decimals=4) != 0):
            rank_num += 1
    return rank_num


class Mattools():
    '''
    Some basic tools of matrix , contain " rref , rank , columnspace and nullspace"
    '''

    def __init__(self, mat):
        self.mat = copy.deepcopy(mat)
        self.nrow, self.ncol = self.mat.shape
        self.stables = []
        self.pivots = []

    def rref(self):
        '''
        Returns the reduced row echelon form of a Matrix.

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,0,0],[3,0,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.rank()

        array([[ 1.        ,  0.        , -1.33333333],
       [-0.        ,  1.        ,  0.3030303 ],
       [ 0.        ,  0.        ,  0.        ]])

        See also
        ========

        def rref(matrix)
        '''
        return rref(self.mat)

    def rank(self):
        '''
        Return the rank of the matrix

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,0,0],[3,0,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.rank()

        2

        See also
        ========

        def rank(matrix)
        '''
        return rank(self.mat)

    def columnspace(self):
        '''
        columnspace is also 值空间

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,0,0],[3,0,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.columnspace()

        array([[ 0.,  0.],
       [ 3.,  0.],
       [ 4., 11.]])
        '''
        if self.rank() == self.ncol:
            return self.mat
        pivots = []
        result_matrix = copy.deepcopy(self.mat)
        rref_state = rref(result_matrix)
        for i in range(rref_state.shape[0]):
            if np.any(np.round(rref_state[i, :], decimals=4) != 0):
                pivots.append(i)
        columspace_mat = np.zeros([self.nrow, len(pivots)])
        for i in range(len(pivots)):
            pivots_row = rref_state[pivots[i], :]
            for j in range(pivots_row.size):
                if pivots_row[j] == 1:
                    columspace_mat[:, i] = result_matrix[:, j]
                    continue
        return columspace_mat

    def nullspace(self):
        '''
        For nullspace of matrix

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[0,0,0],[3,0,-4],[4,11,-2]])
        >>>test_object = Matrix(a)
        >>>test_object.nullspace()

        array([[ 1.33333333],
       [-0.3030303 ],
       [ 1.        ]])
        '''
        if self.rank() == self.ncol:
            return np.zeros([self.ncol, 1])
        vars = [i for i in range(self.ncol)]
        result_matrix = copy.deepcopy(self.mat)
        rref_state = rref(result_matrix)
        for i in range(rref_state.shape[0]):
            if np.any(np.round(rref_state[i, :], decimals=4) != 0):
                self.pivots.append(i)
        for i in range(len(self.pivots)):
            pivots_row = rref_state[self.pivots[i], :]
            for j in range(pivots_row.size):
                if pivots_row[j] == 1:
                    self.stables.append(j)
                    vars.remove(j)
                    break
        nullspace_mat = np.zeros([self.ncol, self.ncol - len(self.pivots)])
        for i in range(len(self.pivots)):
            pivots_row = rref_state[self.pivots[i], :]
            for k in range(len(vars)):
                nullspace_mat[self.stables[i], k] = -pivots_row[vars[k]]
        for i in range(len(vars)):
            nullspace_mat[vars[i], i] = 1
        return nullspace_mat


if __name__ == "__main__":
    pass
