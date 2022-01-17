from factorization.matrixtools import Mattools
import numpy as np
import copy


# 定义类 功能均在类内实现
class PLU(Mattools):
    def __init__(self,mat):
        super().__init__(mat)
        self.mat = copy.deepcopy(mat)  # 深层拷贝，防止初始矩阵内存被重复占用，下同
        self.mat = self.mat.astype(np.float)
        self.mat_stable = copy.deepcopy(mat)
        self.mat_R = mat.shape[0]
        self.mat_C = mat.shape[1]
        self.transtime = 0
        self.mat_P = np.array([i for i in range(self.mat_R)])
        self.mat_P_result = np.zeros([mat.shape[0], mat.shape[1]], dtype=np.float32)
        self.mat_L = np.eye(mat.shape[0], dtype=np.float32)
        self.mat_U = np.zeros([mat.shape[0], mat.shape[1]], dtype=np.float32)
        self.mat_inversion = np.zeros([mat.shape[0], mat.shape[1]], dtype=np.float32)

    # 更新矩阵，始终将该次计算的最大主元移到最上方
    def update_mat(self, row):
        mat_this = abs(self.mat)
        max_index = np.argmax(mat_this[row:, row]) + row
        if row == max_index:
            return self.mat
        else:
            self.transtime += 1
            change = copy.deepcopy(self.mat[max_index, :])
            self.mat[max_index, :] = self.mat[row, :]
            self.mat[row, :] = change
            change_P = copy.deepcopy(self.mat_P[max_index])
            self.mat_P[max_index] = self.mat_P[row]
            self.mat_P[row] = change_P
        return self.mat

    # 分解矩阵，当detail=Trur时输出分解过程
    def PLU_factorization(self, P = True,detail = False):
        '''
        For PLU_factorization of a matrix,when run the program,you should choose a matrix which is
        nonsingular, only then you can see result
        :parm P:bool,when Ture,PLU factorization,when False,LU factorization
        :parm detail: bool,when True,output the process of running

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[1,2,-3,4],[4,8,12,-8],[2,3,2,1],[-3,-1,1,-4]])
        >>>test_object = Matrix(a)
        >>>test_object.PLU_factorization()

        array([[0., 1., 0., 0.],
        [0., 0., 0., 1.],
        [1., 0., 0., 0.],
        [0., 0., 1., 0.]], dtype=float32),
         array([[ 1.        ,  0.        ,  0.        ,  0.        ],
                [-0.75      ,  1.        ,  0.        ,  0.        ],
                [ 0.25      ,  0.        ,  1.        ,  0.        ],
                [ 0.5       , -0.2       ,  0.33333334,  1.        ]],
               dtype=float32),
         array([[  4.,   8.,  12.,  -8.],
                [  0.,   5.,  10., -10.],
                [  0.,   0.,  -6.,   6.],
                [  0.,   0.,   0.,   1.]], dtype=float32)

        See also
        ========

        class PLU(Mattools)
        '''
        if self.ncol != self.nrow and self.rank() != self.nrow:
            error = "This matrix con't PLU factorization! \n Please check its number of rows and colunms to ensure that they are equal! \n More you should check its rank to ensure it's nonsingular!"
            return error
        for i in range(self.mat_R - 1):
            if P:
                self.mat = self.update_mat(i)
                if detail:
                    print("第{}次更新".format(i), "\n", self.mat, "\n", self.mat_P)
            for j in range(self.mat_R - i - 1):
                if self.mat[i + j + 1, i] == 0:
                    self.mat = self.mat
                else:
                    ratio = self.mat[i + j + 1, i] / self.mat[i, i]
                    self.mat[i + j + 1, i + 1:] = self.mat[i + j + 1, i + 1:] - ratio * self.mat[i, i + 1:]
                    self.mat[i + j + 1, i] = ratio
            if detail:
                print("第{}次计算".format(i), "\n", self.mat)
        
        # 拆分矩阵,得到结果
        for i in range(self.mat_R):
            self.mat_P_result[i, self.mat_P[i]] = 1
            for j in range(self.mat_C):
                if i > j:
                    self.mat_L[i, j] = self.mat[i, j]
                else:
                    self.mat_U[i, j] = self.mat[i, j]
        if detail:
            print("P:", "\n", self.mat_P_result)
            print("L:", "\n", self.mat_L)
            print("U:", "\n", self.mat_U)
        return self.mat_P_result, self.mat_L, self.mat_U

    # 进行矩阵求逆，当detail=True时输出运行过程，以便进行bug定位
    # 求逆过程为逐行求逆 类比于Ax=b，x为A逆的每一列，b为单位阵I的每一列
    def get_inversion(self, detail=False):
        '''
        For the inversion of a nonsingular matrix ,when run the program,you should choose a matrix which is
        nonsingular, only then you can see result

        :parm detail: bool,when True,output the process of running

        Example
        ========

        >>>from factorization import Matrix
        >>>import numpy as np
        >>>a = np.array([[1,2,-3,4],[4,8,12,-8],[2,3,2,1],[-3,-1,1,-4]])
        >>>test_object = Matrix(a)
        >>>test_object.get_inversion()

        array([[ 0.16666669,  0.25833336, -1.        , -0.6       ],
       [ 0.3333333 ,  0.06666665,  0.        ,  0.2       ],
       [-0.5       , -0.22499998,  1.        ,  0.2       ],
       [-0.33333334, -0.26666665,  1.        ,  0.2       ]],
        dtype=float32)

        See also
        ========

        class PLU(Mattools)
        '''
        _,_,_ = self.PLU_factorization()
        B = np.eye(self.mat_R)
        for i in range(self.mat_R):
            sub_B = np.dot(self.mat_P_result, B[:, i])
            y = np.zeros([self.mat_R, 1], dtype=np.float32)
            y[0] = sub_B[0]
            for m in range(self.mat_R - 1):
                sum_y = 0
                for n in range(m + 1):
                    sum_y += self.mat_L[m + 1, n] * y[n]
                y[m + 1] = sub_B[m + 1] - sum_y

            if detail:
                print("此时y为：", y, "\n")
                print("此时的乘积Ly为", "\n", np.dot(self.mat_L, y))

            x = np.zeros([self.mat_R, 1], dtype=np.float32)
            x[self.mat_R - 1] = y[self.mat_R - 1] / self.mat_U[self.mat_R - 1, self.mat_R - 1]
            for p in range(self.mat_R - 1):
                sum_x = 0
                for q in range(p + 1):
                    sum_x += self.mat_U[-p - 2, -q - 1] * x[-q - 1]
                x[-p - 2] = (y[-p - 2] - sum_x) / self.mat_U[-p - 2, -p - 2]

            if detail:
                print("此时x为：", x, "\n")
                print("此时的乘积Ux为", "\n", np.dot(self.mat_U, x))
                print("此时的乘积Ax为", "\n", np.dot(self.mat_stable, x))

            self.mat_inversion[:, i] = x.reshape(self.mat_R)
        return self.mat_inversion

    def linear_equation(self,b):
        '''
        For liner system's result.if the matrix is nonsingular, we use the method of PLU factorization
        to solve More detail ,see Matrix.py

        :parm b: the output of the linear system

        See also
        ========

        Matrix.py -> def linear_equation(matrix, b)
        '''
        self.PLU_factorization()
        sub_B = np.dot(self.mat_P_result, b)
        y = np.zeros([self.mat_R, 1], dtype=np.float32)
        y[0] = sub_B[0]
        for m in range(self.mat_R - 1):
            sum_y = 0
            for n in range(m + 1):
                sum_y += self.mat_L[m + 1, n] * y[n]
            y[m + 1] = sub_B[m + 1] - sum_y

        x = np.zeros([self.mat_R, 1], dtype=np.float32)
        x[self.mat_R - 1] = y[self.mat_R - 1] / self.mat_U[self.mat_R - 1, self.mat_R - 1]
        for p in range(self.mat_R - 1):
            sum_x = 0
            for q in range(p + 1):
                sum_x += self.mat_U[-p - 2, -q - 1] * x[-q - 1]
            x[-p - 2] = (y[-p - 2] - sum_x) / self.mat_U[-p - 2, -p - 2]
        self.mat_inversion = x.reshape(self.mat_R)
        return self.mat_inversion



if __name__ == "__main__":
    pass
