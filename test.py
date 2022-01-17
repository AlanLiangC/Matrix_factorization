import numpy as np
from factorization import Matrix
'''
矩阵分析与应用大作业-梁奥-202128014728021

用于项目的测试与验证

统一测试矩阵：
a = np.array([[1,2,-3,4],
              [4,8,12,-8],
              [2,3,2,1],
              [-3,-1,1,-4]])
'''

# 测试程序
def test():
    '''
    测试函数，本函数将测试项目包含的所有功能
    '''
    a = np.array([[1,2,-3,4],
              [4,8,12,-8],
              [2,3,2,1],
              [-3,-1,1,-4]])
    # 实例化对象
    test_object = Matrix(a)
    
    # 行阶梯表示
    print("矩阵的行阶梯表示为：\n",test_object.rref(),"\n")

    # 矩阵的秩
    print("矩阵的秩为：\n",test_object.rank(),"\n")

    # 零空间
    print("矩阵的零空间为：\n",test_object.nullspace(),"\n")

    # 值空间
    print("矩阵的值空间为：\n",test_object.columnspace(),"\n")

    # PLU分解
    print("矩阵的PLU分解为：\n",test_object.PLU_factorization(),"\n")

    # 矩阵的逆
    print("矩阵的逆为：\n",test_object.get_inversion(),"\n")

    # Gram-Schmidt正交化
    print("矩阵的Gram-Schmidt正交化为：\n",test_object.Gram_Schmidt_orthogonalization(),"\n")

    # Householder reduction
    print("矩阵的Householder约简为：\n",test_object.Householder_reduction(),"\n")

    # Givens reduction
    print("矩阵的Givens约简为：\n",test_object.Givens_reduction(),"\n")

    # URV分解
    print("矩阵的URV分解为：\n",test_object.URV_factorization(),"\n")

    # 矩阵的行列式
    print("矩阵的行列式为：\n",test_object.det(),"\n")

    # 线性系统求解，另b=[1,1,1]
    b = np.array([1,1,1,1])
    print("线性系统的解为：\n",test_object.linear_equation(b),"\n")

'''
下面是项目功能的正确性验证，分别写成不同的函数

调用不同的函数以实现不同功能的验证
'''

# 矩阵的秩的验证
def rank_verification():
    '''
    随机生成5000个5x4的0/1矩阵求秩，
    用项目中的函数与numpy库中功能函数的做对比
    '''
    error = 0
    for i in range(5000):
    #     随机生成0/1矩阵
        mat = np.random.randint(0,2,(5,4))
    #     实例化对象
        test = Matrix(mat)
    #     调用rank()函数
        rank1 = test.rank()
    #     调用numpy中的函数
        rank2 = np.linalg.matrix_rank(mat)
        if rank1 != rank2:
            print("出现错误，错误的矩阵是\n{}".format(mat))
            error += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# PLU分解验证
def PLU_verification():
    '''
    随机生成5000个4x4满秩矩阵进行PLU分解，
    并计算是否PA=LU，是则认为程序正确
    '''
    num = 0
    error = 0
    while(num < 5000):
    #     随机生成矩阵
        mat = np.random.randint(0,11,(4,4))
    #     实例化对象
        test = Matrix(mat)
        if test.rank() != 4:
            continue
    #     PLU_factorization()函数
        P,L,U = test.PLU_factorization()
    #     设置四位有效数字的精度
        if np.any(np.round(P@mat,decimals=4) != np.round(L@U,decimals=4)):
            print("P@mat = {}\n,L@U = {}".format(P@mat,L@U))
            print("Error!!!")
            error += 1
        num += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# 矩阵的逆验证
def inv_verification():
    '''
    随机生成5000个4x4满秩矩阵进行求逆，
    并计算是否AA−1=I，是则认为程序正确
    '''
    num = 0
    error = 0
    while(num < 5000):
    #     随机生成矩阵
        mat = np.random.randint(0,11,(4,4))
    #     实例化对象
        test = Matrix(mat)
        if test.rank() != 4:
            continue
    #     get_inversion()函数
        inv = test.get_inversion()
    #     设置三位有效数字的精度
        if np.any(np.round(mat@inv,decimals=3) != np.eye(4)):
            print("mat = {}\n,int = {}\n{}".format(mat,inv,mat@inv))
            print("Error!!!")
            error += 1
        num += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# Gram-Schmidt验证
def QR_verification():
    '''
    随机生成5000个5x3矩阵进行QR分解，
    并计算是否A=QR，是则认为程序正确
    '''
    num = 0
    error = 0
    while(num < 5000):
    #     随机生成矩阵
        mat = np.random.randint(0,11,(5,3))
    #     实例化对象
        test = Matrix(mat)
        if test.rank() != 3:
            continue
    #     Gram_Schmidt_orthogonalization()函数
        Q,R = test.Gram_Schmidt_orthogonalization()
    #     设置四位有效数字的精度
        if np.any(np.round(Q@R,decimals=4) != mat):
            print("Q@R = {}\n,mat = {}".format(Q@R,mat))
            print("Error!!!")
            error += 1
        num += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# Householder reduction验证
def Householder_verification():
    num = 0
    error = 0
    while(num < 5000):
    #     随机生成矩阵
        mat = np.random.randint(0,11,(5,3))
    #     实例化对象
        test = Matrix(mat)
    #     Householder_reduction()函数
        P,R = test.Householder_reduction()
    #     设置四位有效数字的精度
        if np.any(np.round(P@P.T,decimals=4) != np.eye(5)):
            print("P@P.T = {}\n,mat = {}".format(P@P.T,mat))
            print("Error!!!")
            error += 1
        num += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# Givens reduction验证
def Givens_verification():
    '''
    随机生成5000个5x3矩阵进行Givens约简，
    并计算是否PPT=I，是则认为程序正确
    '''
    num = 0
    error = 0
    while(num < 5000):
    #     随机生成矩阵
        mat = np.random.randint(0,11,(5,3))
    #     实例化对象
        test = Matrix(mat)
    #     Givens_reduction()函数
        P,R = test.Givens_reduction()
    #     设置四位有效数字的精度
        if np.any(np.round(P@P.T,decimals=4) != np.eye(5)):
            print("P@P.T = {}\n,mat = {}".format(P@P.T,mat))
            print("Error!!!")
            error += 1
        num += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# URV 分解验证
def URV_verification():
    '''
    随机生成5000个5x3矩阵进行URV分解，
    并计算是否UUT=I,RRT=I,A=URV，是则认为程序正确
    '''
    num = 0
    error = 0
    while(num < 5000):
    #     随机生成矩阵
        mat = np.random.randint(0,11,(5,3))
    #     实例化对象
        test = Matrix(mat)
    #     Givens_reduction()函数
        U,R,V = test.URV_factorization()
    #     设置四位有效数字的精度
        if np.any(np.round(U@U.T,decimals=3) != np.eye(5)) and np.any(np.round(V@V.T,decimals=3) != np.eye(3)) and np.any(np.round(U@R@V.T,decimals=3) != mat) :
            print("mat = {}".format(mat))
            print("Error!!!")
            error += 1
        num += 1
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))

# 矩阵的行列式验证
def det_verification():
    '''
    随机生成5000个5x5的矩阵求行列式，用项目中的函数与numpy库中功能函数的做对比，
    展示前十个矩阵的计算结果
    '''
    error = 0
    for i in range(5000):
    #     随机生成矩阵
        mat = np.random.randint(-10,10,(5,5))
    #     实例化对象
        test = Matrix(mat)
    #     调用det()函数
        det1 = test.det()
    #     调用numpy中的函数
        det2 = np.linalg.det(mat)
        if np.round(det1) != np.round(det2):
            print("出现错误，错误的矩阵是\n{}".format(mat))
            error += 1
        if i < 10:
            print("===================={}===================".format(i))
            print("本项目函数的计算结果：{}".format(det1))
            print(" numpy 库的计算结果：{}".format(det2))
    if error == 0:
        print("未出现错误，测试样本{}个".format(5000))


if __name__ == "__main__":
    # 功能测试
    test()
    # 准确性验证
    rank_verification()
    inv_verification()
    QR_verification()
    Householder_verification()
    Givens_verification()
    URV_verification()
    det_verification()