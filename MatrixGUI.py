# -*- coding: utf-8 -*-
from factorization import Matrix
import numpy as np
import copy
from PyQt5 import QtCore, QtGui, QtWidgets
import sys


class Ui_LiangAo_Matrix_test(object):
    def setupUi(self, LiangAo_Matrix_test):
        LiangAo_Matrix_test.setObjectName("LiangAo_Matrix_test")
        LiangAo_Matrix_test.resize(900, 700)
        self.Matrix_input = QtWidgets.QGroupBox(LiangAo_Matrix_test)
        self.Matrix_input.setGeometry(QtCore.QRect(10, 10, 261, 161))
        self.Matrix_input.setObjectName("Matrix_input")
        self.Mat_Row_Edit = QtWidgets.QLineEdit(self.Matrix_input)
        self.Mat_Row_Edit.setGeometry(QtCore.QRect(80, 20, 171, 20))
        self.Mat_Row_Edit.setObjectName("Mat_Row_Edit")
        self.Row = QtWidgets.QLabel(self.Matrix_input)
        self.Row.setGeometry(QtCore.QRect(30, 20, 41, 21))
        self.Row.setObjectName("Row")
        self.Col = QtWidgets.QLabel(self.Matrix_input)
        self.Col.setGeometry(QtCore.QRect(30, 50, 41, 21))
        self.Col.setObjectName("Col")
        self.Mat_Col_edit = QtWidgets.QLineEdit(self.Matrix_input)
        self.Mat_Col_edit.setGeometry(QtCore.QRect(80, 50, 171, 20))
        self.Mat_Col_edit.setObjectName("Mat_Col_edit")
        self.Mat_edit = QtWidgets.QTextEdit(self.Matrix_input)
        self.Mat_edit.setGeometry(QtCore.QRect(80, 80, 171, 71))
        self.Mat_edit.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.Mat_edit.setObjectName("Mat_edit")
        self.pushButton = QtWidgets.QPushButton(self.Matrix_input)
        self.pushButton.setGeometry(QtCore.QRect(20, 80, 51, 71))
        self.pushButton.setObjectName("pushButton")
        self.MatrixTools = QtWidgets.QGroupBox(LiangAo_Matrix_test)
        self.MatrixTools.setGeometry(QtCore.QRect(10, 180, 261, 171))
        self.MatrixTools.setObjectName("MatrixTools")
        self.pushButton_2 = QtWidgets.QPushButton(self.MatrixTools)
        self.pushButton_2.setGeometry(QtCore.QRect(20, 20, 111, 41))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_3 = QtWidgets.QPushButton(self.MatrixTools)
        self.pushButton_3.setGeometry(QtCore.QRect(140, 20, 111, 41))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_4 = QtWidgets.QPushButton(self.MatrixTools)
        self.pushButton_4.setGeometry(QtCore.QRect(20, 70, 111, 41))
        self.pushButton_4.setObjectName("pushButton_4")
        self.pushButton_5 = QtWidgets.QPushButton(self.MatrixTools)
        self.pushButton_5.setGeometry(QtCore.QRect(140, 70, 111, 41))
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_6 = QtWidgets.QPushButton(self.MatrixTools)
        self.pushButton_6.setGeometry(QtCore.QRect(20, 120, 111, 41))
        self.pushButton_6.setObjectName("pushButton_6")
        self.pushButton_7 = QtWidgets.QPushButton(self.MatrixTools)
        self.pushButton_7.setGeometry(QtCore.QRect(140, 120, 111, 41))
        self.pushButton_7.setObjectName("pushButton_7")
        self.groupBox = QtWidgets.QGroupBox(LiangAo_Matrix_test)
        self.groupBox.setGeometry(QtCore.QRect(10, 360, 261, 171))
        self.groupBox.setObjectName("groupBox")
        self.pushButton_13 = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_13.setGeometry(QtCore.QRect(20, 70, 111, 41))
        self.pushButton_13.setObjectName("pushButton_13")
        self.pushButton_8 = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_8.setGeometry(QtCore.QRect(140, 70, 111, 41))
        self.pushButton_8.setObjectName("pushButton_8")
        self.pushButton_9 = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_9.setGeometry(QtCore.QRect(20, 120, 231, 41))
        self.pushButton_9.setObjectName("pushButton_9")
        self.pushButton_11 = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_11.setGeometry(QtCore.QRect(20, 20, 111, 41))
        self.pushButton_11.setObjectName("pushButton_11")
        self.pushButton_12 = QtWidgets.QPushButton(self.groupBox)
        self.pushButton_12.setGeometry(QtCore.QRect(140, 20, 111, 41))
        self.pushButton_12.setObjectName("pushButton_12")
        self.groupBox_2 = QtWidgets.QGroupBox(LiangAo_Matrix_test)
        self.groupBox_2.setGeometry(QtCore.QRect(280, 10, 601, 631))
        self.groupBox_2.setObjectName("groupBox_2")
        self.show_window = QtWidgets.QTextEdit(self.groupBox_2)
        self.show_window.setGeometry(QtCore.QRect(10, 20, 581, 551))
        self.show_window.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.show_window.setObjectName("show_window")
        self.pushButton_10 = QtWidgets.QPushButton(self.groupBox_2)
        self.pushButton_10.setGeometry(QtCore.QRect(10, 580, 581, 41))
        self.pushButton_10.setObjectName("pushButton_10")
        self.Matrix_input_2 = QtWidgets.QGroupBox(LiangAo_Matrix_test)
        self.Matrix_input_2.setGeometry(QtCore.QRect(10, 540, 261, 101))
        self.Matrix_input_2.setObjectName("Matrix_input_2")
        self.pushButton_14 = QtWidgets.QPushButton(self.Matrix_input_2)
        self.pushButton_14.setGeometry(QtCore.QRect(20, 50, 231, 41))
        self.pushButton_14.setObjectName("pushButton_14")
        self.b_edit = QtWidgets.QLineEdit(self.Matrix_input_2)
        self.b_edit.setGeometry(QtCore.QRect(80, 20, 171, 20))
        self.b_edit.setObjectName("b_edit")
        self.Col_3 = QtWidgets.QLabel(self.Matrix_input_2)
        self.Col_3.setGeometry(QtCore.QRect(40, 20, 31, 21))
        self.Col_3.setObjectName("Col_3")
        self.textBrowser = QtWidgets.QTextBrowser(LiangAo_Matrix_test)
        self.textBrowser.setGeometry(QtCore.QRect(10, 640, 861, 51))
        self.textBrowser.setObjectName("textBrowser")

        self.retranslateUi(LiangAo_Matrix_test)
        self.pushButton.clicked.connect(LiangAo_Matrix_test.show_matrix)
        self.pushButton_2.clicked.connect(LiangAo_Matrix_test.show_row_echelon)
        self.pushButton_3.clicked.connect(LiangAo_Matrix_test.show_rref)
        self.pushButton_4.clicked.connect(LiangAo_Matrix_test.show_rank)
        self.pushButton_5.clicked.connect(LiangAo_Matrix_test.show_det)
        self.pushButton_6.clicked.connect(LiangAo_Matrix_test.show_nullspace)
        self.pushButton_7.clicked.connect(LiangAo_Matrix_test.show_columnspace)
        self.pushButton_11.clicked.connect(LiangAo_Matrix_test.show_PLU)
        self.pushButton_12.clicked.connect(LiangAo_Matrix_test.show_gramschmidt)
        self.pushButton_13.clicked.connect(LiangAo_Matrix_test.show_householder)
        self.pushButton_8.clicked.connect(LiangAo_Matrix_test.show_givens)
        self.pushButton_9.clicked.connect(LiangAo_Matrix_test.show_urv)
        self.pushButton_14.clicked.connect(LiangAo_Matrix_test.show_A)
        self.pushButton_10.clicked.connect(self.show_window.clear)
        QtCore.QMetaObject.connectSlotsByName(LiangAo_Matrix_test)

    def retranslateUi(self, LiangAo_Matrix_test):
        _translate = QtCore.QCoreApplication.translate
        LiangAo_Matrix_test.setWindowTitle(_translate("LiangAo_Matrix_test", "Ui_LiangAo_Matrix_test"))
        self.Matrix_input.setTitle(_translate("LiangAo_Matrix_test", "输入矩阵"))
        self.Row.setText(_translate("LiangAo_Matrix_test", "行数："))
        self.Col.setText(_translate("LiangAo_Matrix_test", "列数："))
        self.pushButton.setText(_translate("LiangAo_Matrix_test", "确认"))
        self.MatrixTools.setTitle(_translate("LiangAo_Matrix_test", "矩阵属性"))
        self.pushButton_2.setText(_translate("LiangAo_Matrix_test", "高斯约简"))
        self.pushButton_3.setText(_translate("LiangAo_Matrix_test", "行阶梯"))
        self.pushButton_4.setText(_translate("LiangAo_Matrix_test", "秩"))
        self.pushButton_5.setText(_translate("LiangAo_Matrix_test", "行列式"))
        self.pushButton_6.setText(_translate("LiangAo_Matrix_test", "零空间"))
        self.pushButton_7.setText(_translate("LiangAo_Matrix_test", "值空间"))
        self.groupBox.setTitle(_translate("LiangAo_Matrix_test", "矩阵分解"))
        self.pushButton_13.setText(_translate("LiangAo_Matrix_test", "Householder"))
        self.pushButton_8.setText(_translate("LiangAo_Matrix_test", "Givens"))
        self.pushButton_9.setText(_translate("LiangAo_Matrix_test", "URV"))
        self.pushButton_11.setText(_translate("LiangAo_Matrix_test", "PLU分解"))
        self.pushButton_12.setText(_translate("LiangAo_Matrix_test", "Gram_Schmidt"))
        self.groupBox_2.setTitle(_translate("LiangAo_Matrix_test", "显示窗口"))
        self.pushButton_10.setText(_translate("LiangAo_Matrix_test", "清空窗口"))
        self.Matrix_input_2.setTitle(_translate("LiangAo_Matrix_test", "线性系统求解"))
        self.pushButton_14.setText(_translate("LiangAo_Matrix_test", "确认"))
        self.Col_3.setText(_translate("LiangAo_Matrix_test", "b ："))
        self.textBrowser.setHtml(_translate("LiangAo_Matrix_test",
                                            "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
                                            "<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
                                            "p, li { white-space: pre-wrap; }\n"
                                            "</style></head><body style=\" font-family:\'SimSun\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
                                            "<p align=\"center\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:28pt; font-weight:600; color:#0a2bff;\">矩阵分解及测试界面</span></p></body></html>"))


class MyPyQT_Form(QtWidgets.QWidget, Ui_LiangAo_Matrix_test):
    '''
    A GUI for the project

    Example
    ========
    >>>python MatrixGUI.py

    '''
    def __init__(self):
        super(MyPyQT_Form, self).__init__()
        self.setupUi(self)
        self.mat = np.ndarray([])

    def show_matrix(self):
        nrow = int(self.Mat_Row_Edit.text())
        ncol = int(self.Mat_Col_edit.text())
        mat_from_str = np.fromstring(self.Mat_edit.toPlainText(), dtype=np.float, sep=",")
        mat = mat_from_str.reshape(nrow, ncol)
        self.mat = copy.deepcopy(mat)
        self.show_window.append("已添加 {} 行 {} 列矩阵\n{}".format(nrow, ncol, np.str(mat)))

    def show_row_echelon(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        self.show_window.append("矩阵的高斯约简为：")
        self.show_window.append("{}".format(test.rref()))

    def show_rref(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        self.show_window.append("矩阵的行阶梯表示为：")
        self.show_window.append("{}".format(test.rref()))

    def show_rank(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        self.show_window.append("矩阵的秩为：{}".format(test.rank()))

    def show_det(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        try:
            self.show_window.append("矩阵的行列式为：{}".format(test.det()))
        except:
            self.show_window.append("The Matrix is not a square matrix , There no det")

    def show_nullspace(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        self.show_window.append("矩阵的零空间为：")
        self.show_window.append("{}".format(test.nullspace()))

    def show_columnspace(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        self.show_window.append("矩阵的值空间为：")
        self.show_window.append("{}".format(test.columnspace()))

    def show_PLU(self):
        test = Matrix(self.mat)
        try:
            P, L, U = test.PLU_factorization()
            self.show_window.append("\n")
            self.show_window.append("矩阵的PLU分解为：")
            self.show_window.append("P：\n{}\n，L：\n{}\n，U：\n{}".format(P, L, U))
        except:
            error = "This matrix con't PLU factorization! \n Please check its number of rows and colunms to ensure that they are equal! \n More you should check its rank to ensure it's nonsingular!"
            self.show_window.append("\n")
            self.show_window.append(error)

    def show_gramschmidt(self):
        test = Matrix(self.mat)
        Q, R = test.Gram_Schmidt_orthogonalization()
        self.show_window.append("\n")
        self.show_window.append("矩阵的Gram_Schmidt正交化结果为：")
        self.show_window.append("Q：\n{}\n，R：\n{}".format(Q, R))

    def show_householder(self):
        test = Matrix(self.mat)
        Q, R = test.Householder_reduction()
        self.show_window.append("\n")
        self.show_window.append("矩阵的Householder约简结果为：")
        self.show_window.append("Q：\n{}\n，R：\n{}".format(Q, R))

    def show_givens(self):
        test = Matrix(self.mat)
        Q, R = test.Givens_reduction()
        self.show_window.append("\n")
        self.show_window.append("矩阵的Givens约简结果为：")
        self.show_window.append("Q：\n{}\n，R：\n{}".format(Q, R))

    def show_urv(self):
        test = Matrix(self.mat)
        U, R, V = test.URV_factorization()
        self.show_window.append("\n")
        self.show_window.append("矩阵的URV分解为：")
        self.show_window.append("U：\n{}\n，R：\n{}\n，V：\n{}".format(U, R, V))

    def show_A(self):
        test = Matrix(self.mat)
        self.show_window.append("\n")
        b = np.fromstring(self.b_edit.text(), dtype=np.float, sep=",")
        self.show_window.append("线性系统的b为\n{}".format(b))
        self.show_window.append("结果为{}".format(test.linear_equation(b)))


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    my_pyqt_form = MyPyQT_Form()
    my_pyqt_form.show()
    sys.exit(app.exec_())
