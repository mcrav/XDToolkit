# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'test.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.setFolderBut = QtWidgets.QPushButton(self.centralwidget)
        self.setFolderBut.setObjectName("setFolderBut")
        self.verticalLayout.addWidget(self.setFolderBut)
        self.line_2 = QtWidgets.QFrame(self.centralwidget)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.verticalLayout.addWidget(self.line_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.atomsDrop = QtWidgets.QComboBox(self.centralwidget)
        self.atomsDrop.setObjectName("atomsDrop")
        self.horizontalLayout.addWidget(self.atomsDrop)
        self.neebInput = QtWidgets.QLineEdit(self.centralwidget)
        self.neebInput.setObjectName("neebInput")
        self.horizontalLayout.addWidget(self.neebInput)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.addNeebsBut = QtWidgets.QPushButton(self.centralwidget)
        self.addNeebsBut.setObjectName("addNeebsBut")
        self.verticalLayout.addWidget(self.addNeebsBut)
        self.addedNeebsLab = QtWidgets.QLabel(self.centralwidget)
        self.addedNeebsLab.setText("")
        self.addedNeebsLab.setObjectName("addedNeebsLab")
        self.verticalLayout.addWidget(self.addedNeebsLab)
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)
        self.testBut = QtWidgets.QPushButton(self.centralwidget)
        self.testBut.setMaximumSize(QtCore.QSize(300, 16777214))
        self.testBut.setObjectName("testBut")
        self.verticalLayout.addWidget(self.testBut)
        self.testAllBut = QtWidgets.QPushButton(self.centralwidget)
        self.testAllBut.setObjectName("testAllBut")
        self.verticalLayout.addWidget(self.testAllBut)
        self.testLab = QtWidgets.QLabel(self.centralwidget)
        self.testLab.setText("")
        self.testLab.setObjectName("testLab")
        self.verticalLayout.addWidget(self.testLab)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 19))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "XD Toolkit Testing"))
        self.setFolderBut.setText(_translate("MainWindow", "Select folder to test"))
        self.addNeebsBut.setText(_translate("MainWindow", "Add correct neighbours"))
        self.testBut.setText(_translate("MainWindow", "Test"))
        self.testAllBut.setText(_translate("MainWindow", "Test everything"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

