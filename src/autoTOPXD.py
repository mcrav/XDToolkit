# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'autoTOPXD.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_autoTOPXD(object):
    def setupUi(self, autoTOPXD):
        autoTOPXD.setObjectName("autoTOPXD")
        autoTOPXD.setWindowModality(QtCore.Qt.ApplicationModal)
        autoTOPXD.resize(500, 450)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../res/flatearth.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        autoTOPXD.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(autoTOPXD)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(10, 10, 10, 10)
        self.verticalLayout.setSpacing(10)
        self.verticalLayout.setObjectName("verticalLayout")
        self.statusLab = QtWidgets.QLabel(autoTOPXD)
        font = QtGui.QFont()
        font.setItalic(True)
        self.statusLab.setFont(font)
        self.statusLab.setText("")
        self.statusLab.setObjectName("statusLab")
        self.verticalLayout.addWidget(self.statusLab)
        self.progressBar = QtWidgets.QProgressBar(autoTOPXD)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout.addWidget(self.progressBar)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.cancelBut = QtWidgets.QPushButton(autoTOPXD)
        self.cancelBut.setMinimumSize(QtCore.QSize(200, 0))
        self.cancelBut.setMaximumSize(QtCore.QSize(200, 16777215))
        self.cancelBut.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.cancelBut.setObjectName("cancelBut")
        self.horizontalLayout.addWidget(self.cancelBut)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.gridLayout.addLayout(self.verticalLayout, 1, 0, 1, 1)

        self.retranslateUi(autoTOPXD)
        QtCore.QMetaObject.connectSlotsByName(autoTOPXD)

    def retranslateUi(self, autoTOPXD):
        _translate = QtCore.QCoreApplication.translate
        autoTOPXD.setWindowTitle(_translate("autoTOPXD", "Auto TOPXD"))
        self.progressBar.setFormat(_translate("autoTOPXD", "%p%"))
        self.cancelBut.setText(_translate("autoTOPXD", "Cancel"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    autoTOPXD = QtWidgets.QWidget()
    ui = Ui_autoTOPXD()
    ui.setupUi(autoTOPXD)
    autoTOPXD.show()
    sys.exit(app.exec_())

