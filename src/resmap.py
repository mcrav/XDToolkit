# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'resmap.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_resmap(object):
    def setupUi(self, resmap):
        resmap.setObjectName("resmap")
        resmap.resize(400, 300)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(231, 143, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(231, 143, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(239, 235, 231))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        resmap.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans Mono")
        font.setPointSize(10)
        resmap.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("flatearth.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        resmap.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(resmap)
        self.gridLayout.setObjectName("gridLayout")
        self.resmapLayout = QtWidgets.QVBoxLayout()
        self.resmapLayout.setContentsMargins(10, 10, 10, 10)
        self.resmapLayout.setObjectName("resmapLayout")
        self.gridLayout.addLayout(self.resmapLayout, 0, 0, 1, 1)

        self.retranslateUi(resmap)
        QtCore.QMetaObject.connectSlotsByName(resmap)

    def retranslateUi(self, resmap):
        _translate = QtCore.QCoreApplication.translate
        resmap.setWindowTitle(_translate("resmap", "Quickplot residual density map"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    resmap = QtWidgets.QWidget()
    ui = Ui_resmap()
    ui.setupUi(resmap)
    resmap.show()
    sys.exit(app.exec_())

