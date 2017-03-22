# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'about.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_aboutBox(object):
    def setupUi(self, aboutBox):
        aboutBox.setObjectName("aboutBox")
        aboutBox.setWindowModality(QtCore.Qt.WindowModal)
        aboutBox.resize(767, 250)
        font = QtGui.QFont()
        font.setFamily("Consolas")
        font.setPointSize(10)
        aboutBox.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("flatearth.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        aboutBox.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(aboutBox)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(aboutBox)
        self.label.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)

        self.retranslateUi(aboutBox)
        QtCore.QMetaObject.connectSlotsByName(aboutBox)

    def retranslateUi(self, aboutBox):
        _translate = QtCore.QCoreApplication.translate
        aboutBox.setWindowTitle(_translate("aboutBox", "About XD Toolkit"))
        self.label.setText(_translate("aboutBox", "<html><head/><body><p><span style=\" font-weight:600;\">XD Toolkit</span><br/>Version 0.2.0</p><p><span style=\" font-weight:600;\">License</span><br/>GNU GPL v3.0</p><p><span style=\" font-weight:600;\">Creator:</span> Matthew Craven</p><p><span style=\" font-weight:600;\">Quickplot:</span> Mads Ry Vogel Jørgensen</p><p><span style=\" font-weight:600;\">XD 2016:</span> XD2016 - A Computer Program Package for Multipole Refinement, Topological Analysis of Charge Densities and Evaluation of Intermolecular Energies from Experimental and Theoretical Structure Factors <br/>Volkov, A.; Macchi, P.; Farrugia, L. J.; Gatti, C.; Mallinson, P.; Richter, T.; Koritsanszky, T. (2016) </p></body></html>"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    aboutBox = QtWidgets.QWidget()
    ui = Ui_aboutBox()
    ui.setupUi(aboutBox)
    aboutBox.show()
    sys.exit(app.exec_())

