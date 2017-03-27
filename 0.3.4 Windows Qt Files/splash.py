# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'splash.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_splash(object):
    def setupUi(self, splash):
        splash.setObjectName("splash")
        splash.resize(400, 300)
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans Mono")
        font.setPointSize(10)
        splash.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("res/flatearth.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        splash.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(splash)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.statusLab = QtWidgets.QLabel(splash)
        self.statusLab.setText("")
        self.statusLab.setObjectName("statusLab")
        self.verticalLayout.addWidget(self.statusLab)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(splash)
        QtCore.QMetaObject.connectSlotsByName(splash)

    def retranslateUi(self, splash):
        _translate = QtCore.QCoreApplication.translate
        splash.setWindowTitle(_translate("splash", "Initializing"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    splash = QtWidgets.QDialog()
    ui = Ui_splash()
    ui.setupUi(splash)
    splash.show()
    sys.exit(app.exec_())

