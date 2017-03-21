# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'wizard.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_wizard(object):
    def setupUi(self, wizard):
        wizard.setObjectName("wizard")
        wizard.setWindowModality(QtCore.Qt.ApplicationModal)
        wizard.resize(400, 357)
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans Mono")
        wizard.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("flatearth.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        wizard.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(wizard)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.wizStatusLab = QtWidgets.QLabel(wizard)
        self.wizStatusLab.setObjectName("wizStatusLab")
        self.verticalLayout.addWidget(self.wizStatusLab)
        self.line = QtWidgets.QFrame(wizard)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)
        self.wizResLab = QtWidgets.QLabel(wizard)
        self.wizResLab.setObjectName("wizResLab")
        self.verticalLayout.addWidget(self.wizResLab)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(wizard)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        self.retranslateUi(wizard)
        self.buttonBox.accepted.connect(wizard.accept)
        self.buttonBox.rejected.connect(wizard.reject)
        QtCore.QMetaObject.connectSlotsByName(wizard)

    def retranslateUi(self, wizard):
        _translate = QtCore.QCoreApplication.translate
        wizard.setWindowTitle(_translate("wizard", "XD Wizard"))
        self.wizStatusLab.setText(_translate("wizard", "Starting wizard"))
        self.wizResLab.setText(_translate("wizard", "Results"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    wizard = QtWidgets.QDialog()
    ui = Ui_wizard()
    ui.setupUi(wizard)
    wizard.show()
    sys.exit(app.exec_())

