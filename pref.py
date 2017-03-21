# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'pref.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_pref(object):
    def setupUi(self, pref):
        pref.setObjectName("pref")
        pref.setWindowModality(QtCore.Qt.NonModal)
        pref.resize(698, 333)
        font = QtGui.QFont()
        font.setFamily("Consolas")
        font.setPointSize(10)
        pref.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("flatearth.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        pref.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(pref)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setContentsMargins(10, 10, 10, 10)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.chooseXDPathBut = QtWidgets.QPushButton(pref)
        self.chooseXDPathBut.setMaximumSize(QtCore.QSize(300, 16777215))
        self.chooseXDPathBut.setObjectName("chooseXDPathBut")
        self.verticalLayout_3.addWidget(self.chooseXDPathBut)
        self.settingsXDLab = QtWidgets.QLabel(pref)
        self.settingsXDLab.setObjectName("settingsXDLab")
        self.verticalLayout_3.addWidget(self.settingsXDLab)
        self.chooseMoleCoolPath = QtWidgets.QPushButton(pref)
        self.chooseMoleCoolPath.setMaximumSize(QtCore.QSize(300, 16777215))
        self.chooseMoleCoolPath.setObjectName("chooseMoleCoolPath")
        self.verticalLayout_3.addWidget(self.chooseMoleCoolPath)
        self.settingsMoleCoolLab = QtWidgets.QLabel(pref)
        self.settingsMoleCoolLab.setObjectName("settingsMoleCoolLab")
        self.verticalLayout_3.addWidget(self.settingsMoleCoolLab)
        self.chooseMercPathBut = QtWidgets.QPushButton(pref)
        self.chooseMercPathBut.setMaximumSize(QtCore.QSize(300, 16777215))
        self.chooseMercPathBut.setObjectName("chooseMercPathBut")
        self.verticalLayout_3.addWidget(self.chooseMercPathBut)
        self.chooseMercPathLab = QtWidgets.QLabel(pref)
        self.chooseMercPathLab.setObjectName("chooseMercPathLab")
        self.verticalLayout_3.addWidget(self.chooseMercPathLab)
        self.chooseTextPathBut = QtWidgets.QPushButton(pref)
        self.chooseTextPathBut.setMaximumSize(QtCore.QSize(300, 16777215))
        self.chooseTextPathBut.setObjectName("chooseTextPathBut")
        self.verticalLayout_3.addWidget(self.chooseTextPathBut)
        self.chooseTextPathLab = QtWidgets.QLabel(pref)
        self.chooseTextPathLab.setObjectName("chooseTextPathLab")
        self.verticalLayout_3.addWidget(self.chooseTextPathLab)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem)
        self.gridLayout.addLayout(self.verticalLayout_3, 0, 0, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(pref)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 2, 0, 1, 1)
        self.anonDatBox = QtWidgets.QCheckBox(pref)
        self.anonDatBox.setObjectName("anonDatBox")
        self.gridLayout.addWidget(self.anonDatBox, 1, 0, 1, 1)

        self.retranslateUi(pref)
        self.buttonBox.accepted.connect(pref.accept)
        self.buttonBox.rejected.connect(pref.reject)
        QtCore.QMetaObject.connectSlotsByName(pref)

    def retranslateUi(self, pref):
        _translate = QtCore.QCoreApplication.translate
        pref.setWindowTitle(_translate("pref", "Preferences"))
        self.chooseXDPathBut.setText(_translate("pref", "Select XD location"))
        self.settingsXDLab.setText(_translate("pref", "Current path:"))
        self.chooseMoleCoolPath.setText(_translate("pref", "Select MoleCoolQt location"))
        self.settingsMoleCoolLab.setText(_translate("pref", "Current path: "))
        self.chooseMercPathBut.setText(_translate("pref", "Select Mercury location"))
        self.chooseMercPathLab.setText(_translate("pref", "Current path:"))
        self.chooseTextPathBut.setText(_translate("pref", "Select text editor location"))
        self.chooseTextPathLab.setText(_translate("pref", "Current path:"))
        self.anonDatBox.setText(_translate("pref", "Send anonymous usage data"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    pref = QtWidgets.QDialog()
    ui = Ui_pref()
    ui.setupUi(pref)
    pref.show()
    sys.exit(app.exec_())

