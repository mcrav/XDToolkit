# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'checkneebs.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_checkneebs(object):
    def setupUi(self, checkneebs):
        checkneebs.setObjectName("checkneebs")
        checkneebs.resize(858, 563)
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans Mono")
        font.setPointSize(9)
        checkneebs.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../res/flatearth.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        checkneebs.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(checkneebs)
        self.gridLayout.setObjectName("gridLayout")
        self.scrollArea = QtWidgets.QScrollArea(checkneebs)
        self.scrollArea.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents_2 = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_2.setGeometry(QtCore.QRect(0, 0, 838, 543))
        self.scrollAreaWidgetContents_2.setObjectName("scrollAreaWidgetContents_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.scrollAreaWidgetContents_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(10, 10, 10, 10)
        self.verticalLayout.setSpacing(10)
        self.verticalLayout.setObjectName("verticalLayout")
        self.infoLab = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        self.infoLab.setWordWrap(True)
        self.infoLab.setObjectName("infoLab")
        self.verticalLayout.addWidget(self.infoLab)
        self.line = QtWidgets.QFrame(self.scrollAreaWidgetContents_2)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout.setSpacing(30)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.tableLab = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        self.tableLab.setText("")
        self.tableLab.setObjectName("tableLab")
        self.horizontalLayout.addWidget(self.tableLab)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.tableLab2 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        self.tableLab2.setText("")
        self.tableLab2.setWordWrap(True)
        self.tableLab2.setObjectName("tableLab2")
        self.horizontalLayout.addWidget(self.tableLab2)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.gridLayout_2.addLayout(self.verticalLayout, 1, 0, 1, 1)


    def retranslateUi(self, checkneebs):
        _translate = QtCore.QCoreApplication.translate
        checkneebs.setWindowTitle(_translate("checkneebs", "Check Neighbours"))
        self.infoLab.setText(_translate("checkneebs", "If nearest neighbours are incorrect for any atom, you must add the local coordinate system for this atom manually in the \"Tools\" tab.\n"
"\n"
"Also, if an H atom in missing/wrongly-included in the nearest neighbours, you must add or delete the reset bond instruction for this atom in the \"RESET BOND\" tab."))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    checkneebs = QtWidgets.QWidget()
    ui = Ui_checkneebs()
    ui.setupUi(checkneebs)
    checkneebs.show()
    sys.exit(app.exec_())

