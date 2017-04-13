# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'sendBug.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_sendBug(object):
    def setupUi(self, sendBug):
        sendBug.setObjectName("sendBug")
        sendBug.resize(597, 300)
        font = QtGui.QFont()
        font.setFamily("Bitstream Vera Sans Mono")
        font.setPointSize(10)
        sendBug.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("flatearth.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        sendBug.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(sendBug)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(10, 10, 10, 10)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(sendBug)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.bugReportInput = QtWidgets.QPlainTextEdit(sendBug)
        self.bugReportInput.setObjectName("bugReportInput")
        self.horizontalLayout_2.addWidget(self.bugReportInput)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(sendBug)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.emailAdInput = QtWidgets.QLineEdit(sendBug)
        self.emailAdInput.setObjectName("emailAdInput")
        self.horizontalLayout.addWidget(self.emailAdInput)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.incFilesBox = QtWidgets.QCheckBox(sendBug)
        self.incFilesBox.setObjectName("incFilesBox")
        self.verticalLayout.addWidget(self.incFilesBox)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(sendBug)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        self.retranslateUi(sendBug)
        self.buttonBox.accepted.connect(sendBug.accept)
        self.buttonBox.rejected.connect(sendBug.reject)
        QtCore.QMetaObject.connectSlotsByName(sendBug)

    def retranslateUi(self, sendBug):
        _translate = QtCore.QCoreApplication.translate
        sendBug.setWindowTitle(_translate("sendBug", "Submit Bug Report"))
        self.label_2.setText(_translate("sendBug", "Message:"))
        self.label.setText(_translate("sendBug", "  Email:"))
        self.emailAdInput.setPlaceholderText(_translate("sendBug", "Only include if you want a response"))
        self.incFilesBox.setText(_translate("sendBug", "Send shelx.ins, shelx.hkl files for debugging."))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    sendBug = QtWidgets.QDialog()
    ui = Ui_sendBug()
    ui.setupUi(sendBug)
    sendBug.show()
    sys.exit(app.exec_())

