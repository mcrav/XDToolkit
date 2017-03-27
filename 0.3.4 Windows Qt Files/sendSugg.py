# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'sendSugg.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_sendSugg(object):
    def setupUi(self, sendSugg):
        sendSugg.setObjectName("sendSugg")
        sendSugg.resize(597, 300)
        font = QtGui.QFont()
        font.setFamily("Consolas")
        font.setPointSize(10)
        sendSugg.setFont(font)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("res/flatearth.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        sendSugg.setWindowIcon(icon)
        self.gridLayout = QtWidgets.QGridLayout(sendSugg)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(10, 10, 10, 10)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(sendSugg)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.suggInput = QtWidgets.QPlainTextEdit(sendSugg)
        self.suggInput.setObjectName("suggInput")
        self.horizontalLayout_2.addWidget(self.suggInput)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(sendSugg)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.emailAdInput = QtWidgets.QLineEdit(sendSugg)
        self.emailAdInput.setText("")
        self.emailAdInput.setObjectName("emailAdInput")
        self.horizontalLayout.addWidget(self.emailAdInput)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(sendSugg)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        self.retranslateUi(sendSugg)
        self.buttonBox.accepted.connect(sendSugg.accept)
        self.buttonBox.rejected.connect(sendSugg.reject)
        QtCore.QMetaObject.connectSlotsByName(sendSugg)

    def retranslateUi(self, sendSugg):
        _translate = QtCore.QCoreApplication.translate
        sendSugg.setWindowTitle(_translate("sendSugg", "Submit Suggestion"))
        self.label_2.setText(_translate("sendSugg", "Message:"))
        self.label.setText(_translate("sendSugg", "  Email:"))
        self.emailAdInput.setPlaceholderText(_translate("sendSugg", "Only include if you are happy to be contacted"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    sendSugg = QtWidgets.QDialog()
    ui = Ui_sendSugg()
    ui.setupUi(sendSugg)
    sendSugg.show()
    sys.exit(app.exec_())

