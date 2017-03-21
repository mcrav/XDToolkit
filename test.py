# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:50:17 2017

@author: Bruker
"""
from PyQt5.QtWidgets import QWidget, QMessageBox, QPushButton, QApplication, QDialog, QLineEdit, QFileDialog
from PyQt5.QtCore import QCoreApplication
import sys
import os

def test():
    filename = QFileDialog.getOpenFileName(None, 'Test Dialog', os.getcwd(), 'All Files(*.*)')
 
def main():
     app = QApplication(sys.argv)

     test()

     sys.exit(app.exec_())

if __name__ == "__main__":
     main()