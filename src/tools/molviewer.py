import sys
import math

from PyQt4.QtOpenGL import QGLWidget
from PyQt4.QtCore import pyqtSignal, QPoint, QSize, Qt
from PyQt4.QtGui import (QColor, QApplication, QHBoxLayout,  QSlider,
                         QWidget)


import OpenGL.GL as gl
import OpenGL.GLU as glu
import OpenGL.GLUT as glut
from numpy import array
import numpy as np
rawAtomPos = {'O(1),asym': (array([ 0.658581,  0.407033,  0.427257]), [('x', 'y', 'z')]), 'O(1),1': (array([-0.158581,  0.907033, -0.427257]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'O(1),3': (array([ 1.158581,  0.092967,  0.427257]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(1),4': (array([-0.658581, -0.407033, -0.427257]), [('-x', '-y', '-z')]), 'O(2),asym': (array([ 0.828859,  0.511055,  0.751749]), [('x', 'y', 'z')]), 'O(2),1': (array([-0.328859,  1.011055, -0.751749]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'O(2),3': (array([ 1.328859, -0.011055,  0.751749]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(2),4': (array([-0.828859, -0.511055, -0.751749]), [('-x', '-y', '-z')]), 'O(3),asym': (array([ 0.928844,  0.175498,  0.721162]), [('x', 'y', 'z')]), 'O(3),1': (array([-0.428844,  0.675498, -0.721162]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'O(3),3': (array([ 1.428844,  0.324502,  0.721162]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(3),4': (array([-0.928844, -0.175498, -0.721162]), [('-x', '-y', '-z')]), 'N(1),asym': (array([ 0.653265,  0.173421,  0.743706]), [('x', 'y', 'z')]), 'N(1),1': (array([-0.153265,  0.673421, -0.743706]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'N(1),3': (array([ 1.153265,  0.326579,  0.743706]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'N(1),4': (array([-0.653265, -0.173421, -0.743706]), [('-x', '-y', '-z')]), 'C(1),asym': (array([ 0.746436,  0.410457,  0.666236]), [('x', 'y', 'z')]), 'C(1),1': (array([-0.246436,  0.910457, -0.666236]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'C(1),3': (array([ 1.246436,  0.089543,  0.666236]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(1),4': (array([-0.746436, -0.410457, -0.666236]), [('-x', '-y', '-z')]), 'C(2),asym': (array([ 0.754494,  0.283584,  0.878906]), [('x', 'y', 'z')]), 'C(2),1': (array([-0.254494,  0.783584, -0.878906]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'C(2),3': (array([ 1.254494,  0.216416,  0.878906]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(2),4': (array([-0.754494, -0.283584, -0.878906]), [('-x', '-y', '-z')]), 'C(3),asym': (array([ 0.889469,  0.216379,  0.96978 ]), [('x', 'y', 'z')]), 'C(3),1': (array([-0.389469,  0.716379, -0.96978 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'C(3),3': (array([ 1.389469,  0.283621,  0.96978 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(3),4': (array([-0.889469, -0.216379, -0.96978 ]), [('-x', '-y', '-z')]), 'H(2),asym': (array([ 0.733868,  0.319768,  1.047163]), [('x', 'y', 'z')]), 'H(2),1': (array([-0.233868,  0.819768, -1.047163]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(2),3': (array([ 1.233868,  0.180232,  1.047163]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(2),4': (array([-0.733868, -0.319768, -1.047163]), [('-x', '-y', '-z')]), 'H(4),asym': (array([ 0.905169,  0.081724,  0.67962 ]), [('x', 'y', 'z')]), 'H(4),1': (array([-0.405169,  0.581724, -0.67962 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(4),3': (array([ 1.405169,  0.418276,  0.67962 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(4),4': (array([-0.905169, -0.081724, -0.67962 ]), [('-x', '-y', '-z')]), 'H(11),asym': (array([ 0.572607,  0.218384,  0.699938]), [('x', 'y', 'z')]), 'H(11),1': (array([-0.072607,  0.718384, -0.699938]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(11),3': (array([ 1.072607,  0.281616,  0.699938]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(11),4': (array([-0.572607, -0.218384, -0.699938]), [('-x', '-y', '-z')]), 'H(12),asym': (array([ 0.652058,  0.108138,  0.876204]), [('x', 'y', 'z')]), 'H(12),1': (array([-0.152058,  0.608138, -0.876204]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(12),3': (array([ 1.152058,  0.391862,  0.876204]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(12),4': (array([-0.652058, -0.108138, -0.876204]), [('-x', '-y', '-z')]), 'H(13),asym': (array([ 0.665591,  0.131178,  0.58416 ]), [('x', 'y', 'z')]), 'H(13),1': (array([-0.165591,  0.631178, -0.58416 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(13),3': (array([ 1.165591,  0.368822,  0.58416 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(13),4': (array([-0.665591, -0.131178, -0.58416 ]), [('-x', '-y', '-z')]), 'H(31),asym': (array([ 0.950043,  0.289619,  1.084379]), [('x', 'y', 'z')]), 'H(31),1': (array([-0.450043,  0.789619, -1.084379]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(31),3': (array([ 1.450043,  0.210381,  1.084379]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(31),4': (array([-0.950043, -0.289619, -1.084379]), [('-x', '-y', '-z')]), 'H(32),asym': (array([ 0.891395,  0.133951,  1.099092]), [('x', 'y', 'z')]), 'H(32),1': (array([-0.391395,  0.633951, -1.099092]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(32),3': (array([ 1.391395,  0.366049,  1.099092]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(32),4': (array([-0.891395, -0.133951, -1.099092]), [('-x', '-y', '-z')])}
atomPos = {}

for atom, info in rawAtomPos.items():
    if atom.split(',')[1]=='asym':
        atomPos[atom] = info[0]
    
class Window(QWidget):

    def __init__(self):
        super(Window, self).__init__()

        self.glWidget = GLWidget()



        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.glWidget)

        self.setLayout(mainLayout)


        
        glut.glutInit()
        self.setWindowTitle("Hello GL")

    def createSlider(self):
        slider = QSlider(Qt.Vertical)

        slider.setRange(0, 360 * 16)
        slider.setSingleStep(16)
        slider.setPageStep(15 * 16)
        slider.setTickInterval(15 * 16)
        slider.setTickPosition(QSlider.TicksRight)

        return slider


class GLWidget(QGLWidget):
    xRotationChanged = pyqtSignal(int)
    yRotationChanged = pyqtSignal(int)
    zRotationChanged = pyqtSignal(int)

    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)

        self.object = 0
        self.xRot = 0
        self.yRot = 0
        self.zRot = 0
        
        self.zoom = 0

        self.lastPos = QPoint()

        self.trolltechGreen = QColor.fromCmykF(0.40, 0.0, 1.0, 0.0)
        self.trolltechPurple = QColor.fromCmykF(0.39, 0.39, 0.0, 0.0)

    def getOpenglInfo(self):
        info = """
            Vendor: {0}
            Renderer: {1}
            OpenGL Version: {2}
            Shader Version: {3}
        """.format(
            gl.glGetString(gl.GL_VENDOR),
            gl.glGetString(gl.GL_RENDERER),
            gl.glGetString(gl.GL_VERSION),
            gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)
        )

        return info

    def minimumSizeHint(self):
        return QSize(50, 50)

    def sizeHint(self):
        return QSize(400, 400)

    def setXRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.xRot:
            self.xRot = angle
            self.xRotationChanged.emit(angle)
            self.update()

    def setYRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.yRot:
            self.yRot = angle
            self.yRotationChanged.emit(angle)
            self.update()

    def setZRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.zRot:
            self.zRot = angle
            self.zRotationChanged.emit(angle)
            self.update()

    def initializeGL(self):
        print(self.getOpenglInfo())

        self.setClearColor(self.trolltechPurple.darker())
        self.object = self.makeObject()

        mat_specular =  (1.0, 1.0, 1.0, 1.0 )
        mat_shininess = 50.0
        light_position = ( 1.0, 1.0, 1.0, 0.0 )
        gl.glClearColor = (0.0, 0.0, 0.0, 0.0)
        gl.glShadeModel = gl.GL_SMOOTH
    
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_SPECULAR, mat_specular)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_SHININESS, mat_shininess)
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, light_position)
    
        gl.glEnable(gl.GL_LIGHTING)
        gl.glEnable(gl.GL_LIGHT0)
        gl.glEnable(gl.GL_DEPTH_TEST)


    def paintGL(self):
        gl.glClear(
            gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        gl.glLoadIdentity()
        #glu.gluLookAt(0, 0, self.zoom, 0, 0, 0, 0, 1, 0)
        gl.glTranslated(0.0, 0.0, -10.0)
        gl.glRotated(self.xRot / 16.0, 1.0, 0.0, 0.0)
        gl.glRotated(self.yRot / 16.0, 0.0, 1.0, 0.0)
        gl.glRotated(self.zRot / 16.0, 0.0, 0.0, 1.0)
        gl.glCallList(self.object)

    def resizeGL(self, width, height):
        side = min(width, height)
        if side < 0:
            return

        gl.glViewport((width - side) // 2, (height - side) // 2, side,
                           side)

        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        gl.glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
        gl.glMatrixMode(gl.GL_MODELVIEW)

    def mousePressEvent(self, event):
        self.lastPos = event.pos()
    
    def wheelEvent(self, event):
        print('mousewheel')
        self.zoom+=1

    def mouseMoveEvent(self, event):
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if event.buttons() & Qt.LeftButton:
            self.setXRotation(self.xRot + 8 * dy)
            self.setYRotation(self.yRot + 8 * dx)
        elif event.buttons() & Qt.RightButton:
            self.setXRotation(self.xRot + 8 * dy)
            self.setZRotation(self.zRot + 8 * dx)

        self.lastPos = event.pos()
        
    def makeObject(self):
   
        genList = gl.glGenLists(1)
        gl.glNewList(genList, gl.GL_COMPILE)


        for atom, pos in atomPos.items():
          
            if atom.startswith('C'):
                color = [0.5, 0.5, 0.5, 1.0]
            elif atom.startswith('O'):
                color = [1.0, 0.0, 0.0, 1.0]
            elif atom.startswith('N'):
                color = [0.0, 0.0, 1.0, 1.0]
            elif atom.startswith('H'):
                color = [1.0,1.0,1.0,1.0]
            self.drawAtom(color,*pos)
            print(pos)

          

        gl.glEndList()
        return genList

    def drawAtom(self, color, x, y, z):
        self.setColor(self.trolltechGreen)
        gl.glTranslatef(x,y,z)
        gl.glPushMatrix()
        
        
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_DIFFUSE, color)
        quadric = glu.gluNewQuadric()
        glu.gluQuadricDrawStyle(quadric, glu.GLU_FILL)
        glu.gluQuadricTexture
        glu.gluQuadricNormals(quadric, glu.GLU_SMOOTH)
        glu.gluSphere(quadric,0.05,20,20)
        gl.glPopMatrix()
        gl.glTranslatef(-x,-y,-z)
        
    def normalizeAngle(self, angle):
        while angle < 0:
            angle += 360 * 16
        while angle > 360 * 16:
            angle -= 360 * 16
        return angle

    def setClearColor(self, c):
        gl.glClearColor(c.redF(), c.greenF(), c.blueF(), c.alphaF())

    def setColor(self, c):
        gl.glColor4f(c.redF(), c.greenF(), c.blueF(), c.alphaF())


if __name__ == '__main__':

    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())