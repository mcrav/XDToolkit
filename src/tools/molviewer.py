import sys
import math

from PyQt4.QtOpenGL import QGLWidget
from PyQt4.QtCore import pyqtSignal, QPoint, QSize, Qt
from PyQt4.QtGui import (QColor, QApplication, QHBoxLayout,  QSlider,
                         QWidget)
from traceback import print_exception

import OpenGL.GL as gl
import OpenGL.GLU as glu
import OpenGL.GLUT as glut
from numpy import array
import numpy as np

def fractional2cartesian(pos, a, b, c, alpha, beta, gamma, volume=None):
    if not volume:
        volume = a*b*c*np.sqrt(1-(np.cos(alpha))**2 - (np.cos(beta))**2 - (np.cos(gamma))**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    return np.dot(np.array([[a, b*np.cos(gamma), c*np.cos(beta)],
                            [0, b*np.sin(gamma), c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))],
                            [0,0,volume/(a*b*np.sin(gamma))]]), np.array(pos))
#rawAtomPos = {'DUM,asym':(array([0,0,0]),'xdfes'), 'O(1),asym': (array([ 0.658581,  0.407033,  0.427257]), [('x', 'y', 'z')]), 'O(1),1': (array([-0.158581,  0.907033, -0.427257]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'O(1),3': (array([ 1.158581,  0.092967,  0.427257]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(1),4': (array([-0.658581, -0.407033, -0.427257]), [('-x', '-y', '-z')]), 'O(2),asym': (array([ 0.828859,  0.511055,  0.751749]), [('x', 'y', 'z')]), 'O(2),1': (array([-0.328859,  1.011055, -0.751749]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'O(2),3': (array([ 1.328859, -0.011055,  0.751749]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(2),4': (array([-0.828859, -0.511055, -0.751749]), [('-x', '-y', '-z')]), 'O(3),asym': (array([ 0.928844,  0.175498,  0.721162]), [('x', 'y', 'z')]), 'O(3),1': (array([-0.428844,  0.675498, -0.721162]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'O(3),3': (array([ 1.428844,  0.324502,  0.721162]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'O(3),4': (array([-0.928844, -0.175498, -0.721162]), [('-x', '-y', '-z')]), 'N(1),asym': (array([ 0.653265,  0.173421,  0.743706]), [('x', 'y', 'z')]), 'N(1),1': (array([-0.153265,  0.673421, -0.743706]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'N(1),3': (array([ 1.153265,  0.326579,  0.743706]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'N(1),4': (array([-0.653265, -0.173421, -0.743706]), [('-x', '-y', '-z')]), 'C(1),asym': (array([ 0.746436,  0.410457,  0.666236]), [('x', 'y', 'z')]), 'C(1),1': (array([-0.246436,  0.910457, -0.666236]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'C(1),3': (array([ 1.246436,  0.089543,  0.666236]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(1),4': (array([-0.746436, -0.410457, -0.666236]), [('-x', '-y', '-z')]), 'C(2),asym': (array([ 0.754494,  0.283584,  0.878906]), [('x', 'y', 'z')]), 'C(2),1': (array([-0.254494,  0.783584, -0.878906]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'C(2),3': (array([ 1.254494,  0.216416,  0.878906]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(2),4': (array([-0.754494, -0.283584, -0.878906]), [('-x', '-y', '-z')]), 'C(3),asym': (array([ 0.889469,  0.216379,  0.96978 ]), [('x', 'y', 'z')]), 'C(3),1': (array([-0.389469,  0.716379, -0.96978 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'C(3),3': (array([ 1.389469,  0.283621,  0.96978 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'C(3),4': (array([-0.889469, -0.216379, -0.96978 ]), [('-x', '-y', '-z')]), 'H(2),asym': (array([ 0.733868,  0.319768,  1.047163]), [('x', 'y', 'z')]), 'H(2),1': (array([-0.233868,  0.819768, -1.047163]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(2),3': (array([ 1.233868,  0.180232,  1.047163]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(2),4': (array([-0.733868, -0.319768, -1.047163]), [('-x', '-y', '-z')]), 'H(4),asym': (array([ 0.905169,  0.081724,  0.67962 ]), [('x', 'y', 'z')]), 'H(4),1': (array([-0.405169,  0.581724, -0.67962 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(4),3': (array([ 1.405169,  0.418276,  0.67962 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(4),4': (array([-0.905169, -0.081724, -0.67962 ]), [('-x', '-y', '-z')]), 'H(11),asym': (array([ 0.572607,  0.218384,  0.699938]), [('x', 'y', 'z')]), 'H(11),1': (array([-0.072607,  0.718384, -0.699938]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(11),3': (array([ 1.072607,  0.281616,  0.699938]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(11),4': (array([-0.572607, -0.218384, -0.699938]), [('-x', '-y', '-z')]), 'H(12),asym': (array([ 0.652058,  0.108138,  0.876204]), [('x', 'y', 'z')]), 'H(12),1': (array([-0.152058,  0.608138, -0.876204]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(12),3': (array([ 1.152058,  0.391862,  0.876204]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(12),4': (array([-0.652058, -0.108138, -0.876204]), [('-x', '-y', '-z')]), 'H(13),asym': (array([ 0.665591,  0.131178,  0.58416 ]), [('x', 'y', 'z')]), 'H(13),1': (array([-0.165591,  0.631178, -0.58416 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(13),3': (array([ 1.165591,  0.368822,  0.58416 ]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(13),4': (array([-0.665591, -0.131178, -0.58416 ]), [('-x', '-y', '-z')]), 'H(31),asym': (array([ 0.950043,  0.289619,  1.084379]), [('x', 'y', 'z')]), 'H(31),1': (array([-0.450043,  0.789619, -1.084379]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(31),3': (array([ 1.450043,  0.210381,  1.084379]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(31),4': (array([-0.950043, -0.289619, -1.084379]), [('-x', '-y', '-z')]), 'H(32),asym': (array([ 0.891395,  0.133951,  1.099092]), [('x', 'y', 'z')]), 'H(32),1': (array([-0.391395,  0.633951, -1.099092]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z')]), 'H(32),3': (array([ 1.391395,  0.366049,  1.099092]), [['0.5-x', '-0.5+y', '-z'], ('x', 'y', 'z'), ('-x', '-y', '-z')]), 'H(32),4': (array([-0.891395, -0.133951, -1.099092]), [('-x', '-y', '-z')])}
atomPos = {'DUM,asym': array([0,0,0]), 'N(1),asym': array([ 9.87384816,  1.42304306,  0.19836874]), 'O(1),asym': array([ 7.45878825,  3.34439748,  0.11396228]), 'C(1),asym': array([ 9.53169452,  3.371745  ,  0.17770516]), 'H(2),asym': array([ 11.5874222 ,   2.62500682,   0.27930984]), 'H(11),asym': array([ 8.61056972,  1.79278644,  0.1866945 ]), 'H(4),asym': array([ 12.66047246,   0.66950482,   0.18127508]), 'H(31),asym': array([ 14.21879193,   2.37705685,   0.28923646]), 'C(2),asym': array([ 11.15595792,   2.32813496,   0.23443064]), 'H(32),asym': array([ 14.29831546,   1.09741193,   0.29316086]), 'O(2),asym': array([ 10.41347615,   4.19837819,   0.20051405]), 'H(12),asym': array([ 10.76199483,   0.8859724 ,   0.23370994]), 'O(3),asym': array([ 12.72747537,   1.44019133,   0.19235558]), 'DUM,asym': array([ 0.,  0.,  0.]), 'H(13),asym': array([ 9.41946547,  1.07633695,  0.15581303]), 'C(3),asym': array([ 13.32097211,   1.77540438,   0.25866947])}
#atomPos = {}
#
#for atom, info in rawAtomPos.items():
#    if atom.split(',')[1]=='asym':
#        #atomPos[atom] = fractional2cartesian(info[0], 10.7764, 9.1947, 4.7788, 90, 106.87, 90)
#        atomPos[atom] = info[0]
#print(atomPos)

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
        self.lastAtomPos = None
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
        gl.glOrtho(0, 10.7764, 0, 9.1947, 0, 4.7788)
        gl.glMatrixMode(gl.GL_MODELVIEW)
    

    def mousePressEvent(self, event):
        self.lastPos = event.pos()
    
    def wheelEvent(self, event):

        gl.glPushMatrix()
        #gl.glScalef(10,10,0)
        #gl.glCallList(self.object)
        
        

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

        centerPos = (0,0,0)
        for atom, pos in atomPos.items():
          
            if atom.startswith('C'):
                color = [0.5, 0.5, 0.5, 1.0]
            elif atom.startswith('O'):
                color = [1.0, 0.0, 0.0, 1.0]
            elif atom.startswith('N'):
                color = [0.0, 0.0, 1.0, 1.0]
            elif atom.startswith('H'):
                continue
                color = [1.0,1.0,1.0,1.0]
            elif atom.startswith('DUM'):
                color = [1.0,1.0,0,1.0]
            if atom.startswith('C(2)'):
                centerPos = pos
            
            self.drawAtom(color, *pos)       
            self.lastAtomPos = pos
        gl.glEndList()
        return genList

    def drawAtom(self, color, x, y, z):
        self.setColor(self.trolltechGreen)
        gl.glTranslatef(x/5,y/5,z/5)
        
        gl.glPushMatrix()
        
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_DIFFUSE, color)
        quadric = glu.gluNewQuadric()
        glu.gluQuadricDrawStyle(quadric, glu.GLU_FILL)
        glu.gluQuadricTexture
        glu.gluQuadricNormals(quadric, glu.GLU_SMOOTH)
        glu.gluSphere(quadric,0.02,10,10)
        gl.glPopMatrix()
        gl.glTranslatef(-x/5,-y/5,-z/5)
      
        
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
        
def customExceptHook(Type, value, traceback):
    print('Exception:')
    print_exception(Type, value, traceback)
    pass

if __name__ == '__main__':
    sys.excepthook = customExceptHook
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())