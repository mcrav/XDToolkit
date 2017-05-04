"""
XDPlotter i3d: A program to plot 3D isosurfaces from XD2006 grd files.

Distributed under the beer license: If you find this program useful - you can 
buy me a beer (or two)...

Mads Joergensen, 2013, Aarhus University

Version tracking: Describe changes and update version number below section. 
0.2:    Higher resolution of the atoms and bonds. Fixed atom bug for atoms with 
        two letter symbols. Possibility to turn off isosurfaces independently. 
        New default bond color. Posibilty to modify atom presets added. Slight 
        modification of GUI layout.
        (October 9th 2013)
0.3     Moved labels away from center of atoms. (January 10th 2014)
"""
version = 0.3

################################################################################

import itertools
import re

import sys

import numpy as np
from mayavi import mlab

from traits.api import HasTraits, Float, Str, Bool, RGBColor, on_trait_change
from traitsui.api import View, Item, HGroup, VGroup
from traitsui.menu import OKCancelButtons

import xd_grd_lib as xd
import atom_dictionary as atomdata


################################################################################
#Very bad colour hack!!!
_NUMERALS = '0123456789abcdefABCDEF'
_HEXDEC = {v: int(v, 16) for v in (x+y for x in _NUMERALS for y in _NUMERALS)}
LOWERCASE, UPPERCASE = 'x', 'X'

def rgb(triplet):
    return float(_HEXDEC[triplet[1:3]])/255.0, float(_HEXDEC[triplet[3:5]])/255.0, float(_HEXDEC[triplet[5:7]])/255.0
################################################################################


def get_version():
    "Version tracking"""
    return "i3d: " + str(version)
    
class param_3d_isosurface(HasTraits):
    """
    3D isosurface object, parameters, GUI setup and plotting
    """
    filename = Str() # Input filename
    dim = Str()
    func = Str()
    # Surface parameters
    show_pos = Bool(True)
    pos_surf = Float(0.30)
    pos_color = RGBColor((0, 0, 1))
    show_neg = Bool(True)
    neg_surf = Float(-0.30)
    neg_color = RGBColor((1, 0, 0))
    opacity = Float(1.0)
    # Atom parameters
    show_atoms = Bool(True)
    atom_size = Float(0.5)
    label_atoms = Bool(True)
    label_color = RGBColor((0, 0, 0))
    label_size = Float(0.5)
    # Bond parameters            
    show_bonds = Bool(True)
#    bond_color = RGBColor((1, 0.6, 0))
    bond_color = RGBColor((0, 0, 0))
    bond_radius = Float(0.04)
    background_color = RGBColor((1, 1, 1)) # Not accesible by GUI
    # Symmetry parameters
    show_symm_atoms = Bool(True)
    show_symm_bonds = Bool(True)
    label_symm_atoms = Bool(False)
    crop_mol = Bool(True) # Crop atom, bonds and labels outside crop_range*datarange?
    crop_range = Float(100.0)

################################################################################


    # Main panel 
    main = View(
                VGroup(
                    Item('filename', label = 'GRD file', style = 'readonly'),
                    HGroup(
                        Item('dim', label = 'Dimension', style = 'readonly'),
                        Item('func', label = '| Function', style = 'readonly'),
                        ),
                    VGroup(
                        HGroup(
                            Item('show_pos', label = 'Show positive isosurface?'),
                            Item('pos_surf', label = '+Isosurface value:', width = -70, enabled_when = 'show_pos'),
                            Item('pos_color', label = '+Isosurface color:', width = -100, enabled_when = 'show_pos')  
                            ),
                        HGroup(
                            Item('show_neg', label = 'Show negative isosurface?'),
                            Item('neg_surf', label = '-Isourface value:', width = -70, enabled_when = 'show_neg'),
                            Item('neg_color', label = '-Isosurface color:', width = -100, enabled_when = 'show_neg')    
                            )
                        ),
                    VGroup(
                        HGroup(
                            Item('opacity', label = 'Isosurface opacity:', width = -50),
                            Item(label = '                                                        ')
                            ),
                        HGroup(
                            Item('show_atoms', label = 'Show atoms?'),
                            Item('atom_size', label = 'Atom_size:', width = -50, enabled_when = 'show_atoms')
                            )
                        ),
                    HGroup(
                        Item('label_atoms', label = 'Label atoms?', enabled_when = 'show_atoms'),
                        Item('label_color', label = 'Label color:', width = -100, enabled_when = 'label_atoms'),
                        Item('label_size', label = 'Label size:', width = -50, enabled_when = 'label_atoms')
                        ),
                    HGroup(
                        Item('show_bonds', label = 'Show bonds?'),
                        Item('bond_color', label = 'Bond color:', width = -100, enabled_when = 'show_bonds'),
                        Item('bond_radius', label = 'Bond radius:', width = -50, enabled_when = 'show_bonds')
                        ),
                    VGroup(
                        HGroup(
                            Item('show_symm_atoms', label = '    Show sym.gen. atoms?', enabled_when = 'show_atoms'),
                            Item('label_symm_atoms', label = '    Label sym.gen. atoms?', enabled_when = 'label_atoms'),
                            Item('show_symm_bonds', label = '    Show sym.gen. bonds?', enabled_when = 'show_bonds')
                            ),
                        HGroup(
                            Item('crop_mol', label = 'Crop molecule at'),
                            Item('crop_range', label = '% of data range.', width = -50, enabled_when = 'crop_mol'),
                            Item(label = '                                                        '),
                            show_left = False
                            ),  
                        ),
                    ),
                buttons = OKCancelButtons, resizable=True, title   = "XDPlotter i3d: " + str(version)                           
            )

################################################################################                
                                                
    @on_trait_change('show_atoms')
    def _show_atoms_changed(self):
        """If show_atoms is deselected, deselect labels and sym. gen. atoms"""
        if self.show_atoms == False:
            self.label_atoms = False
            self.show_symm_atoms = False
            self.label_symm_atoms = False
            
    @on_trait_change('label_atoms')
    def _label_atoms_changed(self):
        """ If show_atoms is deselected, deselect symm. gen. atoms"""
        if self.label_atoms == False:
            self.label_symm_atoms = False

    @on_trait_change('show_bonds')
    def _show_bonds_changed(self):
        """ If show_bonds is deselected, deselect symm. gen. bonds"""
        if self.show_bonds == False:
            self.show_symm_bonds = False
                                            
################################################################################

    def plot(self):
        """
        Plot the isosurfaces
        """
        # Atom colors and covalent radii
        a_color = atomdata.get_atom_color()
        cov_r = atomdata.get_covalent_radii()
        a_color, cov_r = atomdata.change_atom_properties(a_color, cov_r)
        # Dimensions for plot, and atoms
        dim, func, x, y, z, atoms, data = xd.read_xdgrd(self.filename)
        # Clean atom list and correct for offset (origin)
        atoms = xd.clean_atoms(atoms, x[1], y[1], z[1])
        if self.crop_mol:
            atoms = xd.crop_atoms3d(atoms, self.crop_range, x, y, z)
        # Calculate coordinate system
        xgrid, ygrid, zgrid = xd.plot_area(x, y, z)
        # Plot the contours
#        mlab.options.offscreen = False
        plt = mlab.figure(bgcolor = self.background_color, size=(600, 600))
        if self.show_pos:
            mlab.contour3d(xgrid, ygrid, zgrid, data, contours = \
                [self.pos_surf], opacity = self.opacity, color = self.pos_color)
        if self.show_neg:
            mlab.contour3d(xgrid, ygrid, zgrid, data, contours = \
                [self.neg_surf], opacity = self.opacity, color = self.neg_color)

        if self.show_bonds:
            # Itterate over all pairs of atoms in the file
            for pair in itertools.combinations(atoms, 2):
                x1 = float(pair[0][1])
                x2 = float(pair[1][1])
                y1 = float(pair[0][2])
                y2 = float(pair[1][2])
                z1 = float(pair[0][3])
                z2 = float(pair[1][3])
                dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
                a_type1 = re.sub('.*?(_)', '', pair[0][0]).split('(')[0].capitalize()
                a_type2 = re.sub('.*?(_)', '', pair[1][0]).split('(')[0].capitalize()   
                # Plot if distance is smaller than sum of covalent radii         
                if dist <= 2.2:
                    if self.show_symm_bonds:
                        mlab.plot3d([x1,x2], [y1, y2], [z1, z2], \
                                tube_radius = self.bond_radius, color = \
                                self.bond_color, tube_sides = 24)
                    if not self.show_symm_bonds:
                        if pair[0][0][0] != 'X' and pair[1][0][0] != 'X':
                            mlab.plot3d([x1,x2], [y1, y2], [z1, z2], \
                                tube_radius = self.bond_radius, color = \
                                self.bond_color, tube_sides = 24)              
                                                
        if self.show_atoms:
            for atom in atoms:
                if self.show_symm_atoms:
                    a_type = re.sub('.*?(_)', '', atom[0]).split('(')[0].capitalize()
                    mlab.points3d(float(atom[1]), float(atom[2]), \
                                  float(atom[3]), scale_factor = \
                                  self.atom_size, color = rgb(a_color.get(a_type, \
                                  '#ffffff')), resolution = 24)
                if not self.show_symm_atoms:
                        if atom[0][0] != 'X':
                            a_type = re.sub('.*?(_)', '', atom[0]).split('(')[0].capitalize()
                            mlab.points3d(float(atom[1]), float(atom[2]), \
                                          float(atom[3]), scale_factor = \
                                          self.atom_size, color = \
                                          rgb(a_color.get(a_type, '#ffffff')), \
                                          resolution = 24)
            
        if self.label_atoms:    
            for atom in atoms:
                if self.label_symm_atoms:
                    mlab.text3d(float(atom[1])+self.atom_size/3.0, \
                                float(atom[2])+self.atom_size/3.0, \
                                float(atom[3])+self.atom_size/3.0,\
                                atom[0], scale = self.label_size, color = \
                                self.label_color, orient_to_camera = True)
                if not self.label_symm_atoms:
                        if atom[0][0] != 'X':
                            mlab.text3d(float(atom[1])+self.atom_size/3.0, \
                                        float(atom[2])+self.atom_size/3.0, \
                                        float(atom[3])+self.atom_size/3.0, \
                                        atom[0], scale = self.label_size, \
                                        color = self.label_color, \
                                        orient_to_camera = True)                                           
        mlab.show()
        
################################################################################            
            
if __name__ == '__main__':                                     
    param = param_3d_isosurface()
    #param.set(filename = sys.argv[-1])
    param.trait_set(filename = '/home/matt/Work/cosph_cov (copy)/xd_d2rho.grd')
    param.trait_set(dim = '3')
    param.trait_set(func = 'FOU')
    param.configure_traits(view = 'main', kind = 'modal')
    param.plot()
