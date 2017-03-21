
# -*- coding: utf-8 -*-
"""
Library of atom colors and covalent radii used by XDPlotter (h3d, i3d, c2d): A 
program to plot 2D contour maps, 3D height fields or 3D isosurfaces from XD2006 
grd files.
Mads Ry Jørgensen, 2015, Aarhus University
Version tracking: Describe changes and update version number below section.
0.2     Added change_atom_properties() that reads a file if pressent and updates
        the dictionaries accordingly
0.3     Changing all colors from tuples to HEX codes
"""
version = 0.3

################################################################################

import os
import sys

################################################################################

def get_version():
    "Version tracking"""
    return "atom_dictionary: " + str(version)

def change_atom_properties(a_color, cov_r):
    """
    Reads the file 'change_atom_properties.txt' if present in the curent 
    folder and adds the new atomic data to the two atom dictionaries. If the 
    file is not available a file with syntax examples is saved for easy 
    reference.
    """
    if os.path.isfile('change_atom_properties.txt'):
        new_color = {} # Initialize new dictionaries
        new_radii ={}
        atom_changes = open('change_atom_properties.txt','r')
        line = atom_changes.readline()
        while line[0] == '#' or line == '\n': # Read header
            line = atom_changes.readline()
        if line[0:11] == "start_color": # Start color definitions
            entry = atom_changes.readline()
            while entry[0:9] != "end_color" and entry != '': # Read color definitions
                symbol = entry.split('=')[0].strip(' ') # Extract atomic symbol
                color = '#'+entry.split('=')[1].strip(' ')[0:6] # Color as a string
                new_color[symbol] = color
                entry = atom_changes.readline()
            if len(new_color) > 0:
                print ("The following atoms have non standard colors: ", new_color)
            a_color.update(new_color) # Update new values to original dictionary
        line = atom_changes.readline()      
        while line[0] == '#' or line == '\n':
            line = atom_changes.readline()
        if line[0:11] == "start_radii": # Start radii definitions
            entry = atom_changes.readline()
            while entry[0:9] != "end_radii" and entry != '': # Read radii definitions
                symbol = entry.split('=')[0].strip(' ')
                radius = eval(entry.split('=')[1])
                if type(radius) == float: # Add only if radius is a float
                    new_radii[symbol] = radius
                entry = atom_changes.readline()
            if len(new_radii) > 0:
                print ("The following atoms have non standard radii: ", new_radii)
            cov_r.update(new_radii) # Update new values to original dictionary
        atom_changes.close()
    else: # If file is not available create an example file
        example = open('change_atom_properties.txt','w')
        text = """# Edit this file to change atomic color and/or covalent radii for XDPlotter.
# All lines starting with '#' are treated as comments.
#
# Colors: Section start with 'start_color' and ends with 'end_color'. 
# One line per entry specifying first the atomic symbol as a string: e.g. Fe
# equal to (=) color as an HEX code (e.g. ff0000 for red)
# DO NOT USE '#' or ''s!
# Fe = ff0000
#
# Radii: Section start with 'start_radius' and ends with 'end_radius'.
# One line per entry specifying first the atomic symbol as a string: e.g. Fe
# equal to (=) covalent radius as a float (in Aa), e.g. 1.25:
# Fe = 1.25
#
# No empty lines between start and end statements!
start_color
end_color
start_radii
end_radii
"""
        example.write(text) # Write file
        example.close()
    sys.stdout.flush()
    return a_color, cov_r # Return (updated) dictionaries

################################################################################

def get_covalent_radii(): # Dictionary to store covalent radii
    """
    Returns a dictionary with atom symbols as keys (strings) and covalent
    radius as a float.
    """
    cov_r = {}
    cov_r['H'] = 0.38
    cov_r['He'] = 0.32
    cov_r['Li'] = 1.34
    cov_r['Be'] = 0.9
    cov_r['B'] = 0.82
    cov_r['C'] = 0.77
    cov_r['N'] = 0.75
    cov_r['O'] = 0.73
    cov_r['F'] = 0.71
    cov_r['Ne'] = 0.69
    cov_r['Na'] = 1.54
    cov_r['Mg'] = 1.3
    cov_r['Al'] = 1.18
    cov_r['Si'] = 1.11
    cov_r['P'] = 1.06
    cov_r['S'] = 1.02
    cov_r['Cl'] = 0.99
    cov_r['Ar'] = 0.97
    cov_r['K'] = 1.96
    cov_r['Ca'] = 1.74
    cov_r['Sc'] = 1.44
    cov_r['Ti'] = 1.36
    cov_r['V'] = 1.25
    cov_r['Cr'] = 1.27
    cov_r['Mn'] = 1.39
    cov_r['Fe'] = 1.25
    cov_r['Co'] = 1.26
    cov_r['Ni'] = 1.21
    cov_r['Cu'] = 1.38
    cov_r['Zn'] = 1.31
    cov_r['Ga'] = 1.26
    cov_r['Ge'] = 1.22
    cov_r['As'] = 1.19
    cov_r['Se'] = 1.16
    cov_r['Br'] = 1.14
    cov_r['Kr'] = 1.1
    cov_r['Rb'] = 2.11
    cov_r['Sr'] = 1.92
    cov_r['Y'] = 1.62
    cov_r['Zr'] = 1.48
    cov_r['Nb'] = 1.37
    cov_r['Mo'] = 1.45
    cov_r['Tc'] = 1.56
    cov_r['Ru'] = 1.26
    cov_r['Rh'] = 1.35
    cov_r['Pd'] = 1.31
    cov_r['Ag'] = 1.53
    cov_r['Cd'] = 1.48
    cov_r['In'] = 1.44
    cov_r['Sn'] = 1.41
    cov_r['Sb'] = 1.38
    cov_r['Te'] = 1.35
    cov_r['I'] = 1.33
    cov_r['Xe'] = 1.3
    cov_r['Cs'] = 2.25
    cov_r['Ba'] = 1.98
    cov_r['La'] = 1.69
    cov_r['Lu'] = 1.6
    cov_r['Hf'] = 1.5
    cov_r['Ta'] = 1.38
    cov_r['W'] = 1.46
    cov_r['Re'] = 1.59
    cov_r['Os'] = 1.28
    cov_r['Ir'] = 1.37
    cov_r['Pt'] = 1.28
    cov_r['Au'] = 1.44
    cov_r['Hg'] = 1.49
    cov_r['Tl'] = 1.48
    cov_r['Pb'] = 1.47
    cov_r['Bi'] = 1.46
    cov_r['Rn'] = 1.45
    return cov_r
    
################################################################################

def get_atom_color(): # Dictionary to store atom color
    """
    Returns a dictionary with atom symbols as keys (strings) and color as an 
    rgb tuple.
    """
    a_color = {}
    a_color['H'] = '#ffffff'
    a_color['He'] = '#ff1a99'
    a_color['Li'] = '#999999'
    a_color['Be'] = '#999999'
    a_color['B'] = '#268c26'
    a_color['C'] = '#4c4c4c'
    a_color['N'] = '#0000ff'
    a_color['O'] = '#ff0000'
    a_color['F'] = '#00ff00'
    a_color['Ne'] = '#ff1a99'
    a_color['Na'] = '#999999'
    a_color['Mg'] = '#999999'
    a_color['Al'] = '#999999'        
    a_color['Si'] = '#738c99'    
    a_color['P'] = '#800000'
    a_color['S'] = '#ffff00' 
    a_color['Cl'] = '#cc1a8c'
    a_color['Ar'] = '#ff1a99'   
    a_color['K'] = '#999999'      
    a_color['Ca'] = '#999999'         
    a_color['Sc'] = '#999999'
    a_color['Ti'] = '#999999'
    a_color['V'] = '#999999'
    a_color['Cr'] = '#999999'
    a_color['Mn'] = '#8c400d'                                  
    a_color['Fe'] = '#ff4c00'
    a_color['Co'] = '#1a1a73'
    a_color['Ni'] = '#008000'           
    a_color['Cu'] = '#00ffff'
    a_color['Zn'] = '#e6e6e6'        
    a_color['Ga'] = '#999999'
    a_color['Ge'] = '#999999'
    a_color['As'] = '#999999'                   
    a_color['Se'] = '#999999'
    a_color['Br'] = '#8c0000'                                                                
    a_color['Kr'] = '#ff1a99'
    a_color['Rb'] = '#999999'    
    a_color['Sr'] = '#999999'
    a_color['Y'] = '#999999'
    a_color['Zr'] = '#999999'
    a_color['Nb'] = '#999999'
    a_color['Mo'] = '#999999'
    a_color['Tc'] = '#999999'
    a_color['Ru'] = '#999999'
    a_color['Rh'] = '#999999'
    a_color['Pd'] = '#999999'
    a_color['Ag'] = '#999999'
    a_color['Cd'] = '#999999'
    a_color['In'] = '#999999'
    a_color['Sn'] = '#999999'   
    a_color['Sb'] = '#999999'
    a_color['Te'] = '#999999'           
    a_color['I'] = '#800080'
    a_color['Xe'] = '#ff1a99'
    a_color['Cs'] = '#999999'
    a_color['Ba'] = '#999999'
    a_color['La'] = '#999999'
    a_color['Lu'] = '#999999'
    a_color['Hf'] = '#999999'
    a_color['Ta'] = '#999999'
    a_color['W'] = '#999999'
    a_color['Re'] = '#999999'
    a_color['Os'] = '#999999'
    a_color['Ir'] = '#999999'
    a_color['Pt'] = '#999999'
    a_color['Au'] = '#ffd900'
    a_color['Hg'] = '#999999'
    a_color['Tl'] = '#999999'
    a_color['Pb'] = '#999999'
    a_color['Bi'] = '#999999'
    a_color['Rn'] = '#999999'

    return a_color
    
################################################################################