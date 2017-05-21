# -*- coding: utf-8 -*-
"""
QuickPlot: A program to plot 2D contour maps from XD2006 grd files.
The script needs libraries, 'xd_grd_lib.py' and 'atom_dictionay.py' that should
be put in the python path or in the same folder as this file. The modules
matplotlib and numpy are also required.
The script has been tested with:
OSX: python:2.7.6, numpy:1.8.1, matplotlib:1.4.2
The script can be run from the command line with no arguments, one argument or 
two arguments:
- 0 args: If 'xd_fou.grd' is found it is plotted using default parameters. A 
    parameter file 'qp.par' is written.
- 1 arg: The first argument is the name of the grd-file to be plotted. The plot
    will be made using default parmetes and a parameter file 'qp.par' is 
    written.
- 2 args: First arg must be the grd-file and the second arg the parameter file 
    to be used.
    
If a 3D file is supplied the program exits.
The parameter file can be edited to change the look of the plots. Remember to
save the file under a new name as the program will overwrite it if no parameter
file is specified.
A file 'change_atom_properties.txt' is written and can be used to configure 
non-standard cavalent radii and atom colors. Do NOT change the name of this 
file, it will not be overwritten.
The arguments can be files in the cwd or the absolute path, i.e. under Linux 
it is possible to make bash functions with different parameter files e.g.:
	function defplot(){
    		python /ABSPATH/quickplot.py $1 /ABSPATH/def.par
	}
or alaternatively without arguments:   
	function defplot(){
    		python /ABSPATH/quickplot.py xd_def.grd /ABSPATH/def.par
	}
In Windows it is possible to create a batch file located in the search path. Here is an example (quickplot.bat):
	@echo off
	python C:\QuickPlot\quickplot.py %*
	pause
This will pass optional arguments to python. It is also possible to have 
customized parameter files (def_plot.bat):
	@echo off
	python C:\QuickPlot\quickplot.py %1 C:\QuickPlot\par_files\def.par
	pause
Mads Ry Jørgensen, 2015, Aarhus University
Version tracking: Describe changes and update version number below section. 
Colors changed from tuples to single characters
"""
version = 0.2

################################################################################
import os
import sys

import itertools
import re
import configparser

import numpy as np
import matplotlib.pyplot as plt

import xd_grd_lib as xd
import atom_dictionary as atomdata
################################################################################

def center_text(text, width):
    """
    Print text to screen centered in the specified width
    """
    n = int((width - len(text))/2)
    print (n*' ' + text)
    
def print_version():
    """
    Print version information to screen.
    """
    width = 80
    print (width*'#')
    print ('\n')
    center_text("QuickPlot version " + str(version), width)
    center_text("A program to plot XD2006 grd files.", width)
    center_text("Mads Ry Jørgensen, 2015, Aarhus University.\n", width)
    center_text(xd.get_version()+', '+atomdata.get_version()+'\n',width)
    print (width*'#' + "\n")
    sys.stdout.flush()
    
def parse_int_list(par_list):
    """
    Formates the raw input from ConfigParser (e.g. '[1, 2, 4, 8]') to 
    [1, 2, 4, 8] 
    """
    return [int(i) for i in par_list.strip('[').strip(']').split(',')]

def parse_float_tuble(par_tup):
    """
    Formates the raw input from ConfigParser (e.g. '(0, 0, 1)') to 
    (0.0, 0.0, 1.0) 
    """
    return tuple([float(i) for i in par_tup.strip('(').strip(')').split(',')])

def create_qp_par():
    """
    Creates a qp.par file with default parameters 
    """
    qp_par = open('qp.par','w')
    text = """# Feel free to remove items to have a shorter file, but do NOT
# remove the section headers (with [])!
# Remember to rename file if edited!
[contours]
#Linear contours (FOU, DEF)
use_lin_contour = True
pos_lim = 2.0
neg_lim = -2.0
step = 0.05
#Log contours (D2R) [base]*10**[exponent]
base = [1, 2, 4, 8]
exponent = [-2, -1, 0, 1, 2, 3, 4]
# Show zero contour?
zero_cont = False
[lines]
#Colors and contour line styles:
#b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
pos_color = b
pos_line = solid
neg_color = r
neg_line = dashed
zero_color = k 
zero_line = dotted
cont_line_width = 0.8
[atoms]
atom_size = 10
# Swich on/off symmetry generated atoms
show_symm_atoms = True
# Cut-off for atoms out of plane:
atom_cut = 0.2
[bonds]
show_bonds = False
bond_color = k
bond_thickness = 2
#Bonds between symmetry generated atoms
show_symm_bonds = False
[labels]
label_atoms = True
# Atoms has to be plotted to show label 
label_symm_atoms = False
label_color = k
label_size = 15
label_x_offset = 0.1
label_y_offset = 0.1
[save]
# Save file as: png, eps, pdf
save_as = png
"""
    qp_par.write(text)
    qp_par.close()
    print (os.path.join(os.getcwd(),"qp.par") + " created.\n")

################################################################################
def makeResMap(filename, saveFilename = None):
    
    # Set Default parameters  
    fig = plt.Figure()
    fig.add_subplot(1,1,1)
    
    #[contours]
    #Linear contours (FOU, DEF)
    use_lin_contour = True
    pos_lim = 2.0
    neg_lim = -2.0
    step = 0.05
    #Log contours (D2R) [base]*10**[exponent]
    base = [1, 2, 4, 8]
    exponent = [-2, -1, 0, 1, 2, 3, 4]
    # Show zero contour?
    zero_cont = False
    
    #[lines]
    #Colors and contour line styles:
    pos_color = 'b'
    pos_line = 'solid'
    neg_color = 'r'
    neg_line = 'dashed'
    zero_color = 'k' 
    zero_line = 'dotted'
    cont_line_width = 0.8
    
    #[atoms]
    atom_size = 10
    # Swich on/off symmetry generated atoms
    show_symm_atoms = True
    # Cut-off for atoms out of plane:
    atom_cut = 0.2
    
    #[bonds]
    show_bonds = False
    bond_color = 'k'
    bond_thickness = 2
    #Bonds between symmetry generated atoms
    show_symm_bonds = False
    
    #[labels]
    label_atoms = True
    # Atoms has to be plotted to show label 
    label_symm_atoms = False 
    label_color = 'k'
    label_size = 15
    label_x_offset = 0.1
    label_y_offset = 0.1
    
    #[save]
    # Save file as: 'png', 'eps', 'pdf'
    save_as = 'png'
    
    ################################################################################
    print_version()
    
    qp_par = None
    create_qp_par()
#    if len(sys.argv) <= 1:
#        if os.path.isfile('xd_fou.grd'):
#            print ("No grd file given. Will use: 'xd_fou.grd'\n")
#            filename = 'xd_fou.grd'
#            print ("No parameter file given.\nWill create qp.par and use standard parameters.\n")
#            create_qp_par()
#        else:
#            print ("No grd file given.\nPlease specify grd file and the optional parameter file!\n")
#            sys.exit(0)
#    
#    if len(sys.argv) == 2:
#        if os.path.isfile(sys.argv[1]):
#            print (sys.argv[1] + " found! No parameter file given.\nWill create qp.par and use standard parameters.\n")
#            filename = sys.argv[1]
#            qp_par = None
#            create_qp_par()
#        else:
#            print (sys.argv[1] + " not found!\nPlease specify grd file and the optional parameter file!\n")
#            sys.exit(0)
#    
#    if len(sys.argv) >= 3:
#        if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[1]):
#            print (sys.argv[1] + " and " + sys.argv[2] + " found!\n")
#            filename = sys.argv[1]
#            qp_par = sys.argv[2]
#        elif os.path.isfile(sys.argv[1]) and not os.path.isfile(sys.argv[2]):
#            print (sys.argv[1] + "found!\n")
#            print (sys.argv[2] + " not found!\nWill create qp.par and use standard parameters.\n")
#            create_qp_par()
#        else:
#            print (sys.argv[1] + " not found!\nPlease specify grd file and the optional parameter file!\n")
#            sys.exit(0)
    
    ################################################################################
    # Retrive atom colors and covalent radii
    a_color = atomdata.get_atom_color()
    cov_r = atomdata.get_covalent_radii()
    a_color, cov_r = atomdata.change_atom_properties(a_color, cov_r)
    
    # Update parameters based on parameter file
    if qp_par:
        print ("Reading parameter values from " + qp_par + "...")
        config = configparser.RawConfigParser()
        config.read(qp_par)
        #contours
        if config.has_option('contours','use_lin_contour'):
            use_lin_contour = config.getboolean('contours','use_lin_contour')
        if config.has_option('contours','pos_lim'):
            pos_lim = config.getfloat('contours','pos_lim')
        if config.has_option('contours','neg_lim'):
            neg_lim = config.getfloat('contours','neg_lim')
        if config.has_option('contours','step'):
            step = config.getfloat('contours','step')
        if config.has_option('contours','base'):
            base = parse_int_list(config.get('contours','base'))
        if config.has_option('contours','exponent'):
            exponent = parse_int_list(config.get('contours','exponent'))
        if config.has_option('contours','zero_cont'):
            zero_cont = config.getboolean('contours','zero_cont')
        #lines
        if config.has_option('lines','pos_color'):
            pos_color = config.get('lines','pos_color')
        if config.has_option('lines','pos_line'):
            pos_line = config.get('lines','pos_line')
        if config.has_option('lines','neg_color'):
            neg_color = config.get('lines','neg_color')
        if config.has_option('lines','neg_line'):
            neg_line = config.get('lines','neg_line')
        if config.has_option('lines','zero_color'):
            zero_color = config.get('lines','zero_color')
        if config.has_option('lines','zero_line'):
            zero_line = config.get('lines','zero_line')
        if config.has_option('lines','cont_line_width'):
            cont_line_width = config.getfloat('lines','cont_line_width')
        #atoms
        if config.has_option('atoms','atom_size'):
            atom_size = config.getfloat('atoms','atom_size')
        if config.has_option('atoms','show_symm_atoms'):
            show_symm_atoms = config.getboolean('atoms','show_symm_atoms')
        if config.has_option('atoms','atom_cut'):
            atom_cut = config.getfloat('atoms','atom_cut')
        #bonds
        if config.has_option('bonds','show_bonds'):
            show_bonds = config.getboolean('bonds','show_bonds')
        if config.has_option('bonds','bond_color'):
            bond_color = config.get('bonds','bond_color')
        if config.has_option('bonds','bond_thickness'):
            bond_thickness = config.getfloat('bonds','bond_thickness')
        if config.has_option('bonds','show_symm_bonds'):
            show_symm_bonds = config.getboolean('bonds','show_symm_bonds')
        #labels
        if config.has_option('labels','label_atoms'):
            label_atoms = config.getboolean('labels','label_atoms')
        if config.has_option('labels','label_symm_atoms'):
            label_symm_atoms = config.getboolean('labels','label_symm_atoms')
        if config.has_option('labels','label_color'):
            label_color = config.get('labels','label_color')
        if config.has_option('labels','label_size'):
            label_size = config.getfloat('labels','label_size')
        if config.has_option('labels','label_x_offset'):
            label_x_offset = config.getfloat('labels','label_x_offset')
        if config.has_option('labels','label_y_offset'):
            label_y_offset = config.getfloat('labels','label_y_offset')
        #save
        if config.has_option('save','save_as'):
            save_as = config.get('save','save_as')
    
    ################################################################################
    plt.subplot(111)
    fig = plt.figure(facecolor = '#ddd9d5', figsize=(9,7), dpi=300)             
    # Dimensions for plot, and atoms
    dim, func, x, y, z, atoms, data = xd.read_xdgrd(filename)
    
    atoms = xd.clean_atoms(atoms, x[1], y[1], z[1])

    # Check dimensionality of plot
    if dim != 2:
        print ("Grid is not 2 dimensional. Please specify a 2D grid!")
        sys.exit()
    
    # Contours

    if filename == 'xd_fou.grd':
        pos_contours, neg_contours = xd.linear_contour(\
                                        step, pos_lim, neg_lim)
    else:
        pos_contours, neg_contours = xd.log_contour(\
                                        base, exponent)

    xgrid, ygrid = xd.plot_area(x, y, z)

    if (xgrid.shape != data.shape):
        xgrid = xgrid[:-1, :-1]
        ygrid = ygrid[:-1, :-1]
#    xgrid = xgrid[1:, 1:]
#    ygrid = ygrid[1:, 1:]
    
    # Plot contours
    plt.contour(xgrid, ygrid, data, levels = pos_contours, colors = \
                pos_color, linestyles = pos_line, linewidths = \
                cont_line_width)
    plt.contour(xgrid, ygrid, data, levels = neg_contours, colors = \
                neg_color, linestyles = neg_line, linewidths = \
                cont_line_width)
    if zero_cont:
        plt.contour(xgrid, ygrid, data, levels = [0], colors = \
                    zero_color, linestyles = zero_line, \
                    linewidths = cont_line_width)
    
    plt.axes().set_aspect('equal')
    
    if show_bonds:                
        # Itterate over all pairs of atoms in the file
        for pair in itertools.combinations(atoms, 2):
            l1, x1, y1, z1 = pair[0]
            l1, x2, y2, z2 = pair[1]
            dist = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            a_type1 = re.sub('.*?(_)', '', pair[0][0]).split('(')[0]
            a_type2 = re.sub('.*?(_)', '', pair[1][0]).split('(')[0]
            if abs(float(pair[0][3])) <= atom_cut: # Include only atoms near plane
                    if abs(float(pair[1][3])) <= atom_cut:
                    # Plot if distance is smaller than sum of covalent radii         
                        if dist <= cov_r.get(a_type1, 0) + \
                                    cov_r.get(a_type2, 0):
                            if show_symm_bonds: # Plot all bonds
                                plt.plot([x1,x2], [y1, y2], linewidth = \
                                            bond_thickness, color = \
                                            bond_color)
                            # Plot only bonds in asym unit
                            if not show_symm_bonds: 
                                if pair[0][0][0] != 'X' and pair[1][0][0] \
                                    != 'X':
                                    plt.plot([x1,x2], [y1, y2], linewidth =\
                                                bond_thickness, color = \
                                                bond_color)
    
    # Plot atoms
    
    for atom in atoms:
        #JUST THERE FOR CUSTOM LABELS BECAUSE OF PROJECT WHERE SHELX.INS LABELS DON'T MATCH XD.MAS LABELS
#        if atom[0]=='S(007)':
#            atom[0] = 'S(4)'
#        elif atom[0]=='CO(01)':
#            atom[0]='CO(1)'
#        elif atom[0]=='C(01A)':
#            atom[0]='C(19)'
            
    
        if abs(float(atom[3])) <= atom_cut:
            if show_symm_atoms: # Show all atoms
                a_type = re.sub('.*?(_)', '', atom[0]).split('(')[0]
                plt.plot(atom[1], atom[2], marker='o', \
                        mec = (0,0,0), mew = bond_thickness, \
                        mfc = a_color.get(a_type, (0, 0, 0)), ms = \
                                atom_size)
                if label_atoms:
                    if label_symm_atoms: # Label all atoms
                        plt.text(atom[1]+label_x_offset, atom[2] + \
                                    label_y_offset, atom[0], fontsize = \
                                    label_size, color = label_color,\
                                    clip_on=True)
                    if not label_symm_atoms: # Only label asym unit
                        if atom[0][0] != 'X':
                            plt.text(atom[1]+label_x_offset, atom[2] +\
                                        label_y_offset, atom[0], fontsize\
                                        = label_size, color = \
                                        label_color, clip_on=True)
            if not show_symm_atoms: # Show only asym unit
                    if atom[0][0] != 'X':
                        a_type = re.sub('.*?(_)', '', atom[0]).split('(')[0]
                        plt.plot(atom[1], atom[2], marker='o', mec = \
                                    (0,0,0), mew = bond_thickness, mfc = \
                                    a_color.get(a_type, (0, 0, 0)), ms = \
                                    atom_size)
                        if label_atoms: # Label all atoms (in asym unit)
                            plt.text(atom[1]+label_x_offset, atom[2] +\
                                        label_y_offset, atom[0], fontsize\
                                        = label_size, color = \
                                        label_color, clip_on=True)
                                                
    plt.axis([x[3], x[4], y[3], y[4]])
    plt.xticks([])
    plt.yticks([])   
    if saveFilename:   
        print('saving') 
        print(saveFilename)
        print(save_as)
        
        plt.savefig(saveFilename, bbox_inches='tight', pad_inches=0, dpi = 300)
        print('saved')
    #print ('%s_%s%s%s.%s saved in %s/' % (func, atoms[0][0], atoms[1][0], \
     #     atoms[2][0], save_as, os.getcwd()))
    plt.show()

    return fig
    
#makeResMap('/home/matt/Work/Co(SPh)4 Final Files/xd_d2rho.grd')
