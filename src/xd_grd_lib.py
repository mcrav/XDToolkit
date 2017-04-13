

"""
Functions used by XDPlotter (h3d, i3d, c2d): A program to plot 2D contour maps,
3D height fields or 3D isosurfaces from XD2006 grd files.
Mads Joergensen, 2013, Aarhus University
Version tracking: Describe changes and update version number below section. 
0.2     Changed crop_atoms3d() to requires a another parameter: crop_range.
        (October 4th 2013)
0.3     A bug occured in the calculation of the mgrid - some values lead to 
        a grid beeing 1 value larger in each dimension. I think this was due to 
        a numerical issue and by rounding the increment to 3 decimals it seems 
        to work (January 6th 2014)
0.4     Changed clean_atom function to change atom names from e.g FE(1) to 
        Fe(1) and some symmetry generated names. All two letter atomic symbols 
        will be changed. Fixed issue with missing function in some grd files 
        from e.g. ADDGRID. Simplified the expression to calculate coordinates 
        for plotting. For 3D four the data are apparently listed differently 
        than from xdprop - fix for this has been implemented(January 27th 2014)
"""
version = 0.4

################################################################################

import numpy as np
import copy

################################################################################

def get_version():
    "Version tracking"""
    return "xd_grd_lib: " + str(version)

def read_xdgrd(file):
    """
    Read grd file from XD2006
    Returns dimension, function, number of points in xyz, origin and dimensions,
    min and max, atoms (label, x, y, z) and an numpy array with the data
    """
    grd_file = open(file, 'r')

    # Read header and save dimension of file and function
    line = grd_file.readline()
    dim = int(line[0])
    line = grd_file.readline()
    try:
        func = line.split()[-1]
    except IndexError: # If no function is listed, e.g. from ADDGRID
        func = 'NONE'
    line = grd_file.readline()
    while line[0:6] != '! Grid':
        line = grd_file.readline()
    # Read dimensions of grid
    nx, ny, nz = grd_file.readline().split()
    xo, yo, zo = grd_file.readline().split()
    xdim, ydim, zdim = grd_file.readline().split()
    x = (int(nx), float(xo), float(xdim), float(xo)-float(xdim)/2, \
         float(xo)+float(xdim)/2)
    y = (int(ny), float(yo), float(ydim), float(yo)-float(ydim)/2, \
         float(yo)+float(ydim)/2)
    z = (int(nz), float(zo), float(zdim), float(zo)-float(zdim)/2, \
         float(zo)+float(zdim)/2)
    # Store number of atoms
    line = grd_file.readline()
    n_atoms = int(grd_file.readline().split()[0])
    # Read atoms and collect in list of lists
    atoms = []
    for i in range(n_atoms):
        atoms.append(grd_file.readline().split()[0:4])
    # Read until data begins
    while line[0:8] != '! Values':
        line = grd_file.readline()
    # Read rest of file into a long list of points
    data = []
    lines = grd_file.readlines()
    for line in lines:
        for point in line.split():
            data.append(float(point))
    # Convert to np array and reshape
    if dim == 2:
        data = np.array(data, dtype = 'float32')
        data = data.reshape(y[0], x[0])
#ORIGINAL        data = data.reshape(x[0], y[0])
        data = np.swapaxes(data, 0, 1)
    elif dim == 3 and func == 'FOU':
        data = np.array(data, dtype = 'float32')
        data = data.reshape(x[0], y[0], z[0])
        # DO NOT SWAP AXES!
    else:
        data = np.array(data, dtype = 'float32')
        data = data.reshape(x[0], y[0], z[0])
        # Swap x and z, more intuitive with x, y, z
        data = np.swapaxes(data, 0, 2)

    return dim, func, x, y, z, atoms, data

def clean_atoms(atoms, xo, yo, zo):
    """
    Converts the x, y and z cooridnates from strings to floats and corrects for
    the origin offset. Furthermore changes two letter atomic symbols from XX to 
    Xx.
    """
    for atom in atoms:
        if atom[0].find('___') > 0: #If APPLY symm used in XDFOUR e.g. C(1)___1 to X1_C(1)
            atom[0] = 'X'+atom[0].split('___')[1]+'_'+atom[0].split('___')[0]
        elif atom[0].find('__') > 0:
            atom[0] = 'X'+atom[0].split('__')[1]+'_'+atom[0].split('__')[0]
            # Fix e.g. FE(4B) to Fe(4B)
        if len(atom[0].split('_')) == 1 and len(atom[0].split('(')[0]) > 1:
            temp = atom[0].split('(')
            temp[0] = temp[0].capitalize()
            atom[0] = '('.join(temp)
        # Fix symm generated e.g. X1_FE(4B) to X1_Fe(4B)
        elif len(atom[0].split('_')) > 1 and len(atom[0].split('(')[0]) > 1:
            temp = atom[0].split('_')
            temp[1] = temp[1].split('(')
            temp[1][1] = temp[1][1].capitalize()
            temp[1] = '('.join(temp[1])
            atom[0] = '_'.join(temp)
        atom[1] = float(atom[1])+xo
        atom[2] = float(atom[2])+yo
        atom[3] = float(atom[3])+zo
    return atoms 
            
def plot_area(x, y, z):
    """
    Sets up a 2D or 3D grid in the right dimension.
    x, y, z are tuples containing number of points, origin and dimensions.
    """
    nx, xo, xdim, xmin, xmax = x
    nx = int(nx)
    xo = float(xo)
    xdim = float(xdim)
    ny, yo, ydim, ymin, ymax = y
    ny = int(ny)
    yo = float(yo)
    ydim = float(ydim)
    nz, zo, zdim, zmin, zmax = z
    nz = int(nz)
    zo = float(zo)
    zdim = float(zdim)
    
    if nz == 1:
        coord = np.mgrid[xmin+xdim/(2*nx): xmax+xdim/(2*nx): xdim/nx,\
                   ymin+ydim/(2*ny): ymax+ydim/(2*ny): ydim/ny]
    else:
        coord = np.mgrid[xmin+xdim/(2*nx): xmax+xdim/(2*nx): xdim/nx,\
                   ymin+ydim/(2*ny): ymax+ydim/(2*ny): ydim/ny,\
                   zmin+zdim/(2*nz): zmax+zdim/(2*nz): zdim/nz]

    return coord

def linear_contour(step, pos_lim, neg_lim):
    """
    Return two lists of equidistant contour levels: positive and negative.
    """
    pos_contours = np.arange(step, pos_lim+0.1*step, step)
    neg_contours = np.arange(neg_lim, -1*step+0.1*step, step)
    
    return pos_contours, neg_contours

def log_contour(a, b):
    """
    Returns log contours [a]*10^[b] for elements in the lists a and b.
    """
    pos_contours = []
    neg_contours = []
    for i in range(len(a)):
        for j in range(len(b)):
            pos_contours.append(float(a[i])*(10**float(b[j])))
            neg_contours.append(-1*float(a[i])*(10**float(b[j])))
    pos_contours.sort()
    neg_contours.sort(reverse = True)
    pos_contours = np.array(pos_contours)
    neg_contours = np.array(neg_contours)
    
    return pos_contours, neg_contours

def crop_atoms3d(atoms, crop_range, x, y, z):
    """
    Create a new atoms list containing only atoms within the data range times 
    crop_range. Crop_range supplied as a percentage. 
    Should ONLY be used for 3D plots
    """
    cropped_atoms = []
    cr = float(crop_range)/100.0
    for atom in atoms:
        if(atom[1] > x[3]*cr and atom[1] < x[4]*cr and
           atom[2] > y[3]*cr and atom[2] < y[4]*cr and
           atom[3] > z[3]*cr and atom[3] < z[4]*cr):
            cropped_atoms.append(atom)
    return cropped_atoms