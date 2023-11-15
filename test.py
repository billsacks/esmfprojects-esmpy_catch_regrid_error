#!/usr/bin/env python

import esmpy
import numpy as np

# Copied from ESMPy's grid_utilities.py:
def grid_create_from_coordinates(xcoords, ycoords, xcorners=False, ycorners=False, corners=False, domask=False, doarea=False, ctk=esmpy.TypeKind.R8):
    """
    Create a 2 dimensional Grid using the bounds of the x and y coordiantes.
    :param xcoords: The 1st dimension or 'x' coordinates at cell centers, as a Python list or numpy Array
    :param ycoords: The 2nd dimension or 'y' coordinates at cell centers, as a Python list or numpy Array
    :param xcorners: The 1st dimension or 'x' coordinates at cell corners, as a Python list or numpy Array
    :param ycorners: The 2nd dimension or 'y' coordinates at cell corners, as a Python list or numpy Array
    :param domask: boolean to determine whether to set an arbitrary mask or not
    :param doarea: boolean to determine whether to set an arbitrary area values or not
    :param ctk: the coordinate typekind
    :return: grid
    """
    [x, y] = [0, 1]

    # create a grid given the number of grid cells in each dimension, the center stagger location is allocated, the
    # Cartesian coordinate system and type of the coordinates are specified
    max_index = np.array([len(xcoords), len(ycoords)])
    grid = esmpy.Grid(max_index, staggerloc=[esmpy.StaggerLoc.CENTER], coord_sys=esmpy.CoordSys.CART, coord_typekind=ctk)

    # set the grid coordinates using numpy arrays, parallel case is handled using grid bounds
    gridXCenter = grid.get_coords(x)
    x_par = xcoords[grid.lower_bounds[esmpy.StaggerLoc.CENTER][x]:grid.upper_bounds[esmpy.StaggerLoc.CENTER][x]]
    gridXCenter[...] = x_par.reshape((x_par.size, 1))

    gridYCenter = grid.get_coords(y)
    y_par = ycoords[grid.lower_bounds[esmpy.StaggerLoc.CENTER][y]:grid.upper_bounds[esmpy.StaggerLoc.CENTER][y]]
    gridYCenter[...] = y_par.reshape((1, y_par.size))

    # create grid corners in a slightly different manner to account for the bounds format common in CF-like files
    if corners:
        grid.add_coords([esmpy.StaggerLoc.CORNER])
        lbx = grid.lower_bounds[esmpy.StaggerLoc.CORNER][x]
        ubx = grid.upper_bounds[esmpy.StaggerLoc.CORNER][x]
        lby = grid.lower_bounds[esmpy.StaggerLoc.CORNER][y]
        uby = grid.upper_bounds[esmpy.StaggerLoc.CORNER][y]

        gridXCorner = grid.get_coords(x, staggerloc=esmpy.StaggerLoc.CORNER)
        for i0 in range(ubx - lbx - 1):
            gridXCorner[i0, :] = xcorners[i0+lbx, 0]
        gridXCorner[i0 + 1, :] = xcorners[i0+lbx, 1]

        gridYCorner = grid.get_coords(y, staggerloc=esmpy.StaggerLoc.CORNER)
        for i1 in range(uby - lby - 1):
            gridYCorner[:, i1] = ycorners[i1+lby, 0]
        gridYCorner[:, i1 + 1] = ycorners[i1+lby, 1]

    # add an arbitrary mask
    if domask:
        mask = grid.add_item(esmpy.GridItem.MASK)
        mask[:] = 1
        mask[np.where((1.75 <= gridXCenter.any() < 2.25) &
                      (1.75 <= gridYCenter.any() < 2.25))] = 0

    # add arbitrary areas values
    if doarea:
        area = grid.add_item(esmpy.GridItem.AREA)
        area[:] = 5.0

    return grid


weight_file = 'weights.nc'
xcoords1 = np.array([1,2,3,4,5])
ycoords1 = np.array([1,2,3,4,5])
grid1 = grid_create_from_coordinates(xcoords1, ycoords1)

xcoords2 = np.array([1.5, 2.5, 3.5, 4.5])
ycoords2 = np.array([1.5, 2.5, 3.5, 4.5])
grid2 = grid_create_from_coordinates(xcoords2, ycoords2)

field1 = esmpy.Field(grid1, staggerloc=esmpy.StaggerLoc.CENTER)
field1.data[:] = 1.0

field2 = esmpy.Field(grid2, staggerloc=esmpy.StaggerLoc.CENTER)

regrid1 = esmpy.Regrid(field1, field2,
                       regrid_method=esmpy.RegridMethod.BILINEAR, unmapped_action=esmpy.UnmappedAction.IGNORE,
                       filename=weight_file)

# This should work:
regrid_again = esmpy.RegridFromFile(field1, field2, weight_file)

xcoords3 = np.array([1,2,3])
ycoords3 = np.array([1,2,3])
grid3 = grid_create_from_coordinates(xcoords3, ycoords3)

xcoords4 = np.array([1.5, 2.5])
ycoords4 = np.array([1.5, 2.5])
grid4 = grid_create_from_coordinates(xcoords4, ycoords4)

field3 = esmpy.Field(grid3, staggerloc=esmpy.StaggerLoc.CENTER)
field3.data[:] = 1.0

field4 = esmpy.Field(grid4, staggerloc=esmpy.StaggerLoc.CENTER)

# This shouldn't work:
try:
    regrid_should_fail = esmpy.RegridFromFile(field3, field4, weight_file)
except ValueError:
    print("Do something here to handle the error")
