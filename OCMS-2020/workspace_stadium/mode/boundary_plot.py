#!/usr/bin/env python
# -*- Coding: utf-8 -*-

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.


from logging import getLogger, basicConfig, DEBUG
logger = getLogger(__name__)

import sys
import numpy
import matplotlib.pyplot
from functools import reduce

SYM_CODE_NONE    = 0  # The cavity does not have any mirror symmetry
SYM_CODE_X       = 1  # The cavity is symmetric with respect to the x-axis
SYM_CODE_X_AND_Y = 2  # The cavity is symmetric with respect to both x- and y-axis
SYM_CODE_C4v     = 4  # The cavity is symmetric with respect to x=y and the y-axis

# Specify the ratio of the length of normal vectors per max(width,height),
# where width=(xmax-xmin), height=(ymax-ymin) of boundary line(s) on the right-top graph
NORMAL_VECTOR_LENGTH_RATIO = 0.1


# The container class for storing data of a (single) boundary
class Boundary:
    def __init__(self, beginning_point_index):
        self.beginning_point_index = beginning_point_index
        self.points = []
        self.normal_vectors = []
        self.ds = []
        self.k = []
        self.data_types = [ int, float, float, float, float, float, float ]


def load_boundary_data(boundary_file_name):

    # Contains data of each boundary (an instance of the Boundary class) at each element.
    # (boundary/ries may be composed of multiple lines, e.g. "annular" model.)
    boundaries = []

    # The first boundary
    boundary = Boundary(0)
    boundaries.append(boundary)

    symmetry_code = SYM_CODE_NONE
    next_point_index = 0

    with open(boundary_file_name, 'r') as boundary_lines:
        for line in boundary_lines:

            # header lines
            if line.strip().startswith('#'):

                # parse header
                line = line.replace(' ', '').replace('\t', '').replace('\r', '').replace('\n', '')
                header_name  = line.split('=')[0]
                header_value = line.split('=')[1]

                # read header and store value
                if header_name == "#SYM":
                    symmetry_code = int(header_value)

            # blank-lines (as separators between multiple boundaries)
            elif len(line.split()) == 0:

                # Create an new instance of Boundary
                boundary = Boundary(next_point_index)
                boundaries.append(boundary)

            # data lines
            else:
                # read data and store them in the current instance of Boundary
                data = [ boundary.data_types[i](t) for i, t in enumerate(line.split()) ]
                logger.debug(data)
                boundary.points.append(numpy.array([data[1], data[2]]))
                boundary.normal_vectors.append(numpy.array([data[3], data[4]]))
                boundary.ds.append(data[5])
                boundary.k.append(data[6])
                next_point_index += 1

    return boundaries, symmetry_code


def load_domain_data(domain_file_name):
    points = []
    data_types = [ float, float ]
    with open(domain_file_name, 'r') as boundary:
        for line in boundary:
            data = [ data_types[i](t) for i, t in enumerate(line.split()) ]
            logger.debug(data)
            points.append(numpy.array([data[0], data[1]]))
    return points


def plot_quadrant_boundary(xs, ys, axes=None, color='blue', linewidth=2, marker=None, markercolor='red'):
    if not axes:
        axes = matplotlib.pyplot.axes()

    xs = numpy.array(xs)
    xs = numpy.concatenate([ xs, -xs[::-1], -xs, +xs[::-1] ])

    ys = numpy.array(ys)
    ys = numpy.concatenate([ ys, +ys[::-1], -ys, -ys[::-1] ])

    #axes.plot(numpy.insert(xs, -1, xs[0]), numpy.insert(ys, -1, ys[0]), color=color, linewidth=linewidth, marker=marker, markeredgecolor=markercolor, markerfacecolor=markercolor)

    # add an end-point with the same coordinates as the beginning-point to close the line.
    xs = numpy.append(xs, xs[0])
    ys = numpy.append(ys, ys[0])

    # draw a line
    axes.plot(xs, ys, color=color, linewidth=linewidth, marker=marker, markeredgecolor=markercolor, markerfacecolor=markercolor)


def plot_x_axis_mirror_boundary(xs, ys, axes=None, color='blue', linewidth=2, marker=None, markercolor='red'):
    if not axes:
        axes = matplotlib.pyplot.axes()

    # supplement mirror data,
    # where numpy.concatenate merges arrays, and zs[::-1] is the reversed-order array of zs.
    xs = numpy.array(xs)
    xs = numpy.concatenate([ xs, xs[::-1] ])
    ys = numpy.array(ys)
    ys = numpy.concatenate([ ys, -ys[::-1] ])

    # add an end-point with the same coordinates as the beginning-point to close the line.
    xs = numpy.append(xs, xs[0])
    ys = numpy.append(ys, ys[0])

    # draw a line
    axes.plot(xs, ys, color=color, linewidth=linewidth, marker=marker, markeredgecolor=markercolor, markerfacecolor=markercolor)


def plot_asymmetric_boundary(xs, ys, axes=None, color='blue', linewidth=2, marker=None, markercolor='red'):
    if not axes:
        axes = matplotlib.pyplot.axes()

    xs = numpy.array(xs)
    ys = numpy.array(ys)

    # add an end-point with the same coordinates as the beginning-point to close the line.
    xs = numpy.append(xs, xs[0])
    ys = numpy.append(ys, ys[0])

    # draw a line
    axes.plot(xs, ys, color=color, linewidth=linewidth, marker=marker, markeredgecolor=markercolor, markerfacecolor=markercolor)

def plot_C4v_boundary(xs, ys, axes=None, color='blue', linewidth=2, marker=None, markercolor='red'):
    if not axes:
        axes = matplotlib.pyplot.axes()

    xs = numpy.array(xs)
    ys = numpy.array(ys)

    xs_ext = numpy.concatenate([ xs, ys[::-1], -ys, -xs[::-1], -xs, -ys[::-1], ys, xs[::-1] ])
    ys_ext = numpy.concatenate([ ys, xs[::-1], xs, ys[::-1], -ys, -xs[::-1], -xs, -ys[::-1] ])

    # add an end-point with the same coordinates as the beginning-point to close the line.
    xs = numpy.append(xs_ext, xs[0])
    ys = numpy.append(ys_ext, ys[0])

    # draw a line
    axes.plot(xs, ys, color=color, linewidth=linewidth, marker=marker, markeredgecolor=markercolor, markerfacecolor=markercolor)

def plot_boundary(points, axes=None, color='blue', linewidth=2, marker=None, markercolor='red', symmetry_code=SYM_CODE_NONE):
    # split points list to two lists of x and y coordinates
    xs, ys = zip(*points)

    if not axes:
        axes = matplotlib.pyplot.axes()

    if symmetry_code == SYM_CODE_C4v :
        logger.debug('C4v')
        plot_C4v_boundary(xs, ys, axes=axes, color=color, linewidth=linewidth, marker=marker, markercolor=markercolor)
    elif symmetry_code == SYM_CODE_X_AND_Y :
        logger.debug('quadrant')
        plot_quadrant_boundary(xs, ys, axes=axes, color=color, linewidth=linewidth, marker=marker, markercolor=markercolor)
    elif symmetry_code == SYM_CODE_X :
        logger.debug('x_axis_mirror')
        plot_x_axis_mirror_boundary(xs, ys, axes=axes, color=color, linewidth=linewidth, marker=marker, markercolor=markercolor)
    else:
        plot_asymmetric_boundary(xs, ys, axes=axes, color=color, linewidth=linewidth, marker=marker, markercolor=markercolor)

    if axes:
        axes.set_xlabel('x')
        axes.set_ylabel('y')
    else:
        matplotlib.pyplot.xlabel('x')
        matplotlib.pyplot.ylabel('y')


def plot_domain(points, axes=None, color='yellow', marker='.', size=1):
    # split points list to two lists of x and y coordinates
    xs, ys = zip(*points)

    if not axes:
        axes = matplotlib.pyplot.axes()

    axes.scatter(xs, ys, c=color, marker=marker, s=size)


def plot_normal_vectors(points, normal_vectors, length, axes=None):
    # split points list to two lists of x and y coordinates
    xs, ys = zip(*points)
    logger.debug('xs = {xs}'.format(xs=xs))
    logger.debug('ys = {ys}'.format(ys=ys))

    a = length
    vectors = sum([ (lambda p1, p2: [ [ p1[0], p2[0] ], [ p1[1], p2[1] ] ])(p, p + 0.9 * a * n) for p, n in zip(points, normal_vectors) ], [])
    logger.debug('vectors = {vectors}'.format(vectors=vectors))

    if not axes:
        axes = matplotlib.pyplot.axes()

    axes.plot(xs, ys, marker='x', markeredgecolor='red')
    axes.plot(*vectors, color='blue', linewidth=0)

    for p, n in zip(points, normal_vectors):
        an = a * n
        axes.arrow(p[0], p[1], an[0], an[1], color='blue', linewidth=2)

    if axes:
        axes.set_xlabel('x')
        axes.set_ylabel('y')
    else:
        matplotlib.pyplot.xlabel('x')
        matplotlib.pyplot.ylabel('y')


def plot_ds(ds, beginning_point_index, axes=None):
    if not axes:
        axes = matplotlib.pyplot.axes()

    end_point_index = beginning_point_index + len(ds) - 1
    indices = numpy.arange(beginning_point_index, end_point_index + 1)

    axes.plot(indices, ds, '.', markeredgecolor='red', markerfacecolor='red')

    if axes:
        axes.set_xlabel('n')
        axes.set_ylabel('ds')
    else:
        matplotlib.pyplot.xlabel('n')
        matplotlib.pyplot.ylabel('ds')


def plot_curvatures(k, beginning_point_index=None, axes=None):
    if not axes:
        axes = matplotlib.pyplot.axes()

    end_point_index = beginning_point_index + len(k) - 1
    indices = numpy.arange(beginning_point_index, end_point_index + 1)

    axes.plot(indices, k, '.', markeredgecolor='red', markerfacecolor='red')

    if axes:
        axes.set_xlabel('n')
        axes.set_ylabel('curvature')
    else:
        matplotlib.pyplot.xlabel('n')
        matplotlib.pyplot.ylabel('curvature')


# get xmax,xmin,ymax,ymin of all points of boundary lines 
# (does not contain mirror projected points)
def get_ranges_of_all_points_of_boundaries(boundaries):
    xmax = -sys.float_info.max
    xmin = sys.float_info.max
    ymax = -sys.float_info.max
    ymin = sys.float_info.max
    xmaxSwapped = False
    xminSwapped = False
    ymaxSwapped = False
    yminSwapped = False

    if len(boundaries) == 0 : 
        raise ValueError("Threre is no boundary line")

    for boundary in boundaries:
        xs, ys = zip(*boundary.points)
        xmax_part = numpy.amax(xs)
        xmin_part = numpy.amin(xs)
        ymax_part = numpy.amax(ys)
        ymin_part = numpy.amin(ys)
        if xmax_part > xmax :
            xmax = xmax_part
            xmaxSwapped = True
        if xmin_part < xmin :
            xmin = xmin_part
            xminSwapped = True
        if ymax_part > ymax :
            ymax = ymax_part
            ymaxSwapped = True
        if ymin_part < ymin :
            ymin = ymin_part
            yminSwapped = True

    if xmax<xmin or ymax<ymin or not (xmaxSwapped and xminSwapped and ymaxSwapped and yminSwapped):
        raise ValueError("Failed to find max/min coordinate values. Data might contain invalid values.")

    return xmax, xmin, ymax, ymin


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Tool for visualizing the information of a cavity shape and its boundary.")
    # parser.add_argument("boundary", help="load cavity boundary data")
    parser.add_argument("--boundary", default=None, action="store", required=True, help="load cavity boundary data")
    parser.add_argument("--domain", default=None, action="store", help="load cavity domain data")
    parser.add_argument('--debug', default=False, action='store_true', help='output debug log')
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(DEBUG)

    # load boundary data
    #points, normal_vectors, ds, k = load_data(args.boundary)
    #boundary_points, normal_vectors, ds, k = load_boundary_data(args.boundary)
    boundaries, symmetry_code = load_boundary_data(args.boundary)

    # get max/min values of x and y axes of the boundary line on the right-top graph
    xmax, xmin, ymax, ymin = get_ranges_of_all_points_of_boundaries(boundaries)

    # auto adjustment of the length of normal vectors
    width = xmax - xmin
    height = ymax - ymin
    normal_vector_length = max(width, height) * NORMAL_VECTOR_LENGTH_RATIO

    matplotlib.pyplot.figure(figsize=(6.5, 6))

    ax1 = matplotlib.pyplot.subplot(2, 2, 1)
    ax2 = matplotlib.pyplot.subplot(2, 2, 2)
    ax3 = matplotlib.pyplot.subplot(2, 2, 3)
    ax4 = matplotlib.pyplot.subplot(2, 2, 4)


    for boundary in boundaries:
        plot_boundary(boundary.points, axes=ax1, symmetry_code=symmetry_code)
        plot_normal_vectors(boundary.points, boundary.normal_vectors, normal_vector_length, axes=ax2)
        plot_ds(boundary.ds, boundary.beginning_point_index, axes=ax3)
        plot_curvatures(boundary.k, boundary.beginning_point_index, axes=ax4)

    ax1.set_aspect('equal', 'datalim')
    ax1.grid()
    ax2.set_aspect('equal', 'datalim')
    ax2.grid()


    # load and plot domain data
    if args.domain != None:
        domain_points = load_domain_data(args.domain)
        plot_domain(domain_points, axes=ax1)

    matplotlib.pyplot.tight_layout()

    matplotlib.pyplot.show()


try:
    if __name__ == '__main__':
        basicConfig(
            format='[%(asctime)s] %(module)s.%(funcName)s %(levelname)s -> %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        sys.exit(main(sys.argv))
except KeyboardInterrupt:
    sys.stderr.write("\nkeyboard interrupted.\n\n")
    sys.exit(1)
