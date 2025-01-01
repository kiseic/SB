#!/usr/bin/env python
# -*- Coding: utf-8 -*-

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.


import sys
import numpy
import matplotlib.pyplot


det_params = {
    'ixmax' : int,
    'iymax' : int
    }


def search_local_minimals(z):
    imax, jmax = z.shape
    res = []
    for j in range(jmax):
        for i in range(imax):
            zp = _i, _j = search_local_minimal(z, i, j)
            # ignore local minimals on boundary of search region of wavenumber space
            if _i > 0 and _i < imax-1 and _j > 0 and _j < jmax-1:
            # if True:
                res.append(zp)
    return list(set(res))


def search_local_minimal(z, i, j):
    imax, jmax = z.shape
    cur = nxt = (i, j)
    rels = numpy.array([ (+1, 0), (0, +1), (-1, 0), (0, -1) ])
    curz = nxtz = z[i, j]
    while(True):
        for i, j in rels + cur:
            if i >= 0 and i < imax and j >= 0 and j < jmax:
                if z[i, j] < nxtz:
                    nxtz = z[i, j]
                    nxt = (i, j)
        if curz > nxtz:
            curz = nxtz
            cur = nxt
        else:
            return nxt


def get_det_params(params, expr):
    expr = expr[1:].strip().split()
    for i in range(len(expr))[::3]:
        if len(expr) >= i+3 and expr[i] in det_params.keys() and expr[i+1] == '=':
            params[expr[i]] = det_params[expr[i]](expr[i+2])


def load_det_data(filename):
    params = {}
    x = [] ; y = [] ; z = []
    with open(filename, 'r') as mtv:
        num_lines = 0
        for line in mtv:
            line = line.strip()
            if line[0] == '#':
                get_det_params(params, line)
            else:
                _x, _y, _z = line.split()
                x.append(float(_x))
                y.append(float(_y))
                z.append(float(_z))
                num_lines += 1
    imax = params['ixmax']
    jmax = params['iymax']
    x = numpy.array(x).reshape(jmax, imax)
    y = numpy.array(y).reshape(jmax, imax)
    z = numpy.array(z).reshape(jmax, imax)
    return x, y, z


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Tool for finding local minimums of a determinant value distribution in the complex wave number space.")
    parser.add_argument('filename', help='determinant value data file')
    #parser.add_argument('--show', default=False, action='store_true', help='show plot figure')
    parser.add_argument('--nodisplay', default=False, action='store_true', help='don\'t display plot figure')
    parser.add_argument('--savefig', default=False, action='store_true', help='save figure image as png')
    #parser.add_argument('--zeropoints', default=False, action='store_true', help='output local minimums found')
    parser.add_argument('--nozeropoints', default=False, action='store_true', help='don\'t output local minimums found')
    parser.add_argument('--colormap', action='store', type=str, default='viridis_r', help='set colormap(default is \'viridis_r\')')
    parser.add_argument('--annotate', default=False, action='store_true', help='show annotations of local minimums found')
    parser.add_argument('--labelsize', action='store', type=float, default=18, help='set font size of axes label (default is 18)')
    parser.add_argument('--ticksize', action='store', type=float, default=14, help='set font size of tick labels (default is 14)')
    parser.add_argument('--textsize', action='store', type=float, default=9, help='set font size of text (default is 9)')
    parser.add_argument('--textcolor', action='store', type=str, default='black', help='set text color(default is \'black\')')
    args = parser.parse_args()

    x, y, z = load_det_data(args.filename)

    ilocmins = search_local_minimals(z)
    locmins = []
    for i, j in ilocmins:
        locmins.append((x[i, j], y[i, j]))

    #if args.zeropoints:
    if not args.nozeropoints:
        for i, j in ilocmins:
            print(x[i, j], y[i, j], z[i, j])

    if (not args.nodisplay) or args.savefig:
        xmin = x[0, 0]
        xmax = x[-1, -1]
        ymin = y[0, 0]
        ymax = y[-1, -1]
        extent=[xmin, xmax, ymin, ymax]

        fig, axes = matplotlib.pyplot.subplots()

        # interpolation = 'nearest'
        # interpolation = 'bilinear'
        interpolation = 'bicubic'
        matplotlib.pyplot.imshow(z, extent=extent, origin='lower', interpolation=interpolation, zorder=1, cmap=args.colormap)

        if len(locmins) == 0:
            print('! caution: no local minimums are found in this area.')
        else:
            matplotlib.pyplot.scatter(*zip(*locmins), marker='+', color='black', zorder=2)

        if args.annotate:
            for i, j in ilocmins:
                _x, _y, _z = x[i, j], y[i, j], z[i, j]
                matplotlib.pyplot.annotate(str('({x}, {y})'.format(x=_x, y=_y)), xy=(_x, _y), xytext=(5,5), textcoords='offset pixels', color=args.textcolor, fontsize=args.textsize)

        jmax, imax = z.shape
        dmx = (xmax - xmin) / (imax - 1)
        dmy = (ymax - ymin) / (jmax - 1)
        axes.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=dmx))
        axes.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=dmy))
        matplotlib.pyplot.grid(b=True, which='both', alpha=0.3)
        matplotlib.pyplot.xlabel('Re k', fontsize=args.labelsize)
        matplotlib.pyplot.ylabel('Im k', fontsize=args.labelsize)
        matplotlib.pyplot.tick_params(labelsize=args.ticksize)

        matplotlib.pyplot.tight_layout()

        #if args.show:
        if not args.nodisplay:
            matplotlib.pyplot.show()
        if args.savefig:
            print('save image file: {filename} ....'.format(filename=args.filename+'.png'))
            matplotlib.pyplot.savefig(args.filename+'.png', transparent=True)


try:
    if __name__ == '__main__':
        sys.exit(main(sys.argv))
except KeyboardInterrupt:
    sys.stderr.write("\nkeyboard interrupted.\n\n")
    sys.exit(1)
