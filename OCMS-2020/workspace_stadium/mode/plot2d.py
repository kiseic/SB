#!/usr/bin/env python
# -*- Coding: utf-8 -*-

# Copyright(C) 2019,2020 Telecognix Corporation. All rights reserved.
# See OCMS-License-OCMS-2020-Basic.txt/OCMS-License-OCMS-2020-Extension-Tools.txt in OCMS top directory.
# Contact ocms@telecognix.com for further information.


# Tomohiro Miyasaka (miyasaka@telecognix.com)

from logging import getLogger, basicConfig, DEBUG
logger = getLogger(__name__)

import sys
import re
import bz2

import numpy
import matplotlib.pyplot
import mpl_toolkits.axes_grid1

import boundary_plot


CONTENT_HUSIMI = 'husimi'
CONTENT_WAVEFUNCTION = 'wavefunction'


symbols = { '$', '%', '#' }

# header settings for "mtv" format files
pattern_mtv_header = r'\s*([a-zA-Z]+)\s*=\s*("[^"]*"|[+-.0-9eEdD]+|[tT][rR][uU][eE]|[fF][aA][lL][sS][eE])'
mtv_params = {
    'nx' : int, 'ny' : int, 'nsteps' : int,
    'xmin' : float , 'xmax' : float, 'ymin' : float , 'ymax' : float, 
    'cmin' : float, 'cmax' : float, 'xyratio' : float,
    'toplabel' : str, 'subtitle' : str,
    'contfill' : bool,
    'content' : str
}

# header settings for "read" and "complex" format files
pattern_gen2_format_header = r'\s*([a-zA-Z]+)\s*=\s*("[^"]*"|[+-.0-9eEdD]+|[tT][rR][uU][eE]|[fF][aA][lL][sS][eE])'
gen2_format_params = {
    'nx' : int, 'ny' : int, 'nsteps' : int,
    'xmin' : float , 'xmax' : float, 'ymin' : float , 'ymax' : float, 
    'vmin' : float, 'vmax' : float, 'xyratio' : float,
    'toplabel' : str, 'subtitle' : str,
    'contfill' : bool,
    'content' : str
}

content_dataformat_dict = { 'wavefunction':'complex', 'husimi':'real' }


def parse_mtv_header(line):
    res = re.findall(pattern_mtv_header, line)
    ret = []
    for var, val in res:
        var = var.lower()
        if mtv_params[var] is str:
            val = val[1:-1]
        else:
            val = mtv_params[var](val)
        ret.append((var, val))
    return ret


#unused (> read_data)
def read_mtv(file, params, data):
    for line in file:
        line = line.strip()
        if line[0] == '#':
            pass
        elif line[0] == '$':
            pass
        elif line[0] == '%':
            results = parse_mtv_header(line)
            for var, val in results:
                params[var] = val
        else:
            data.append(float(line))
    return params, data


#unused (> load_data)
def load_mtv(filename):
    params = {}
    data = []
    if re.match('^.*\.bz2', filename): # bz2ed text
        with bz2.open(filename, 'rt') as file:
            params, data = read_mtv(file, params, data)
    else: # raw text
        with open(filename, 'r') as file:
            params, data = read_mtv(file, params, data)
    return numpy.array(data).reshape((params['ny'], params['nx'])), params


def parse_gen2_format_header(line):
    res = re.findall(pattern_gen2_format_header, line)
    ret = []
    for var, val in res:
        var = var.lower()

        if not var in gen2_format_params:
            continue

        if gen2_format_params[var] is str:
            val = val[1:-1]
        else:
            val = gen2_format_params[var](val)

        ret.append((var, val))

    return ret


def read_data(file, params, data, format_name):

    # Flags to reduce comparison costs of branchings in the loop
    is_auto_format    = (format_name == "auto")
    is_real_format    = (format_name == "real")
    is_complex_format = (format_name == "complex")
    is_mtv_format     = (format_name == "mtv")
    is_gen2_formats = is_auto_format or is_complex_format or is_real_format

    for line in file:
        line = line.strip()

        # Read header lines
        if line[0] == '$':
            pass
        elif line[0] == '#':
            if is_gen2_formats:
                results = parse_gen2_format_header(line)
                for var, val in results:
                    params[var] = val
        elif line[0] == '%':
            if is_mtv_format:
                results = parse_mtv_header(line)
                for var, val in results:
                    params[var] = val

        # Read data lines
        else:
            # "auto" format: Determines the actual data format at the first data line (after parsing header lines).
            if is_auto_format:

                # The default data format is 'real' when 'auto' is specified.
                is_auto_format    = False
                is_complex_format = False
                is_real_format    = True

                # If 'content' is specified in headers, set suitable format for the content automatically.
                if 'content' in params:
                    if params['content'] in content_dataformat_dict:
                        suitable_format_for_content = content_dataformat_dict[ params['content'] ]
                        is_real_format    = (suitable_format_for_content == "real")
                        is_complex_format = (suitable_format_for_content == "complex")

            # "real" or "mtv" formats: Reads single-column data.
            if is_real_format or is_mtv_format:
                intensity = float(line)
                data.append(intensity)

            # "complex" format: Reads double-column data.
            elif is_complex_format:
                columns = line.split()
                re = float(columns[0])
                im = float(columns[1])
                intensity = re*re + im*im
                data.append(intensity)
            else:
                sys.stderr.write("Unsupported format: " + format_name)
                return None

    return params, data


def load_data(filename, format_name):
    params = {}
    data = []
    if re.match('^.*\.bz2', filename): # bz2ed text
        with bz2.open(filename, 'rt') as file:
            params, data = read_data(file, params, data, format_name)
    else: # raw text
        with open(filename, 'r') as file:
            params, data = read_data(file, params, data, format_name)

    return numpy.array(data).reshape((params['ny'], params['nx'])), params


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Tool for plotting a wave function and a Husimi distribution.")
    parser.add_argument('filename', help='wave function or Husimi distribution data file')
    parser.add_argument('--savefig', default=False, action='store_true', help='save figure image as png')
    parser.add_argument('--colormap', action='store', type=str, default='jet', help='set colormap(default is \'jet\')')
    parser.add_argument('--colorbarformat', default="{:.1E}", type=str, action='store', help='set format of colorbar ticks')
    parser.add_argument('--nocolorbar', default=False, action='store_true', help='don\'t display colorbar')
    parser.add_argument('--nocolorbarticks', default=False, action='store_true', help='don\'t display ticks of colorbar')
    parser.add_argument('--rightmargin', action='store', type=float, metavar="[0.0 to 1.0]", help='set ratio of right-margin')
    parser.add_argument('--vscale', action='store', type=str, default='linear', metavar="[linear/log]", help='set the v-axis scale (default is \'linear\')')
    parser.add_argument('--vmin', action='store', type=float, default=None)
    parser.add_argument('--vmax', action='store', type=float, default=None)
    parser.add_argument('--boundary', action='store', type=str, default=None, help='draw cavity boundary')
    parser.add_argument('--color', action='store', type=str, default='white', help='set color of boundary line (default is \'white\')')
    parser.add_argument('--linewidth', action='store', type=float, default=2.0, help='set width of boundary line (default is 2)')
    parser.add_argument('--labelsize', action='store', type=float, default=18, help='set font size of axes label (default is 18)')
    parser.add_argument('--ticksize', action='store', type=float, default=14, help='set font size of tick labels (default is 14)')
    parser.add_argument('--textsize', action='store', type=float, default=11, help='set font size of text (default is 11)')
    parser.add_argument('--dataformat', action='store', type=str, default='auto', metavar="[auto/real/complex/mtv]", help='set format of data file (default is \'auto\')')
    parser.add_argument('--debug', default=False, action='store_true', help='output debug log')
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(DEBUG)

    wfunc, params = load_data(args.filename, args.dataformat)

    extent = [ params['xmin'], params['xmax'], params['ymin'], params['ymax'] ]
    fig = matplotlib.pyplot.figure(figsize=(8,6))
    axes = fig.subplots(1, 1)

    vmin = None ; vmax = None

    # get the value of "vmax" from the header if it is defined.
    if 'vmax' in params.keys():    # for "real" and "complex" formats
        vmax = params['vmax']
    elif 'cmax' in params.keys():  # for "mtv" format
        vmax = params['cmax']

    # get the value of "vmin" from the header if it is defined.
    if 'vmin' in params.keys():    # for "real" and "complex" formats
        vmin = params['vmin']
    elif 'cmin' in params.keys():  # for "mtv" format
        vmin = params['cmin']

    # if "vmin" is not defined in the header, determin the value of it by scanning data
    if not vmin:
        vmin = sys.float_info.max
        for ix in range(params['nx']):
            for iy in range(params['ny']):
                if wfunc[iy,ix] < vmin:
                    vmin = wfunc[iy,ix]


    if 'xyratio' in params.keys():
        aspect = params['xyratio'] * (params['xmax'] - params['xmin']) / (params['ymax'] - params['ymin'])
    else:
        aspect = 'equal'

    if args.vmin:
        vmin = args.vmin
    if args.vmax:
        vmax = args.vmax

    if args.vscale.lower() == 'log':
        vscale = matplotlib.colors.LogNorm(1e-5, vmax)
    else:
        vscale = matplotlib.colors.Normalize(1e-5, vmax)


    # plot data
    image = matplotlib.pyplot.imshow(wfunc, vmin=vmin, vmax=vmax, norm=vscale, extent=extent, origin='lower', aspect=aspect, cmap=args.colormap)


    matplotlib.pyplot.xlabel('x', fontsize=args.labelsize)
    matplotlib.pyplot.ylabel('y', fontsize=args.labelsize)
    matplotlib.pyplot.title(params['toplabel'], pad=16, fontsize=args.textsize)
    matplotlib.pyplot.text(0.5, 1.01, params['subtitle'], horizontalalignment='center', verticalalignment='bottom', family='monospace', transform=axes.transAxes, fontsize=args.textsize)

    if args.boundary:
        boundaries, symmetry_code = boundary_plot.load_boundary_data(args.boundary)
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        logger.debug('xlim = {xlim}, ylim = {ylim}'.format(xlim=xlim, ylim=ylim))
        for boundary in boundaries:
            logger.debug('points = {}'.format(boundary.points))
            boundary_plot.plot_boundary(boundary.points, axes=axes, color='white', linewidth=args.linewidth, symmetry_code=symmetry_code)
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)
        axes.tick_params(labelsize=args.ticksize)

    with_colorbar = not args.nocolorbar # temporary variable for readability
    if with_colorbar:
        # For husimi-plot: simplify the configuration of the colorbar to avoid "overlength" behaviour
        if 'content' in params and params['content'].lower() == CONTENT_HUSIMI:
            colorbar = matplotlib.pyplot.colorbar()
        else:
            divider = mpl_toolkits.axes_grid1.make_axes_locatable(axes)
            cax = divider.append_axes("right", size="3%", pad=0.05)
            colorbar = matplotlib.pyplot.colorbar(image, cax=cax)
        colorbar.ax.set_xlabel("")
        colorbar.ax.set_ylabel("")

    matplotlib.pyplot.tight_layout()

    # Modification of ticklabels of the colorbar
    # (should be executed after when contents of tickslabels of the colorber are determined)
    if with_colorbar and args.colorbarformat:
        colorbar_ticks_list = colorbar.ax.get_yticklabels()
        colorbar_formatted_ticks_list = [ args.colorbarformat.format(tick._y) for tick in colorbar_ticks_list]
        colorbar.ax.set_yticklabels(colorbar_formatted_ticks_list)

    if args.rightmargin:
        matplotlib.pyplot.subplots_adjust(right = 1.0 - args.rightmargin)

    if args.nocolorbarticks:
        colorbar.ax.tick_params(right=False, labelright=False)

    if args.savefig:
        print('save image file: {filename} ....'.format(filename=args.filename+'.png'))
        matplotlib.pyplot.savefig(args.filename+'.png', transparent=True)
    else:
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
