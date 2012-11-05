#!/usr/bin/env python
#
# Copyright (C) 2012 Matthias Vogelgesang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import argparse
import phantoms
import numpy as np
from libtiff import TIFF


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate phantoms and \
            sinograms from ellipses model.')

    parser.add_argument('-d', '--predefined', type=str, metavar='DATASET', required=True)
    parser.add_argument('-i', '--input-file', type=str, metavar='FILE')
    parser.add_argument('--width', type=int, required=True)
    parser.add_argument('--height', type=int, required=True)
    parser.add_argument('-a', '--angle-step', type=float)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--sinogram', action='store_true')
    group.add_argument('-p', '--phantom', action='store_true')

    args = parser.parse_args()

    try:
        ds = phantoms.predefined[args.predefined]
    except:
        print "Data set `%s' is not defined." % args.predefined
        sys.exit(1)

    if args.sinogram:
        sinogram = phantoms.sinogram(args.width, np.pi / args.height, args.height, ds)
        tif = TIFF.open('sinogram-%ix%i.tif' % (args.width, args.height), mode='w')
        tif.write_image(sinogram.astype(np.float32))
    if args.phantom:
        phantom = phantoms.slice(args.width, args.height, ds)
        tif = TIFF.open('phantom-%ix%i.tif' % (args.width, args.height), mode='w')
        tif.write_image(phantom.astype(np.float32))
