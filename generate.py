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

import argparse
import phantoms
import numpy as np
from libtiff import TIFF


if __name__ == '__main__':
    width, height = (256, 256)
    phantom = phantoms.slice(width, height, phantoms.shepp_logan)
    # sinogram = phantom_sinogram(256, np.pi / 256, 256,
    #         [(-0.3, -0.4, 0.5, 0.3, 45.0, 1.0)])
    #sinogram = phantom_sinogram(256, np.pi / 256, 256, sh_data)

    # tif = TIFF.open('sinogram-direct.tif', mode='w')
    # tif.write_image(sinogram.astype(np.float32))

    tif = TIFF.open('phantom.tif', mode='w')
    tif.write_image(phantom.astype(np.float32))
