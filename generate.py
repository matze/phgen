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
