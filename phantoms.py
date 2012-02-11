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
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as np


# Data taken from "Principles of Computerized Tomographic Imaging" by Kak and
# Slaney (1999). It is in the same format as table 3.1: x1, y1, A, B, rotation
# and refractive index
shepp_logan = [
        (0.0, 0.0, 0.92, 0.69, 90.0, 2.0),
        (0.0, -0.0184, 0.874, 0.6624, 90.0, -0.98),
        (0.22, 0.0, 0.31, 0.11, 72.0, -0.02),
        (-0.22, 0.0, 0.41, 0.16, 108.0, -0.02),
        (0.0, 0.35, 0.25, 0.21, 90.0, 0.01),
        (0.0, 0.1, 0.046, 0.046, 0.0, 0.01),
        (0.0, -0.1, 0.046, 0.046, 0.0, 0.01),
        (-0.08, -0.605, 0.046, 0.023, 0.0, 0.01),
        (0.0, -0.605, 0.023, 0.023, 0.0, 0.01),
        (0.06, -0.605, 0.046, 0.023, 90.0, 0.01),
        ]

predefined = { 'shepp-logan' : shepp_logan }

def _rotation_matrix(A, B, angle, deg=True):
    if deg:
        angle = angle * np.pi / 180.0

    Q = np.mat([[np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]])
    X = np.mat([[1. / A**2, 0.0], [0.0, 1. / B**2]])
    return Q.T * X * Q


def slice(width, height, ellipses, deg=True):
    phantom = np.zeros((height, width), dtype=np.float)
    v = np.mat([[0.0], [0.0]])

    for x1, y1, A, B, angle, phi in ellipses:
        R = _rotation_matrix(A, B, angle, deg)

        for x in xrange(width):
            for y in xrange(height):
                rx = (x - width / 2) / float(width / 2)
                ry = (y - height / 2) / float(height / 2)

                v[0, 0] = rx - x1
                v[1, 0] = ry - y1
                thr = v.T * R * v
                if thr <= 1.0:
                    phantom[y, x] += phi

    return phantom


def sinogram(width, angle_step, num_angles, ellipses, deg=True):
    sinogram = np.zeros((num_angles, width))
    thetas = np.linspace(0.0, num_angles * angle_step, num_angles)
    T = np.linspace(-1.0, 1.0, width)

    for x1, y1, A, B, angle, phi in ellipses:
        s = np.sqrt(x1**2 + y1**2)
        gamma = np.pi / 2. if x1 == 0 else np.arctan(y1 / x1)
        if y1 < 0.0:
            gamma = -gamma
        d = 1.0 if x1 <= 0 else - 1.0
        alpha = angle * np.pi / 180.0 if deg else angle
        refr = 2 * phi * A * B

        for x, t in enumerate(T):
            for y, theta in enumerate(thetas):
                theta_ = theta - alpha
                a2 = A**2 * np.cos(theta_)**2 + B**2 * np.sin(theta_)**2
                t_ = t + d * s * np.cos(gamma - theta)
                if np.abs(t_) <= np.sqrt(a2):
                    sinogram[y, x] += refr * np.sqrt(a2 - t_**2) / a2

    return sinogram
