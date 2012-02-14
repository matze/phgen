/* Copyright (C) 2012 Matthias Vogelgesang
 *
 * This program is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <tiffio.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

static const unsigned int num_shepp_logan = 10;

static double shepp_logan[] = {
        0.0, 0.0, 0.92, 0.69, 90.0, 2.0,
        0.0, -0.0184, 0.874, 0.6624, 90.0, -0.98,
        0.22, 0.0, 0.31, 0.11, 72.0, -0.02,
        -0.22, 0.0, 0.41, 0.16, 108.0, -0.02,
        0.0, 0.35, 0.25, 0.21, 90.0, 0.01,
        0.0, 0.1, 0.046, 0.046, 0.0, 0.01,
        0.0, -0.1, 0.046, 0.046, 0.0, 0.01,
        -0.08, -0.605, 0.046, 0.023, 0.0, 0.01,
        0.0, -0.605, 0.023, 0.023, 0.0, 0.01,
        0.06, -0.605, 0.046, 0.023, 90.0, 0.01
};
        

static float *compute_slice(size_t width, size_t height, const double *ellipses, size_t num_ellipses, bool is_deg)
{
    float *phantom = (float *) malloc(sizeof(float) * width * height);

    if (phantom == NULL)
        return NULL;

    /* This could be incorrect if double is not IEEE compliant */
    memset(phantom, 0, sizeof(float) * width * height);

    const double w2 = width / 2.0;
    const double h2 = height / 2.0;

    for (int i = 0; i < num_ellipses; i++) {
        const double x1 = *(ellipses++);
        const double y1 = *(ellipses++);
        const double A2 = *(ellipses) * (*ellipses);
        ellipses++;
        const double B2 = *(ellipses) * (*ellipses);
        ellipses++;
        const double alpha = is_deg ? M_PI * (*ellipses) / 180.0 : (*ellipses);
        ellipses++;
        const double phi = *(ellipses++);

        const double ca = cos(alpha);
        const double sa = sin(alpha);

        const double R[2][2] = {
            { ca*ca / A2 + sa*sa / B2, -ca*sa / A2 + ca*sa / B2 },
            { -sa*ca / A2 + ca*sa / B2, sa*sa / A2 + ca*ca / B2 }
        };

#pragma omp parallel for
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                double vx = (x - w2) / w2 - x1;
                double vy = (y - h2) / h2 - y1;
                double thr = vx * (vx * R[0][0] + vy * R[1][0]) + vy * (vx * R[0][1] + vy * R[1][1]);

                if (thr <= 1.0)
                    phantom[y * width + x] += phi;
            } 
        }
    }

    return phantom;
}

static float *compute_sinogram(size_t width, double angle_step, size_t num_angles,
        const double *ellipses, size_t num_ellipses, bool is_deg)
{
    float *sinogram = (float *) malloc(sizeof(float) * width * num_angles);

    if (sinogram == NULL)
        return NULL;
    
    /* This could be incorrect if double is not IEEE compliant */
    memset(sinogram, 0, sizeof(float) * width * num_angles);

    const double t_step = 2.0 / (width - 1); /* Including end-point +1.0 */

    for (int i = 0; i < num_ellipses; i++) {
        const double x1 = *(ellipses++);
        const double y1 = *(ellipses++);
        const double A = *(ellipses++);
        const double B = *(ellipses++);
        const double alpha = is_deg ? M_PI * (*ellipses) / 180.0 : (*ellipses);
        ellipses++;
        const double phi = *(ellipses++);

        double gamma = x1 == 0.0 ? M_PI / 2.0 : atan(y1 / x1);
        gamma = y1 < 0.0 ? -gamma : gamma;

        const double s = sqrt(x1*x1 + y1*y1);
        const double r = 2.0 * phi * A * B;
        const double d = x1 <= 0.0 ? 1.0 : -1.0;

#pragma omp parallel for
        for (int y = 0; y < num_angles; y++) {
            const double theta = y * angle_step;
            const double theta_ = theta - alpha;
            const double c_theta = cos(theta_);
            const double s_theta = sin(theta_);
            const double a2 = A * A * c_theta * c_theta + B * B * s_theta * s_theta;
            const double r_a2 = r / a2;

            double t = -1.0 + d * s * cos(gamma - theta);

            for (int x = 0;  x < width; x++) {
                if (fabs(t) <= sqrt(a2))
                    sinogram[y * width + x] += r_a2 * sqrt(a2 - t*t);
                
                t += t_step;
            }
        }
    }
    
    return sinogram;
}

static void write_tiff(float *image, const char *name, size_t width, size_t height)
{
    TIFF *tif = TIFFOpen(name, "w");

    if (tif == NULL)
        return;

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    const size_t rows_per_strip = TIFFDefaultStripSize(tif, 0);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rows_per_strip);

    for (int y = 0; y < height; y++, image += width)
        TIFFWriteScanline(tif, image , y, 0);

    TIFFClose(tif);
}

static void usage(void)
{
    fprintf(stderr, 
            "usage: generate --width WIDTH --height HEIGHT [Options] [-]\n" \
            "Options:\n" \
            "  -s, --sinogram\tGenerate sinogram\n" \
            "  -p, --phantom\t\tGenerate slice phantom\n" \
            "  -a, --angle-step S\tAngle between two projections\n" \
            "  -w, --width WIDTH\tWidth of phantom or projection\n" \
            "  -h, --height HEIGHT\tHeight of phantom or number of projections\n" \
            "  -\t\t\tRead ellipses data from stdin\n");
    exit(EXIT_SUCCESS);
}

double *read_ellipses(size_t *num_ellipses)
{
    const char *delimiter = " ,;\t\n";
    const size_t num_bytes = 128;
    char buffer[num_bytes];

    double *ellipses = NULL;
    const int num_params = 6;
    size_t n = 0;
    int current_param = 0;
    double params[num_params];

    while (fgets(buffer, num_bytes, stdin) != NULL) {
        char *pch = strtok(buffer, delimiter);

        while (pch != NULL) {
            params[current_param++] = atof(pch);

            if (current_param == num_params) {
                current_param = 0;
                n++;
                ellipses = (double *) realloc(ellipses, sizeof(double) * num_params * n);

                for (int i = 0; i < num_params; i++)
                    ellipses[(n-1)*num_params + i] = params[i]; 
            }

            pch = strtok(NULL, delimiter);
        }
    }

    *num_ellipses = n;
    return ellipses;
}

int main(int argc, const char **argv)
{
    int getopt_ret, option_index;

    static struct option long_options[] = {
        { "sinogram", no_argument, 0, 's' },
        { "phantom",  no_argument, 0, 'p' },
        { "width", required_argument, 0, 'w' },
        { "height", required_argument, 0, 'h' },
        { "help", no_argument, 0, '?' },
        { "angle-step", required_argument, 0, 'a' },
        {0, 0, 0, 0}       
    };

    bool read_from_stdin = strcmp(argv[argc-1], "-") == 0;
    bool generate_phantom = false, generate_sinogram = false;
    double angle_step = -1.0;
    int width = -1, height = -1;

    while (1) {
        getopt_ret = getopt_long(argc, (char *const *) argv, "spw:h:?:a:n:", long_options,
                &option_index);

        if (getopt_ret == -1) 
            break;

        switch (getopt_ret) {
            case 's':
                    generate_sinogram = true;
                    break;
            case 'p':
                    generate_phantom = true;
                    break;
            case 'w':
                    width = atoi(optarg);
                    break;
            case 'h':
                    height = atoi(optarg);
                    break;
            case '?':
                    usage();
            default:
                    break;
        }
    }

    size_t num_ellipses = num_shepp_logan;
    double *ellipses = shepp_logan;

    if (read_from_stdin)
        ellipses = read_ellipses(&num_ellipses);

    if ((width < 0) || (height < 0) || (!generate_phantom && !generate_sinogram))
        usage();

    if (generate_phantom) {
        float *phantom = compute_slice(width, height, ellipses, num_ellipses, true);
        write_tiff(phantom, "phantom.tif", width, height);
        free(phantom);
    }

    if (generate_sinogram) {
        angle_step = angle_step < 0.0 ? M_PI / height : angle_step;
        float *sinogram = compute_sinogram(width, angle_step, height, ellipses, num_ellipses, true);
        write_tiff(sinogram, "sinogram.tif", width, height);
        free(sinogram);
    }
    return 0;
}
