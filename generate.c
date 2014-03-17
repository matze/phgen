/* Copyright (C) 2012 Matthias Vogelgesang
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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
 
static const double M_PI = 3.1415926535897932384626433832795;
 
static const unsigned N_SHEPP_LOGAN_ELLIPSES = 10;
static const unsigned N_ELLIPSE_PARAMS = 6;
 
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
 
typedef struct {
    bool read_from_stdin;
    bool generate_phantom;
    bool generate_sinogram;
    bool split;
    int width;
    int height;
    double angle_step;
    int n_frames;
    int n_vectors;
    int n_ellipsesize;
    int n_alpha;
    int n_phi;
    double *vectors;
    double *ellipsesize;
    double *alpha;
    double *phi;
} Options;
 
static float *
compute_slice (size_t width,
               size_t height,
               const double *ellipses,
               size_t num_ellipses,
               bool is_deg)
{
    float *phantom = (float *) malloc (sizeof (float) * width * height);
 
    if (phantom == NULL)
        return NULL;
 
    /* This could be incorrect if double is not IEEE compliant */
    memset (phantom, 0, sizeof (float) * width * height);
 
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
 
        const double ca = cos (alpha);
        const double sa = sin (alpha);
 
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
 
static float *
compute_sinogram (size_t width,
                  double angle_step,
                  size_t num_angles,
                  const double *ellipses,
                  size_t num_ellipses,
                  bool is_deg)
{
    float *sinogram = (float *) malloc (sizeof (float) * width * num_angles);
 
    if (sinogram == NULL)
        return NULL;
 
    /* This could be incorrect if double is not IEEE compliant */
    memset (sinogram, 0, sizeof (float) * width * num_angles);
 
    const double t_step = 2.0 / (width - 1); /* Including end-point +1.0 */
 
    for (int i = 0; i < num_ellipses; i++) {
        const double x1 = *(ellipses++);
        const double y1 = *(ellipses++);
        const double A = *(ellipses++);
        const double B = *(ellipses++);
        const double alpha = is_deg ? M_PI * (*ellipses) / 180.0 : (*ellipses);
        ellipses++;
        const double phi = *(ellipses++);
        const double gamma = atan2 (y1, x1);
 
        const double s = sqrt (x1*x1 + y1*y1);
        const double r = 2.0 * phi * A * B;
 
#pragma omp parallel for
        for (int y = 0; y < num_angles; y++) {
            const double theta = y * angle_step;
            const double theta_ = theta - alpha;
            const double c_theta = cos (theta_);
            const double s_theta = sin (theta_);
            const double a2 = A * A * c_theta * c_theta + B * B * s_theta * s_theta;
            const double r_a2 = r / a2;
 
            double t = -1.0 - s * cos (gamma - theta);
 
            for (int x = 0; x < width; x++) {
                if (fabs (t) <= sqrt (a2))
                    sinogram[y * width + x] += r_a2 * sqrt (a2 - t*t);
 
                t += t_step;
            }
        }
    }
 
    return sinogram;
}
 
static void
write_tiff (float *image,
            const char *name,
            size_t width,
            size_t height,
            int y)
{
    TIFF *tif = TIFFOpen (name, "w");
 
    if (tif == NULL)
        return;
 
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
 
    const size_t rows_per_strip = TIFFDefaultStripSize (tif, 0);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, rows_per_strip);
 
    for (; y < height; y++, image += width)
        TIFFWriteScanline (tif, image , y, 0);
 
    TIFFClose (tif);
}
 
static void
usage (void)
{
    fprintf (stderr,
             "usage: generate -s/-p -w WIDTH -h HEIGHT [Options] [-]\n" \
             "Options:\n" \
             " -s, --sinogram\tGenerate sinogram\n" \
             " -p, --phantom\t\tGenerate slice phantom\n" \
             "     --split\t\tSplit sinogram into single lines\n" \
             " -a, --angle-step S\tAngle between two projections\n" \
             " -w, --width WIDTH\tWidth of phantom or projection\n" \
             " -h, --height HEIGHT\tHeight of phantom or number of projections\n" \
             " -n, --num-frames\tNumber of frames to generate\n" \
             " -t, --translate VEC\tString containing vector elements (x1 y1 [x2 y2 ...])\n" \
             " -r, --resize AB\tString containing vector elements (x1 y1 [x2 y2 ...])\n" \
             " -e, --resize alpha\tAlpha angle of phantom or projection\n" \
             " -f, --resize phi\tPhi angle of phantom or projection\n" \
             " -\t\t\tRead ellipses data from stdin\n");
    exit (EXIT_SUCCESS);
}
 
double *
read_ellipses (size_t *num_ellipses)
{
    const char *delimiter = " ,;\t\n";
    const size_t num_bytes = 128;
    char buffer[num_bytes];
 
    double *ellipses = NULL;
    size_t n = 0;
    int current_param = 0;
    double params[N_ELLIPSE_PARAMS];
 
    while (fgets (buffer, num_bytes, stdin) != NULL) {
        char *pch = strtok (buffer, delimiter);
 
        while (pch != NULL) {
            params[current_param++] = atof (pch);
 
            if (current_param == N_ELLIPSE_PARAMS) {
                current_param = 0;
                n++;
                ellipses = (double *) realloc (ellipses, sizeof (double) * N_ELLIPSE_PARAMS * n);
 
                for (int i = 0; i < N_ELLIPSE_PARAMS; i++)
                    ellipses[(n-1)*N_ELLIPSE_PARAMS + i] = params[i];
            }
 
            pch = strtok (NULL, delimiter);
        }
    }
 
    *num_ellipses = n;
    return ellipses;
}
 
static double *
parse_vectors (char *s, int *n)
{
    static const char *delimiter = " ";
    char *token;
    double *array;
 
    token = strtok (s, delimiter);
    array = NULL;
 
    while (token != NULL) {
        double x;
 
        x = atof (token);
        token = strtok (NULL, delimiter);
 
        if (token != NULL) {
            double y;
 
            *n += 1;
            y = atof (token);
 
            array = realloc (array, *n * 2 * sizeof (double));
            array[(*n - 1) * 2] = x;
            array[(*n - 1) * 2 + 1] = y;
            token = strtok (NULL, delimiter);
        }
    }
    return array;
}

static double *
parse_angle (char *s, int *n)
{
    static const char *delimiter = " ";
    char *token;
    double *array;
 
    token = strtok (s, delimiter);
    array = NULL;
 
    while (token != NULL) {
        double x;
        
        *n += 1;
        x = atof (token);

        array = realloc (array, *n * 1 * sizeof (double));
        array[*n - 1] = x;
        token = strtok (NULL, delimiter);
    }
    return array;
} 
 
static Options *
options_read (int argc,
              const char **argv
)
{
    Options *opts;
    int getopt_ret;
    int option_index;
 
    enum {
        OPT_HELP = '?',
        OPT_ANGLE_STEP = 'a',
        OPT_SINOGRAM = 's',
        OPT_HEIGHT = 'h',
        OPT_PHANTOM = 'p',
        OPT_WIDTH = 'w',
        OPT_NUM_FRAMES = 'n',
        OPT_TRANSLATE = 't',
        OPT_RESIZE = 'r',
        OPT_ALPHA = 'e',
        OPT_PHI = 'f',
        OPT_SPLIT = 1,
    };
 
    static struct option long_options[] = {
        { "sinogram", no_argument, 0, OPT_SINOGRAM },
        { "phantom", no_argument, 0, OPT_PHANTOM },
        { "split", no_argument, 0, OPT_SPLIT },
        { "width", required_argument, 0, OPT_WIDTH },
        { "height", required_argument, 0, OPT_HEIGHT },
        { "help", no_argument, 0, OPT_HELP },
        { "angle-step", required_argument, 0, OPT_ANGLE_STEP },
        { "num-frames", required_argument, 0, OPT_NUM_FRAMES },
        { "translate", required_argument, 0, OPT_TRANSLATE },
        { "resize", required_argument, 0, OPT_RESIZE },
        { "alpha", required_argument, 0, OPT_ALPHA },
        { "phi", required_argument, 0, OPT_PHI },
        {0, 0, 0, 0}
    };
 
    opts = malloc (sizeof (Options));
 
    opts->read_from_stdin = strcmp (argv[argc-1], "-") == 0;
    opts->n_frames = 1;
    opts->n_vectors = 0;
    opts->n_ellipsesize = 0;
    opts->n_alpha = 0;
    opts->n_phi = 0;
    opts->vectors = NULL;
    opts->ellipsesize = NULL;
    opts->alpha = NULL;
    opts->phi = NULL;
    opts->width = opts->height = opts->angle_step = -1;
 
    while (1) {
        getopt_ret = getopt_long (argc, (char *const *) argv, "spw:h:?:a:n:t:r:e:f:", long_options,
                &option_index);
 
        if (getopt_ret == -1)
            break;
 
        switch (getopt_ret) {
            case OPT_SINOGRAM:
                opts->generate_sinogram = true;
                break;
            case OPT_PHANTOM:
                opts->generate_phantom = true;
                break;
            case OPT_SPLIT:
                opts->split = true;
                break;
            case OPT_WIDTH:
                opts->width = atoi (optarg);
                break;
            case OPT_HEIGHT:
                opts->height = atoi (optarg);
                break;
            case OPT_NUM_FRAMES:
                opts->n_frames = atoi (optarg);
                break;
            case OPT_TRANSLATE:
                opts->vectors = parse_vectors (optarg, &opts->n_vectors);
                break;
            case OPT_RESIZE:
                opts->ellipsesize = parse_vectors (optarg, &opts->n_ellipsesize);
                break;
            case OPT_ALPHA:
                opts->alpha = parse_angle (optarg, &opts->n_alpha);
                break;
            case OPT_PHI:
                opts->phi = parse_angle (optarg, &opts->n_phi);
                break;
            case OPT_HELP:
                usage ();
            default:
                break;
        }
    }
 
    return opts;
}
 
static void
options_free (Options *opts)
{
    free (opts->vectors);
    free (opts->ellipsesize);
    free (opts->alpha);
    free (opts->phi);
    free (opts);
}
 
 
static void
animate_ellipses (double *ellipses,
                  size_t n_ellipses,
                  Options *opts)
{
    size_t max_ellipses;
    size_t num_ellipsesize;
    size_t num_alpha;
    size_t num_phi;
     
    max_ellipses = opts->n_vectors > n_ellipses ? n_ellipses : opts->n_vectors;
    for (int i = 0; i < max_ellipses; i++) {
        int index = i * N_ELLIPSE_PARAMS;
        ellipses[index] += opts->vectors[i*2];
        ellipses[index+1] += opts->vectors[i*2+1];
    }
 
    num_ellipsesize = opts->n_ellipsesize;
    for (int i = 0; i < num_ellipsesize; i++) {
        int index = i * N_ELLIPSE_PARAMS;
        ellipses[index+2] += opts->ellipsesize[i*2];
        ellipses[index+3] += opts->ellipsesize[i*2+1];
    }

    num_alpha = opts->n_alpha;
    for (int i = 0; i < num_alpha; i++) {
        int index = i * N_ELLIPSE_PARAMS;
        ellipses[index+4] += opts->alpha[i];
    }

    num_phi = opts->n_phi;
    for (int i = 0; i < num_phi; i++) {
        int index = i * N_ELLIPSE_PARAMS;
        ellipses[index+5] += opts->phi[i];
    }
}
 
int
main (int argc,
      const char **argv)
{
    Options *opts;
    size_t n_ellipses;
    double *ellipses;
 
    opts = options_read (argc, argv);
 
    if ((opts->width < 0) || (opts->height < 0) || (!opts->generate_phantom && !opts->generate_sinogram))
        usage ();
 
    if (opts->read_from_stdin) {
        ellipses = read_ellipses (&n_ellipses);
    }
    else {
        ellipses = shepp_logan;
        n_ellipses = N_SHEPP_LOGAN_ELLIPSES;
    }
 
    for (int i = 0; i < opts->n_frames; i++) {
        if (opts->generate_phantom) {
            char filename[FILENAME_MAX];
            float *phantom;
 
            phantom = compute_slice (opts->width, opts->height, ellipses, n_ellipses, true);
 
            snprintf (filename, FILENAME_MAX, "phantom-%05i.tif", i);
            write_tiff (phantom, filename, opts->width, opts->height, 0);
            free (phantom);
        }
 
        if (opts->generate_sinogram) {
            char filename[FILENAME_MAX];
            double angle_step;
            float *sinogram;
 
            angle_step = opts->angle_step < 0.0 ? M_PI / opts->height : opts->angle_step;
            sinogram = compute_sinogram (opts->width, angle_step, opts->height, ellipses, n_ellipses, true);

            if (opts->split) {
                for (int j = 0; j < opts->height; j++) {
                    snprintf (filename, FILENAME_MAX, "sinogram-%05i-%05i.tif", i, j);
                    write_tiff (sinogram, filename, opts->width, j + 1, j);
                }
            }
            else {
                snprintf (filename, FILENAME_MAX, "sinogram-%05i.tif", i);
                write_tiff (sinogram, filename, opts->width, opts->height, 0);
            }
            free (sinogram);
        }
 
        animate_ellipses (ellipses, n_ellipses, opts);
    }
 
    options_free (opts);
    return 0;
}
