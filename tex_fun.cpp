/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include <cstdlib>

#define PI (float) 3.14159265358979323846

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

int tex_fun(float u, float v, GzColor color)
{
    unsigned char		pixel[3];
    unsigned char     dummy;
    char  		foo[8];
    int   		i, j;
    FILE* fd;

    if (reset) {          /* open and load texture file */
        fd = fopen("texture", "rb");
        if (fd == NULL) {
            fprintf(stderr, "texture file not found\n");
            exit(-1);
        }
        fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
        image = (GzColor*)malloc(sizeof(GzColor) * (xs + 1) * (ys + 1));
        if (image == NULL) {
            fprintf(stderr, "malloc for texture image failed\n");
            exit(-1);
        }

        for (i = 0; i < xs * ys; i++) {	/* create array of GzColor values */
            fread(pixel, sizeof(pixel), 1, fd);
            image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
            image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
            image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
        }

        reset = 0;          /* init is done */
        fclose(fd);
    }

    /* bounds-test u,v to make sure nothing will overflow image array bounds */
    u = min(max(u, 0.0f), 1.0f);
    v = min(max(v, 0.0f), 1.0f);

    /* Scale (u, v) to image size */
    float scaledU = u * (xs - 1);
    float scaledV = v * (ys - 1);

    /* Get the integer pixel positions */
    int u0 = (int)floor(scaledU);
    int u1 = min(u0 + 1, xs - 1);
    int v0 = (int)floor(scaledV);
    int v1 = min(v0 + 1, ys - 1);

    /* Calculate the interpolation weights */
    float s = scaledU - u0;
    float t = scaledV - v0;

    /* Get the four corner colors */
    GzColor color00, color01, color10, color11;

    for (int k = 0; k < 3; ++k) {
        color00[k] = image[u0 + v0 * xs][k];
        color01[k] = image[u0 + v1 * xs][k];
        color10[k] = image[u1 + v0 * xs][k];
        color11[k] = image[u1 + v1 * xs][k];
    }

    /* Perform bilinear interpolation */
    for (int k = 0; k < 3; ++k) {
        color[k] = (1 - s) * (1 - t) * color00[k] +
            (1 - s) * t * color01[k] +
            s * (1 - t) * color10[k] +
            s * t * color11[k];
    }

    return GZ_SUCCESS;
}
int ptex_fun(float u, float v, GzColor color)
{
    /* bounds-test u,v to make sure nothing will overflow image array bounds */
    u = min(max(u, 0.0f), 1.0f);
    v = min(max(v, 0.0f), 1.0f);

    /* Define parameters for the spiral pattern */
    float spiral_frequency = 10.0f;  /* Frequency of the spiral */
    float spiral_width = 0.1f;       /* Width of the spiral */

    /* Calculate angle and radius for the spiral pattern */
    float angle = atan2(v - 0.5f, u - 0.5f) / (2.0f * PI); /* Normalize angle between -0.5 and 0.5 */
    float radius = sqrt((u - 0.5f) * (u - 0.5f) + (v - 0.5f) * (v - 0.5f)); /* Distance from center (0.5, 0.5) */

    /* Create a spiral pattern using a sine wave */
    float spiral = 0.5f + 0.5f * sin(spiral_frequency * angle - radius / spiral_width * 2.0f * PI);

    /* Assign color based on the spiral pattern value */
    color[RED] = spiral;
    color[GREEN] = 1.0f - spiral;
    color[BLUE] = spiral;

    return GZ_SUCCESS;
}


/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

