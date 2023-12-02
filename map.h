// Laura Toma
//
#ifndef __map_h
#define __map_h

#include "pixel_buffer.h"
#include "grid.h"


#include <stdlib.h>
#include <string.h>
#include <stdio.h> 
#include <time.h>
#include <math.h>
#include <assert.h>

#define PI 3.1415
#define FLAT_ASPECT 0
#define Z_FACTOR 5



typedef float Color[3];


// Laura Toma
//
#include "pixel_buffer.h"
#include "grid.h"
#include "map.h"


#include <stdlib.h>
#include <string.h>
#include <stdio.h> 
#include <time.h>
#include <math.h>




////////////////////////////////////////////////////////////
//returns the distance between (x,y) and (a,b)
float dist (u16 x, u16 y, u16 a, u16 b);



////////////////////////////////////////////////////////////
/*
a b c
d e f
g h i 
return nodata if e is nodata, or on th eboundary, or has a nodata neighbor
else: return (c + 2f + i - a - 2d -g)/8
 */
float dzdx(const Grid* grid, int r, int c);


////////////////////////////////////////////////////////////
/*
a b c
d e f
g h i 
return nodata if e is nodata, or on th eboundary, or has a nodata neighbor
else: return (g + 2h + i - a - 2b -c)/8
 */
float dzdy(const Grid* grid, int r, int c);

////////////////////////////////////////////////////////////
float slope_rad(const Grid* grid, float dzdx, float dzdy);

////////////////////////////////////////////////////////////
float aspect_rad(const Grid* grid, float dzdx, float dzdy);


//populate with values hillshade_grid
void generate_hillshade(const Grid* grid, Grid* hillshade_grid);

//populate with values slope_grid
void generate_slope(const Grid* grid, Grid* slope_grid);

//populate with values aspect_grid
void generate_aspect(const Grid* grid, Grid* aspect_grid);



////////////////////////////////////////////////////////////
//set the pixel buffer to correpsond to a grayscale image of grid
void grid_grayscale_to_pixelbuffer(const Grid* grid, PixelBuffer pb);


////////////////////////////////////////////////////////////
//helper function 
//return the index i of the interval values[i], values[i+1] that contains h 
// nb_intervals: nb intervals
//values: an array of nb_intervals + 1 values
int find_interval(float h, const Grid* grid, float* values, int nb_intervals);


////////////////////////////////////////////////////////////
//helper function 
//int_colors and values are arrays of size NB_COLORS_INTERVALS +1
void set_intervals_colors_with_lowground(const Grid* grid, Color* int_colors, float* values);



////////////////////////////////////////////////////////////
//set the pixel buffer to correpsond to color discrete intervals image
//of grid
void grid_discretecolorintervals_withlowarea_to_pixelbuffer(const Grid* grid, PixelBuffer pb);




///////////////////////////////////////////////////////////
//set the pixel buffer to correpsond to color discrete intervals image
//of grid
void grid_discretecolorintervals_to_pixelbuffer(const Grid* grid, PixelBuffer pb);






////////////////////////////////////////////////////////////
//set the pixels in the pixel buffer to correpsond to color gradient
//intervals image of grid
void grid_gradientcolorintervals_to_pixelbuffer(const Grid* grid, PixelBuffer pb);



////////////////////////////////////////////////////////////
//set the pixels in the pixel buffer to correspond to color gradient
//intervals image of grid, with the first interval set from grid->min_value : LOW_AREA_THRESHOLD
void grid_gradientcolorintervals_withlowarea_to_pixelbuffer(const Grid* grid, PixelBuffer pb);



////////////////////////////////////////////////////////////
//this grid contains values 0..255, put  those in a pixelbuffer
void grid_to_pixelbuffer(const Grid* grid, PixelBuffer pb);




///////////////////////////////////////////////////////////
// grid: each point in the grid is labeled with the SLR that floods
// it. color each value with a diff color, first BLUE colors, then
// random
void grid_flood_to_pixelbuffer(const Grid* grid, PixelBuffer pb, float rise_incr, float rise_end);



////////////////////////////////////////////////////////////
// overlay pb2 on top of pb1, with given transparency alpha. when
// alpha=1, show pb1; when alpha=0, show pb2. write the output in pb1
void pixelbuffer_overlay(PixelBuffer pb1, PixelBuffer pb2, float alpha);



///////////////////////////////////////////////////////////
// flow_grid: grid of flow accumulation  values
// nintervals: number of height intervals 
// color gradient of blue. height intervals are not equal but logarithmic 
void grid_flowaccu_to_pixelbuffer(const Grid* flow_grid, PixelBuffer pb);

//grayscale, but flat areas shown in red
void grid_flowdir_to_pixelbuffer(const Grid* fd_grid, PixelBuffer pb); 


//viewshed bitmap 
void grid_vis_to_pixelbuffer(const Grid* vis_grid, int vr, int vc, PixelBuffer pb);
  

////////////////////////////////////////////////////////////
// CONTOURS 
////////////////////////////////////////////////////////////
/*
  PARAM i, j:  row and col in grid
  PARAM h: height of (i,j) in the grid 
  PARAM grid: a grid
  PARAM iso: height of isoline 
  return 1 if this cell is at height h, or if  height h passes between this cell and one of its neighbors
*/ 
int cell_on_isoline(float iso, int i, int j, float h, const Grid* grid);


//helper function 
//populate countour_heights and contour_values with values
void setup_contours(const Grid* grid, int nb_contour_intervals, float* contour_heights, Color* contour_colors);


////////////////////////////////////////////////////////////
//overlay contours on this pixelbuffer
void overlay_contours_on_pixelbuffer(PixelBuffer pb, const Grid* grid, int nb_contour_intervals, float* contour_heights, Color* contour_colors);


////////////////////////////////////////////////////////////
//overlay contours on this pixelbuffercrtea two bitmaps, contours over
//grayscale and contours over hillshade
void grid_contours_to_pixelbuffer(const Grid* grid, int nb_contour_intervals, const Grid* hillshade_grid,  PixelBuffer pb);


#endif
