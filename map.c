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


//to compare float values 
#define EPSILON 10 

//for maps that use color intervals
const int NB_COLOR_INTERVALS = 7;
//Color NODATA_COLOR = {1, 1, 1};  // white
Color NODATA_COLOR = {1, 1, 0};  // yellow


Color LOW_GROUND_COLOR1 =  {.68, .85, .90}; //light blue rgb(173,216,230) 
Color LOW_GROUND_COLOR2 =  {.64, .76, .68}; //cambridge blue 
float LOW_AREA_THRESHOLD = 1; //under this height, color with LOW_GROUND_COLOR

//colors used by interval color maps 
Color INTERVAL_COLOR[NB_COLOR_INTERVALS + 1] = {
   {0.9, .9, .6}, // yellow sand
   {.56, .83, 0}, //sheen green rgb(143,212,0)
  
   {.1, .5, 0},  //green
   {0, .44, .24}, //dartmouth green gb(0,112,60)
   
   {.13, .26, .12}, //myrtle rgb(33,66,30)
   
   // {.51, .41, .32}, //brown 
   //{1, .5, .3},  //orange
   {1, .2, .2}, // fire -orange

   //{0.7, 0.7, .4},   //yellow 
   //{.63, .32, .18}, //brown wood
   {0.3, 0.3, 0.3}, //dark gray 
   {0.1, .1, .1} //black 
};

const int NB_BLUE_COLORS = 4; 
Color BLUE[NB_BLUE_COLORS] = {
  {.68, .85, .90}, //light blue rgb(173,216,230) 
  {.31, .53, .97}, //royal rgb(79,134,247)
  {.15, .23, .89}, //palatinate rgb(39,59,226)
  //{.64, .76, .68} //cambridge blue
  {.4, .6, .68}
} ;


//for flow direction
Color FLAT_DIR_COLOR = {1, 0, 0};//red
Color PIT_DIR_COLOR = {1, .5, 0}; //orange? 


//for viewshed
Color VIS_COLOR = {0, 1, 0}; //green
Color INVIS_COLOR = {.6, .6, .6}; //gray
float VIEWPOINT_RADIUS =  .01;
Color VIEWPOINT_COLOR = {1, .2, .3}; 


//for hillshading
const float SUN_zenith_deg = 45;
const float SUN_azimuth_deg = 315;


//for contours 
//const int DEFAULT_NB_CONTOUR_INTERVALS = 5; 
const float FIRST_CONTOUR_COLOR[] = { .8, 0, 0}; 
const float LAST_CONTOUR_COLOR[] = {.1, 0, .9}; 


//for flooding
Color DANGER_ZONE_COLOR = {1, .2, .2};






float  dist(u16 x1, u16 y1, u16 x2, u16 y2) {

  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)); 
} 


////////////////////////////////////////////////////////////
/*
a b c
d e f
g h i 

if e is nodata :           return nodata 
if any neighbor is nodata: return nodata
if (r,c) on the byoundary: return nodata

else return (c + 2f + i - a - 2d -g)/(8*cell)

 */
float dzdx(const Grid* grid, int r, int c) {

  assert (inside_grid(grid, r, c)); 
  if (grid_is_nodata(grid, r, c)) return grid->nodata_value;

  //if on the boundary, set it to 0 
  if (r == 0 || r==grid->nrows-1 || c==0 || c==grid->ncols-1) return grid->nodata_value;
  
  float res =0;
  for (int l= -1; l<= 1; l++) {
    for (int j=-1; j<= 1; j++) {

      if (!inside_grid(grid, r+l, c+j)) continue;
      if (l==0 && j==0) continue;

      if (grid_is_nodata(grid, r+l, c+j)) return grid->nodata_value;
      
      if (l==-1 && j==-1) res -= grid_get(grid, r+l, c+j); //a
      if (l== 0 && j==-1) res -= 2*grid_get(grid, r+l, c+j); //2d
      if (l== 1 && j==-1) res -= grid_get(grid, r+l, c+j); //g

      if (l==-1 && j== 1) res += grid_get(grid, r+l, c+j); //c
      if (l== 0 && j== 1) res += 2*grid_get(grid, r+l, c+j); //2f
      if (l== 1 && j== 1) res += grid_get(grid, r+l, c+j); //i
    } 
  }
  return res/(8*grid->cellsize); 
}


////////////////////////////////////////////////////////////
/*
a b c
d e f
g h i 

if e is nodata :           return nodata 
if any neighbor is nodata: return nodata
if (r,c) on the byoundary: return nodata

else return  (g + 2h + i)- (a + 2b + c) )/(8*cell)

NOTE WE GO DY = SOUTH - NORTH TO GO CLOCKWISE, BECAUSE  SUN AZIMUTH IS MEASURED CLOCKWISE WITH 315 BEING NW
if we do dy  = north - south then we need to set sun azimuth as 135 instead of 315
 */
float dzdy(const Grid* grid, int r, int c) {

  assert (inside_grid(grid, r, c)); 
  if (grid_is_nodata(grid, r, c)) return grid->nodata_value;

  //if on the boundary, set it to 0 
  if (r == 0 || r==grid->nrows-1 || c==0 || c==grid->ncols-1) return grid->nodata_value;
  
  float res =0;
  for (int l= -1; l<= 1; l++) {
    for (int j=-1; j<= 1; j++) {
      if (!inside_grid(grid, r+l, c+j)) continue;
      if (l==0 && j==0) continue;

      //if any neighbor is nodata, return 0 
      if (grid_is_nodata(grid, r+l, c+j)) return grid->nodata_value;


      if (l== 1 && j==-1) res += grid_get(grid, r+l, c+j);//g
      if (l== 1 && j== 0) res += 2*grid_get(grid, r+l, c+j); //2h
      if (l== 1 && j== 1) res += grid_get(grid, r+l, c+j);////i
      
      if (l==-1 && j==-1) res -= grid_get(grid, r+l, c+j); //a
      if (l==-1 && j== 0) res -= 2*grid_get(grid, r+l, c+j);  //2b
      if (l==-1 && j== 1) res -= grid_get(grid, r+l, c+j); //c
    } 
  }
  return res/(8*grid->cellsize); 
}


////////////////////////////////////////////////////////////
float slope_rad(const Grid* grid, float dzdx, float dzdy) {

  if ((dzdx == grid->nodata_value) || (dzdy == grid->nodata_value)) return grid->nodata_value;
  
  float slope = atan(Z_FACTOR* sqrt(dzdx*dzdx + dzdy*dzdy));
  //printf("%.1f ", slope);
  return slope; 
} 

////////////////////////////////////////////////////////////
//The direction the steepest downslope direction is facing is the
//aspect. Aspect in radians is defined in the range of 0 to 2pi, with
//0 toward east. The aspect is determined under the rules in the
//following algorithm:
float aspect_rad(const Grid* grid, float dzdx, float dzdy) { 

  if ((dzdx == grid->nodata_value) || (dzdy == grid->nodata_value)) return grid->nodata_value;

  float aspect = 0;
  if (dzdx != 0) { 
    aspect = atan2(dzdy, dzdx); 
    if (aspect< 0) aspect  += 2 * PI; 
  } else  { 
    if (dzdy > 0)  aspect = PI / 2;
    else if (dzdy < 0) aspect = 2 * PI - PI / 2;
    else aspect = FLAT_ASPECT;  //slope zero, flat 
  }
  return aspect; 
}



////////////////////////////////////////////////////////////
void generate_hillshade(const Grid* grid,  Grid* hillshade_grid) {
  
  //for hillshading 
  //sun zenith in radians relative to vertical 
  float sun_zenith_rad = (90- SUN_zenith_deg) * PI/180;
 
  //sun azimuth in radians; change the azimuth angle from its
  //geographic unit (compass direction) to a mathematic unit (right
  //angle). The direction of the illumination source, azimuth, is
  //specified in degrees. The hillshade formula requires this angle to
  //be in units of radians. First, change the azimuth angle from its
  //geographic unit (compass direction) to a mathematic unit (right
  //angle). Next, convert the azimuth angle to radians.
  float sun_azimuth_math =  360 - SUN_azimuth_deg + 90;
  if (sun_azimuth_math > 360) sun_azimuth_math -= 360; 
  float sun_azimuth_rad = sun_azimuth_math *PI/180; 

  // printf("sun zenith rad = %.1f\n", sun_zenith_rad);
  //printf("sun azimuth rad = %.1f\n", sun_azimuth_rad);

  for (int i=0; i< grid->nrows; i++) {
    for (int j=0; j< grid->ncols; j++) {
      
      float dx = dzdx(grid, i, j);
      float dy = dzdy(grid, i, j);
      float terrain_slope_rad = slope_rad(grid, dx, dy);
      float terrain_aspect_rad = aspect_rad(grid, dx, dy); 
      float hillshade; 
      if (grid_is_nodata(grid, i, j) ||
	   (terrain_slope_rad == grid->nodata_value) ||
	  (terrain_aspect_rad == grid->nodata_value)) { 
	hillshade = grid->nodata_value;
      } else {
	//hillshade value at (i,j) 
	hillshade = cos(sun_zenith_rad) * cos(terrain_slope_rad) +
	  sin(sun_zenith_rad) * sin(terrain_slope_rad) * cos(sun_azimuth_rad - terrain_aspect_rad);
	
	//can this be bigger than 1?  no 
	//if (hillshade > 255) assert(0); 
	//printf("shade =%.1f ", hillshade);
	//can this be negative?  should not 
	if (hillshade<0)  hillshade = 0; 
	
	hillshade *= 255.0;
	
      } 
      grid_set(hillshade_grid, i, j, hillshade); 
    } 
  }
  //printf("some stats:\n");
  //printf("\thillshade range: [%.1f, %.1f]\n", hillshade_grid->min_value, hillshade_grid->max_value); 
}

////////////////////////////////////////////////////////////
void generate_slope(const Grid* grid, Grid* slope_grid) {
  

  for (int i=0; i< slope_grid->nrows; i++) {
    for (int j=0; j< slope_grid->ncols; j++) {
      float dx = dzdx(grid, i, j);
      float dy = dzdy(grid, i, j);
      float terrain_slope_rad = slope_rad(grid, dx, dy);
      grid_set(slope_grid, i, j, terrain_slope_rad);
    } 
  }
  printf("some stats:\n");
  printf("\tslope range: [%.1f, %.1f]\n", slope_grid->min_value, slope_grid->max_value);
}

void generate_aspect(const Grid* grid, Grid* aspect_grid) {
  

  for (int i=0; i< aspect_grid->nrows; i++) {
    for (int j=0; j< aspect_grid->ncols; j++) {
      float dx = dzdx(grid, i, j);
      float dy = dzdy(grid, i, j);
      float terrain_aspect_rad = aspect_rad(grid, dx, dy); 
      grid_set(aspect_grid, i, j, terrain_aspect_rad);
    } 
  }
  printf("some stats:\n");
  printf("\taspect range: [%.1f, %.1f]\n", aspect_grid->min_value, aspect_grid->max_value);
}







static const u8 max_rgb_value = 0xffu;

////////////////////////////////////////////////////////////
//set the pixel buffer to correpsond to a grayscale image of grid
void grid_grayscale_to_pixelbuffer(const Grid* grid, PixelBuffer pb) { 

    //write values in pb
  float color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(grid, y, x); //height of this pixel 
      if (grid_is_nodata(grid, y, x)) {
	color[0]= NODATA_COLOR[0]; color[1]= NODATA_COLOR[1];color[2]= NODATA_COLOR[2];  
      } else {
	//interpolate the color
	// hmin is set to {0,0,0} and hmax is set to {1,1,1}
	float c = (h - grid->min_value)/ (grid->max_value - grid->min_value);

	//hmin to {0,0,0} and hmax to {1,1,1} 
	color[0] = c; 
	color[1] = c; 
	color[2] = c; 
      }
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
}


////////////////////////////////////////////////////////////
//helper function 
//return the index i of the interval values[i], values[i+1] that contains h 
// nb_intervals: nb intervals
//values: an array of nb_intervals + 1 values
int find_interval(float h, const Grid* grid, float* values, int nb_intervals) {

  //assert((h  == grid->nodata_value) || ((h >= grid->min_value) && (h <= grid->max_value)));
  assert(h  != grid->nodata_value);
  assert((h >= grid->min_value) && (h <= grid->max_value)); 
  
  // float dh = (grid->max_value - grid->min_value)/NB_COLOR_INTERVALS; 
  for (int i=1; i<= nb_intervals; i++) {
    if (h <= values[i])  return i-1; 
  }
  //perhaps rounding error?
  if (h < values[nb_intervals] + EPSILON) return nb_intervals-1;
  
  //if we are here it's an error 
  printf("find_interval: h=%f, nb. intervals = %d, min=%.2f, max=%f\n",
	 h, nb_intervals, values[0], values[nb_intervals]);
  assert(0); 
} 



////////////////////////////////////////////////////////////
//helper function 
//int_colors and values are arrays of size NB_COLORS_INTERVALS +1
void set_intervals_colors_with_lowground(const Grid* grid, Color* int_colors, float* values) { 

  //copy INTERVAL_COLORS in int_colors 
  for (int i=0; i< NB_COLOR_INTERVALS+1; i++) {
    for (int j=0; j<3; j++)
      int_colors[i][j] = INTERVAL_COLOR[i][j]; 
  }
  
  //set values 
  if (grid->min_value < LOW_AREA_THRESHOLD) {
    //first interval 
    values[0] = grid->min_value;
    values[1] = LOW_AREA_THRESHOLD;
    //set first two colors to  low area colors 
    for (int i=0; i<3; i++) int_colors[0][i] = LOW_GROUND_COLOR1[i];
    for (int i=0; i<3; i++) int_colors[1][i] = LOW_GROUND_COLOR2[i];
    
    float dh = (grid->max_value - LOW_AREA_THRESHOLD) / (NB_COLOR_INTERVALS-1); 
    for (int i=2; i< NB_COLOR_INTERVALS; i++) {
      values[i] = values[i-1] + dh; 
    }
  } else {
    //there is no low area, divide the intervals equally
    float dh = (grid->max_value - grid->min_value) / NB_COLOR_INTERVALS; 
    values[0] = grid->min_value;
    for (int i=1; i< NB_COLOR_INTERVALS; i++) {
      values[i] = values[i-1] + dh; 
    }
  }
  //avoid  rounding errors 
  values[NB_COLOR_INTERVALS] = grid->max_value;

  //print
  //printf("set_interval_colors_with_lowground: \n");
  //printf("\tgrid range [%.1f, %.1f], intervals: [", grid->min_value, grid->max_value );
  //for (int i=0; i<= NB_COLOR_INTERVALS; i++)
  //  printf("%.2f ", values[i]);
  //printf("]\n"); 
} 



////////////////////////////////////////////////////////////
//set the pixel buffer to correpsond to color discrete intervals image
//of grid
void grid_discretecolorintervals_withlowarea_to_pixelbuffer(const Grid* grid, PixelBuffer pb) { 

  //set up the colors used, as a copy of the global INTERVAL COLORS
  //because we'll modify them if there is low area
  Color * int_colors = (Color*) malloc((NB_COLOR_INTERVALS +1) * sizeof(Color));
  // set up the intervals
  float* values = (float*) malloc((NB_COLOR_INTERVALS+1)*sizeof(float));

  printf("discrete_colorinterval_with_lowground_to_pixelbuffer: \n");
  set_intervals_colors_with_lowground(grid, int_colors, values); 
  printf("\tgrid range [%.1f, %.1f], intervals: [", grid->min_value, grid->max_value );
  for (int i=0; i<= NB_COLOR_INTERVALS; i++)
    printf("%.2f ", values[i]);
  printf("]\n");
  
  float color[3];
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(grid, y, x); //height of this pixel 

      if (grid_is_nodata(grid, y, x)) {
	//NODATA COLOR 
	color[0] = NODATA_COLOR[0]; color[1] = NODATA_COLOR[1]; color[2] = NODATA_COLOR[2]; 
      } else {
	
	//find its interval 
	int k = find_interval(h, grid, values, NB_COLOR_INTERVALS);
	assert(k >= 0 && k <= NB_COLOR_INTERVALS); 
	
	//set its color to the interval color 
	for (int i=0; i<3; i++) color[i] = int_colors[k][i]; 
	
      }
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
}


///////////////////////////////////////////////////////////
//set the pixel buffer to correpsond to color discrete intervals image
//of grid
void grid_discretecolorintervals_to_pixelbuffer(const Grid* grid, PixelBuffer pb) { 

  // set up the intervals
  float* values = (float*)malloc((NB_COLOR_INTERVALS+1)*sizeof(float)); 
  float dh = (grid->max_value - grid->min_value) / NB_COLOR_INTERVALS; 
  values[0] = grid->min_value;
  for (int i=1; i< NB_COLOR_INTERVALS; i++) {
    values[i] = values[i-1] + dh; 
  }
  //avoid  rounding errors 
  values[NB_COLOR_INTERVALS] = grid->max_value;

   //print
  printf("discrete_colorintervals_to_pixelbuffer: \n");
  printf("\tgrid range [%.1f, %.1f], intervals: [", grid->min_value, grid->max_value );
 //  printf("color intervals: [");
  for (int i=0; i<= NB_COLOR_INTERVALS; i++)
    printf("%.2f ", values[i]);
  printf("]\n");

  
  float color[3];
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(grid, y, x); //height of this pixel 

      if (grid_is_nodata(grid, y, x)) {
	color[0] = NODATA_COLOR[0]; color[1] = NODATA_COLOR[1]; color[2] = NODATA_COLOR[2]; //nodata 
      } else {
	
	//find its  interval
	int k = find_interval(h, grid, values, NB_COLOR_INTERVALS);
	assert(k >= 0 && k <= NB_COLOR_INTERVALS);
		
	//set its color to the interval color 
	for (int i=0; i<3; i++) color[i] = INTERVAL_COLOR[k][i]; 
	
      }
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
}


///////////////////////////////////////////////////////////
// grid:  each point in the grid is labeled with the SLR that floods
// it. color each value  with a diff color, first use BLUE[], then random
void grid_flood_to_pixelbuffer(const Grid* grid, PixelBuffer pb, float rise_incr, float rise_end) {

  //nb levels of slr. the values in the grid will be incr, 2*incr, ... 
  int levels = (rise_end - rise_incr)/ rise_incr + 1; 
  //printf("nb levels = %d\n", levels ); 

  Color* mycolors = (Color*)malloc(levels*sizeof(Color)); 
  for (int i = 0; i < NB_BLUE_COLORS; i++) {
    mycolors[i][0] = BLUE[i][0]; 
    mycolors[i][1] =  BLUE[i][1]; 
    mycolors[i][2] =  BLUE[i][2]; 
  } 
  for (int i = NB_BLUE_COLORS; i <= levels; i++) {
    mycolors[i][0] = (float)rand() / ((float)RAND_MAX);
    mycolors[i][1] = (float)rand() / ((float)RAND_MAX);
    mycolors[i][2] = .7* (float)rand() / ((float)RAND_MAX);
    
  } 
  
  float color[3];
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(grid, y, x); //flood level  of this pixel 
          
      if (grid_is_nodata(grid, y, x)) {
	color[0] = NODATA_COLOR[0]; color[1] = NODATA_COLOR[1]; color[2] = NODATA_COLOR[2]; //nodata 
      } else {
	int k; 
	for (k =0; k<= levels; k++) {
	  if (h<= k*rise_incr) break; 
	}
	k--; 
	assert(k <= levels); 

	  //use color k 
	for (int i=0; i<3; i++) color[i] = mycolors[k][i]; 
      }//else 
    
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);

    } //for x
  }//for y 
}





////////////////////////////////////////////////////////////
//set the pixels in the pixel buffer to correpsond to color gradient
//intervals image of grid
void grid_gradientcolorintervals_to_pixelbuffer(const Grid* grid, PixelBuffer pb) { 

  // set up the intervals
  float* values = (float*) malloc((NB_COLOR_INTERVALS+1)*sizeof(float)); 
  float dh = (grid->max_value - grid->min_value) / NB_COLOR_INTERVALS; 
  values[0] = grid->min_value;
  for (int i=1; i<= NB_COLOR_INTERVALS; i++) {
    values[i] = values[i-1] + dh; 
  }
   //print
   printf("gradient_colorintervals_to_pixelbuffer: \n");
   printf("\tgrid range [%.1f, %.1f], intervals: [", grid->min_value, grid->max_value );
   // printf("color intervals: [");
  for (int i=0; i<= NB_COLOR_INTERVALS; i++)
    printf("%.2f ", values[i]);
  printf("]\n");
  
    
  float color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(grid, y, x); //height of this pixel
      
      if (grid_is_nodata(grid, y, x)) {
	for (int i=0; i<3; i++) color[i] = NODATA_COLOR[i]; //nodata
      } else {
	
	//find the interval that this color is in 
	int k = find_interval(h, grid, values, NB_COLOR_INTERVALS);
	assert(k >= 0 && k <= NB_COLOR_INTERVALS);
	
	//interpolate between INTERVAL_COLOR[x] and INTERVAL_COLOR[x+1]
	float h1 = values[k]; 
	float h2 =  values[k+1]; 
	float c = (h - h1)/ (h2-h1);
	//printf("h1=%.1f, h2=%.1f, h=%.1f, c=%.1f\n", h1, h2, h, c); 
	
	//set its color to the interval color 
	for (int i=0; i<3; i++)
	  color[i] = INTERVAL_COLOR[k][i] + c * (INTERVAL_COLOR[k+1][i] - INTERVAL_COLOR[k][i]); 
      }
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
}



////////////////////////////////////////////////////////////
//set the pixels in the pixel buffer to correspond to color gradient
//intervals image of grid, with the first interval set from grid->min_value : LOW_AREA_THRESHOLD
void grid_gradientcolorintervals_withlowarea_to_pixelbuffer(const Grid* grid, PixelBuffer pb) { 

  //set up the colors used, as a copy of the global INTERVAL COLORS
  //because we'll modify them if there is low area
  Color * int_colors = (Color*) malloc((NB_COLOR_INTERVALS +1) * sizeof(Color));
  // set up the intervals
  float* values = (float*) malloc((NB_COLOR_INTERVALS+1)*sizeof(float));

  set_intervals_colors_with_lowground(grid, int_colors, values); 
  printf("gradient_colorinterval_with_lowground_to_pixelbuffer: \n");
  printf("\tgrid range [%.1f, %.1f], intervals: [", grid->min_value, grid->max_value );
  for (int i=0; i<= NB_COLOR_INTERVALS; i++)
    printf("%.2f ", values[i]);
  printf("]\n");
  
  
  float color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(grid, y, x); //height of this pixel
      
      if (grid_is_nodata(grid, y, x)) {
	for (int i=0; i<3; i++) color[i] = NODATA_COLOR[i]; //nodata 
      } else {
	//find the interval that this color is in 
	int k = find_interval(h, grid, values, NB_COLOR_INTERVALS);
	assert(k >= 0 && k <= NB_COLOR_INTERVALS);
	
	//interpolate between INTERVAL_COLOR[x] and INTERVAL_COLOR[x+1]
	float h1 = values[k]; 
	float h2 =  values[k+1]; 
	float c = (h - h1)/ (h2-h1);
	//printf("h1=%.1f, h2=%.1f, h=%.1f, c=%.1f\n", h1, h2, h, c); 
	
	//set its color to the interval color 
	for (int i=0; i<3; i++)
	  color[i] = int_colors[k][i] + c * (int_colors[k+1][i] - int_colors[k][i]); 
      }
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }

}



////////////////////////////////////////////////////////////
//this grid contains values 0..255, put  those in a pixelbuffer
void grid_to_pixelbuffer(const Grid* grid, PixelBuffer pb) { 
  //write values in pb
  u16 rgb_color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float shade = grid_get(grid, y, x); //shade  of this pixel 
      if (grid_is_nodata(grid, y, x)) {
	//set this pixel to nodata color 
	for (int i=0; i<3; i++) rgb_color[i] = (u16)255*NODATA_COLOR[i]; 
      } else {
	//set this pixel to its value 
	rgb_color[0] = (u16) shade; rgb_color[1] = (u16) shade; rgb_color[2] = (u16) shade; 
      }
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);

    }//for x
  }//for y 
}



////////////////////////////////////////////////////////////
//grayscale, but flat areas shown in red
void grid_flowdir_to_pixelbuffer(const Grid* fd_grid, PixelBuffer pb) { 
  //write values in pb
  float color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      int  dir = grid_get(fd_grid, y, x); //value of this pixel 
      if (grid_is_nodata(fd_grid, y, x)) {
	color[0]= NODATA_COLOR[0]; color[1]= NODATA_COLOR[1];color[2]= NODATA_COLOR[2];  
      } else if (dir == FLAT_DIR) {
	for (int i=0; i<3; i++) color[i] = FLAT_DIR_COLOR[i]; 
      } else if (dir==PIT_DIR) {
	for (int i=0; i<3; i++) color[i] = PIT_DIR_COLOR[i]; 
      } else { 
      	//interpolate the color
	// hmin is set to {0,0,0} and hmax is set to {1,1,1}
	float c = (dir - fd_grid->min_value)/ (fd_grid->max_value - fd_grid->min_value);
	
	//hmin to {0,0,0} and hmax to {1,1,1} 
	color[0] = c; 
	color[1] = c; 
	color[2] = c; 
      }
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
}




///////////////////////////////////////////////////////////
// flow_grid: grid of flow accumulation  values
// nintervals: number of height and color intervals 
// color gradient of blue. height intervals are not equal but logarithmic 
void grid_flowaccu_to_pixelbuffer(const Grid* flow_grid, PixelBuffer pb) {

  int nintervals = 4;
  
  // set up the values intervals
  float* values = (float*) malloc((nintervals+1)*sizeof(float)); 

  //float range = (flow_grid->max_value - flow_grid->min_value);
  float maxfa = flow_grid->max_value; 
  float logn = log(maxfa);
  //double n = (double)flow_grid->nrows * flow_grid->ncols; 
  printf("grid_flow_accu_to_pixelbuffer:\n\tmax_fa=%.1f, logn=%.1f\n", maxfa, logn); 
  values[0] = flow_grid->min_value;
  values[1] = 1.1;
  values[2] = 5*log(maxfa); // 2 for flat areas; 5 for steep mountainous area
  values[3] = 2*values[2];  //
  values[4] = maxfa; 
  
   //print  
  printf("\tgrid range [%.1f, %.1f]", flow_grid->min_value,flow_grid->max_value );
  printf(" flow intervals: [");
  for (int i=0; i<= nintervals; i++)
    printf("%.2f ", values[i]);
  printf("]\n");

  //set up colors of intervals 
  Color* mycolors = (Color*)malloc((nintervals+1)*sizeof(Color));

  //first  interval is for points on ridges 
  //mycolors[0][0] = 1;  mycolors[0][1] = 1;  mycolors[0][2] = 1;//white
  //mycolors[0][0] = 1;  mycolors[0][1] = .6;  mycolors[0][2] = .6;//peach  (255, 153, 153)
  //mycolors[0][0] = 0;  mycolors[0][1] = .4;  mycolors[0][2] = 0;//green  (255, 153, 153)
  // mycolors[0][0] = 102/255.0;  mycolors[0][1] = 204.0/255;  mycolors[0][2] = 0;//green  (255, 153, 153)
   mycolors[0][0] = 210/255.0;  mycolors[0][1] = 180/255.0;  mycolors[0][2] = 140/255.0;//tan
  
  //second interval is for mediumvery low  values of flow 
  //  mycolors[1][0] = .8; mycolors[1][1] = 1; mycolors[1][2] = 1;
  mycolors[1][0] = .8; mycolors[1][1] = .8; mycolors[1][2] = .8;
    
  //third color is for smaller rivers  medium blue 
  mycolors[2][0] = 0; mycolors[2][1] = .5; mycolors[2][2] = 1;
  
  //fourth color is for larger flow in [values[3]: values[4]]  darker blue 
  mycolors[3][0] = 0; mycolors[3][1] = 0; mycolors[3][2] = .7;
  
  
  
  float color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      float h = grid_get(flow_grid, y, x); //height of this pixel
      
      if (grid_is_nodata(flow_grid, y, x)) {
	for (int i=0; i<3; i++) color[i] = NODATA_COLOR[i]; //nodata
      } else {
	
	int k; 
	if (h <= values[1]) k = 0;
	else if (h <= values[2]) k = 1;
	else if (h <= values[3]) k = 2;
	else  k = 3; 
	
	//set its color to the interval color 
	for (int i=0; i<3; i++) color[i] = mycolors[k][i]; 
      } 
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }//for x
  }//for y 
}


///////////////////////////////////////////////////////////
void grid_vis_to_pixelbuffer(const Grid* vis_grid, int vr, int vc, PixelBuffer pb) {

  float rad = (float)VIEWPOINT_RADIUS * (float)vis_grid->nrows; 
  printf("rad = %.1f\n", rad);
  
  //write values in pb
  float color[3]; 
  for (u16 y=0; y < pb.height; y++) {
    for (u16 x=0; x< pb.width; x++) {
      
      //int  vis = grid_get(vis_grid, y, x); //value of this pixel 
      if (grid_is_nodata(vis_grid, y, x)) {
	for (int i=0; i<3; i++) color[i]= NODATA_COLOR[i]; 
      } else if (dist((u16)vr, (u16)vc, (u16)y, (u16)x) < rad) { 
	  for (int i=0; i<3; i++) color[i] = VIEWPOINT_COLOR[i];
      } else { 
	float vis = grid_get(vis_grid, y, x);
	if (vis == 1) {
	  //visible
	  for (int i=0; i<3; i++) color[i] = VIS_COLOR[i]; 
	} else {
	  //invisible
	  for (int i=0; i<3; i++) color[i] = INVIS_COLOR[i]; 
	}
      } 
      
      //bring it from [0,1] to [0, 255]
      u16 rgb_color[3] = {
	(u16) (color[0] * max_rgb_value),
	(u16) (color[1] * max_rgb_value),
	(u16) (color[2] * max_rgb_value)}; 
      
      write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
} 


////////////////////////////////////////////////////////////
// overlay pb2 on top of pb1, with given transparency alpha. when
// alpha=1, show pb1; when alpha=0, show pb2. write the output in pb1
void pixelbuffer_overlay(PixelBuffer pb1, PixelBuffer pb2, float alpha) {

  assert(alpha >=0 && alpha <=1);

  //shoudl also assert that they have same size
  assert(pb1.height == pb2.height && pb1.width == pb2.width);
  
   for (u16 y=0; y < pb1.height; y++) {
    for (u16 x=0; x< pb1.width; x++) {
      
      const u32 pixel1 = pb1.pixels[y * pb1.width + x];
      const u32 pixel2 = pb2.pixels[y * pb2.width + x];
      
      //get colors 
      //static const u32 alpha_mask = 0xffu << 24u;
      //const u32 pixel = alpha_mask | (u32) (b << 16u) | (u32) (g << 8u) | r;
      u16 color1[3], color2[3];
      u32 red_mask = 0xff;
      u32 green_mask = 0xff << 8u;
      u32 blue_mask = 0xff << 16u; 
      color1[0] = pixel1 & red_mask; 
      color1[1] = (pixel1 & green_mask) >> 8u; 
      color1[2] = (pixel1 & blue_mask) >> 16u; 
      //printf("red mask = %ux green mask = %ux blue mask = %ux\n", red_mask, green_mask, blue_mask);
      //printf("pixel1=%x, r=%x g=%x b=%x\n", pixel1, color1[0], color1[1], color1[2] ); 

      color2[0] = pixel2 & red_mask; 
      color2[1] = (pixel2 & green_mask) >> 8u; 
      color2[2] = (pixel2 & blue_mask) >> 16u; 

     
      u16 rgb_color[3];  
      for (int i=0; i< 3; i++) {
	rgb_color[i] = color2[i] + (u16) (alpha * (float) (color1[i] - color2[i])); 
      }

        write_pixel_to_buffer(&pb1, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
    }
  }
  

} 



////////////////////////////////////////////////////////////
/*
  PARAM i, j:  row and col in grid
  PARAM h: height of (i,j) in the grid 
  PARAM grid: a grid
  PARAM iso: height of isoline 
  return 1 if this cell is at height h, or if  height h passes between this cell and one of its neighbors
*/ 
int cell_on_isoline(float iso, int i, int j, float h, const Grid* grid) {

  //if this is nodata, return 0 
  if (grid_is_nodata(grid, i, j))  return 0;
  
  //if cell is at this height, then it is on the isoline 
  if (h == iso) return 1;

  //if the isoline goes between this cell and one if its neighbors, then return true
  for (int l = -1; l<= 1; l++) {
    for (int k = -1; k <= 1; k++)  {
      if (l==0 && k==0) continue;
      if (!inside_grid(grid, i+l, j+k)) continue;
      float hneighb = grid_get(grid,i+l, j+k);
      if (h < iso && iso < hneighb) return 1;

      //this creates a double line which is too think for small sets but  good for large sets
      if (hneighb < iso && iso < h) return 1; 
    }//for k
  }//for l 
  
  return 0; 
} 



////////////////////////////////////////////////////////////
//helper function 
//populate countour_heights and contour_values with values
void setup_contours(const Grid* grid, int nb_contour_intervals, float* contour_heights, Color* contour_colors) { 
  
  printf("contour lines:  hmin = %.1f, hmax = %.1f, countour intervals=%d\n", \
	 grid->min_value, grid->max_value,  nb_contour_intervals);
 
  
  printf("\tcontour lines at heights: "); 
  for (int i=0; i< nb_contour_intervals; i++) {
    contour_heights[i] =  grid->min_value + (grid->max_value - grid->min_value)/nb_contour_intervals/2; 
    contour_heights[i] += 	i*(grid->max_value - grid->min_value)/nb_contour_intervals; 
    printf("%.2f ", contour_heights[i]);
  }
  printf("\n");

  //set up contour colors: the colors of the contour lines
  //interpolate between first color and last color, based on the
  //height
  
  float eps[3];
  for (int i=0; i< 3; i++) 
    eps[i] = (float)(LAST_CONTOUR_COLOR[i] -  FIRST_CONTOUR_COLOR[i])/(nb_contour_intervals-1);
  for (int l=0; l< nb_contour_intervals; l++) {
    for (int i=0; i< 3; i++) 
      contour_colors[l][i] = FIRST_CONTOUR_COLOR[i] + l* eps[i] ;
    //printf("contour color[%d] = %.1f, %.1f, %.1f  ", l,
    //contour_colors[l][0], contour_colors[l][1], contour_colors[l][2] );
    //printf("\n"); 
  }
}




////////////////////////////////////////////////////////////
//overlay contours on this pixelbuffer
void overlay_contours_on_pixelbuffer(PixelBuffer pb, const Grid* grid, int nb_contour_intervals, float* contour_heights, Color* contour_colors) { 
  
    float color[3]; 
    for (u16 y=0; y < pb.height; y++) {
      for (u16 x=0; x< pb.width; x++) {
	
	float h = grid_get(grid, y, x); //height of this pixel 
	//checking if point (x, y, h) is on isoline i then set its color to isoline color 
	for (int i=0; i< nb_contour_intervals; i++) {
	  if (cell_on_isoline(contour_heights[i], y, x, h, grid)) {
	    color[0] = contour_colors[i][0];
	    color[1] = contour_colors[i][1];
	    color[2] = contour_colors[i][2];
	    
	    //bring it from [0,1] to [0, 255]
	    u16 rgb_color[3] = {
	      (u16) (color[0] * max_rgb_value),
	      (u16) (color[1] * max_rgb_value),
	      (u16) (color[2] * max_rgb_value)}; 
	    
	    write_pixel_to_buffer(&pb, x, y, rgb_color[0], rgb_color[1], rgb_color[2]);
	  }
	} //for int i  
      }///for x
    }//for y
}




void grid_contours_to_pixelbuffer(const Grid* grid, int nb_contour_intervals, const Grid* hillshade_grid,  PixelBuffer pb) {

  float* contour_heights = NULL;
  Color* contour_colors  = NULL;
  contour_heights = (float*) malloc(nb_contour_intervals * sizeof(float));
  assert(contour_heights);
  contour_colors = (Color*) malloc(nb_contour_intervals *sizeof(Color));
  assert(contour_colors);
  //populate with values 
  setup_contours(grid, nb_contour_intervals, contour_heights, contour_colors); 

   //set buffer for grayscale grid 
  grid_grayscale_to_pixelbuffer(grid, pb);
  //overlay contours 
  overlay_contours_on_pixelbuffer(pb, grid, nb_contour_intervals, contour_heights, contour_colors);
  printf("writing map.contours.bmp\n"); 
  save_pixel_buffer_to_file(&pb, "map.contours.bmp");
  
    
  //set buffer for hillshade grid 
  grid_to_pixelbuffer(hillshade_grid, pb);
  //overlay contours 
  overlay_contours_on_pixelbuffer(pb, grid, nb_contour_intervals, contour_heights, contour_colors);
  printf("writing map.contours-hilshade.bmp\n"); 
  save_pixel_buffer_to_file(&pb, "map.contours-hillshade.bmp");
  
}

