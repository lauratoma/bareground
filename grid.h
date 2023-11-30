// Laura Toma
//
#ifndef __grid_h
#define __grid_h

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


typedef struct _grid {
  unsigned int     ncols;
  unsigned int     nrows;
  float* data;
  
  float   xllcorner;
  float   yllcorner;
  float   cellsize;
  float   nodata_value;


  //computed 
  float   min_value;
  float   max_value;
  int nb_nodata_values; 
} Grid;


typedef struct gridpoint {
    unsigned int r, c;
} GridPoint;



#define FLAT_DIR -1
#define PIT_DIR  -2


// create a grid from file 
Grid* grid_read_from_file(char* filename); 

Grid* grid_clone(Grid* grid); 

  
//create a grid which is a clone of the grid passed as parameter
Grid* grid_init_from(Grid* grid);

//create a grid  of these specs
Grid* grid_init_from_specs(int nrows, int ncols,
			   float xllcorner, float yllcorner, float cellsize, float nodata_value,
			   float initial_value); 
//print basic info
void grid_print_stats(const Grid* g, char* name);

//write grid to file, arcascii format 
void  grid_write_to_file(char* filename, Grid* grid);

//free space 
void  grid_free(Grid* grid);


//return true if (r,c) is inside the grid , including the boundary 
static inline int inside_grid(const Grid* grid, int r, int c) {
    return r >= 0 && r < grid->nrows && c >= 0 && c < grid->ncols;
}

//return true if (r,c) is on the boundary 
static inline int on_boundary_grid(const Grid* grid, int r, int c) {
  return r == 0 || r == grid->nrows-1 || c == 0 || c == grid->ncols-1;
}

//return grid[r][c]
static inline float grid_get(const Grid* grid, int r, int c) {
  assert(inside_grid(grid, r, c));
  return grid->data[r * grid->ncols + c];
}
//return average value 
float grid_get_avg_value(const Grid* grid);

//return nb nodata values
int grid_count_nodata(const Grid* grid); 

  
//return true if grid(r,c) is nodata 
static inline bool  grid_is_nodata(const Grid* grid, int r, int c) {
  assert(inside_grid(grid, r, c)); 
  return (int)grid->nodata_value == (int)grid_get(grid, r, c);
}

//set grid[r][c]=val
static inline void grid_set(Grid* grid, int r, int c, const float val) {
  assert(inside_grid(grid, r, c)); 
  grid->data[r * grid->ncols + c] = val;
  if (val != grid->nodata_value) {
    grid->min_value = fmin(val, grid->min_value);
    grid->max_value = fmax(val, grid->max_value);
  } else {
    grid->nb_nodata_values++; 
  } 
}

//set grid(r,c) = nodata 
static inline void  grid_set_nodata(Grid* grid, int r, int c) {
  assert(inside_grid(grid, r, c)); 
  grid->data[r * grid->ncols + c] = grid->nodata_value;
  grid->nb_nodata_values++; 
} 

//interpolate nodata values from nearest neighbors 
void grid_interpolate_nodata(Grid* grid);

//recompute min, max, nodata 
void grid_reset_stats(Grid* grid); 

#endif
