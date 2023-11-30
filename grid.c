// Laura Toma
//
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "grid.h"



// Returns an empty grid object
static Grid* grid_init() {
  Grid* grid;
  grid = (Grid*) malloc(sizeof(Grid));
  assert(grid);
  grid->data = NULL;
  grid->nrows = grid->ncols = 0; 
  grid->min_value = INT_MAX;
  grid->max_value = -INT_MAX;
  grid->nb_nodata_values = 0; 
  return grid;
}


//helper function 
//read header from file into the grid object; grid must be allocated 
static void grid_read_header(FILE* in_file, Grid* grid) {

  assert(grid);
  int res = fscanf(in_file, "ncols %d\n",        &grid->ncols);
  assert(res > 0); 
  res = fscanf(in_file, "nrows %d\n",        &grid->nrows);
  assert(res > 0); 
  res = fscanf(in_file, "xllcorner %f\n",    &grid->xllcorner);
  assert(res > 0); 
  res = fscanf(in_file, "yllcorner %f\n",    &grid->yllcorner);
  assert(res > 0); 
  res = fscanf(in_file, "cellsize %f\n",     &grid->cellsize);
  assert(res > 0); 
  res = fscanf(in_file, "NODATA_value %f\n", &grid->nodata_value);
  assert(res > 0); 
}


// helper function 
// allocate space for the grid data.
static void grid_malloc_data(Grid* grid) {
  assert(grid);
  assert(grid->data == NULL); 
  if (!grid->data) {
    grid->data = (float*)malloc(grid->nrows * grid->ncols * sizeof(float));
    assert(grid->data);
  }
}
// helper function 
// allocate space for the grid data.
static void grid_malloc_data_initialize(Grid* grid, float initial_value) {
  assert(grid); 
  assert(grid->data == NULL); 
  if (!grid->data) {
    grid->data = (float*)malloc(grid->nrows * grid->ncols*sizeof(float));
    assert(grid->data);
    for (int i=0; i< grid->nrows * grid->ncols; i++)
      grid->data[i] = initial_value; 

    if (initial_value != grid->nodata_value) {
      grid->min_value = initial_value;
      grid->max_value = initial_value; 
    }
  }
}

//helper function 
static void grid_read_data(FILE* in_file, Grid* grid) {
  float val;
  for (int r = 0; r < grid->nrows; r++) {
    for (int c = 0; c < grid->ncols; c++) {
      int res = fscanf(in_file, "%f ", &val);
      assert(res > 0); 
      grid_set(grid, r, c, val);
    }//for c
  }//for r
}


//helper function 
// Read grid from file and return it 
static Grid* grid_read(FILE* in_file) {
  Grid* grid = grid_init();
  grid_read_header(in_file, grid);
  grid_malloc_data(grid);
  grid_read_data(in_file, grid);
  return grid;
}


/* ******************************************************************** */
//create a grid and fill it with data read from file 
Grid* grid_read_from_file(char* filename) {

  FILE* infile = fopen(filename, "r");
  if (!infile) {
    printf("cannot open file %s\n", filename);
    exit(1); 
  } 
  Grid* grid = grid_read(infile);
  fclose(infile);
  return grid; 
} 
/* ******************************************************************** */





/* ******************************************************************** */
// print stats about the grid 
void grid_print_stats(const Grid* g, char* name) {

  //compute the average
  float avg = grid_get_avg_value(g); 
  int nnodata = grid_count_nodata(g); 
  printf("grid %s (%p): \tn=%ld [rows=%d,cols=%d], range=[%.2f, %.2f], avg value=%.1f "
	 "NODATA=%d (%.1f%%)\n",
	 name, g, 
	 (long int)g->nrows*(long int)g->ncols, g->nrows, g->ncols, 
	 g->min_value, g->max_value,  avg,
	 nnodata, 100.0*((double)nnodata)/((double)g->nrows*g->ncols));
	 //100.0*((double)g->nb_nodata_values)/((double)g->nrows*g->ncols));
}


//return nb nodata values
int grid_count_nodata(const Grid* grid) { 
  assert(grid);
  int nnodata=0; 
  for (int r=0; r<grid->nrows; r++) {
    for (int c=0; c<grid->ncols; c++) {
      if (grid_is_nodata(grid, r, c)) {
	nnodata++;
      }
    }
  }
  return nnodata;  
}
    
//return average value 
float grid_get_avg_value(const Grid* grid) {

  //compute the average
  float sum=0;
  int nvalues = 0; 
  for (int i=0; i<grid->nrows*grid->ncols; i++) {
    if (grid->data[i]!= grid->nodata_value) {
      sum+=grid->data[i];
      nvalues++; 
    }//int i 
  }
  
  //return sum /(float)(grid->nrows*grid->ncols);
  return sum /(float)(nvalues); 
}


/* ******************************************************************** */


  
//helper function 
static void grid_write_header(FILE* out_file, Grid* grid) {
  fprintf(out_file, "ncols %d\n",        grid->ncols);
  fprintf(out_file, "nrows %d\n",        grid->nrows);
  fprintf(out_file, "xllcorner %f\n",    grid->xllcorner);
  fprintf(out_file, "yllcorner %f\n",    grid->yllcorner);
  fprintf(out_file, "cellsize %f\n",     grid->cellsize);
  fprintf(out_file, "NODATA_value %f\n", grid->nodata_value);
}


// helper function 
// Write the complete asc file for a grid.
static void grid_write(FILE* out_file, Grid* grid) {

  grid_write_header(out_file, grid);
  for (int r = 0; r < grid->nrows; r++) {
    for (int c = 0; c < grid->ncols; c++) {
      fprintf(out_file, "%f ", grid_get(grid, r, c));
    }
    fprintf(out_file, "\n");
  }
}

 

/* ******************************************************************** */
// Write the complete asc file for a grid.
void grid_write_to_file(char* filename, Grid* grid) {

  FILE* outfile = fopen(filename, "w+");
  if (!outfile) {
    printf("cannot open file %s\n", filename);
    exit(1); 
  } 
  grid_write(outfile, grid);
  fclose(outfile); 
} 
/* ******************************************************************** */



// helper function  
static void grid_copy_header(Grid* grid, Grid* new_grid) {
  assert (grid && new_grid); 
  new_grid->ncols =        grid->ncols;
  new_grid->nrows =        grid->nrows;
  new_grid->xllcorner =    grid->xllcorner;
  new_grid->yllcorner =    grid->yllcorner;
  new_grid->cellsize =     grid->cellsize;
  new_grid->nodata_value = grid->nodata_value;
}



/* ******************************************************************** */
// Initialize a grid based on an existing grid, copying e.g. its dimensions.
Grid* grid_init_from(Grid* grid) {
  //Grid* new_grid = grid_init();
  assert(grid && grid->data); 
  Grid* new_grid = (Grid*) malloc(sizeof(Grid));
  assert(new_grid); 
  //grid_copy_header(grid, new_grid);
  new_grid->nrows = grid->nrows;
  new_grid->ncols = grid->ncols;
  new_grid->xllcorner = grid->xllcorner;
  new_grid->yllcorner = grid->yllcorner;
  new_grid->cellsize = grid->cellsize;
  new_grid->nodata_value = grid->nodata_value;
  
  new_grid->min_value = INT_MAX;
  new_grid->max_value = -INT_MAX;
  new_grid->nb_nodata_values = 0;
  
  // grid_malloc_data(new_grid);
  printf("grid_init_from: allocating %d, %d\n", grid->nrows, grid->ncols); 
  new_grid->data = (float*)malloc(grid->nrows * grid->ncols * sizeof(float));
  assert(new_grid->data); 
  return new_grid;
}
/* ******************************************************************** */



/* ******************************************************************** */
// Initialize a grid based on an existing grid, copying e.g. its dimensions.
Grid* grid_clone(Grid* grid) {
  //Grid* new_grid = grid_init();
  assert(grid && grid->data); 
  Grid* new_grid = (Grid*) malloc(sizeof(Grid));
  assert(new_grid); 
  //grid_copy_header(grid, new_grid);
  new_grid->nrows = grid->nrows;
  new_grid->ncols = grid->ncols;
  new_grid->xllcorner = grid->xllcorner;
  new_grid->yllcorner = grid->yllcorner;
  new_grid->cellsize = grid->cellsize;
  new_grid->nodata_value = grid->nodata_value;
  
  new_grid->min_value = INT_MAX;
  new_grid->max_value = -INT_MAX;
  new_grid->nb_nodata_values = 0;
  
  // grid_malloc_data(new_grid);
  printf("grid_init_from: allocating %d, %d\n", grid->nrows, grid->ncols); 
  new_grid->data = (float*)malloc(grid->nrows * grid->ncols * sizeof(float));
  assert(new_grid->data);
  for (int i=0; i< grid->nrows * grid->ncols; i++) { 
    new_grid->data[i] = grid->data[i];
  }
  return new_grid;
}
/* ******************************************************************** */




//recompute min, max, nodata 
void grid_reset_stats(Grid* grid) {

  printf("reset_stats: start\n");
  assert(grid && grid->data);
  grid->min_value = INT_MAX;
  grid->max_value = -INT_MAX;
  grid->nb_nodata_values = 0;
  
  for (int i=0; i< grid->nrows * grid->ncols; i++) { 
    float h = grid->data[i]; 
    //assert(h == grid->nodata_value || h >= grid->min_value);
    if (h!= grid->nodata_value && h < 10) {
      printf("ERROR reset_stats: i=%d h=%f\n", i, h);
      assert(0); 
      exit(1); 
    } 
    
    if (h==grid->nodata_value) {
      grid->nb_nodata_values++;
    } else {
      
      if (h < grid->min_value) grid->min_value = h; 
      if (h > grid->max_value) grid->max_value = h; 
    }
  }//for i
  printf("\tmin set to %f, max set to %f\n", grid->min_value, grid->max_value);
  printf("reset_stats: end\n");
} 


/* ******************************************************************** */
//create a grid of these specs 
Grid* grid_init_from_specs(int nrows, int ncols,
			   float   xllcorner, float yllcorner, float cellsize,
			   float nodata_value,  float initial_value) {

  //Grid* grid = grid_init();
  Grid* grid = (Grid*) malloc(sizeof(Grid));
  grid->nrows = nrows;
  grid->ncols = ncols;
  grid->xllcorner = xllcorner;
  grid->yllcorner = yllcorner;
  grid->cellsize = cellsize;
  grid->nodata_value = nodata_value;

 
  //initialize  to initial_value 
  //grid_malloc_data_initialize(grid, initial_value);
  grid->data = (float*)malloc(grid->nrows * grid->ncols * sizeof(float));
  assert(grid->data);
  for (int i=0; i< grid->nrows * grid->ncols; i++) {
    grid->data[i] = initial_value; 
  }
  
  grid->min_value = initial_value;
  grid->max_value = initial_value; 
    
  return grid;
}
/* ******************************************************************** */



/* ******************************************************************** */
// Free a grid and its associated malloced data;
void grid_free(Grid* grid) {
  assert(grid);
  if (grid->data) free(grid->data);
  free(grid);
}
/* ******************************************************************** */



//interpolate nodata values from nearest neighbors 
void grid_interpolate_nodata(Grid* grid) { 

  int nnodata=0; 
  for (int r=0; r<grid->nrows; r++) {
    for (int c=0; c<grid->ncols; c++) {
      if (!grid_is_nodata(grid, r, c)) continue;
      
      //(r,c) is nodata 
      //get average of neighbors that have data 
      float sum = 0, max = grid->nodata_value, min = grid->max_value;
      int nneighb = 0; 
      for (int i=-1; i<=1; i++) {
	for (int j=-1; j<=1; j++) {
	  if (i==0 && j==0) continue; 
	  if (inside_grid(grid, r+i, c+j) && !grid_is_nodata(grid, r+i, c+j)) {
	    float neighb = grid_get(grid, r+i, c+j);
	    sum += neighb; 
	    if (max < neighb) max = neighb;
	    if (min > neighb) min = neighb;
	    nneighb++;
	  }
	} //for j 
      }//for i
      
      //set this point to neighb average
      if (nneighb > 0) {
	grid_set(grid, r, c, sum/nneighb);
	//grid_set(grid, r, c, max);
	//grid_set(grid, r, c, min);
	nnodata++;
      }
     
      
    } //for c
  } //for r
  grid_reset_stats(grid); 
  printf("grid_interpolate: interpolated %d nodata cells\n", nnodata); 
}
