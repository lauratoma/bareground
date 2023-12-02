
#include "lidar.hpp"
#include "grid.h"
#include "map.h"
#include "pixel_buffer.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


#include <vector>
using namespace std; 




/************************************************************/
int main(int argc, char** argv) {

  
  //read number of points from user
  if (argc!=2) {
    printf("usage: %s lidarfile.txt\n", argv[0]);
    exit(1); 
  }

  

  //the lidar points
  lidar_point_cloud* lpoints;
  
  //the DSM grid
  Grid* dsm_grid; 
  
  //the DTM grid
  Grid* bare_grid; 
  
  //read the lidar points from file 
  lpoints = read_lidar_from_file(argv[1]); 
  //lidar_print(lpoints);

  
  //CREATE A GRID FROM THE LAST RETURNS 
  dsm_grid = lidar_to_dsm(lpoints); 
  //grid_print_stats(dsm_grid, "dsm");  
  //grid_write_to_file("dsm.asc", dsm_grid);


  // commented out because dsm_grid is NULL 
  /*
  // BITMAPS 
  const PixelBuffer pb = init_pixel_buffer(dsm_grid->ncols, dsm_grid->nrows);
  
  Grid * hillshade_grid = grid_init_from(dsm_grid);
  generate_hillshade(dsm_grid, hillshade_grid);
  grid_to_pixelbuffer(hillshade_grid, pb);
  grid_free(hillshade_grid);
  printf("writing map.dsm.hillshade.bmp\n");  fflush(stdout); 
  save_pixel_buffer_to_file(&pb, "map.dsm.hillshade.bmp");
  
  grid_grayscale_to_pixelbuffer(dsm_grid, pb);
  printf("writing map.dsm.grayscale.bmp\n"); fflush(stdout); 
  save_pixel_buffer_to_file(&pb, "map.dsm.grayscale.bmp");
  */
  

  // CREATE THE BARE GROUND GRID BY ELIMINATING NOISE, BUILDINGS AND VEGETATION 
  bare_grid = lidar_to_bareground(lpoints, dsm_grid);


  //clean up
  free(lpoints);
  //these fail now because they are NULL
  // deinit_pixel_buffer(&pb);
  // grid_free(dsm_grid);
  //grid_free(bare_grid);
  

}
