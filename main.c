
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
  Grid* dtm_grid; 
  
  //this populates the global that holds the points
  lpoints = read_lidar_from_file(argv[1]); 
  //lidar_print(lpoints);
  
  //create the DSM 
  dsm_grid = lidar_to_dsm(lpoints); 
  grid_print_stats(dsm_grid, "dsm");  
  //grid_write_to_file("dsm.asc", dsm_grid);


  //bitmaps
  // FOR THE BITMAPS 
  const PixelBuffer pb = init_pixel_buffer(dsm_grid->ncols, dsm_grid->nrows);
  
  Grid * hillshade_grid = grid_init_from(dsm_grid);
  generate_hillshade(dsm_grid, hillshade_grid);
  grid_print_stats(hillshade_grid, "hillshade");
  grid_to_pixelbuffer(hillshade_grid, pb);
  grid_free(hillshade_grid);
  printf("writing map.dsm.hillshade.bmp\n");  fflush(stdout); 
  save_pixel_buffer_to_file(&pb, "map.dsm.hillshade.bmp");
  
  grid_grayscale_to_pixelbuffer(dsm_grid, pb);
  printf("writing map.dsm.grayscale.bmp\n"); fflush(stdout); 
  save_pixel_buffer_to_file(&pb, "map.dsm.grayscale.bmp");

  
  //want the bare ground DTM grid
  //erode 1
  dtm_grid= grid_erode(dsm_grid);
  grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.erode1.grayscale.bmp");

  //erode 2
  dtm_grid = grid_erode(dtm_grid);
  grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.erode2.grayscale.bmp");
   
  
  //erode 3
  dtm_grid = grid_erode(dtm_grid);
  grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.erode3.grayscale.bmp");

  //dilate 
  dtm_grid = grid_dilate(dtm_grid);
  grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.dilate1.grayscale.bmp");
  
  
  //free the pixel buffer
  deinit_pixel_buffer(&pb);
  grid_free(dsm);
  grid_free(dtm);
  
}
