
#include "lidar.hpp"
#include "grid.h"
#include "pixel_buffer.h"
#include "map.h"


#include <assert.h>

#include <vector>
using namespace std; 


const float  DEFAULT_NODATA_VALUE = -9999;


//used in deciding teh resolution of the grid. If the average number
//of points in a cell is larger than this, the cell is made smaller
const float DSM_AVG_COUNTS_THRESHOLD = 2; 




//adds point p  to  lp 
void lidar_add_point(lidar_point_cloud* lp, lidar_point p) {

  //make sure its a valid pointer 
  assert(lp);

  //add the point
  lp->data.push_back(p);
  
  //update bounding box
  if (lidar_size(lp) == 1) {
    lp->minx = lp->maxx = p.x; 
    lp->miny = lp->maxy = p.y; 
    lp->minz =lp-> maxz = p.z; 
  } else {
    if (lp->minx > p.x) lp->minx=p.x; 
    if (lp->maxx < p.x) lp->maxx = p.x; 
    if (lp->miny > p.y) lp->miny=p.y; 
    if (lp->maxy < p.y) lp->maxy = p.y; 
    if (lp->minz > p.z) lp->minz=p.z; 
    if (lp->maxz < p.z) lp->maxz = p.z; 
  }
} 



/*
  reads lidar points from file and  populates points
  NOTE: file.txt must be obtained from file.las with 'pdal translate'
  reads the points from file in global array points
*/
lidar_point_cloud* read_lidar_from_file(char* fname) {

  lidar_point_cloud* points = (lidar_point_cloud*)malloc(sizeof(lidar_point_cloud));
    
  FILE* file = fopen(fname, "r"); 
  if (!file) {
    printf("read_lidar:from_file: cannot open file %s\n",  fname);
    exit(1); 
  }

  //read first line
  //"X","Y","Z","Intensity","ReturnNumber","NumberOfReturns",
  //"ScanDirectionFlag","EdgeOfFlightLine","Classification","ScanAngleRank",
  //"UserData","PointSourceId","GpsTime"
  char line[1000];
  if (fgets(line,sizeof line,file)!= NULL) /* read a line from a file */ {
    printf("%s",line); //print the file contents on stdout.
  } else {
    printf("read_lidar_fom_file: cannot read from file\n");
    exit(1); 
  } 
  
  lidar_point p;
  float retnb, nbret, classif; //dir, edge, angle, user, pid, gps;
  while (1) {
    /*
      if (fscanf(file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
	       &p.x, &p.y, &p.z,
	       &p.intensity, &retnb, &nbret,
	       &dir, &edge, &classif,
	       &angle, &user, &pid,
	       &gps
	       ) < 13)
    */
    if (fscanf(file, "%f,%f,%f,%f,%f,%f", &p.x, &p.y, &p.z,&retnb, &nbret, &classif) < 6)    { 
      //either an error or all done 
      break; 
    }
    //we got another point
    p.return_number = (int) retnb;
    p.nb_of_returns = (int)nbret;
    p.code=(int)classif;
    p.mycode=0; //everything unclassified
    
    //printf("point: x=%f, y=%f, z=%f, intensity=%f return_nb=%d nb_returns=%d classif=%d\n",
    //  p.x, p.y, p.z, p.intensity, p.return_number, p.nb_of_returns, p.code); 
	
    //insert the point 
    lidar_add_point(points, p); 
    
  } //while 

  //done reading points 
  fclose(file); 
  
  //print info about the points that were read 
  printf("read total %d points\n", (int)lidar_size(points)); 
  printf("\tbounding box:  x=[%.2f, %.2f], y=[%.2f,%.2f], z=[%.2f,%.2f]\n",
	 points->minx, points->maxx, points->miny, points->maxy, points->minz, points->maxz); 

  
  return points; 
}


//print all the points 
void lidar_print(const lidar_point_cloud* lp) {
  
   for (int i=0; i< lidar_size(lp); i++) {
     
     //curent lidar point 
     lidar_point curp = lidar_get(lp, i);
     printf("i=%d: (x=%.2f, y=%.2f, z=%.2f, %d, %d)\n",
	    i, curp.x, curp.y, curp.z, curp.return_number, curp.nb_of_returns); 
   } 
}



//lp: the lidar points 
//return: a grid of last returns 
Grid* lidar_to_dsm(const lidar_point_cloud* lp) {


  return NULL;
} 




//creates and returns a new grid, where each cell is eroded
Grid* grid_erode(const Grid* grid) {

  printf("erode\n"); fflush(stdout); 
 
  return NULL; 
}



//creates and returns a new grid, where each cell is dilated
Grid* grid_dilate(const Grid* grid) {

  printf("dilate:\n");
  
  return NULL; 
}
 


//lp: the lidar points
//dsm: the grid of last returns
//return: a grid corresponding to the bare ground
Grid* lidar_to_bareground(const lidar_point_cloud* lp, const Grid* dsm_grid) {

  assert(lp && dsm_grid);

  Grid* dtm_grid; 
    
  //FOR THE BITMAPS 
  const PixelBuffer pb = init_pixel_buffer(dsm_grid->ncols, dsm_grid->nrows);
  Grid* tmpgrid;
   
  //erode 1
  dtm_grid= grid_erode(dsm_grid);
  //grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.erode1.grayscale.bmp");
  
  //erode 2
  tmpgrid = grid_erode(dtm_grid);
  grid_free(dtm_grid);
  dtm_grid = tmpgrid;
  //grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.erode2.grayscale.bmp");
   
  
  //erode 3
  tmpgrid = grid_erode(dtm_grid);
  grid_free(dtm_grid);
  dtm_grid = tmpgrid;
  //grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.erode3.grayscale.bmp");

  //dilate 
  tmpgrid = grid_dilate(dtm_grid);
  grid_free(dtm_grid);
  dtm_grid = tmpgrid;
  // grid_grayscale_to_pixelbuffer(dtm_grid, pb);
  save_pixel_buffer_to_file(&pb, "map.dilate1.grayscale.bmp");


  //do more stuff 



  
  //cleanup
  deinit_pixel_buffer(&pb);

  return dtm_grid;
} 







/* ************************************************************ 
lidar classification codes

0 bever classified 
1 unassigned 
2 ground 
3 low vegetation 
4 medium vegetation 
5 high vegetation
6 building 
7 low point (noise)
8 model key-point (mass point)
9 water
10 railroad
11 road
12 overlap 
13 wire-guard (shield)
14 wire-conductor (phase)
15 transmission tower
17 bridge
18 hight point (noise)
19-255 reserved for asprs definition
*/

/* for every point p, it sets p.mycode to one of the codes above */
void classify(lidar_point_cloud* points) {

  //float minheight, max_height; 
  for (int i=0; i< lidar_size(points); i++) {

    lidar_point p = lidar_get(points, i); 

    //vegetation: points with > 1 return, and not last return  
    if ((p.nb_of_returns>1) && (p.return_number != p.nb_of_returns)) {
      p.mycode=4;
      lidar_set(points, i, p); 
    }

    //under the vegetation is ground
    if((p.nb_of_returns > 1) && (p.return_number == p.nb_of_returns)) {
      p.mycode=2;
      lidar_set(points, i, p);
    }
  } 
}//classify 




