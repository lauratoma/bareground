#ifndef __CLASSIFY_HPP
#define  __CLASSIFY_HPP


#include "grid.h"

#include <vector>
using namespace std; 





typedef struct _lidar_point {
  float x,y,z; 
  float intensity;

  int return_number; //the number of this return
  int nb_of_returns; //how many returns this pulse has

  int code;  //classification code read from file 
  
  int mycode; //classification code that we'll assign to this point

  //there's more info available for a point (such as
  //scanDirectionF;ag, EdgeOfFlightLine, etc) but we don't store it

} lidar_point;


typedef struct _lidar_data {

  vector<lidar_point> data; 

  //bounding box
  float  minx, maxx, miny, maxy, minz, maxz; 
  
} lidar_point_cloud; 



//returns size (= nb points) 
static inline int lidar_size(const lidar_point_cloud* lp) {
  assert(lp); 
  return lp->data.size();
}

//returns the i-th point
static inline lidar_point lidar_get(const lidar_point_cloud* lp, int i) {

  assert(lp); 
  assert(i>= 0 && i < lidar_size(lp)); 
  return  lp->data.at(i);  // (*(lp->data))[i]; 
}

//sets the i-th point
static inline void lidar_set(lidar_point_cloud* lp, int i, lidar_point p) {

  assert(lp); 
  assert(i>= 0 && i < lidar_size(lp)); 
  lp->data.at(i) = p;  
}

//print all points 
void lidar_print(const lidar_point_cloud* lp);


/*
  reads lidar points from file and  populates points
  NOTE: file.txt must be obtained from file.las with 'pdal translate'
  reads the points from file in global array points
*/
lidar_point_cloud* read_lidar_from_file(char* fname); 


//adds point p  to  points
void lidar_add_point(lidar_point_cloud* lp, lidar_point p); 

//lp: the lidar points 
//return: a grid of last returns 
Grid* lidar_to_dsm(const lidar_point_cloud* lp); 


//lp: the lidar points
//dsm: the grid of last returns
//return: a grid corresponding to the bare ground
Grid* lidar_to_bareground(const lidar_point_cloud* lp,  const Grid* dsm);

//creates and returns a new grid, where each cell is dilated
Grid* grid_dilate(const Grid* grid);

//creates and returns a new grid, where each cell is eroded
Grid* grid_erode(const Grid* grid);


/* ************************************************************ 
lidar classification codes

0 never classified 
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
void classify(lidar_point_cloud & points);

#endif 
