#ifndef MYIMAGE_H
#define MYIMAGE_H

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
#include <string>
#include "mySuperpixel.h"
#include "myBundler.h"
#include "DataInterface.h"
#include "parameters.h"

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::random;
using namespace mrpt::math;
using namespace std;

#ifndef vec2dd
#define vec2dd vector<vector<double> >
#endif
#ifndef vec2di
#define vec2di vector<vector<int> >
#endif
#ifndef vec2db
#define vec2db vector<vector<bool> >
#endif
/* The number of iterations run by the clustering algorithm. */
#define NR_ITERATIONS 10

class myImage {

	private:

	void makeStrings(bool createSlicFolder, string fullpath); //creates strings used for saving and displaying
	//void assignSuperpixels(myBundler &bundler, vector<mySuperpixel> &allSuperpixels); //assign 3d points
	void assignSuperpixels(vector<Point2D> &measurements, vector<int> &camidx, vector<Point3D> &point_data, vector<int> &ptidx, vector<mySuperpixel> &allSuperpixels);
	int step, ns; //Step size per cluster, distance (nc) and colour (ns) parameters
  	double compute_dist(int ci, CvPoint pixel, CvScalar colour); //get distance between a center and a pixel
  	CvPoint find_local_minimum(CvPoint center, IplImage* labImage); //get pixel with lowest gradient in a 3x3 nbhd
  	void generate_superpixels(vector<mySuperpixel> &allSuperpixels, int num_sprpxl);
  	void create_connectivity();
  	void display_contours(CvScalar colour, vector<mySuperpixel> &allSuperpixels);
  	void colour_with_cluster_means(vector<mySuperpixel> &allSuperpixels);
  	void createSuperpixelContours(vector<mySuperpixel> &allSuperpixels);

	public:

	myImage(vector<Point2D> &measurements, vector<int> &camidx, vector<Point3D> &point_data, vector<int> &ptidx, int imIndex, string fullpath, string fullpath_LowRes, vector<mySuperpixel> &allSuperpixels, int num_sprpxl);
	~myImage();
	string filename;
  	string imagePath_original;
	string imagePath_labspace;
	string imagePath_contours;
	string imagePath_clustermeans;
	string imagePath_lap;
	int imageIndex;
	int w;
	int h;
	int w_LowRes;
	int h_LowRes;
	vector<int> superpixelIndices; //ives index for allSuperpixels
	vec2dd centers; //The LAB and xy values of the centers.
  	vector<int> center_counts; //The number of occurences of each center
  	void display_center_grid(CvScalar colour, IplImage* labImage);
	void seeSuperpixelByID(int ID);

};

#endif


