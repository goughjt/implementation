#ifndef mySuperpixel_H
#define mySuperpixel_H

#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
#include <string>
#include <mrpt/math/ransac.h>
#include <mrpt/math/ops_vectors.h>
#include <mrpt/math/CMatrix.h>
#include <mrpt/math/CMatrixD.h>
#include <mrpt/math/ops_matrices.h>
#include <mrpt/random.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::random;
using namespace mrpt::math;
using namespace std;

class mySuperpixel {
	public:
	
	void see();
	
	mySuperpixel(int ind, int i, int sprpxl, vector<double> centerPnt, string imPath_original);
	~mySuperpixel();
	
	void hull();
        
	int imageIndex;
	string imagePath_original;
	int sprpxlIndex; //index of sprpxls in image allImages[imageIndex]
	int index; //index in allSprpxls
	
	int acceptedPlaneIndex;
	
	
	vector<int> neigh;
	vector<int> neighC;
	vector<int> neighB;
	vector<int> neighBC;
	vector<int> gradients;
	
	double coverage;
	
	vector<double> center; //5D vector, 3 colour, 2 pixel coord
	int hasPlane; //1 for true, 0 for false, 2 for without hard
	mrpt::math::TPlane initPlane;
	mrpt::math::TPlane perturbedPlane;

	vector<mrpt::math::TPoint3D> allPnts; //3D vector of point position
	vector<mrpt::math::TPoint3D> inliers;
	vector<mrpt::math::TPoint3D> perturbedInliers;
	
	vector<mrpt::math::TPoint3D> initialVertices;
	vector<mrpt::math::TPoint3D> finalVertices;
	
	vector<mrpt::math::TPoint3D> vertexNormals;
	
	vector<double> pntVwx; //point x coord
	vector<double> pntVwy; //point y coord
	vector<cv::Point> contours;
	vector<cv::Point> hullPnts;
	vector<int> colour;

	void clear_data();

};

#endif


