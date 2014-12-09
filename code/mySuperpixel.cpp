
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <float.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "mySuperpixel.h"
#include <mrpt/math/ransac.h>
#include <mrpt/math/ops_vectors.h>
#include <mrpt/math/CMatrix.h>
#include <mrpt/math/CMatrixD.h>
#include <mrpt/math/ops_matrices.h>
#include <mrpt/random.h>
#include <iostream>
#include <fstream>

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::random;
using namespace mrpt::math;
using namespace std;

void mySuperpixel::see(){
	IplImage* image = cvLoadImage(imagePath_original.c_str(), 1);
	ostringstream convert;
	convert << "See sprpxl " << sprpxlIndex << " in image " << imageIndex;
	string showTitle = convert.str();
	for (int j = 0; j < contours.size(); j++) {
		cvSet2D(image, contours[j].y, contours[j].x, CV_RGB(255,0,0));
	}
	cvShowImage(showTitle.c_str(), image);
	cvWaitKey(0);
	cvReleaseImage(&image);
}

mySuperpixel::mySuperpixel(int ind, int i, int sprpxl, vector<double> centerPnt, string imPath_original) {
	index = ind;
	imageIndex=i;
	sprpxlIndex=sprpxl;
	center=centerPnt;
	hasPlane=0;
	imagePath_original=imPath_original;
	acceptedPlaneIndex=0;
}

mySuperpixel::~mySuperpixel() {
    clear_data();
}

void mySuperpixel::clear_data() {
	//initPlane.forget()
}


