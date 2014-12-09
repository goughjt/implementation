#ifndef MYRANSAC_H
#define MYRANSAC_H

#include <mrpt/math/ransac.h>
#include<mrpt/math/lightweight_geom_data.h>
#include "mySuperpixel.h"


//void TestRANSAC(mrpt::math::CMatrixDouble dataIn, mrpt::math::TPlane &initPlane);
void myRANSAC(mySuperpixel* superpixel, double distanceThreshold);
void ransacInliers(mySuperpixel* superpixel, double distanceThreshold);
bool fitToAllInliers(mySuperpixel* superpixel);
bool fitToPerturbedInliers(mySuperpixel* superpixel);
void ransacPerturbedInliers(mySuperpixel* superpixel, double distanceThreshold);

#endif
