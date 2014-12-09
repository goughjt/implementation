
#ifndef labellingCost_H
#define labellingCost_H


#include "mySuperpixel.h"
#include "myImage.h"
#include "parameters.h"
#include <mrpt/math/ransac.h>
#include <mrpt/math/ops_vectors.h>
#include <mrpt/math/CMatrix.h>
#include <mrpt/math/CMatrixD.h>
#include <mrpt/math/ops_matrices.h>
#include <mrpt/random.h>

using namespace std;


double PI =3.14;
double alpha = 0.1;

double betweenCost(mySuperpixel &s1, mySuperpixel &s2, int nbrIndex){
	double sigmaNum = 2;
	double sigmaColour = 0.05; //shaping parameter
	double cij = exp(-sqrt(pow(s1.colour[0]-s2.colour[0],2)+pow(s1.colour[1]-s2.colour[1],2)+pow(s1.colour[2]-s2.colour[2],2))/(256*sigmaColour));
	double gamma = alpha*0.2;
	double wij = 1-exp(-(s1.neighBC[nbrIndex])/sigmaNum);
	return gamma*wij*cij;
}

double withinCost(mySuperpixel &s1, mySuperpixel &s2, int nbrIndex){
	double beta = alpha;
	double sigmaColour = 0.05; //shaping parameter
	double sigmaBoundary = 0.1;
	double sigmaGradient = 0.05;
	//consider this neighbour
	double cij = exp(-sqrt(pow(s1.colour[0]-s2.colour[0],2)+pow(s1.colour[1]-s2.colour[1],2)+pow(s1.colour[2]-s2.colour[2],2))/(256*sigmaColour));
	double gij = exp(abs(s1.gradients[nbrIndex]/s1.neighC[nbrIndex])/(256*sigmaGradient));
	double wij = 1-exp(-(s1.neighC[nbrIndex]/min(s1.contours.size(),s2.contours.size()))/sigmaBoundary);
	if((alpha*cij+beta*gij)*wij<0){
		return 0;
	}else{
		return (alpha*cij+beta*gij)*wij;
	}
}

double unaryCost(mySuperpixel &s, CameraT camera_data, TPlane plane){
		int sigmaFit =3;
		int sigmaRays =5;
		double delta=PI/36.0;
		float Dunary[3];
		float campnt[3];
		camera_data.GetCameraCenter(campnt);
		TPoint3D cameraPnt(campnt[0], campnt[1], campnt[2]);
		double totalDist=0;
		int totalRays=0;
		for(int j = 0; j<s.allPnts.size(); j++){
			double d = plane.distance(s.allPnts[j]);		
			totalDist +=exp(-pow(d,2)/(2*pow(*distThreshold,2)));
			if(plane.distance(cameraPnt)<cameraPnt.distanceTo(s.allPnts[j]) && plane.distance(s.allPnts[j])>*distThreshold){
				totalRays++;
			}
		}
		Dunary[0]+=exp(totalDist/sigmaFit);
		Dunary[1]+=1-exp(totalRays/sigmaRays);
		double largestAngle=0	;
		for(int j = 0; j<s.initialVertices.size(); j++){
			double angley=std::acos((s.vertexNormals[j].x*plane.coefs[0]+s.vertexNormals[j].y*plane.coefs[1]+s.vertexNormals[j].z*plane.coefs[2])/pow((s.vertexNormals[j].x*s.vertexNormals[j].x+s.vertexNormals[j].y*s.vertexNormals[j].y + s.vertexNormals[j].z*s.vertexNormals[j].z)*(plane.coefs[0]*plane.coefs[0]+plane.coefs[1]*plane.coefs[1]+plane.coefs[2]*plane.coefs[2]),0.5));
			if(angley>largestAngle){
				largestAngle-=largestAngle;
				largestAngle+=angley;
			}
		}
		if(largestAngle>((PI/2)-delta)){
			Dunary[2]+= 0.5 + 0.5*cos(PI*(largestAngle-(PI/2))/delta);
		}
		return Dunary[0] + Dunary[1] + Dunary[2];
}

double lCost(vector<mySuperpixel> &allSuperpixels, vector<CameraT> &camera_data){

	float Dunary=0;
	float Vw;
	float Vb;
	
	for(int i =0; i<allSuperpixels.size(); i++){
		mySuperpixel s = allSuperpixels[i];
		TPlane plane = s.initPlane;
		Dunary+=unaryCost(s,camera_data[s.imageIndex],plane);
	}
	
	for(int i =0; i<allSuperpixels.size(); i++){
		mySuperpixel s1 = allSuperpixels[i];
		for(int j=0; j<s1.neigh.size(); j++){
			if(s1.neigh[j]>i ){//otherwise the sprpxl has already been calculated
				mySuperpixel s2 = allSuperpixels[s1.neigh[j]];
				Vw+=withinCost(s1, s2, j);
			}
		}
		for(int j=0; j<s1.neighB.size(); j++){
			//if(s.neighB[j]>i ){ //not needed due to way neighb initialized?
			mySuperpixel s2 = allSuperpixels[s1.neighB[j]];
			Vb+=betweenCost(s1, s2, j);
		}
	}
	
	return Dunary +Vw +Vb;

}


#endif

