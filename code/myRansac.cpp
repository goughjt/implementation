#include <mrpt/math/ops_vectors.h>
#include <mrpt/math/CMatrix.h>
#include <mrpt/math/CMatrixD.h>
#include <mrpt/math/ops_matrices.h>
#include <mrpt/random.h>
#include <mrpt/math/ransac.h>
#include <mrpt/gui/CDisplayWindow3D.h>
#include <mrpt/random.h>
#include <mrpt/utils/CTicTac.h>
#include <mrpt/poses/CPose3D.h>
#include <mrpt/opengl/CGridPlaneXY.h>
#include <mrpt/opengl/CPointCloud.h>
#include <mrpt/opengl/CTexturedPlane.h>
#include "myRansac.h"
#include "RedSVD.h"

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::gui;
using namespace mrpt::math;
using namespace mrpt::random;
using namespace mrpt::poses;
using namespace std;

size_t nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

void  ransac3Dplane_fit(
	//const CMatrixDouble  &allData,
	vector<TPoint3D> &allData,
	const vector_size_t  &useIndices,
	vector< CMatrixDouble > &fitModels )
{
	//cout<<"welcome to fit func"<<endl;
	ASSERT_(useIndices.size()==3);
	//cout<<"assertion made"<<endl;
	try
	{
	//cout<<"going to make plane"<<endl;
		TPlane  plane( allData[0],allData[1],allData[2] );
	//cout<<"plane made"<<endl;
		fitModels.resize(1);
	//cout<<"fitModels resized"<<endl;
		CMatrixDouble &M = fitModels[0];
	//cout<<"M pointer"<<endl;

		M.setSize(1,4);
	//cout<<"M size"<<endl;
		for (size_t i=0;i<4;i++)
			M(0,i)=plane.coefs[i];
	//cout<<"M to plane coef "<<endl;
	}
	catch(exception &)
	{
		//cout<<p1.x<<p1.y<<p1.z<<endl;
		//cout<<p2.x<<p2.y<<p2.z<<endl;
		//cout<<p3.x<<p3.y<<p3.z<<endl;
		fitModels.clear();
		return;
	}



}

void ransac3Dplane_distance(
	//const CMatrixDouble &allData,
	vector<TPoint3D> &allData,
	const vector< CMatrixDouble > & testModels,
	const double distanceThreshold,
	unsigned int & out_bestModelIndex,
	vector_size_t & out_inlierIndices,
	double & totalDist)
{
	ASSERT_( testModels.size()==1 )
	out_bestModelIndex = 0;
	const CMatrixDouble &M = testModels[0];

	ASSERT_( size(M,1)==1 && size(M,2)==4 )

	TPlane  plane;
	plane.coefs[0] = M(0,0);
	plane.coefs[1] = M(0,1);
	plane.coefs[2] = M(0,2);
	plane.coefs[3] = M(0,3);

	const size_t N = allData.size();
	out_inlierIndices.clear();
	out_inlierIndices.reserve(100);
	for (size_t i=0;i<N;i++)
	{
		const double d = plane.distance(allData[i]);		
		totalDist = totalDist + exp(-pow(d,2)/(2*pow(distanceThreshold,2)));
		if (d<distanceThreshold)
			out_inlierIndices.push_back(i);
	}
}



/** Return "true" if the selected points are a degenerate (invalid) case.
  */
bool ransac3Dplane_degenerate(
	//const CMatrixDouble &allData,
	vector<TPoint3D> &allData,
	const mrpt::vector_size_t &useIndices )
{
	return false;
}

void ransacPerturbedInliers(mySuperpixel* s, double distanceThreshold){
	CMatrixDouble best_model;
	vector_size_t best_inliers;
	double leastTotDist=10000;
	CTicTac	tictac;
	size_t TIMES=100;
	for (size_t iters=0;iters<TIMES;iters++){
		math::RANSAC::execute(
			s->perturbedInliers,
			ransac3Dplane_fit,
			ransac3Dplane_distance,
			ransac3Dplane_degenerate,
			distanceThreshold,
			3,  // Minimum set of points
			best_inliers,
			best_model,
			leastTotDist,
			false // iters==0   // Verbose
			);		
	}
	if(size(best_model,1)==1 && size(best_model,2)==4 && best_inliers.size()>2){
		//double normi = pow(best_model(0,0)*best_model(0,0)+best_model(0,1)*best_model(0,1)+best_model(0,2)*best_model(0,2),0.5);
		s->perturbedPlane.coefs[0]=best_model(0,0);
		s->perturbedPlane.coefs[1]=best_model(0,1);
		s->perturbedPlane.coefs[2]=best_model(0,2);
		s->perturbedPlane.coefs[3]=best_model(0,3);
		s->hasPlane=1;
	}else{
		s->hasPlane=0;
	}
}



void ransacInliers(mySuperpixel* s, double distanceThreshold){
	CMatrixDouble best_model;
	vector_size_t best_inliers;
	double leastTotDist=10000;
	CTicTac	tictac;
	const size_t TIMES=100;
	for (size_t iters=0;iters<TIMES;iters++){
		math::RANSAC::execute(
			s->inliers,
			ransac3Dplane_fit,
			ransac3Dplane_distance,
			ransac3Dplane_degenerate,
			distanceThreshold,
			3,  // Minimum set of points
			best_inliers,
			best_model,
			leastTotDist,
			false // iters==0   // Verbose
			);		
	}
	if(size(best_model,1)==1 && size(best_model,2)==4 && best_inliers.size()>2){
		s->initPlane.coefs[0]=best_model(0,0);
		s->initPlane.coefs[1]=best_model(0,1);
		s->initPlane.coefs[2]=best_model(0,2);
		s->initPlane.coefs[3]=best_model(0,3);
		s->hasPlane=1;
	}else{
		s->hasPlane=0;
	}
}

void myRANSAC(mySuperpixel* s, double distanceThreshold)
{
	if(s->allPnts.size()==3){
		try{
			s->initPlane = TPlane(s->allPnts[0], s->allPnts[1], s->allPnts[2]);
			s->inliers.push_back(s->allPnts[0]);
			s->inliers.push_back(s->allPnts[1]);
			s->inliers.push_back(s->allPnts[2]);
			s->hasPlane=1;
		}
		catch(exception &)
		{
			s->hasPlane=0;
		}
	}else{
		CMatrixDouble best_model;
		vector_size_t best_inliers;
		double leastTotDist=10000;
		CTicTac	tictac;
		const size_t TIMES=100;
		for (size_t iters=0;iters<TIMES;iters++){
			math::RANSAC::execute(
				s->allPnts,
				ransac3Dplane_fit,
				ransac3Dplane_distance,
				ransac3Dplane_degenerate,
				//DIST_THRESHOLD,
				distanceThreshold,
				3,  // Minimum set of points
				best_inliers,
				best_model,
				leastTotDist,
				false // iters==0   // Verbose
				);
			
		}
			//cout << "Computation time: " << tictac.Tac()*1000.0/TIMES << " ms" << endl;
		if(size(best_model,1)==1 && size(best_model,2)==4 && best_inliers.size()>2){
			//cout << "RANSAC finished: Best model: " << best_model << endl;
			//cout << "Best inliers: " << best_inliers << endl;
			s->initPlane.coefs[0]=best_model(0,0);
			s->initPlane.coefs[1]=best_model(0,1);
			s->initPlane.coefs[2]=best_model(0,2);
			s->initPlane.coefs[3]=best_model(0,3);
			s->hasPlane=1;
			s->inliers.clear();
			for(int i = 0; i < best_inliers.size(); i++){
				//TPoint3D p(s->allPnts[best_inliers[i]]),s->allPnts(1, best_inliers[i]), s->allPnts(2, best_inliers[i]));
				s->inliers.push_back(s->allPnts[best_inliers[i]]);
			}
		}else{
			s->hasPlane=0;
			//cout<<"NO PLANE!!!"<<endl;
		}
	}

}


bool fitToAllInliers(mySuperpixel* s){

	//cout<<"start"<<endl;
	int sizey = (s->inliers).size();
	//int sizey = 100;
	if(sizey<3){cout<<"YUP"<<endl; return false;}
	mrpt::math::CMatrixDouble mx(sizey, 4);
	//cout<<"mx rows: "<<size(mx,1)<<endl;
	//cout<<"mx cols: "<<size(mx,2)<<endl;
	for(int i =0; i<sizey; i++){
		mx(i,0) = s->inliers[i].x;
		mx(i, 1) = s->inliers[i].y;
		mx(i, 2) = s->inliers[i].z;
/*		mx(i,0) = i;
		mx(i, 1) = pow(-1,i)*i*i;
		mx(i, 2) = 0;*/
		mx(i, 3) = 1;

	}
	
	/*
	cout<<"BEFORE"<<endl;
	cout<<mx<<endl;
	cout<<mx(1,2)<<endl;
	*/
	
	RedSVD::RedSVD<CMatrixDouble> redd(mx);

	mrpt::math::CMatrixDouble mxu = redd.matrixU();
	mrpt::math::CMatrixDouble mxv = redd.matrixV();
	mrpt::math::CMatrixDouble sVals = redd.singularValues();
	int sizey2 = size(mxv,2)-1;
	
	/*cout<<"mxu"<<endl;
	cout<<mxu<<endl;
	cout<<"mxv"<<endl;
	cout<<mxv<<endl;
	cout<<mxv(0,sizey2)<<endl;
	/*cout<<"sizey2: "<<sizey2<<endl;
	cout<<"sVals"<<endl;
	cout<<sVals<<endl;*/

	/*for(int k=0; k< s->allPnts.size(); k++){
		double heya=0;
		heya+=s->allPnts[k].x*s->initPlane.coefs[0];
		heya+=s->allPnts[k].y*s->initPlane.coefs[1];
		heya+=s->allPnts[k].z*s->initPlane.coefs[2];
		heya+=s->initPlane.coefs[3];
		cout<<"eqn 1: "<<heya<<endl;
	}

	for(int k=0; k< s->allPnts.size(); k++){
		double heya=0;
		heya+=s->allPnts[k].x*s->initPlane.coefs[0];
		heya+=s->allPnts[k].y*s->initPlane.coefs[1];
		heya+=s->allPnts[k].z*s->initPlane.coefs[2];
		heya+=s->initPlane.coefs[3];
		cout<<"eqn 2: "<<heya<<endl;
	}*/
	
	s->initPlane.coefs[0]=mxv(0, sizey2);
	s->initPlane.coefs[1]=mxv(1, sizey2);
	s->initPlane.coefs[2]=mxv(2, sizey2);
	s->initPlane.coefs[3]=mxv(3, sizey2);
	
	if(s->initPlane.coefs[0]==0 &&
	s->initPlane.coefs[1]==0 &&
	s->initPlane.coefs[2]==0 &&
	s->initPlane.coefs[3]==0){
		//cout<<"ZERO FELLAHS"<<endl;
		s->hasPlane=0;
		return false;
	}else{
		return true;
	}
	
	//cout<<"end"<<endl;

}


bool fitToPerturbedInliers(mySuperpixel* s){

	//cout<<"start"<<endl;
	int sizey = (s->perturbedInliers).size();
	//int sizey = 100;
	if(sizey<3){cout<<"YUP"<<endl; return false;}
	mrpt::math::CMatrixDouble mx(sizey, 4);
	//cout<<"mx rows: "<<size(mx,1)<<endl;
	//cout<<"mx cols: "<<size(mx,2)<<endl;
	for(int i =0; i<sizey; i++){
		mx(i,0) = s->perturbedInliers[i].x;
		mx(i, 1) = s->perturbedInliers[i].y;
		mx(i, 2) = s->perturbedInliers[i].z;
/*		mx(i,0) = i;
		mx(i, 1) = pow(-1,i)*i*i;
		mx(i, 2) = 0;*/
		mx(i, 3) = 1;

	}
	
	/*
	cout<<"BEFORE"<<endl;
	cout<<mx<<endl;
	cout<<mx(1,2)<<endl;
	*/
	
	RedSVD::RedSVD<CMatrixDouble> redd(mx);

	mrpt::math::CMatrixDouble mxu = redd.matrixU();
	mrpt::math::CMatrixDouble mxv = redd.matrixV();
	mrpt::math::CMatrixDouble sVals = redd.singularValues();
	int sizey2 = size(mxv,2)-1;
	
	/*cout<<"mxu"<<endl;
	cout<<mxu<<endl;
	cout<<"mxv"<<endl;
	cout<<mxv<<endl;
	cout<<mxv(0,sizey2)<<endl;
	/*cout<<"sizey2: "<<sizey2<<endl;
	cout<<"sVals"<<endl;
	cout<<sVals<<endl;*/

	/*for(int k=0; k< s->allPnts.size(); k++){
		double heya=0;
		heya+=s->allPnts[k].x*s->initPlane.coefs[0];
		heya+=s->allPnts[k].y*s->initPlane.coefs[1];
		heya+=s->allPnts[k].z*s->initPlane.coefs[2];
		heya+=s->initPlane.coefs[3];
		cout<<"eqn 1: "<<heya<<endl;
	}

	for(int k=0; k< s->allPnts.size(); k++){
		double heya=0;
		heya+=s->allPnts[k].x*s->initPlane.coefs[0];
		heya+=s->allPnts[k].y*s->initPlane.coefs[1];
		heya+=s->allPnts[k].z*s->initPlane.coefs[2];
		heya+=s->initPlane.coefs[3];
		cout<<"eqn 2: "<<heya<<endl;
	}*/
	
	s->perturbedPlane.coefs[0]=mxv(0, sizey2);
	s->perturbedPlane.coefs[1]=mxv(1, sizey2);
	s->perturbedPlane.coefs[2]=mxv(2, sizey2);
	s->perturbedPlane.coefs[3]=mxv(3, sizey2);
	
	if(s->perturbedPlane.coefs[0]==0 &&
	s->perturbedPlane.coefs[1]==0 &&
	s->perturbedPlane.coefs[2]==0 &&
	s->perturbedPlane.coefs[3]==0){
		//cout<<"ZERO FELLAHS"<<endl;
		s->hasPlane=0;
		return false;
	}else{
		return true;
	}
	
	//cout<<"end"<<endl;

}

