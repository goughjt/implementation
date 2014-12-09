
//#include "graph.h"
#include "graph.cpp"
#include "maxflow.cpp"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <float.h>
#include <fstream>
#include <iostream>
#include <mrpt/math/ransac.h>
#include <mrpt/math/ops_vectors.h>
#include <mrpt/math/CMatrix.h>
#include <mrpt/math/CMatrixD.h>
#include <mrpt/math/ops_matrices.h>
#include <mrpt/random.h>
#include "myRansac.h"
#include "mySuperpixel.h"
#include "myImage.h"
#include "util.h"
#include "labellingCost.h"
#include "mergeSort.h"

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::random;
using namespace mrpt::math;
using namespace std;

vector<float> averageCoords;
vector<CameraT> camera_data;
vector<Point3D> point_data;
vector<Point2D> measurements;
vector<int> ptidx;
vector<int> camidx;
vector<string> names; 
vector<int> ptc;

void see(string filepath, string title="title"){
	IplImage* image = cvLoadImage(filepath.c_str(), 1);
	cvShowImage(title.c_str(), image);
	cvWaitKey(0);
	cvReleaseImage(&image);
}

void seeThis(IplImage * image, string title="title"){
	cvShowImage(title.c_str(), image);
	cvWaitKey(0);
	cvReleaseImage(&image);
}

TPoint3D superpixelPointTo3D(mySuperpixel* s, float x, float y, TPlane plane){
	float e[9];
	float T[3];
	float f = camera_data[s->imageIndex].GetFocalLength();
	camera_data[s->imageIndex].GetInvertedR9T(e, T);
	TPoint3D t1(-(e[0]*T[0]+e[3]*T[1]+e[6]*T[2]), -(e[1]*T[0]+e[4]*T[1]+e[7]*T[2]), -(e[2]*T[0]+e[5]*T[1]+e[8]*T[2]));
	x=-x/mogrifyFactor;
	y=y/mogrifyFactor;
	TPoint3D t2(e[0]*((x/f)-T[0])+e[3]*((y/f)-T[1])+e[6]*(1-T[2]), e[1]*((x/f)-T[0])+e[4]*((y/f)-T[1])+e[7]*(1-T[2]), e[2]*((x/f)-T[0])+e[5]*((y/f)-T[1])+e[8]*(1-T[2]));
	TLine3D myLine(t1, t2);
	double k = -(plane.coefs[3]+plane.coefs[0]*myLine.pBase.x+plane.coefs[1]*myLine.pBase.y+plane.coefs[2]*myLine.pBase.z)/(plane.coefs[0]*myLine.director[0]+plane.coefs[1]*myLine.director[1]+plane.coefs[2]*myLine.director[2]);
	TPoint3D p(myLine.pBase.x+k*myLine.director[0], myLine.pBase.y+k*myLine.director[1], myLine.pBase.z+k*myLine.director[2]);
	return p;
}



TPoint3D superpixelPointTo3D(mySuperpixel* s, float x, float y, ofstream& osPF, int j, TPlane plane){
	float e[9];
	float T[3];
	float f = camera_data[s->imageIndex].GetFocalLength();
	camera_data[s->imageIndex].GetInvertedR9T(e, T);
	TPoint3D t1(-(e[0]*T[0]+e[3]*T[1]+e[6]*T[2]), -(e[1]*T[0]+e[4]*T[1]+e[7]*T[2]), -(e[2]*T[0]+e[5]*T[1]+e[8]*T[2]));
	x=-x/mogrifyFactor;
	y=y/mogrifyFactor;
	TPoint3D t2(e[0]*((x/f)-T[0])+e[3]*((y/f)-T[1])+e[6]*(1-T[2]), e[1]*((x/f)-T[0])+e[4]*((y/f)-T[1])+e[7]*(1-T[2]), e[2]*((x/f)-T[0])+e[5]*((y/f)-T[1])+e[8]*(1-T[2]));
	TLine3D myLine(t1, t2);
	double k = -(plane.coefs[3]+plane.coefs[0]*myLine.pBase.x+plane.coefs[1]*myLine.pBase.y+plane.coefs[2]*myLine.pBase.z)/(plane.coefs[0]*myLine.director[0]+plane.coefs[1]*myLine.director[1]+plane.coefs[2]*myLine.director[2]);
	TPoint3D p(myLine.pBase.x+k*myLine.director[0], myLine.pBase.y+k*myLine.director[1], myLine.pBase.z+k*myLine.director[2]);
	//cout<<plane.contains(p)<<endl;
	osPF<<s->imageIndex<<":"<<x<<":"<<y<<":"<<s->allPnts[0]<<":"<<s->allPnts[1]<<":"<<s->allPnts[2]<<":"<<p.x<<":"<<p.y<<":"<<p.z<<":"<<t1.x<<":"<<t1.y<<":"<<t1.z<<":"<<t2.x<<":"<<t2.y<<":"<<t2.z<<":"<<k<<", "<<plane.coefs[0]<<":"<<plane.coefs[1]<<":"<<plane.coefs[2]<<":"<<plane.coefs[3]<<":\n";
	return p;
}


TPoint3D superpixelPointTo3D_NORMALS(mySuperpixel* s, float x, float y, TPlane plane){
	float e[9];
	float T[3];
	float f = camera_data[s->imageIndex].GetFocalLength();
	camera_data[s->imageIndex].GetInvertedR9T(e, T);
	TPoint3D t1(-(e[0]*T[0]+e[3]*T[1]+e[6]*T[2]), -(e[1]*T[0]+e[4]*T[1]+e[7]*T[2]), -(e[2]*T[0]+e[5]*T[1]+e[8]*T[2]));
	x=-x/mogrifyFactor;
	y=y/mogrifyFactor;
	TPoint3D t2(e[0]*((x/f)-T[0])+e[3]*((y/f)-T[1])+e[6]*(1-T[2]), e[1]*((x/f)-T[0])+e[4]*((y/f)-T[1])+e[7]*(1-T[2]), e[2]*((x/f)-T[0])+e[5]*((y/f)-T[1])+e[8]*(1-T[2]));
	TLine3D myLine(t1, t2);
	TPoint3D pN(myLine.director[0], myLine.director[1], myLine.director[2]);
	s->vertexNormals.push_back(pN);
	double k = -(plane.coefs[3]+plane.coefs[0]*myLine.pBase.x+plane.coefs[1]*myLine.pBase.y+plane.coefs[2]*myLine.pBase.z)/(plane.coefs[0]*myLine.director[0]+plane.coefs[1]*myLine.director[1]+plane.coefs[2]*myLine.director[2]);
	TPoint3D p(myLine.pBase.x+k*myLine.director[0], myLine.pBase.y+k*myLine.director[1], myLine.pBase.z+k*myLine.director[2]);
	return p;
}

int main(int argc, char *argv[]) {


	char* filename = (char*)(folderPath + nvmSuffix).c_str();
	
	ofstream osPnt("/home/james/mySite/aaaZePoints.html");
	ifstream part1("resources/nvmToJS_1.html");
	osPnt << part1.rdbuf();
	part1.close();	 
	LoadNVM(filename,  camera_data, point_data, measurements, ptidx, camidx, names, ptc, osPnt, pointScale, pointSize, averageCoords);
	ifstream part2("resources/nvmToJS_2.html");
	osPnt << part2.rdbuf();
	part2.close();	
	osPnt.close();
	
	while(*nr_sprpxl<maxNrSp){
	time_t start_time = 	time(NULL);
	//cout<<*nr_sprpxl<<" here"<<endl;
	//cout<<"SANITY CHECK: There are "<<point_data.size()<<" points"<<endl;
	//myBundler bundler(folderPath, mogrifyFactor);	
	vector<myImage> allImages;
	vector<mySuperpixel> allSuperpixels;
	allImages.clear();
	allSuperpixels.clear();
	
 	for(int i=0; i<names.size(); i++){
		//printf("Taking image %i of %i , ", i+1, names.size());
		myImage im(measurements, camidx, point_data, ptidx, i, (char*) (folderPath + names[i]).c_str(), (char*) (folderPath_LowRes + names[i]).c_str(), allSuperpixels, *nr_sprpxl);
		//myImage im(bundler, i, (char*) (folderPath + names[i]).c_str(), (char*) (folderPath_LowRes + names[i]).c_str(), allSuperpixels, nr_sprpxl, nc);
		allImages.push_back(im);
		cout<<i+1<<" of "<<names.size()<<" SLICed w/ "<<im.superpixelIndices.size()<<" superpixels"<<endl;
	}
	
	//cout<<"my_init done"<<endl;
	
	int totalWith;
	int totalWithout;
	int totalWithoutHard;
	int totalWithMerged;
	
	cout<<"distThreshold: totalWith: totalWithout: totalWithoutHard: totalWith: totalWithout: totalWithoutHard: totalWithMerged: removed byFilter: removed byMerge: acceptedPlanes size: total2"<<endl;
	
	while(*distThreshold<maxDT){
	//cout<<"here1"<<endl;

	ostringstream dToS;
	dToS << *distThreshold;
	
	string gtyui;
	if(loopDT){
		gtyui = "/home/james/mySite/aaaZePlanes" + dToS.str() + ".html"; 
	}else{
		gtyui = "/home/james/mySite/aaaZePlanes.html"; 
	}
	ofstream os(gtyui.c_str());
	ifstream part1("resources/nvmToJS_1.html");
	os << part1.rdbuf();
	part1.close();	 

	totalWith=0;
	totalWithout=0;
	totalWithoutHard=0;
	totalWithMerged=0;
	
	/*
	*For loop to make TPlane3d for each sprpxl, if possible
	*/
	//cout<<"Creating initPlanes..."<<endl;
	int maxPntsPrSpr = 0;
	int initPlaneC=0;
	for(int i=0; i<allSuperpixels.size(); i++){
		mySuperpixel* s = &allSuperpixels[i];
		s->hasPlane=0;
		if(s->allPnts.size()>maxPntsPrSpr){
			maxPntsPrSpr-=maxPntsPrSpr;
			maxPntsPrSpr+=s->allPnts.size();
		}
		if(s->allPnts.size()>2){
			myRANSAC(s, *distThreshold);
			if(s->hasPlane==1){
				if(!fitToAllInliers(s)){
					totalWithout++;
				}else{
					totalWith++;
				}
			}else{
				totalWithout++;
			}			
		} else {
			totalWithoutHard++;
			s->hasPlane=2;
		}
		if(totalWith+totalWithout+totalWithoutHard!=i+1){
			cout<<"gone wrong on sprpxl "<<i<<" ... SHOULD NOT SEE THIS!!!"<<endl;
			return 0;
		}
		if((int)100.0*i/allSuperpixels.size()>initPlaneC){
			initPlaneC+=10;
			//cout<<"initPlaneC "<<initPlaneC<<"\% complete"<<endl;
		}
	}

	//loop to test how good my point formula is AND sort sprpxls
	ofstream osPF("/home/james/implementation/pointFormulaTest.txt");
	osPF << "im:x:y:X:Y:Z:x:y:z:t1.x:t1.y:t1.z:t2.x:t2.y:t2.z:k:pl1:pl2:pl3:pl4: \n";
	for(int i=0; i<allSuperpixels.size(); i++){
		mySuperpixel* s = &allSuperpixels[i];
		if(s->hasPlane==1){
			for(int j = 0; j<s->pntVwx.size(); j++){
				TPoint3D p = superpixelPointTo3D(s, s->pntVwx[j], s->pntVwy[j], osPF, j, s->initPlane);
			}
		}
	}
	osPF.close();
	//cout<<"point formula test complete"<<endl;
	
		
	//cout<<"The best of the best?"<<endl;
	int bestTot=0;
	int bestSpr=0;
	for(int i=0; i<allSuperpixels.size(); i++){
		mySuperpixel* s = &allSuperpixels[i];
		s->coverage=0;
		if(s->hasPlane==1){
			for(int k=0; k<point_data.size(); k++){
				float xyz[3];
				point_data[k].GetPoint(xyz);
				TPoint3D p(xyz[0], xyz[1], xyz[2]);
				if(s->initPlane.distance(p)<*distThreshold){
					s->coverage++;	
				}
			}
			//s->coverage/=point_data.size();
			if(s->coverage>bestTot){
				bestTot-=bestTot;
				bestTot+=s->coverage;
				bestSpr-=bestSpr;
				bestSpr+=i;
			}
		}
	}
	
	//cout<<"initializing s_i...";
	vector<int> s_i;
	for(int i=0; i<allSuperpixels.size(); i++){
		s_i.push_back(i);
		//mySuperpixel* s = &allSuperpixels[i];
		//cout<<s->imageIndex<<": "<<s->sprpxlIndex<<": "<<s->coverage<<endl;
	}
	//cout<<"sorting...";
	mergeSort(s_i, 0, s_i.size()-1, allSuperpixels); 
	//cout<<"complete"<<endl;
	/*for(int i=0; i<allSuperpixels.size(); i++){
		mySuperpixel* s = &allSuperpixels[s_i[i]];
		//cout<<s->imageIndex<<": "<<s->sprpxlIndex<<": "<<s->coverage<<endl;
	}*/
	/*return 0;
		*/
		
		
	//cout<<"distThreshold: totalWith: totalWithout: totalWithoutHard: totalWith: totalWithout: totalWithoutHard: totalWithMerged: removed byFilter: removed byMerge: acceptedPlanes size: total2"<<endl;
	cout<<*distThreshold<<": "<<totalWith<<": "<<totalWithout<<": "<<totalWithoutHard<<": ";;
	
	//Section 3.1.2 PLANE FILTERING
	
	vector<TPlane> acceptedPlanes;
	acceptedPlanes.clear();
	int Nmc=20;
	int remByFilter=0;
	int remByMerge=0;
	int filterMerge=0;
	
	vector<double> backCP;
	backCP.push_back(0);
	backCP.push_back(0);
	backCP.push_back(0);
	mySuperpixel sBackground(allSuperpixels.size(), -1, 0, backCP, "does not exist");
	//for each sprpxl w/ plane, 
	for(int i=0; i<allSuperpixels.size(); i++){
		mySuperpixel* s = &allSuperpixels[s_i[i]];
		if(s->hasPlane==1){
			//project 2d hull to plane and get inital vertices
			for(int j=0; j<s->hullPnts.size(); j++){
				TPoint3D p = superpixelPointTo3D(s, s->hullPnts[j].x, s->hullPnts[j].y, s->initPlane);
				s->initialVertices.push_back(p);
			}
			double vertexDisplacement=0;
			for(int j=0; j<Nmc; j++){
				//perturb inliers
				s->perturbedInliers.clear();
				for(int k=0; k<s->inliers.size(); k++){
					s->perturbedInliers.push_back(s->inliers[k]);
					//generate two random angles between 0 and 2pi
					double ang1 = M_PI*(rand()%2);
					double ang2 = M_PI*(rand()%2);
					//generate perturbation vector
					double delz = *distThreshold*sin(ang2);
					double dely = *distThreshold*sin(ang2)*sin(ang2);
					double delx = *distThreshold*sin(ang2)*cos(ang2);
					//perturb inliers
					s->perturbedInliers[k].x+=delx;
					s->perturbedInliers[k].y+=dely;
					s->perturbedInliers[k].z+=delz;
				}
				fitToPerturbedInliers(s);
				//project 2d hull to plane and get sum vertex displacement / Nmc
				for(int k=0; k<s->hullPnts.size(); k++){
					TPoint3D p;
					if(k==s->hullPnts.size()-1){
						p = superpixelPointTo3D_NORMALS(s, s->hullPnts[j].x, s->hullPnts[j].y, s->perturbedPlane);
					}else{
						p = superpixelPointTo3D(s, s->hullPnts[j].x, s->hullPnts[j].y, s->perturbedPlane);
					}
					vertexDisplacement+=pow(pow(s->initialVertices[k].x-p.x,2) + pow(s->initialVertices[k].y-p.y,2) + pow(s->initialVertices[k].z-p.z,2), 0.5)/(Nmc*(s->initialVertices.size()));
				}
			}
			if(exp(-vertexDisplacement/(*distThreshold))<qth||  s->hasPlane==0){
				totalWith--; 
				totalWithout++;
				remByFilter++;
				s->hasPlane = 0;
				//cout<<"adding to totalwithou "<<totalWithout+1<<endl;
			}
		}
		//MERGING HERE
		if(acceptedPlanes.size()>0){
			if(s->hasPlane==0 || s->hasPlane==2){
				s->inliers = s->allPnts;
			}
			for(int j=0; j<acceptedPlanes.size(); j++){
				for(int k=0; k<s->inliers.size(); k++){
					//cout<<acceptedPlanes[j].distance(s->inliers[k])<<" vs "<<*distThreshold<<endl;
					if(acceptedPlanes[j].distance(s->inliers[k])>*distThreshold){
						break;
					}else {
						if(k==s->inliers.size()-1){
							s->initPlane = acceptedPlanes[j];
							s->acceptedPlaneIndex = j;
							switch(s->hasPlane){
								case 0: totalWithout--;
				//cout<<"taking from totalwithou "<<totalWithout<<endl;
											break;
								case 1: totalWith--;
											break;
								case 2: totalWithoutHard--;
											break;
								default: ; //a guess!
								}
							s->hasPlane=3;
							totalWithMerged++;
							remByMerge++;
							j=acceptedPlanes.size();
						}
					}
				}
			}
		}
			//cout<<"x "<<s->hasPlane<<endl;
		if(s->hasPlane==1){
			s->acceptedPlaneIndex = acceptedPlanes.size();
			acceptedPlanes.push_back(s->initPlane);
		}
		bool sendToBack=false;
		if(s->hasPlane==0 || s->hasPlane==2){
			if(s->inliers.size()!=0){
				for(int j=0; j<acceptedPlanes.size(); j++){
					for(int k=0; k<s->inliers.size(); k++){
						if(acceptedPlanes[j].distance(s->inliers[k])<*distThreshold){
							break;
						}else if(k==s->inliers.size()-1 && j==acceptedPlanes.size()-1){
							sendToBack=true;
						}
					}
				}
			}
		}
		if(sendToBack){
			for(int l=0; l<s->inliers.size(); l++){
				TPoint3D p(s->inliers[l].x, s->inliers[l].y, s->inliers[l].z);
				sBackground.allPnts.push_back(p);
			}
			//UNCOMMENT FOR BACKGROUND 1 of 2
			//s->acceptedPlaneIndex=acceptedPlanes.size();
			//s->hasPlane=4;
		}
		
		if((int)100.0*i/allSuperpixels.size()>filterMerge){
			//cout<<"filterMerge "<<filterMerge<<"\% complete"<<endl;
			filterMerge+=10;
		}
	}
	
	cout<<totalWith<<": "<<totalWithout<<": "<<totalWithoutHard<<": "<<totalWithMerged<<": "<<remByFilter<<": "<<remByMerge<<": "<<acceptedPlanes.size()<<": "<<totalWith+totalWithout+totalWithoutHard+totalWithMerged<<":";
	
	if(acceptedPlanes.size()!=0){
		myRANSAC(&sBackground, *distThreshold);
		fitToAllInliers(&sBackground);
		//UNCOMMENT FOR BACKGROUND 2 of 2
		//acceptedPlanes.push_back(sBackground.initPlane);
		//cout<<"background plane "<<sBackground.initPlane.coefs[0]<<", "<<sBackground.initPlane.coefs[1]<<sBackground.initPlane.coefs[2]<<", "<<sBackground.initPlane.coefs[3]<<endl;
	}
	
	//cout<<"Graph cuts.."<<endl;
	
	//try graph cuts here
	for(int alpha=0; alpha< acceptedPlanes.size();alpha++){
		TPlane alphaPlane = acceptedPlanes[alpha];
		//source-> stays as is 0    .......... sink -> switches to alpha 1
		typedef Graph<double,double,double> GraphType;
		GraphType *g = new GraphType( allSuperpixels.size(),  allSuperpixels.size()*4); 
		int totNeigh=0;
		for(int i=0; i<allSuperpixels.size(); i++){
			g -> add_node();
			totNeigh++;
		}
		for(int i=0; i<allSuperpixels.size(); i++){
		//cout<<"about to maxflow1 "<<i<<endl;
			mySuperpixel s1 = allSuperpixels[i];
			if(s1.acceptedPlaneIndex==alpha){
				g -> add_tweights( i,  unaryCost(s1, camera_data[s1.imageIndex], alphaPlane), 9999999); //unary weight
			}else{
				g -> add_tweights( i, unaryCost(s1, camera_data[s1.imageIndex], alphaPlane), unaryCost(s1, camera_data[s1.imageIndex], acceptedPlanes[s1.acceptedPlaneIndex]) ); //unary weight
			}
		//cout<<"about to maxflow 2"<<endl;
			for(int j = 0; j< s1.neigh.size(); j++){
				if(s1.neigh[j]>i ){//otherwise the sprpxl has already been calculated
					mySuperpixel s2 = allSuperpixels[s1.neigh[j]];
					if(s1.acceptedPlaneIndex!=s2.acceptedPlaneIndex && (s1.acceptedPlaneIndex==alpha || s2.acceptedPlaneIndex==alpha)){
						g -> add_edge(i, s1.neigh[j], withinCost(s1, s2, j), withinCost(s1, s2, j) );
					}
				}
			}
		//cout<<"about to maxflow3333333333"<<endl;
			for(int j = 0; j< s1.neighB.size(); j++){
				if(s1.neighB[j]!=i ){//just in case
					mySuperpixel s2 = allSuperpixels[s1.neighB[j]];
					if(s1.acceptedPlaneIndex!=s2.acceptedPlaneIndex && (s1.acceptedPlaneIndex==alpha || s2.acceptedPlaneIndex==alpha)){
						g -> add_edge( i, s1.neighB[j], betweenCost(s1, s2, j), betweenCost(s1, s2, j) );
					}
				}
			}
		}
		
		int flow = g -> maxflow();
		int totalSwitched=0;
		int totalStayed=0;
		
		for(int i=0; i<allSuperpixels.size(); i++){
			mySuperpixel s1 = allSuperpixels[i];
			if(s1.acceptedPlaneIndex!=alpha){
				if (g->what_segment(i) == GraphType::SOURCE){
					//This fellah didnt switch to alpha
					totalStayed++;
				}else{
					totalSwitched++;
					s1.acceptedPlaneIndex = alpha;
				}
			}
		}
		//cout<<"sw, st: "<<totalSwitched<<", "<<totalStayed<<endl;
		delete g;
	}
	
	
	//cout<<"Dense reconstruction.."<<endl;
	
	//cout<<maxPntsPrSpr<<endl;
	for(int i=0; i<allSuperpixels.size(); i++){
		mySuperpixel* s = &allSuperpixels[s_i[i]];
		TPlane plane;
		//if(s->acceptedPlaneIndex!=0){continue;}
		//if(s->hasPlane==1 || s->hasPlane==3 || s->hasPlane==4){
		if(s->hasPlane==1 || s->hasPlane==3){
		//if(1){
			TPlane plane = acceptedPlanes[s->acceptedPlaneIndex];
			TPoint3D p = superpixelPointTo3D(s, s->center[3], s->center[4], plane);
			double rotAng = std::acos(1/pow(1+pow(plane.coefs[0]/plane.coefs[2],2)+pow(plane.coefs[1]/plane.coefs[2],2), 0.5));
			os << "var material = new THREE.MeshBasicMaterial( {color: 0x" <<std::hex << ((1 << 24) + (s->colour[2] << 16) + (s->colour[1] << 8) + s->colour[1]) << std::dec <<"}); \n var geometry = new THREE.BoxGeometry( " << geomSize << ", " << geomSize << ", 1); \n	var cube = new THREE.Mesh( geometry, material ); \n scene.add( cube ); var myAxis = new THREE.Vector3("<<-1000*plane.coefs[1]<<","<<1000*plane.coefs[0]<<","<<0<<"); rotateAroundWorldAxis(cube, myAxis, " << rotAng << "); \n cube.position.x += " << (p.x-averageCoords[0])*pointScale << "; \n cube.position.y += " << (p.y-averageCoords[1])*pointScale << "; \n cube.position.z += " << (p.z-averageCoords[2])*pointScale << "; \n ";
			
			/*cout<<"HERE spr: "<<s->sprpxlIndex<<" in im "<<s->imageIndex<<endl;
			cout<<"work? "<<s->colour[0]<<", "<<s->colour[1]<<", "<<s->colour[2]<<endl;
			cout<<"work? "<<plane.coefs[0]<<", "<<plane.coefs[1]<<", "<<plane.coefs[2]<<", "<<plane.coefs[3]<<endl;
			cout<<"p: "<<p.x<<", "<<p.y<<", "<<p.z<<endl;
			cout<<"allpnts: "<<s->allPnts<<endl;
			cout<<"center: "<<s->center[3]<<", "<<s->center[4]<<endl;
			cout<<"allx: "<<s->pntVwx<<endl;
			cout<<"ally: "<<s->pntVwy<<endl;
			cout<<i<<endl;
					cout<<"now for the eqn test, use each allPnt to feed into coefs"<<endl;
					for(int k=0; k< s->allPnts.size(); k++){
						double heya=0;
						heya+=s->allPnts[k].x*plane.coefs[0];
						heya+=s->allPnts[k].y*plane.coefs[1];
						heya+=s->allPnts[k].z*plane.coefs[2];
						heya+=plane.coefs[3];
						cout<<heya<<endl;
					}
					
			s->see();
			break;*/
				
			/*cout<<"NORMALIZED PLANE COEFS"<<endl;
			double normy = pow(pow(plane.coefs[0],2)+pow(plane.coefs[1],2)+pow(plane.coefs[2],2), 0.5);
			cout<<"-("<<plane.coefs[0]/normy<<"*x + "<<plane.coefs[1]/normy<<"*y + "<<plane.coefs[3]/normy<<")/"<<plane.coefs[2]/normy<<endl;
			cout<<"POINTS"<<endl;
			for(int j=0; j< size(s->allPnts,2); j++){
				cout<<s->allPnts(0, j)<<", "<<s->allPnts(1, j)<<", "<<s->allPnts(2, j)<<endl;
				//cout<<s->pntVwx[j]<<", "<<s->pntVwy[j]<<endl;
			}*/
			/*s->see();
			
			cout<<s->contours.size()<<endl;
			for(int j=0; j< s->contours.size(); j++){
				cout<<s->contours[j].x<<", "<<s->contours[j].y<<endl;
			}
			for(int j=0; j< size(s->allPnts,2); j++){
				cout<<s->allPnts(0, j)<<", "<<s->allPnts(1, j)<<", "<<s->allPnts(2, j)<<endl;
				cout<<s->pntVwx[j]<<", "<<s->pntVwy[j]<<endl;
			}
			for(int j=0; j<4; j++){
				cout<<plane.coefs[j]<<endl;
			}
			cout<<"now for the eqn test, use each allPnt to feed into coefs"<<endl;
			for(int j=0; j< size(s->allPnts,2); j++){
				double heya=0;
				heya+=s->allPnts(0, j)*plane.coefs[0];
				heya+=s->allPnts(1, j)*plane.coefs[1];
				heya+=s->allPnts(2, j)*plane.coefs[2];
				heya+=plane.coefs[3];
				cout<<heya<<endl;
			}
			cout<<"now for the distance test, use each allPnt to feed into distance_fxn"<<endl;
			for(int j=0; j< size(s->allPnts,2); j++){
				//TPoint3D p = TPoint3D(s->allPnts(0, j), s->allPnts(1, j), s->allPnts(2, j));
				cout<<plane.distance(TPoint3D(s->allPnts(0, j), s->allPnts(1, j), s->allPnts(2, j)))<<endl;
			}*/
			
		}else{
			
	//cout<<"MISSING"<<endl;
		}
	}
	    
	ifstream part2("resources/nvmToJS_2.html");
	os << part2.rdbuf();
	part2.close();
	
	os.close();
	
	//cout<<"ALMOST END"<<endl;
	//mySuperpixel* s = &allSuperpixels[s_i[bestSpr]];
	//s->see();
	//cout<<"here"<<endl;
	
	//double energy = lCost(allSuperpixels, camera_data);

	//see("/home/james/implementation/guitarImages/Tanglewood_LowRes/Images_SLIC/tanglewood1_contours.JPG");

	if(!loopNrSp){
		time_t end_time = time(NULL);
		cout<<end_time-start_time<<endl;
	}
	if(loopDT){
		(*distThreshold)+=incDT;
	}else{
		*distThreshold+=maxDT;
	}
	
	}

	if(loopNrSp){
		time_t end_time = time(NULL);
		cout<<end_time-start_time<<endl;
		(*nr_sprpxl)+=incNrSp;
		*distThreshold-=*distThreshold;
		*distThreshold+=rememberDisty;
	}else{
		*nr_sprpxl+=maxNrSp;
	}
	
	}
	cout<<"END"<<endl;


	//std::system("firefox ~/mySite/aaaZePlanes.html");


	return 1;

}



