
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <float.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "myBundler.h"
//#include "util.h"

using namespace std;

myBundler::myBundler(const string folderpath, double mogrifyFactor) {
    clear_data();
    readBundler(folderpath, mogrifyFactor);
}

myBundler::~myBundler() {
    clear_data();
}

void myBundler::clear_data() {
	numCameras=0;
	numPoints=0;
	focalLengths.clear();
	rot1.clear();
	rot2.clear();
	rot3.clear();
	translation.clear();
	pntPos.clear(); //a vector of position 3-vectors
	pntClr.clear(); //a vector of colour 3-vectors
	pntVws.clear(); //a vector of view-list n-vectors
	pntKey.clear(); //a vector of key-list n-vectors
	pntVwx.clear(); //a vector of xcoor-list n-vectors
	pntVwy.clear(); //a vector of ycoor-list n-vectors
}


void myBundler::readBundler(const string folderpath,double mogrifyFactor) {
	string s;
	double d;
	int inty;
	
	const string bundlerSuffix = "CMVS/densey.nvm.cmvs/00/bundle.rd.out";
	char* filename = (char*)(folderpath + bundlerSuffix).c_str() ;

	ifstream in(filename);
	getline(in,s);
	in>>s;
	numCameras=atoi(s.c_str());
	in>>s;
	numPoints=atoi(s.c_str());

//	mrpt::math::CMatrixDouble dataIn(3, numPoints);

	for(int i=0;i<numCameras;i++){
	
		vector<double> fl;
		for (int j = 0; j < 3; j++) {
			in>>s;
			d=strtod(s.c_str(),NULL);
			fl.push_back(d);
		}
		focalLengths.push_back(fl);
		
		vector<double> r1;
		for (int j = 0; j < 3; j++) {
			in>>s;
			d=strtod(s.c_str(),NULL);
			r1.push_back(d);
		}
		rot1.push_back(r1);
		
		vector<double> r2;
		for (int j = 0; j < 3; j++) {
			in>>s;
			d=strtod(s.c_str(),NULL);
			r2.push_back(d);
		}
		rot2.push_back(r2);
		
		vector<double> r3;
		for (int j = 0; j < 3; j++) {
			in>>s;
			d=strtod(s.c_str(),NULL);
			r3.push_back(d);
		}
		rot3.push_back(r3);
		
		vector<double> trans;
		for (int j = 0; j < 3; j++) {
			in>>s;
			d=strtod(s.c_str(),NULL);
			trans.push_back(d);
		}
		translation.push_back(trans);
	}

	for(int i=0;i<numPoints;i++){
	
		vector<double> fl;
		for (int j = 0; j < 3; j++) {
			in>>s;
			d=strtod(s.c_str(),NULL); //PntPos
			fl.push_back(d);
		}
		pntPos.push_back(fl);
		
		vector<int> r1;
		for (int j = 0; j < 3; j++) {
			in>>s;
			//d=strtod(s.c_str(),NULL);
			inty=atoi(s.c_str());
			r1.push_back(inty);
		}
		pntClr.push_back(r1);
		
		in>>s;
		pntVwSize.push_back(atoi(s.c_str()));
		vector<int> p1;
		vector<int> p2;
		vector<double> p3;
		vector<double> p4;
		for(int j=0; j<pntVwSize[i]; j++){
			in>>s;
			inty=atoi(s.c_str());
			p1.push_back(inty);
			in>>s;
			inty=atoi(s.c_str());
			p2.push_back(inty);
			in>>s;
			d=strtod(s.c_str(),NULL);
			p3.push_back(d*mogrifyFactor);
			in>>s;
			d=strtod(s.c_str(),NULL);
			p4.push_back(d*mogrifyFactor); //y should be negative. Change all bunder refences Think children of fxns)
		}
		pntVws.push_back(p1);
		pntKey.push_back(p2);
		pntVwx.push_back(p3);
		pntVwy.push_back(p4);
	}

	cout<<"SANITY CHECK: There are "<<numPoints<<" points"<<endl;
	cout<<"SANITY CHECK: There are "<<pntPos.size()<<" points"<<endl;
	//GOT ALL BUNDLER OUTPUT NOW ----------------------------------------

}

