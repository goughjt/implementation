#ifndef MYBUNDLER_H
#define MYBUNDLER_H

#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
#include <string>

#ifndef vec2dd
#define vec2dd vector<vector<double> >
#endif
#ifndef vec2di
#define vec2di vector<vector<int> >
#endif
#ifndef vec2db
#define vec2db vector<vector<bool> >
#endif

using namespace std;

class myBundler {
	private:
		void readBundler(const string folderpath, double mogrifyFactor);
	
	public:
	
	myBundler(const string folderpath, double mogrifyFactor);
	~myBundler();
        
	int numCameras;
	int numPoints;
	
	//FOR CAMERAS
	//Each vec2dd is size numCameras
	vec2dd focalLengths;
	vec2dd rot1;
	vec2dd rot2;
	vec2dd rot3;
	vec2dd translation;
	
	//3D Points
	//Each ved2dd and vec2di
	vec2dd pntPos; //a vector of position 3-vectors
	vec2di pntClr; //a vector of colour 3-vectors
	
	vector<int> pntVwSize;
	vec2di pntVws; //a vector of view-list n-vectors
	vec2di pntKey; //a vector of key-list n-vectors
	vec2dd pntVwx; //a vector of xcoor-list n-vectors in standard Bundler format
	vec2dd pntVwy; //a vector of ycoor-list n-vectors in standard Bundler format

	vec2dd pntVwxScaled; //a vector of xcoor-list n-vectors scaled to [0,1]
	vec2dd pntVwyScaled; //a vector of ycoor-list n-vectors scaled to [0,1]

	void clear_data();

};

#endif


