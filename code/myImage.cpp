#include <cstdlib>
#include <stdlib.h>
#include "myImage.h"
#include <iostream>
//ABOVE HERE DEFO NEEDED

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
#include <string>
#include <fstream>
#include <mrpt/math/CMatrixD.h>
#include "parameters.h"

using namespace mrpt;
using namespace mrpt::utils;
using namespace mrpt::random;
using namespace mrpt::math;
using namespace std;

/* The cluster assignments and distance values for each pixel. */
//vec2di clusters;
//vec2dd distances;
vector< vector< int > > clusters;
vector< vector< double > > distances;

//CameraT c;

void myImage::display_contours(CvScalar colour, vector<mySuperpixel> &allSuperpixels) {
	
	IplImage* image = cvLoadImage(imagePath_original.c_str(), 1);
	IplImage* lapImage = cvCreateImage(cvGetSize(image), 8, 1);
	
	cvCvtColor(image, lapImage, CV_BGR2GRAY);
	cvSmooth(lapImage, lapImage, CV_GAUSSIAN, 3, 0, 0, 0);
	
	cvLaplace(lapImage, lapImage, 3);
	cvSaveImage(imagePath_lap.c_str(), lapImage);
	
	
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	vector<CvPoint> contours;
	vec2db istaken;
	for (int i = 0; i < w; i++) { 
		vector<bool> nb;
		for (int j = 0; j < h; j++) {
			nb.push_back(false);
		}
		istaken.push_back(nb);
	}
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			int nr_p = 0;	//number of nbrs not same sprpxl index 
			vector<int> nbrs; //neighbours, as sprpxl indices
			for (int k = 0; k < 8; k++) {
				int x = i + dx8[k], y = j + dy8[k];
				//check that x and y are not boundary points
				if (x >= 0 && x < w && y >= 0 && y < h) {
					//if this pixel has not already been assigned and if sprpxls different
					if (!istaken[x][y] && clusters[i][j] != clusters[x][y]) {
						vector<int>::iterator findIt = std::find(nbrs.begin(), nbrs.end(), superpixelIndices[clusters[x][y]]);
						int pos =std::distance(nbrs.begin(), findIt);
						if(findIt == nbrs.end()){
							nbrs.push_back(superpixelIndices[clusters[x][y]]);
						}
						nr_p += 1;
					}
				}
			}
			if (nr_p >= 2) {
				contours.push_back(cvPoint(i,j));
				istaken[i][j] = true;
				for (int k = 0; k < nbrs.size(); k++) {
					mySuperpixel* s = &allSuperpixels[superpixelIndices[clusters[i][j]]];
					vector<int>::iterator findIt = std::find(s->neigh.begin(), s->neigh.end(), nbrs[k]);
					int pos = std::distance(s->neigh.begin(), findIt);
					if(findIt == s->neigh.end()){
						(s->gradients).push_back(0);
						(s->neighC).push_back(0);
						(s->neigh).push_back(nbrs[k]);
					}
					s->neighC[pos]++;
					int value =  ((uchar *)(lapImage->imageData + i*lapImage->widthStep))[j*lapImage->nChannels +0];
					s->gradients[pos]+=value;
					
					mySuperpixel* s2 = &allSuperpixels[nbrs[k]];
					vector<int>::iterator findIt2 = std::find(s2->neigh.begin(), s2->neigh.end(), superpixelIndices[clusters[i][j]]);
					int pos2 = std::distance(s2->neigh.begin(), findIt2);
					if(findIt2 == s2->neigh.end()){
						(s2->gradients).push_back(0);
						s2->neighC.push_back(0);
						s2->neigh.push_back(superpixelIndices[clusters[i][j]]);
					}
					s2->neighC[pos2]++;
					s2->gradients[pos2]+=value;
				}
			}
		}
	}
	for (int i = 0; i < (int)contours.size(); i++) {
		//cout<<"paint it red"<<endl;
		cvSet2D(image, contours[i].y, contours[i].x, colour);
	}
	cvSaveImage(imagePath_contours.c_str(), image);
	cvReleaseImage(&image);
	cvReleaseImage(&lapImage);
	
	/*for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]==-1){
				cout<<"loop6: "<<i<<", "<<j<<", "<<clusters[i][j]<<endl;
				cout<<"SHOULD NOT SEE THIS: "<<clusters[9999999][9999999]<<endl;
			}
		}
	}*/
}

void myImage::colour_with_cluster_means(vector<mySuperpixel> &allSuperpixels) {
	IplImage* image = cvLoadImage(imagePath_original.c_str(), 1);
	vector<CvScalar> colours(centers.size());
	//cout<<"here"<<endl;
	//if(imageIndex==6){cout<<"at start:"<<clusters[519][429]<<endl;}
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			//if(imageIndex==6){cout<<"at inner start:"<<i<<", "<<j<<clusters[519][429]<<endl;}
			int index=clusters[i][j];
			if(index!=-1){
				CvScalar colour = cvGet2D(image, j, i);
				//if(imageIndex==6 && i==639){cout<<"at inner mid1:  "<<i<<", "<<j<<", "<<clusters[519][429]<<", "<<index<<", "<<clusters[i][j]<<endl;}
				colours[index].val[0] += colour.val[0];
				//if(imageIndex==6&& i==639){cout<<"at inner mid2:  "<<i<<", "<<j<<", "<<clusters[519][429]<<endl;}
				colours[index].val[1] += colour.val[1];
				colours[index].val[2] += colour.val[2];
			}
		}
	}
	//if(imageIndex==6){cout<<"at mid:"<<clusters[519][429]<<endl;}
	//cout<<"here"<<endl;
	for (int i = 0; i < (int)colours.size(); i++) {
		colours[i].val[0] /= center_counts[i];
		colours[i].val[1] /= center_counts[i];
		colours[i].val[2] /= center_counts[i];
		mySuperpixel* s = &allSuperpixels[superpixelIndices[i]];
		s->colour.push_back(colours[i].val[0]);
		s->colour.push_back(colours[i].val[1]);
		s->colour.push_back(colours[i].val[2]);
		convexHull(cv::Mat(s->contours), s->hullPnts, false);
	}
	//cout<<"here"<<endl;
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			//cout<<"clusters: "<<clusters[i][j]<<endl;
			if(clusters[i][j]>9999999){cout<<"problem: "<<i<<", "<<j<<", "<<clusters.size()<<", "<<clusters[i].size()<<", "<<colours.size()<<", "<<clusters[i][j]<<endl;}
			CvScalar ncolour = colours[clusters[i][j]];
			//cout<<"ncolour: "<<ncolour.val[0]<<", "<<ncolour.val[1]<<", "<<ncolour.val[2]<<endl;
			cvSet2D(image, j, i, ncolour);
		}
	}
	//cout<<"here"<<endl;
   cvSaveImage(imagePath_clustermeans.c_str(), image);
	cvReleaseImage(&image);
	
	
	/*for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]==-1){
				cout<<"loop5: "<<i<<", "<<j<<", "<<clusters[i][j]<<endl;
				cout<<"SHOULD NOT SEE THIS: "<<clusters[9999999][9999999]<<endl;
			}
		}
	}*/
	
}



void myImage::seeSuperpixelByID(int ID){
	IplImage* image = cvLoadImage(imagePath_original.c_str(), 1);
	ostringstream convert;
	convert << ID;
	string showTitle=(string)"See sprpxl " + convert.str();
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	vector<CvPoint> contours;
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]!=ID){
				for (int k = 0; k < 8; k++) {
					int x = i + dx8[k], y = j + dy8[k];
					if (x >= 0 && x < w && y >= 0 && y < h && clusters[x][y] == ID) {
						contours.push_back(cvPoint(i,j));
						continue;
					}
				}
			}
		}
	}
	for (int i = 0; i < (int)contours.size(); i++) {
		cvSet2D(image, contours[i].y, contours[i].x, CV_RGB(255,0,0));
	}
	cvShowImage(showTitle.c_str(), image);
	cvWaitKey(0);
	cvReleaseImage(&image);
}

void myImage::createSuperpixelContours(vector<mySuperpixel> &allSuperpixels){
	IplImage* image = cvLoadImage(imagePath_original.c_str(), 1);
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	vec2db isTaken;
	for (int i = 0; i < w; i++) { 
		vector<bool> nb;
		for (int j = 0; j < h; j++) {
			nb.push_back(false);
		}
		isTaken.push_back(nb);
	}
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			mySuperpixel* s = &allSuperpixels[superpixelIndices[clusters[i][j]]];
			for (int k = 0; k < 8; k++) {
				int x = i + dx8[k], y = j + dy8[k];
				if (x >= 0 && x < w && y >= 0 && y < h) {
					if(!isTaken[x][y] && clusters[x][y] != clusters[i][j]){
						s->contours.push_back(cvPoint(x,y));
						isTaken[x][y]=true;
						continue;
					}
				}
			}
		}
	}
	cvReleaseImage(&image);
	
	/*for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]==-1){
				cout<<"loop3: "<<i<<", "<<j<<", "<<clusters[i][j]<<endl;
				cout<<"SHOULD NOT SEE THIS: "<<clusters[9999999][9999999]<<endl;
			}
		}
	}*/
}

void myImage::makeStrings(bool createSlicFolder, string fullpath_LowRes){
	string folderpath;
	filename=fullpath_LowRes.substr(fullpath_LowRes.rfind("/")+1,-1);
	folderpath=fullpath_LowRes.substr(0,fullpath_LowRes.rfind("/"));
	if(imageIndex==0){std::system(((string)"mkdir "+folderpath + (string)"/Images_SLIC").c_str());}
	std::string slicfolderpath=folderpath + "/Images_SLIC/";
	imagePath_original=fullpath_LowRes;
	imagePath_labspace=slicfolderpath + filename.substr(0,filename.find(".")) + (string)"_labspace" + filename.substr(filename.rfind("."),-1);
	imagePath_contours=slicfolderpath + filename.substr(0,filename.find(".")) + (string)"_contours" + filename.substr(filename.rfind("."),-1);
	imagePath_clustermeans=slicfolderpath + filename.substr(0,filename.find(".")) + (string)"_clusterMeans" + filename.substr(filename.rfind("."),-1);
	imagePath_lap=slicfolderpath + filename.substr(0,filename.find(".")) + (string)"_lap" + filename.substr(filename.rfind("."),-1);
	
}
/*
void myImage::assignSuperpixels(myBundler &bundler, vector<mySuperpixel> &allSuperpixels){
	for(int pnt=0; pnt<bundler.numPoints; pnt++){
		for(int pntView=0; pntView<bundler.pntVws[pnt].size(); pntView++){
			//I have to change the coord system from bundler to image
			if(bundler.pntVws[pnt][pntView]==imageIndex){
				int sprpxlNum = clusters[floor((w/2)+bundler.pntVwx[pnt][pntView])][floor((h/2)-bundler.pntVwy[pnt][pntView])];
				mySuperpixel* superpixel = &allSuperpixels[superpixelIndices[sprpxlNum]];
				superpixel->allPnts.setSize(3,size(superpixel->allPnts,2)+1);
				superpixel->allPnts(0,size(superpixel->allPnts,2)-1)=bundler.pntPos[pnt][0];
				superpixel->allPnts(1,size(superpixel->allPnts,2)-1)=bundler.pntPos[pnt][1];
				superpixel->allPnts(2,size(superpixel->allPnts,2)-1)=bundler.pntPos[pnt][2];
				//superpixel->pntVwx.push_back(floor((w/2)+bundler.pntVwx[pnt][pntView]));
				//superpixel->pntVwy.push_back(floor((h/2)-bundler.pntVwy[pnt][pntView]));
				//Put mogrifyfactor here?
				superpixel->pntVwx.push_back(bundler.pntVwx[pnt][pntView]);
				superpixel->pntVwy.push_back(bundler.pntVwy[pnt][pntView]);
			}
		}
	}
}
*/

void myImage::assignSuperpixels(vector<Point2D> &measurements, vector<int> &camidx, vector<Point3D> &point_data, vector<int> &ptidx, vector<mySuperpixel> &allSuperpixels){
	vector< vector< int > > ptIdvecOfSpID;
	for(int i=0; i<point_data.size(); i++){
		vector<int> veccy;
		ptIdvecOfSpID.push_back(veccy);
	}
	for(int i=0; i<measurements.size(); i++){
		if(camidx[i]==imageIndex){
			float x,y;
			measurements[i].GetPoint2D(x,y);
			
			x=x*mogrifyFactor;
			y=y*mogrifyFactor;
			
			int sprpxlNum = clusters[floor((w/2)+x)][floor((h/2)+y)];
			ptIdvecOfSpID[ptidx[i]].push_back(superpixelIndices[sprpxlNum]);
			mySuperpixel* s = &allSuperpixels[superpixelIndices[sprpxlNum]];
			float xyz[3];
			point_data[ptidx[i]].GetPoint(xyz);
			TPoint3D p(xyz[0], xyz[1], xyz[2]);
			s->allPnts.push_back(p);
			//LOOK TO SEE IF ANYTHING ELSE HAS HAD THIS PTDIDX YET?
			for(int j=0; j<ptIdvecOfSpID[ptidx[i]].size(); j++){
				vector<int>::iterator findIt = std::find(s->neighB.begin(), s->neighB.end(), ptIdvecOfSpID[ptidx[i]][j]);
				int pos = std::distance(s->neighB.begin(), findIt);
				if(pos== s->neighB.size()){
					s->neighBC.push_back(0);
					s->neighB.push_back(ptIdvecOfSpID[ptidx[i]][j]);
				}
				//cout<<s->neighBC.size()<<endl;
				s->neighBC[pos]++;
			}
			s->pntVwx.push_back(x);
			s->pntVwy.push_back(y);
			
		}
			
	}
	
	/*for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]==-1){
				cout<<"loop4: "<<i<<", "<<j<<", "<<clusters[i][j]<<endl;
				cout<<"SHOULD NOT SEE THIS: "<<clusters[9999999][9999999]<<endl;
			}
		}
	}*/
	
}
	
//myImage::myImage(myBundler &bundler, int imIndex, string fullpath, string fullpath_LowRes, vector<mySuperpixel> &allSuperpixels, int nr_superpixels, int nc) {
myImage::myImage(vector<Point2D> &measurements, vector<int> &camidx, vector<Point3D> &point_data, vector<int> &ptidx, int imIndex, string fullpath, string fullpath_LowRes, vector<mySuperpixel> &allSuperpixels, int num_sprpxl) {
	//cout<<"hi"<<endl;
	imageIndex=imIndex;
	//cout<<"hi"<<endl;
	makeStrings(imageIndex==0, fullpath_LowRes);
	//printf(" generating superpixels...\n");
 	generate_superpixels(allSuperpixels, num_sprpxl); //takes v long time
	//printf(" finished superpixels\n");
	//cout<<"hi"<<endl;
 	//create_connectivity();
	////cout<<"hi"<<endl;
	//assignSuperpixels(bundler, allSuperpixels);
	createSuperpixelContours(allSuperpixels);
	//cout<<"hi"<<endl;
	assignSuperpixels(measurements, camidx, point_data, ptidx, allSuperpixels);
	//cout<<"hi"<<endl;
	colour_with_cluster_means(allSuperpixels);
	//cout<<"OUT"<<endl;
	display_contours(CV_RGB(255,0,0), allSuperpixels);
	//cout<<"hi"<<endl;
}

myImage::~myImage() {
}

void myImage::generate_superpixels(vector<mySuperpixel> &allSuperpixels, int num_sprpxl) {
	IplImage* image;
	image = cvLoadImage(imagePath_original.c_str(), 1);
	w = image->width;
	h = image->height;
	step = sqrt((w * h) / (double) num_sprpxl);
	cvCvtColor(image, image, CV_BGR2Lab);
	cvSaveImage(imagePath_labspace.c_str(), image);
	clusters.clear();
	distances.clear();
	for (int i = 0; i < w; i++) { 
		vector<int> cr;
		vector<double> dr;
		for (int j = 0; j < h; j++) {
			cr.push_back(-1);
			dr.push_back(FLT_MAX);
		}
		clusters.push_back(cr);
		distances.push_back(dr);
	}
	for (int i = step; i < w - step/2; i += step) {
		for (int j = step; j < h - step/2; j += step) {
			vector<double> center;
			CvPoint ncPnt = find_local_minimum(cvPoint(i,j), image);
			CvScalar colour = cvGet2D(image, ncPnt.y, ncPnt.x);
			center.push_back(colour.val[0]);
			center.push_back(colour.val[1]);
			center.push_back(colour.val[2]);
			center.push_back(ncPnt.x);
			center.push_back(ncPnt.y);
			centers.push_back(center);
			center_counts.push_back(0);
		}
		//cout<<centers.size()<<endl;

		for (int j = 0; j < w; j++) {
			for (int k = 0;k < h; k++) {
				distances[j][k] = FLT_MAX;
			}
		}
		
		int maxj=0;
		for (int j = 0; j < (int) centers.size(); j++) {
			/* Only compare to pixels in a 2 x step by 2 x step region. */
			for (int k = centers[j][3] - step; k < centers[j][3] + step; k++) {
				for (int l = centers[j][4] - step; l < centers[j][4] + step; l++) {
					if (k >= 0 && k < w && l >= 0 && l < h) {
						CvScalar colour = cvGet2D(image, l, k);
						double d = compute_dist(j, cvPoint(k,l), colour);
						if (d < distances[k][l]) {
							distances[k][l] = d;
							clusters[k][l] = j;
							if(j>9999999){
								cout<<"WATH: "<<k<<", "<<l<<", "<<j<<", "<<i<<", "<<j<<endl;
								clusters[66666][88888]=15;	
							}
							if(j>maxj){maxj=j;}
						}
					}
				}
			}
		}
		for (int j = 0; j < (int) centers.size(); j++) {
			centers[j][0] = centers[j][1] = centers[j][2] = centers[j][3] = centers[j][4] = 0;
			center_counts[j] = 0;
		}
		for (int j = 0; j < w; j++) {
	//		cout<<"j: "<<j<<endl;
			for (int k = 0; k < h; k++) {
				
				int c_id = clusters[j][k];
				if (c_id != -1) {
					CvScalar colour = cvGet2D(image, k, j);
					centers[c_id][0] += colour.val[0];
					centers[c_id][1] += colour.val[1];
					centers[c_id][2] += colour.val[2];
					centers[c_id][3] += j;
					centers[c_id][4] += k;
					center_counts[c_id] ++;
				}
			}
		}
		/* Normalize the clusters. */
		for (int j = 0; j < (int) centers.size(); j++) {
			centers[j][0] /= center_counts[j];
			centers[j][1] /= center_counts[j];
			centers[j][2] /= center_counts[j];
			centers[j][3] /= center_counts[j];
			centers[j][4] /= center_counts[j];
		}
		//cout<<"EM algorithm "<<(int)(100*(i+1)/NR_ITERATIONS)<<"\% complete"<<endl;
	}
	cvReleaseImage(&image);
	//make the superpixels
	for(int i=0; i<centers.size(); i++){
		vector<double> c;
		c.push_back(centers[i][0]);
		c.push_back(centers[i][1]);
		c.push_back(centers[i][2]);
		c.push_back(centers[i][3]-w/2);
		c.push_back(centers[i][4]-h/2);
		mySuperpixel superpixel(allSuperpixels.size(), imageIndex, i, c, imagePath_original);
		allSuperpixels.push_back(superpixel);
		superpixelIndices.push_back(allSuperpixels.size()-1);
	}
	/*for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]==-1){
				cout<<"loop1: "<<i<<", "<<j<<", "<<clusters[i][j]<<endl;
				cout<<"SHOULD NOT SEE THIS: "<<clusters[9999999][9999999]<<endl;
			}
		}
	}*/
	//cout<<"end of generate_sprpxls"<<endl;
}

/*
 * Compute the distance between a cluster center and an individual pixel.
 *
 * Input : The cluster index (int), the pixel (CvPoint), and the Lab values of
 *         the pixel (CvScalar).
 * Output: The distance (double).
 */
double myImage::compute_dist(int ci, CvPoint pixel, CvScalar colour) {
    //double dc = sqrt(pow(centers[ci][0] - colour.val[0], 2) + pow(centers[ci][1]
    //        - colour.val[1], 2) + pow(centers[ci][2] - colour.val[2], 2));
    //double ds = sqrt(pow(centers[ci][3] - pixel.x, 2) + pow(centers[ci][4] - pixel.y, 2));
    
    //return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
    //cout<<"nc: "<<nc<<endl;
    //cout<<"ns: "<<ns<<endl;
    double dc_new = (pow(centers[ci][0] - colour.val[0], 2) + pow(centers[ci][1] - colour.val[1], 2) + pow(centers[ci][2] - colour.val[2], 2))/pow(nc,2);
    double ds_new = (pow(centers[ci][3] - pixel.x, 2) + pow(centers[ci][4] - pixel.y, 2))/pow(step,2);
    
    return sqrt(dc_new + ds_new);
    
    //double w = 1.0 / (pow(ns / nc, 2));
    //return sqrt(dc) + sqrt(ds * w);
}

CvPoint myImage::find_local_minimum(CvPoint center, IplImage* image) {
    double min_grad = FLT_MAX;
    CvPoint loc_min = cvPoint(center.x, center.y);
    
    for (int i = center.x-1; i < center.x+2; i++) {
        for (int j = center.y-1; j < center.y+2; j++) {
            CvScalar c1 = cvGet2D(image, j+1, i);
            CvScalar c2 = cvGet2D(image, j, i+1);
            CvScalar c3 = cvGet2D(image, j, i);
            /* Convert colour values to grayscale values. */
            double i1 = c1.val[0];
            double i2 = c2.val[0];
            double i3 = c3.val[0];
            /*double i1 = c1.val[0] * 0.11 + c1.val[1] * 0.59 + c1.val[2] * 0.3;
            double i2 = c2.val[0] * 0.11 + c2.val[1] * 0.59 + c2.val[2] * 0.3;
            double i3 = c3.val[0] * 0.11 + c3.val[1] * 0.59 + c3.val[2] * 0.3;*/
            
            /* Compute horizontal and vertical gradients and keep track of the
               minimum. */
            if (sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3,2)) < min_grad) {
                min_grad = fabs(i1 - i3) + fabs(i2 - i3);
                loc_min.x = i;	
                loc_min.y = j;
            }
        }
    }
    return loc_min;
}

void myImage::create_connectivity() {
    int label = 0, adjlabel = 0;
    const int lims = (w*h) / ((int)centers.size());
    
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    
    /* Initialize the new cluster matrix. */
    vec2di new_clusters;
    for (int i = 0; i < w; i++) { 
        vector<int> ncl;
        for (int j = 0; j < h; j++) {
            ncl.push_back(-1);
        }
        new_clusters.push_back(ncl);
    }

    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            if (new_clusters[i][j] == -1) {
                vector<CvPoint> elements;
                elements.push_back(cvPoint(i, j));
            
                /* Find an adjacent label, for possible use later. */
                for (int k = 0; k < 4; k++) {
                    int x = elements[0].x + dx4[k], y = elements[0].y + dy4[k];
                    
                    if (x >= 0 && x < w && y >= 0 && y < h) {
                        if (new_clusters[x][y] >= 0) {
                            adjlabel = new_clusters[x][y];
                        }
                    }
                }
                
                int count = 1;
                for (int c = 0; c < count; c++) {
                    for (int k = 0; k < 4; k++) {
                        int x = elements[c].x + dx4[k], y = elements[c].y + dy4[k];
                        
                        if (x >= 0 && x < w && y >= 0 && y < h) {
                            if (new_clusters[x][y] == -1 && clusters[i][j] == clusters[x][y]) {
                                elements.push_back(cvPoint(x, y));
                                new_clusters[x][y] = label;
                                count += 1;
                            }
                        }
                    }
                }
                
                /* Use the earlier found adjacent label if a segment size is
                   smaller than a limit. */
                if (count <= lims >> 2) {
                    for (int c = 0; c < count; c++) {
                        new_clusters[elements[c].x][elements[c].y] = adjlabel;
                    }
                    label -= 1;
                }
                label += 1;
            }
        }
    }

	for (int i = 0; i < w; i++) { 
		vector<int> ncl;
		for (int j = 0; j < h; j++) {
			clusters[i][j]=new_clusters[i][j];
		}
	}
	
	/*for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			if(clusters[i][j]==-1){
				cout<<"loop2: "<<i<<", "<<j<<", "<<clusters[i][j]<<endl;
				cout<<"SHOULD NOT SEE THIS: "<<clusters[9999999][9999999]<<endl;
			}
		}
	}*/
}

void myImage::display_center_grid(CvScalar colour, IplImage* image) {
    for (int i = 0; i < (int) centers.size(); i++) {
        cvCircle(image, cvPoint(centers[i][3], centers[i][4]), 2, colour, 2);
    }
}


