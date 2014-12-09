
#ifndef parameters_H
#define parameters_H

//use 0.05 for tanglewood
//use 0.06 for gibsonA
static double rememberDisty=0.04;
static double disty=0.01;

double * const distThreshold = &disty; // will be always pointing to this var
const bool loopDT = true;
const double maxDT = 0.04;
const double incDT = 0.01;

const double mogrifyFactor = 0.15;

static int nry=400;

int * const nr_sprpxl = &nry;
const bool loopNrSp = false;
const int maxNrSp = 650;
const int incNrSp = 50;

const int nc=20;
const int ns=8;

const double qth = 0.1; //plane equality threshold
const string folderPath = "./guitarImages/Tanglewood/";
const string folderPath_LowRes = "./guitarImages/Tanglewood_LowRes/";
const string nvmSuffix = "CMVS/densey.nvm";
const int pointScale = 300;
const int pointSize=5;

const int geomSize = 20;

//const void increaseDT();

//{
//	*distThreshold+=0.01;
//}


#endif

