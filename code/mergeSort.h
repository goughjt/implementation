#ifndef mergeSort_h
#define mergeSort_h

#include<vector>
#include "mySuperpixel.h"

using namespace std;

void merge(vector<int> &s_i, int first, int mid, int last, vector<mySuperpixel> &sprpxls){
	int temp[s_i.size()];
	int first1=first, last1 = mid;
	int first2=mid+1, last2 = last;
	int index=first1;
	while(first1<=last1 && first2<=last2){
		mySuperpixel s1 = sprpxls[s_i[first1]];
		mySuperpixel s2 = sprpxls[s_i[first2]];
		if(s1.coverage>s2.coverage){
			temp[index]=s_i[first1];
			first1++;
		}else{
			temp[index]=s_i[first2];
			first2++;
		}
		index++;
	}
	while(first1<=last1){
		temp[index]=s_i[first1];
		first1++;
		index++;
	}
	while(first2<=last2){
		temp[index]=s_i[first2];
		first2++;
		index++;
	}
	for(index = first; index<=last; index++){
		s_i[index] = temp[index];
	}
		
}


void mergeSort(vector<int> &s_i, int min, int max, vector<mySuperpixel> &sprpxls){
	if(min<max){
		int mid=(min+max)/2;
		mergeSort(s_i, min, mid, sprpxls);
		mergeSort(s_i, mid+1, max, sprpxls);
		merge(s_i, min, mid, max, sprpxls);
	}
}

#endif
