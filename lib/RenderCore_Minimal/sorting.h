#pragma once
//MergeSort Algorithm from: https://www.geeksforgeeks.org/cpp-program-for-quicksort/
//Sorts the indices on the x, y or z. One method for each to prevent exessive checking




/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
	array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot. Needed for quicksort */
int partitionX(vector<uint> &indices, vector<float3> &centroids, int low, int high);
int partitionY(vector<uint> &indices, vector<float3> &centroids, int low, int high);
int partitionZ(vector<uint> &indices, vector<float3> &centroids, int low, int high);

//sorts the indices based on the list of centroids (centroids has to be copied, and therefore is not a reference or constant)
void quickSortL1byL2X(vector<uint> &indices, vector<float3> centroids, int low, int high);
void quickSortL1byL2Y(vector<uint> &indices, vector<float3> centroids, int low, int high);
void quickSortL1byL2Z(vector<uint> &indices, vector<float3> centroids, int low, int high);

//returns the index to find the index of the triangle for which the centroid on the left is left of the line, 
//and the centroid on the right, is on the right side of the line for that dimension.
//so centroids[indices["returned int"]] will yield the corresponding centroid.
int binarySearchX(const vector<uint> &indices, const vector<float3> &centroids, int p, int r, float num);
int binarySearchY(const vector<uint> &indices, const vector<float3> &centroids, int p, int r, float num);
int binarySearchZ(const vector<uint> &indices, const vector<float3> &centroids, int p, int r, float num);