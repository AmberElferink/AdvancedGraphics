#pragma once
#include "core_settings.h"


//MergeSort Algorithm from: https://www.geeksforgeeks.org/cpp-program-for-quicksort/
//Sorts the indices on the x, y or z. One method for each to prevent exessive checking




/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
	array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partitionX(vector<uint> &indices, vector<float3> &centroids, int low, int high)
{
	int i = (low - 1);  // Index of smaller element 

	float pivot = centroids[high].x;    // pivot 

	for (int j = low; j <= high - 1; j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (centroids[j].x <= pivot)
		{
			i++;    // increment index of smaller element 
			Swap(centroids[i], centroids[j]);
			Swap(indices[i], indices[j]);
		}
	}
	Swap(centroids[i + 1], centroids[high]);
	Swap(indices[i + 1], indices[high]);
	return (i + 1);
}

int partitionY(vector<uint> &indices, vector<float3> &centroids, int low, int high)
{
	int i = (low - 1);  // Index of smaller element 
	float pivot = centroids[high].y;    // pivot 

	for (int j = low; j <= high - 1; j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (centroids[j].y <= pivot)
		{
			i++;    // increment index of smaller element 
			Swap(centroids[i], centroids[j]);
			Swap(indices[i], indices[j]);
		}
	}
	Swap(centroids[i + 1], centroids[high]);
	Swap(indices[i + 1], indices[high]);
	return (i + 1);
}

int partitionZ(vector<uint> &indices, vector<float3> &centroids, int low, int high)
{
	int i = (low - 1);  // Index of smaller element 


	float pivot = centroids[high].z;    // pivot 

	for (int j = low; j <= high - 1; j++)
	{
		// If current element is smaller than or 
		// equal to pivot 
		if (centroids[j].z <= pivot)
		{
			i++;    // increment index of smaller element 
			Swap(centroids[i], centroids[j]);
			Swap(indices[i], indices[j]);
		}
	}

	Swap(centroids[i + 1], centroids[high]);
	Swap(indices[i + 1], indices[high]);
	return (i + 1);
}

void quickSortL1byL2X(vector<uint> &indices, vector<float3> centroids, int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
			at right place */
		int pi = partitionX(indices, centroids, low, high);

		// Separately sort elements before 
		// partition and after partition 
		quickSortL1byL2X(indices, centroids, low, pi - 1);
		quickSortL1byL2X(indices, centroids, pi + 1, high);
	}
}

void quickSortL1byL2Y(vector<uint> &indices, vector<float3> centroids, int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
			at right place */
		int pi = partitionY(indices, centroids, low, high);

		// Separately sort elements before 
		// partition and after partition 
		quickSortL1byL2Y(indices, centroids, low, pi - 1);
		quickSortL1byL2Y(indices, centroids, pi + 1, high);
	}
}

void quickSortL1byL2Z(vector<uint> &indices, vector<float3> centroids, int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
			at right place */
		int pi = partitionZ(indices, centroids, low, high);

		// Separately sort elements before 
		// partition and after partition 
		quickSortL1byL2Z(indices, centroids, low, pi - 1);
		quickSortL1byL2Z(indices, centroids, pi + 1, high);
	}
}


//binary search adapted from: https://www.tutorialspoint.com/binary-search-in-cplusplus
int binarySearchX(const vector<uint> &indices, const vector<float3> &centroids, int p, int r, float num) {
	if (p <= r) {
		int mid = (p + r) / 2;
		float midCoordinateL = centroids[indices[mid]].x;
		float midCoordinateR = centroids[indices[mid + 1]].x;
		if (midCoordinateL <= num && midCoordinateR > num)
			return mid;
		if (midCoordinateL > num)
			return binarySearchX(indices, centroids, p, mid - 1, num);

		return binarySearchX(indices, centroids, mid + 1, r, num);
	}
	return -1;
}

int binarySearchY(const vector<uint> &indices, const vector<float3> &centroids, int p, int r, float num) {
	if (p <= r) {
		int mid = (p + r) / 2;
		float midCoordinateL = centroids[indices[mid]].y;
		float midCoordinateR = centroids[indices[mid + 1]].y;
		if (midCoordinateL <= num && midCoordinateR > num)
			return mid;
		if (midCoordinateL > num)
			return binarySearchY(indices, centroids, p, mid - 1, num);

		return binarySearchY(indices, centroids, mid + 1, r, num);
	}
	return -1;
}

int binarySearchZ(const vector<uint> &indices, const vector<float3> &centroids, int p, int r, float num) {
	if (p <= r) {
		int mid = (p + r) / 2;
		float midCoordinateL = centroids[indices[mid]].z;

		float midCoordinateR = centroids[indices[mid + 1]].z;


		if (midCoordinateL <= num && midCoordinateR > num)
			return mid;
		if (midCoordinateL > num)
			return binarySearchZ(indices, centroids, p, mid - 1, num);
		
		return binarySearchZ(indices, centroids, mid + 1, r, num);
	}
	return -1;
}