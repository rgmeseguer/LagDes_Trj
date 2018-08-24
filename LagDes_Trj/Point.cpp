#include "Point.h"

#include <iostream>				// std::cout, std::fixed

Point::Point() { keys = {}; }
Point::~Point() {}

void Point::open(vector<string>loc) { complete[loc] = false; times[loc] += 1; } //Opens an existing point
void Point::close(vector<string>loc) { complete[loc] = true; }					//Closes an opened point

//Create a new point using its "location" as a key
void Point::addPoint(vector<string>location)				
{
	P[location] = {};										//Initialize the point
	length[location] = 0;									//Length of the point begins at 0
	complete[location] = true;								//The point begins closed
	times[location] = 0;									//Number of times the point is opened
	keys.push_back(location);								//Add its key to the list of keys

}
//Saves the average of the point in a file
void Point::SavePav(ofstream &sfile)
{
	for (size_t i = 0; i < keys.size(); i++)				//Run over all the keys created
	{
		vector<string> key = keys[i];						//Get the key					

#pragma region Average Calculation

		Pm[key] = { 0.,0. };						//Calculate the Average of the Point
		
		for (size_t j = 0; j < P[key].size(); j++)				//Add all the LD values of the Point
		{

			Pm[key][0] += P[key][j][0];
			Pm[key][1] += P[key][j][1];
		}
		Pm[key][0] /= times[key];								//Divide it by the total number of values
		Pm[key][1] /= times[key];								//Divide it by the total number of values

#pragma endregion

		for (size_t j = 0; j < key.size(); j++)					//Saves the point in the file
		{
			sfile << key[j] << ' ';
		}
		sfile << Pm[key][0] << ' ';
		sfile << Pm[key][1] << ' ';
		sfile << Pm[key][0] + Pm[key][1] << ' ';
		sfile << endl;
	}
}
bool Point::isKey(vector<string>key)
{
	for (size_t i = 0; i < keys.size(); i++)
	{
		if (key == keys[i]) { return true; }
	}
	return false;
}
;



