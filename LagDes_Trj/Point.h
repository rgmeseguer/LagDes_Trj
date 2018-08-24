#pragma once
//#include <stdio.h>
#include <fstream>				// std::ofstream, print
#include <vector>				// std::vector
#include <string>				// string
#include <map>					// map

using namespace std;

class Point
{
	//map<vector<string>, double> tau;
	vector<vector<string>> keys;				//List of keys in the map
public:
	Point();
	~Point();
	map<vector<string>, vector<vector<double>>> P;	//Vector of LD values of a point 
	map<vector<string>, vector<double>> Pm;			//Average value of the LD of a point
	map<vector<string>, bool> complete;				//Status of the point
	map<vector<string>, double> length;				//Completion of the point
	map<vector<string>, int> times;					//times a point has been opened


	void addPoint(vector<string>);				//Create a new point using its "location" as a key
	void open(vector<string>);					//Opens an existing point
	void close(vector<string>);					//Closes an opened point
	void SavePav(ofstream&);					//Calculates the average of the points and saves them
	bool isKey(vector<string>key);					//Check the existence of a key

};

