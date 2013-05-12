#include "filter_utils.h"

#include <fstream>
#include <iostream>
#include <cstdlib>

void read_lines(string& fname, vector<double>& out) {
	std::ifstream fhandle(fname.c_str());
	char line[1000];
	bool first = true;
	int index = 0;
	do {
		fhandle.getline(line, 1000);
		if(first) {
			cout<<"Skipping the first line"<<endl;
			first = false;
			continue;
		}
		if(!fhandle.eof())
			out.push_back(atof(line));

	} while(!fhandle.eof());

}