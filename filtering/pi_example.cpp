#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

void read_lines(char * fname, double * out);

int main() {
	const int lines = 1000;
	double z[lines];
	read_lines("C:\\Users\\hildingur\\Dropbox\\vol_filter (2)\\z.csv", z);
	for(int i = 0; i < lines; i++)
		cout<<z[i]<<endl;
	return 0;
}

void read_lines(char * fname, double * out) {
	std::ifstream fhandle(fname);
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
			out[index++] = atof(line);
	} while(!fhandle.eof());

}