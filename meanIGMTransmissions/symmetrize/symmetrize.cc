#include <iostream>
#include <vector>
#include <fstream>
#include <string>
using namespace std;

int main() {

    vector<double> wl;
    vector<double> mu;

    string fname = "./../nonCompleteStats/stats_2.9zSource.dat";
    string oname = "./../meanTransmission_2.9zSource.dat";
    ifstream ifile;
    ifile.open(fname.c_str());

    // Find line length
    string temp;
    int len = 0;
    double temp2;
    getline(ifile, temp);
    getline(ifile, temp);
    while(ifile >> temp2)
        len+=1;
    ifile.close();

    // Read in WL and mean
    ifile.open(fname.c_str());
    for(int i=0; i<len; i++) {
        ifile >> temp2;
        wl.push_back(temp2);
    }
    for(int i=0; i<len; i++) {
        ifile >> temp2;
        mu.push_back(temp2);
    }

    ifile.close();


    // go to wl = 11000A
    double nwl = 0.0;
    while(wl[wl.size()-1] < 11000.e-10) {
        nwl = wl[wl.size()-1]+0.1e-10;
        wl.push_back(nwl);
        mu.push_back(1.0);
    }

    // Output as columns
    ofstream ofile;
    ofile.open(oname.c_str());

    for(int i=0; i<wl.size(); i++)
        ofile << wl[i] << "    " << mu[i] << endl;
    ofile.close();


    return 0;

};
