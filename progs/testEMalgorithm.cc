#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>
#include <algorithm>

// sophya libraries
//#include "histinit.h"
#include "fiosinit.h"

#include "em.h"
#include "mydefrg.h"
#include "geneutils.h"


#define PI 3.141592
/*



*/
void usage(void);
void usage(void) {
	cout << endl<<" Usage: test [...options...]" << endl<<endl;
    cout << " -n : number of distributions to fit "<< endl;
	cout << endl;
    }
    
bool myfunction (int i, int j) {
  return (i==j);
}
    
int main(int narg, char* arg[]) {

	cout << " ==== test.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfile = "testfiles/EMresults.txt";
    string infileroot = "mcode/dist";
    int nModel = 2;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:n:i:")) != -1) {
        switch (c) {
            case 'i' :
	            infileroot = optarg;
	            break;
            case 'o' :
	            outfile = optarg;
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nModel);
                break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

  //-- end command line arguments
  cout << endl;
  
  int rc = 1;  
    try {  // exception handling try bloc at top level
  
    RandomGenerator rg;

	ifstream inp;
	ofstream outp;
	string infile;
	
	outp.open(outfile.c_str());
	
	for (int i=0; i<nModel; i++) {
	
	    
	    if (nModel>1) {
	        stringstream ss;
	        ss << i+1;
	        infile = infileroot + ss.str() + ".txt";
	        }
	    else
	        infile = infileroot;
	    cout <<"     Reading in file "<< infile << endl << endl;
	    inp.open(infile.c_str());

	    vector<double> data;

        int nGaussian = 3;
        int nData = 1000;
        for (int i=0; i<nData; i++) {
            double val;
            inp >> val;
            data.push_back(val);
            }
	    inp.close();
	    
        EMAlgorithm emAlg(data, rg, nGaussian);
        modelGaussian model = emAlg.doExpectationMaximization();
        cout << endl;
        
        cout <<"     Printing unsorted results:"<< endl;
        for (int i=0; i< nGaussian; i++)
            cout << model.mus[i] <<"  ";
        cout << endl;
        for (int i=0; i< nGaussian; i++)
            cout << model.sigmas[i] <<"  ";
        cout << endl;
        for (int i=0; i< nGaussian; i++)
            cout << model.weights[i] <<"  ";
        cout << endl;
        cout << endl;
        
        vector<double> rmeans,rsigs,rweights;
        for (int i=0; i< nGaussian; i++)
            rmeans.push_back(model.mus[i]);
        
        vector<int> isort;
        vector<double> rmeansSorted = sortAndGetIndices(rmeans,isort);  
        for (int i=0; i< nGaussian; i++) {
            rsigs.push_back(model.sigmas[isort[i]]);
            rweights.push_back(model.weights[isort[i]]);
            }
                    
        //cout <<"     Converged in " << model.convNum <<" steps "<<endl;
        //cout <<"     Printing the estimated means of the distribution"<<endl;
        for (int i=0; i< nGaussian; i++){
            //cout << "     "<< rmeansSorted[i] <<"  ";
            outp << rmeansSorted[i] <<"  ";
            }
        //cout << endl;
        //cout <<"     Printing the estimated sigmas of the distribution"<<endl;
        for (int i=0; i< nGaussian; i++){
            //cout << "     "<< rsigs[i] <<"  ";
            outp << rsigs[i] <<"  ";
            
            }
        for (int i=0; i< nGaussian; i++)
            outp << rweights[i] <<"  ";
        outp << endl;
         
        cout <<"     Printing sorted results:"<< endl;
        for (int i=0; i< nGaussian; i++)
            cout << rmeansSorted[i] <<"  ";
        cout << endl;
        for (int i=0; i< nGaussian; i++)
            cout << rsigs[i] <<"  ";
        cout << endl;
        for (int i=0; i< nGaussian; i++)
            cout << rweights[i] <<"  ";
        cout << endl;
        cout << endl;
        cout << endl;

        }

    outp.close();
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " test.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " test.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " test.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of test.cc program  Rc= " << rc << endl;
  return rc;	
}

