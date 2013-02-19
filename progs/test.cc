// testing branching: can i just index the change?
#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>
//#include <algorithm>

// sophya libraries
//#include "histinit.h"
#include "fiosinit.h"
#include "mydefrg.h"

#include "geneutils.h"
#include "matrix.h"
// root
//#include "root_plots.h"
//#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TPrincipal.h"
#include <cmath>
//#include "iomanip.h"
#include "TRandom.h"
#define PI 3.141592
/*



*/
void usage(void);
void usage(void) {
	cout << endl<<" Usage: test [...options...]" << endl<<endl;
		
	cout << endl;
    }
    
bool myfunction (int i, int j) {
  return (i==j);
}

/*template<class T> struct index_cmp {
        index_cmp(const T arr) : arr(arr) {}
        bool operator()(const size_t a, const size_t b) const
            { return arr[a] < arr[b]; }
        const T arr;
        };
template<class T> 
class index_cmp {
    const T arr;
    public:
        index_cmp(const T arr) : arr(arr) {};
        bool operator()(const size_t a, const size_t b) const
            { return arr[a] < arr[b]; }
        
        };*/
    
int main(int narg, char* arg[]) {

	cout << " ==== test.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfile = "testfiles/TESTCODE.txt";
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:")) != -1) 
	{
	switch (c) 
		{
	  case 'o' :
	    outfile = optarg;
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
  
    // testing vector sorting
    
    /*vector<double> myVector;
    myVector.push_back(5.);
    myVector.push_back(1.);
    myVector.push_back(6.);
    
    for (int i=0; i<myVector.size(); i++)
        cout << myVector[i] <<"  ";
    cout << endl;
    
    sort(myVector.begin(),myVector.end());
    for (int i=0; i<myVector.size(); i++)
        cout << myVector[i] <<"  ";
    cout << endl;*/
    
    
    

    vector<double> a;
    a.push_back(5.); a.push_back(1.); a.push_back(7.);
    cout << " Initial a: ";
    for (int i=0; i<a.size(); i++)
        cout << a[i] <<"  ";
    cout << endl;
    
    vector<int> indices;
    vector<double> sorteda = sortAndGetIndices(a, indices);

    cout << " Final a: ";
    for (int i=0; i<sorteda.size(); i++)
        cout << sorteda[i] <<"  ";
    cout << endl;
    cout << " Indices: ";
    for (int i=0; i<indices.size(); i++)
        cout << indices[i] <<"  ";
    cout << endl;
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

