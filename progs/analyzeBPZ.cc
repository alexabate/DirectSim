/**
 * @file  analyzeBPZ.cc
 * @brief Analyzes BPZ output
 *
 * Could add more information here I think
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2012
 * @date 2012
 *
 */
#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya libraries
//#include "histinit.h"
#include "fiosinit.h"
#include "mydefrg.h"

#include "hpoly.h"
#include "genericfunc.h"
#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "mydefrg.h"
#include "sedfilter.h"
#include "simdata.h"


#define PI 3.141592

void usage(void);
void usage(void) {
	cout << endl<<" Usage: analyzeBPZ [...options...]" << endl<<endl;
	cout << "  Input catalog is catalog created by BPZ and optionally, the PROBS file"<<endl;
	cout << "  that contains the final PDF of each galaxy"<< endl<<endl;

	cout << " -i INCAT:             name of BPZ catalog to read in "<<endl;
	cout << " -p PROBS:             name of file containing galaxy PDF's to read in"<<endl;
	cout << " -o OUTFILE:           name of file to output results to (saved to output/)"<<endl;
	cout << " -n NGAL               number of galaxies in the BPZ catalog"<<endl;
	cout << " -b ZMIN,DZ,NZ         redshift bins to study photo-z distribution in ";
	cout << " [DEFAULT=0,0.6,5] "<<endl;
	cout << " -z DZMIN,DZMAX,NDZ    bins to use for dz distribution histogram";
	cout << " [DEFAULT=-0.5,0.5,200] "<<endl;
	cout << " -m MAGCUT             apply magnitude quality cut, only keep";
	cout << " galaxies with MAG<=MAGCUT [DEFAULT MAGCUT=25.3] "<<endl;
	cout << " -c ODDSCUT            apply odds cut, only keep galaxies with";
	cout << " ODDS>=ODDSCUT [DEFAULT ODDSCUT=0] "<<endl;
	cout << " -s ISSPECBIN          bin based on spec-z instead of phot-z "<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

    cout << " ==== analyzeBPZ.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string bpzcat, bpzprobs, outfileroot, outfile;
    
    // Number of galaxies in the catalog
    long nGal = 1000;
    
    // Redshift bins
    double zMin = 0;
    double dz = 0.6;
    int nz = 5;
    
    // DZ bins
    double dzMin =-0.5;
    double dzMax = 0.5;
    int ndz = 200;
    
    // Quality cuts
    double magCut = 25.3;
    double oddsCut = 0.;
    
    // Things that might want to be arguments later ....
    // Number of comments at head of BPZ catalog
    int nComments = 63;
    // Column number's of quantities in BPZ catalog (zero indexed)
    int colZB = 1;   // Bayesian photo-z
    int colTB = 4;   // Best-fit type
    int colODDS = 5; // ODDS parameter
    int colZS = 9;   // Spectroscopic redshift
    int colM0 = 10;  // Magnitude in prior filter
    
    // base binning on spectroscopic redshift instead of photometric?
    bool isSpecBin = false;
    
    // is PROBS file available?
    bool isProbs = false;
    
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hsi:p:o:b:n:m:c:")) != -1) {
	    switch (c) {
	        case 'i' :
	            bpzcat = optarg;
	            break;
	        case 'p' :
	            bpzprobs = optarg;
	            isProbs = true;
	            break;
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
	            break;
	        case 'b' :
	            sscanf(optarg,"%lf,%lf,%d",&zMin,&dz,&nz);
		        break;
		    case 'n' :
		        sscanf(optarg,"%ld",&nGal);
		        break;
		    case 'm' :
		        sscanf(optarg,"%lf",&magCut);
		        break;
		    case 'c' :
		        sscanf(optarg,"%lf",&oddsCut);
		        break;
		    case 's' :
		        isSpecBin = true;
		        break;
	        case 'h' :
		        default :
		        usage(); return -1;
		    }
	    }
	    
	
    //-- end command line arguments
    cout <<"     Reading in photo-z from catalog "<< bpzcat <<endl;
    cout <<"     "<< nGal <<" galaxies to read in from the file "<< endl;
    cout <<"     Only keeping galaxies with mag <= "<< magCut <<" and ODDS >= "<< oddsCut <<endl;
    cout <<"     Writing out results to files beginning "<< outfileroot <<endl;
    cout <<"     Studying photo-z distribution in "<< nz;
    if (isSpecBin)
        cout <<" spec-redshift bins \n";
    else
        cout <<" phot-redshift bins \n";
    cout <<"     with spacing "<< dz <<" starting from z = "<< zMin <<endl;
    if (isProbs)
        cout <<"     Reading in photo-z PDF's from file "<< bpzprobs <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
    ifstream inp2;
	ofstream outp;
	
	// Open up BPZ catalog
	inp.open(bpzcat.c_str());
	
	// Print the header
	cout <<"     Printing header of file "<< bpzcat << endl;
	char a[256];
	for (int i=0; i<nComments; i++) {
	    inp.getline(a,256);
	    cout <<"     "<< a <<endl;
	    }
	cout << endl;
	    
	// want to compute sigmaz, bz and etaz as a function of zs
	// want to histogram up dz
	
	// Set up histograms
    cout <<"     Set up DZ histograms and stats vectors"<<endl;
    vector<Histo*> dzHistograms, zsHistograms;
    vector<double> sigmaz(nz,0.), biasz(nz,0.), etaz(nz,0.);
    vector<long> nGalsBin(nz,0.);
    double bW;
    cout <<"     "<< ndz <<" bins with "<< dzMin <<" < dz < "<< dzMax <<endl;
    for (int i=0; i<nz; i++) {
        dzHistograms.push_back(new Histo());
        dzHistograms[i]->ReSize(dzMin,dzMax,ndz);
        bW = dzHistograms[i]->BinWidth();
        
        zsHistograms.push_back(new Histo());
        zsHistograms[i]->ReSize(0.,4,100);
	    }
    cout <<"     Bin width of dz histogram = "<< bW << endl;
    cout << endl;
    
    TArray<double> pBinSum;
    long bSize = 1e6;
    char b[bSize];
    double zGridMin = 0.01, dzGrid = 0.01;
    int nZGrid = 1000;
    if (isProbs) {
        // Open up photo-z PDF file
        inp2.open(bpzprobs.c_str());
        cout <<"     Printing header of file "<< bpzprobs << endl;
        inp2.getline(b,bSize);
        cout <<"     "<< b << endl;
        cout << endl;
        
        // Set up P(z) summed over all galaxies in the bin
        sa_size_t ndim = 2;
        sa_size_t mydim[ndim];
        
        
        mydim[0] = nz; mydim[1] = nZGrid;
        pBinSum.SetSize(ndim, mydim);
        pBinSum = 0.;
        }
    
    
    
    // Loop over catalog and add to Histogram 
    cout <<"     Loop over catalog and add to Histo "<<endl;
    string deliminator = "  ";
    string deliminator2 = " ";
    long cnt = 0;
    for (sa_size_t i=0; i<nGal; i++){
		
        cout <<"     On galaxy "<< i+1 <<" of "<< nGal;// <<endl;

        // read in row from files
        inp.getline(a,256);     // BPZ catalog file
        vector<double> row = getDataFileRow(a, deliminator);
        if (isProbs)
            inp2.getline(b, bSize);  // BPZ .probs file (final PDF(z) of each galaxy)


        // galaxy data from BPZ catalog file
        double zs = row[colZS];
        double zp = row[colZB];
        double odds = row[colODDS];
        double m0 = row[colM0];
        double delz = (zp - zs)/(1. + zs);
        
        // If galaxy is within the "optimal" sample
        if ( (m0<=magCut) && (odds>=oddsCut) ) {
        			
	        // Find which redshift bin this delz belongs in
	        int index;
	        if (isSpecBin)
                index=(int)floor((zs-zMin)/dz);
            else 
                index=(int)floor((zp-zMin)/dz);

            // if galaxy lies within binning range ...
            if (index<nz && index>=0) {
                dzHistograms[index]->Add(delz);
                
                if (!isSpecBin)
                    zsHistograms[index]->Add(zs);

                // Add galaxy values to correct bin
                nGalsBin[index]++;
                sigmaz[index] += delz*delz;
                biasz[index] += delz;
                if (abs(delz)>0.15)
                    etaz[index]++;
                    
                // count galaxy
                cnt++;
                
                if (isProbs) {
                    // Add probability to correct bin 
                    vector<double> rowProbs = getDataFileRow(b, deliminator);
                    cout <<", size of row of probabilities = "<< rowProbs.size();
                    for (int j=0; j<nZGrid; j++)
                        pBinSum(index,j) += rowProbs[j+1];
                    }
                    
                    
                } // end if galaxy lies in binning range
            else {
                cout <<", rejecting galaxy: z = " << zs;
               if (index<0)
                   cout <<" and WARNING! z<0!"<<endl;
		        }
            cout << endl;
            
            }// end if galaxy is in optimal sample
        
		}// end loop over galaxies
    cout << endl;
    
    // close the files
    inp.close();
    if (isProbs)
        inp2.close();
    
    
    // Summary stats
    cout <<"     "<< 100.*(1.-(double)cnt/nGal) <<" percent of galaxies thrown out"<< endl;
    cout <<endl; 
    cout <<"     Number of galaxies in each bin: ";
    for (int j=0; j<nz; j++)
        cout <<"     "<< nGalsBin[j] <<"  ";
    cout << endl;
    cout <<"     Bias in each bin: ";
    for (int j=0; j<nz; j++){
        biasz[j]/=nGalsBin[j];
        cout <<"     "<< biasz[j] <<"  ";
        }
    cout << endl;
    cout <<"     Sigma in each bin: ";
    for (int j=0; j<nz; j++){
        sigmaz[j] = sqrt(sigmaz[j]/nGalsBin[j]);// - biasz[j]*biasz[j];
        cout <<"     "<< sigmaz[j] <<"  ";
        }
    cout << endl;
    cout <<"     Catastrophic fraction of galaxies in each bin: ";
    for (int j=0; j<nz; j++){
        etaz[j] = etaz[j]/nGalsBin[j];
        cout <<"     "<< etaz[j] <<"  ";
        }
    cout << endl;
    cout << endl;

				
    // Write histograms of DZ to text file
    outfile = outfileroot + "_dzHistograms.txt";
    cout <<"     Writing DZ histograms to file " <<outfile<<endl;
    outp.open(outfile.c_str(), ofstream::out);
			
    for (int i=0; i<ndz; i++) {
			
        outp << dzHistograms[0]->BinCenter(i) <<"  ";
			    
        for(int j=0; j<nz; j++)
            outp << dzHistograms[j]->operator()(i) <<" ";
        outp << endl;
		}
    outp.close();

    // Write histograms of zs to text file
    outfile = outfileroot + "_zsHistograms.txt";
    cout <<"     Writing zs histograms to file " <<outfile<<endl;
    outp.open(outfile.c_str(), ofstream::out);
			
    for (int i=0; i<ndz; i++) {
			
        outp << zsHistograms[0]->BinCenter(i) <<"  ";
			    
        for(int j=0; j<nz; j++)
            outp << zsHistograms[j]->operator()(i) <<" ";
        outp << endl;
		}
    outp.close();
    
    // Write photo-z statistics to text file
    outfile = outfileroot + "_stats.txt";
    cout <<"     Writing stats to file " <<outfile<<endl;
    outp.open(outfile.c_str(), ofstream::out);
		
    for(int j=0; j<nz; j++) {
		
        outp << zMin + j*dz + dz/2. <<"  ";
        outp << sigmaz[j] << "  ";
        outp << biasz[j] << "  ";
        outp << etaz[j] << "  ";
        outp << nGalsBin[j] << "  ";
        outp << endl;
        
	    }
    outp.close();
    
    if (isProbs) {
        // Write sum of probabilities in each bin to a text file
        outfile = outfileroot + "_sumProbs.txt";
        cout <<"     Writing sum of probabilities to file " <<outfile<<endl;
        outp.open(outfile.c_str(), ofstream::out);
	    for (int i=0; i<nZGrid; i++) {
	        outp << zGridMin + i*dzGrid <<"  ";
	        for(int j=0; j<nz; j++)
                outp << pBinSum(j,i) <<"  ";
	        outp << endl;    
            }
        outp.close();
        }
        
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " analyzeBPZ.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " analyzeBPZ.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " analyzeBPZ.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of analyzeBPZ.cc program  Rc= " << rc << endl;
  return rc;	
};
