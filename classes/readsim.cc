#include "readsim.h"

ReadSim::ReadSim(TArray<r_8> drho, int nbadplanes) 
: d_(drho)
// reads in SimLSS output which is a 3D cube of delta values
// delta is overdensity: (rho-rho^bar) / rho^bar
// the way the FFT is performed in SimLSS means there are 1 or 2 extra planes in the first dimension
// these are removed here
{
	drho_ = d_(Range(0, d_.SizeX()-nbadplanes-1), Range::all(), Range::all());
	cout << "    Mass2Gal::Mass2Gal(TArray<r_4> drho,  int nbadplanes=";
	cout << nbadplanes << ")" << endl<<endl;
	
	cout << "    Print TArray properties:"<<endl;
	drho_.Show();
	cout<<endl;

};


void ReadSim::ReadHeader(FitsInOutFile& fin)
// reads in fits file header so the parameters of the simulated cube can be used
{
    cout <<"     ReadSim::ReadHeader()"<<endl;
	// read in values from fits header
	string DXs, DYs, DZs, DKXs, DKYs, DKZs, NXs, NYs, NZs, zrefs, idmids;
	DKXs=fin.KeyValue("DKX");
	DKYs=fin.KeyValue("DKY");
	DKZs=fin.KeyValue("DKZ");
	
	DXs=fin.KeyValue("DX");
	DYs=fin.KeyValue("DY");
	DZs=fin.KeyValue("DZ");
	
	NXs=fin.KeyValue("NX");
	NYs=fin.KeyValue("NY");
	NZs=fin.KeyValue("NZ");
	
	zrefs=fin.KeyValue("ZREF");
	idmids=fin.KeyValue("KZREF");
	
	// Fourier space spacing
	Dkx_=atof(DKXs.c_str());//(double)(drho_.Info()["DKX"];
	Dky_=atof(DKYs.c_str());//(double)(drho_.Info()["DKY"];
	Dkz_=atof(DKZs.c_str());//(double)(drho_.Info()["DKZ"];
	
	// Read space cube cell spacing
	Dx_=atof(DXs.c_str());//(double)(drho_.Info()["DX"];
	Dy_=atof(DYs.c_str());//(double)(drho_.Info()["DY"];
	Dz_=atof(DZs.c_str());//(double)(drho_.Info()["DZ"];
	
	// Number of pixels in each direction
	Nx_=(sa_size_t)atof(NXs.c_str());//(double)(drho_.Info()["NX"];
	Ny_=(sa_size_t)atof(NYs.c_str());//(double)(drho_.Info()["NY"];
	Nz_=(sa_size_t)atof(NZs.c_str());//(double)(drho_.Info()["NZ"];
	
	
	zref_=atof(zrefs.c_str()); // Redshift of the centre of the cube
	// indices of the centre pixel
	idmidz_=(long)ceil(atof(idmids.c_str())); 
	idmidy_=(long)ceil(atof(NYs.c_str())/2);
	idmidx_=(long)ceil(atof(NXs.c_str())/2);
	
	cout << "    Print SimLSS cube properties ...."<<endl;
	cout << "    dk's: (DKX,DKY,DKZ)="<< Dkx_<<","<< Dky_ <<","<< Dkz_<<endl;
	cout << "    dx's: (DX,DY,DZ)="<< Dx_ <<","<< Dy_ <<","<< Dz_ <<endl;
	cout << "    N's:  (NX,NY,NZ)="<< Nx_ <<","<< Ny_ <<","<< Nz_ <<endl;
	cout << "    zref: "<< zref_ <<endl;
	cout << "    index of centre pixel (z,y,x)="<< idmidz_ <<","<< idmidy_;
	cout <<","<< idmidx_ <<endl<<endl;
	cout <<"     EXIT ReadSim::ReadHeader()"<<endl;
};


sa_size_t ReadSim::CleanNegativeMassCells()
// the way SimLSS works means there can be values of delta<-1
// this is equivalent to nonlinear structure formation
// however it is unphysical as it is equivalent to negative mass
// therefore if the cell has delta <-1 it sets delta=-1
{
    sa_size_t nbad = 0;

    for(sa_size_t iz=0; iz<drho_.SizeZ(); iz++) 
        for(sa_size_t iy=0; iy<drho_.SizeY(); iy++) 
	        for(sa_size_t ix=0; ix<drho_.SizeX(); ix++) 
	            if (drho_(ix, iy, iz)<-1.)  {drho_(ix, iy, iz) =-1; nbad++;}

  cout << "    ReadSim::CleanNegativeMassCells() nbad=" << nbad << endl;
  return nbad;
}
