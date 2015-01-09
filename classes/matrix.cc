#include "matrix.h"

namespace SOPHYA {

SVD::SVD(TMatrix<double> &a)
: m(a.NRows()), n(a.NCols()), u(a) , v(n,n), w(n)
// The single argument is A.  The SVD computation is done by
// decompose and the results sorted by reorder
// NOTE matrix A is copied to matrix u in this construtor
{
	
	eps=1e-6;
	decompose();
	reorder();
	tsh=0.5*sqrt(m+n+1.)*w(0)*eps; // default threshold for nonzero 
					       // singular values

};

//	int SVD::rank(double thresh=-1.)
//		// Return the rank of A, after zeroing any singular values smaller than thresh.
//		// If thresh is negative, a defulat value based on estimated roundoff is used
//	{
//		int j,nr=0;
//		tsh = ( thresh>=0. ? thresh : 0.5*sqrt(m+n+1.)*w(0)*eps);
//		for (j=0; j<=n; j++)
//		{
//			if (w(j)>tsh)
//				nr++;
//		}
//		return nr;
//	};

//	int SVD::nullity(double thresh=-1.)
//		// Return the nullity of A, after zeroing any singular values smaller than thresh.
//		// Default value as above
//	{
//		int j,nn=0;
//		tsh = ( thresh>=0. ? thresh : 0.5*sqrt(m+n+1.)*w(0)*eps);
//		for (j=0; j<=n; j++)
//		{
//			if (w(j)<=tsh)
//				nn++;
//		}
//		return nn;
//	};
//
//	TMatrix<double> SVD::range(double thresh=-1.)
//		// Give an orthonormal basis for the range of A as the columns of a returned matrix.
//		// thresh as above
//	{
//		int i,j,nr=0;
//		TMatrix<double> rnge(m,rank(thresh));
//		for (j=0; j<n; j++)
//		{
//			if (w(j) > tsh)
//			{
//				for (i=0; i<m;i++)
//					rnge(i,nr)=u(i,j);
//				nr++;
//			}
//		}
//
//		return rnge;
//	};
//
//	TMatrix<double> SVD::nullspace(double thresh=-1.)
//		// Give an orthonormal basis for the nullspace of A as the columns of a returned matrix.
//		// thresh as above
//	{
//
//		int j,jj,nn=0;
//		TMatrix<double> nullsp(n,nullity(thresh));
//		for (j=0; j<n; j++)
//		{
//			if (w(j) <= tsh)
//			{
//				for (jj=0; jj<n;jj++)
//					nullsp(jj,nn)=v(jj,j);
//				nn++;
//			}
//		}
//
//		return nullsp;
//
//	};

void SVD::solve(TVector<double> &b, TVector<double> &x, double thresh =-1. )
// Solves A · X = B for a vector X using the pseudo inverse A 
// as obtained by SVD. If positive, thresh is the threshold value below which
// singular values are considered as zero.  If thresh is negative, a default 
// based on expected roundoff error is used.
// THIS IS THE SAME AS THE SVBKSB ROUTINE
	{
		int jj,j,i;
		double s;

		if (b.Size() != m || x.Size() != n)
			throw ParmError("SVD:: solve bad sizes");

		TVector<double> tmp(n);

		tsh = ( thresh>=0. ? thresh : 0.5*sqrt(m+n+1.)*w(0)*eps);
		for (j=0;j<n;j++) // Calculate U^T B
		{
			s=0.0;
			if (w(j)>tsh) // Nonzero result only if w_j is nonzero
			{
				for (i=0;i<m;i++)
					s+=u(i,j)*b(i);
				s/=w(j);// This is the divide by w_j
			}
			tmp(j)=s;
		}

		for (j=0;j<n;j++)// Matrix multiply by V to get answer
		{
			s=0.0;
			for (jj=0;jj<n;jj++)
				s+=v(j,jj)*tmp(jj);
			x(j)=s;
		}

};


void SVD::solve(TMatrix<double> &b, TMatrix<double> &x, double thresh=-1.)
// Solves m sets of n equations A · X = B using the pseudoinverse of A. The 
// right hand sides are input as b(0..n-1,0..m-1), while x(0..n-1,0..m-1) 
// returns the solutions. thresh as above
{

	int i,j,m=b.NCols();
	if (b.NRows() !=n || x.NRows() != n || b.NCols() != x.NCols() )
		throw ParmError("SVD::solve bad sizes");
	TVector<double> xx(n);
	for (j=0; j<m;j++)
		{
		for (i=0; i<n; i++)
			xx(i)=b(i,j);
		solve(xx,xx,thresh);
		for (i=0;i<n;i++)
			x(i,j)=xx(i);
		}
};


void SVD::decompose() 
// Given the matrix A stored in u(0..m-1,0..n-1), this routine computes its
// singular value decomposition, A = U · W · V^T and stores the results in the
// matrices u and v, and the vector w.
// Matrix A is decomposed into 
// * "left" eigenvector matrix U (whose columns are the eigenvectors of A*A^T)
// * diagonal matrix W (or I guess this could be a vector) whose elements are  
//   the singular values of A
// * "right" eigenvector matrix V (whose columns are the eigenvectors of A^T*A)
// THIS IS THE SAME AS THE SVDCMP ROUTINE
{
	bool flag;
	int i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;
	TVector<double> rv1(n);   
	g = scale = anorm = 0.0;// Householder reduction to bidiagonal form.
	for (i=0;i<n;i++) {
		l=i+2;
		rv1(i)=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += abs(u(k,i));
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					u(k,i) /= scale;
					s += u(k,i)*u(k,i);
				    }
				f=u(i,i);
				// Hard coded SIGN operation
				double val;
				if (f>=0.0)
					val=abs(sqrt(s));
				else 
					val=-abs(sqrt(s));
				g = -val;
				h=f*g-s;
				u(i,i)=f-g;
					for (j=l-1;j<n;j++) {
						for (s=0.0,k=i;k<m;k++) s += u(k,i)*u(k,j);
						f=s/h;
						for (k=i;k<m;k++) u(k,j)+= f*u(k,i);
					    }
					for (k=i;k<m;k++) u(k,i) *= scale;
				}
			}
			w(i)=scale *g;
			g=s=scale=0.0;
			if (i+1 <= m && i+1 != n) {
				for (k=l-1;k<n;k++) scale += abs(u(i,k));
				if (scale != 0.0) {
					for (k=l-1;k<n;k++) {
						u(i,k) /= scale;
						s += u(i,k)*u(i,k);
					    }
					f=u(i,l-1);
					// Hard coded SIGN operation
					double val;
					if (f>=0.0)
						val=abs(sqrt(s));
					else 
						val=-abs(sqrt(s));
					g = -val;
					h=f*g-s;
					u(i,l-1)=f-g;
					for (k=l-1;k<n;k++) rv1(k)=u(i,k)/h;
					for (j=l-1;j<m;j++) {
						for (s=0.0,k=l-1;k<n;k++) s += u(j,k)*u(i,k);
						for (k=l-1;k<n;k++) u(j,k) += s*rv1(k);
					    }
					for (k=l-1;k<n;k++) u(i,k) *= scale;
				    }
			    }
			// hard coded the MAX operation below
			if ( anorm > (abs(w(i))+abs(rv1(i))) )
				anorm=anorm;
			else
				anorm=(abs(w(i))+abs(rv1(i)));
		}              
		for (i=n-1;i>=0;i--) { //Accumulation of right-hand transformations.
			if (i < n-1) {
				if (g != 0.0) {
					for (j=l;j<n;j++) //Double division to avoid poss 
						v(j,i)=(u(i,j)/u(i,l))/g; // underflow
					for (j=l;j<n;j++) {
						for (s=0.0,k=l;k<n;k++) s += u(i,k)*v(k,j);
						for (k=l;k<n;k++) v(k,j) += s*v(k,i);
						}
					}
				for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
				}
			v(i,i)=1.0;
			g=rv1(i);
			l=i;
			}
		// Hard coded MIN operation
		int istart;		
		if ( m < n )
			istart=m;
		else
			istart=n;

		for (i=istart-1;i>=0;i--) {//  Accumulation of left-hand transformations.
			l=i+1;
			g=w(i);
			for (j=l;j<n;j++) u(i,j)=0.0;
			if (g != 0.0) {
				g=1.0/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<m;k++) s += u(k,i)*u(k,j);
					f=(s/u(i,i))*g;
					for (k=i;k<m;k++) u(k,j) += f*u(k,i);
				}
				for (j=i;j<m;j++) u(j,i) *= g;
			} else for (j=i;j<m;j++) u(j,i)=0.0;
			++u(i,i);
		}    
		for (k=n-1;k>=0;k--) { //Diagonalization of the bidiagonal form: Loop over
			for (its=0;its<30;its++) {// singular values, and over allowed iterations.
				flag=true;                  
				for (l=k;l>=0;l--) { //Test for splitting.
					nm=l-1;
					if (l == 0 || abs(rv1(l)) <= eps*anorm) {
						flag=false;
						break;
						}	
					if (abs(w(nm)) <= eps*anorm) break;
					}
				if (flag) {        
					c=0.0;  //Cancellation of rv1(l), if l > 0.
					s=1.0;
					for (i=l;i<k+1;i++) {
						f=s*rv1(i);
						rv1(i)=c*rv1(i);
						if (abs(f) <= eps*anorm) break;
						g=w(i);
						h=pythag(f,g);
						w(i)=h;
						h=1.0/h;
						c=g*h;
						s = -f*h;
						for (j=0;j<m;j++) {
							y=u(j,nm);
							z=u(j,i);
							u(j,nm)=y*c+z*s;
							u(j,i)=z*c-y*s;
							}
						}
					}	
				z=w(k);             
				if (l == k) { //Convergence.    
					if (z < 0.0) { //Singular value is made nonnegative.
						w(k) = -z;
						for (j=0;j<n;j++) v(j,k) = -v(j,k);
					}
					break;
				}
				if (its == 29) 
					throw ParmError("no convergence in 30 svdcmp iterations"); 
				x=w(l); // Shift from bottom 2-by-2 minor.
				nm=k-1;
				y=w(nm);
				g=rv1(nm);
				h=rv1(k);
				f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
				g=pythag(f,1.0);
				// Hard code SIGN operation
				double val;
				if (f>=0.0)
					val=abs(g);
				else 
					val=-abs(g);
				f=((x-z)*(x+z)+h*((y/(f+val))-h))/x;            
				c=s=1.0; // Next QR transformation:
				for (j=l;j<=nm;j++) {
					i=j+1;
					g=rv1(i);
					y=w(i);
					h=s*g;
					g=c*g;
					z=pythag(f,h);
					rv1(j)=z;
					c=f/z;
					s=h/z;
					f=x*c+g*s;
					g=g*c-x*s;
					h=y*s;
					y *= c;
					for (jj=0;jj<n;jj++) {
						x=v(jj,j);
						z=v(jj,i);
						v(jj,j)=x*c+z*s;
						v(jj,i)=z*c-x*s;
					}
					z=pythag(f,h);   
					w[j]=z; // Rotation can be arbitrary if z = 0.
					if (z) {
						z=1.0/z;
						c=f*z;
						s=h*z;
					}
					f=c*g+s*y;
					x=c*y-s*g;
					for (jj=0;jj<m;jj++) {
						y=u(jj,j);
						z=u(jj,i);
						u(jj,j)=y*c+z*s;
						u(jj,i)=z*c-y*s;
					}
				}
				rv1(l)=0.0;
				rv1(k)=f;
				w(k)=x;
			}
		}
	
};


void SVD::reorder()
// Given the output of decompose, this routine sorts the singular values, and 
// corresponding columns of u and v, by decreasing magnitude. Also, signs of 
// corresponding columns are flipped so as to maximize the number of positive elements.
{
		
	int i,j,k,s,inc=1;
	double sw;
	TVector<double> su(m), sv(n);
	do { inc *= 3; inc++; } while (inc <= n);     //Sort. The method is Shell’s sort.
		do {                                      //(The work is negligible as compared
			inc /= 3;                             //to that already done in decompose.)
			for (i=inc;i<n;i++) {     
				sw = w(i);
				for (k=0;k<m;k++) su(k) = u(k,i);
				for (k=0;k<n;k++) sv(k) = v(k,i);
				j = i;
				while (w(j-inc) < sw) {
					w(j) = w(j-inc);
					for (k=0;k<m;k++) u(k,j) = u(k,j-inc);
					for (k=0;k<n;k++) v(k,j) = v(k,j-inc);
					j -= inc;
					if (j < inc) break;
				    }
				w[j] = sw;
				for (k=0;k<m;k++) u(k,j) = su(k);
				for (k=0;k<n;k++) v(k,j) = sv(k);
			}
		} while (inc > 1);                                      
		for (k=0;k<n;k++) { //Flip signs.
			s=0;
			for (i=0;i<m;i++) if (u(i,k)< 0.) s++;
			for (j=0;j<n;j++) if (v(j,k)< 0.) s++;
			if (s > (m+n)/2) {
				for (i=0;i<m;i++) u(i,k) = -u(i,k);
				for (j=0;j<n;j++) v(j,k) = -v(j,k);
			}
		}
		
};

double SVD::pythag(const double a, const double b)
// Computes (a^2 + b^2)^(1/2) without destructive underflow or overflow.
// THIS IS THE SAME AS THE DPYTHAG ROUTINE
{
		double absa=abs(a), absb=abs(b);
		// Hard code SQR operation
		return (absa > absb ? absa*sqrt(1.0+(absb/absa)*(absb/absa)) :
			(absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb))));
};


TMatrix<double> SVD::matrix_inv_svd()
// Invert matrix, formula is:
// inv(A) = V * inv(W) * U^T
{

		// solve with svbksb for each column of inverse matrix
		TVector<double> BB(n);
		TMatrix<double> Y(m,n);

		// what is this part doing?
		for (int lcv1=0; lcv1<n; lcv1++)
		{
			for (int lcv2=0; lcv2<m ; lcv2++) // set be according to lcv1th column
			{
				if ( lcv2 == lcv1 )
					BB(lcv2)=1.0;
				else
					BB(lcv2)=0.0;

			}
			TVector<double> x(n);
			// same as solve() routine
			//CALL SVBKSB(UU,WW,VV,na,na,nap,nap,BB,Y(1,lcv1))
			solve(BB,x);

			for (int lcv2=0; lcv2<m; lcv2++)
				Y(lcv2,lcv1)=x(lcv2);
		} // Y is not the inverse of AI

		return Y;

};

	
TArray<double> Transpose(TArray<double> matrix)
{
	int nx=matrix.SizeX();
	int ny=matrix.SizeY();
	TArray<double> transpose;
	int ndim=2;
	sa_size_t mydim[ndim];
	mydim[0]=ny; mydim[1]=nx;
	transpose.SetSize(ndim, mydim);
	
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
			transpose(j,i)=matrix(i,j);
			
	return transpose;
};


TMatrix<double> Mult(TMatrix<double> matrix1, TMatrix<double> matrix2)
{
    
    int nx1=matrix1.NRows();
    int nInner=matrix1.NCols();
    int nx2=matrix2.NRows();
	int ny2=matrix2.NCols();

	if ( nInner!=nx2){
	    stringstream ss, ss2;
	    ss << nInner;
	    ss2 << nx2;
	    string emsg = "Inner matrix dimensions don't agree: .." + ss.str() + ")x(";
	    emsg += ss.str() + "...";
		throw ParmError(emsg);
		}
	
	TMatrix<double> mult(nx1,ny2);

	
	//cout <<"size of matrix 1 = "<<matrix1.NRows()<<"x"<<matrix1.NCols()<<endl;
    //cout <<"size of matrix 2 = "<<matrix2.NRows()<<"x"<<matrix2.NCols()<<endl; 
    //cout <<"size of final matrix = "<<nx<<"x"<<ny<<endl;
	for (int i=0; i<nx1; i++)
		for (int j=0; j<ny2; j++) {
			double sum=0;
			for (int k=0; k<nInner; k++)
				sum+=matrix1(i,k)*matrix2(k,j);
			mult(i,j)=sum;
			}
    //cout <<" multiplied!"<<endl;
	return mult;

};


TArray<double> Mult(TArray<double> matrix1, TArray<double> matrix2)
{

    int n1x=matrix1.SizeX();
    int nInner=matrix1.SizeY();
    int n2x=matrix2.SizeX();
	int n2y=matrix2.SizeY();

	if ( nInner!=n2x)
		throw ParmError("Inner matrix dimensions don't agree");
	
	TArray<double> mult;
	int ndim=2;
	
	sa_size_t mydim[ndim];
	mydim[0]=n1x; mydim[1]=n2y;
	mult.SetSize(ndim, mydim);

    for (int i=0; i<n1x; i++) // row number of matrix1
		for (int j=0; j<n2y; j++) { // col number of matrix2
		
			double sum=0;
			for (int k=0; k<nInner; k++) // col number of matrix1/row number of matrix2
				sum+=matrix1(i,k)*matrix2(k,j);
			mult(i,j)=sum;
			
			}
				
	return mult;

};

/*TArray<double> doMatrixMult(TArray<double> mult, int n1x, int n2y, int nInner)
{ 




};*/

/*TArray<double> columnRowMult(vector<double> columnVector,vector<double> rowVector)
{
		
	int nCol=columnVector.size();
	int nRow=rowVector.size();
	
	TArray<double> mult;
	int ndim=2;
	
	sa_size_t mydim[ndim];
	mydim[0]=nCol; mydim[1]=nRow;
	mult.SetSize(ndim, mydim);

	for (int i=0; i<nRow; i++)
		for (int j=0; j<nCol; j++) {
			mult(j,i)=columnVector[j]*rowVector[i];
			}
				
	return mult;

};*/

/*double vectorMagnitude(vector<double> vect)
{
    int n = vect.size();
    double sum = 0;
    for (int i=0; i<n; i++)
        sum += (vect[i]*vect[i]);
    return sum/2.;
}*/

} // end namespace sophya
