#include "em.h"

modelGaussian EMAlgorithm::doExpectationMaximization()
{
    // number of data points
    int n = data_.size();

    TArray<double> R = initialization();
    
    bool isConverged = false;
    int t=0; // iteration number
    model_.llhs.push_back(-1e10); // = -inf to start with
    
    
    while (!isConverged && t<maxIter_) {
        t++;
        //cout << " Iteration number ="<<t<< endl;
        
        // maximization step: basically calculates this iteration's parameter
        // values, updates model_
        maximization(R);
        
        // expectation step
        R = expectation();

        // get new labels
        getLabels(R,n);

        // convergence test
        if (  ( model_.llhs[t]-model_.llhs[t-1] ) < (tol_*abs(model_.llhs[t])) )
            isConverged = true;

        }
        
    int convNum = t-1;
    
    if (isConverged)
        cout << "     Converged in "<< convNum <<" steps.\n";
    else {
        cout << "     Not converged in "<< maxIter_ <<" steps.\n";
        convNum = -1;
        }
        
    model_.convNum = convNum;
    return model_;
};

TArray<double> EMAlgorithm::initialization()
{

    int n = data_.size();
    
    // initialize the model to the number of Gaussians and data points
    initializeModel(n);

    // random initialization
    //cout <<"     Number of Gaussian models = "<< nGaussian_ << endl;
    //cout <<"     Number of data points = "<< n << endl;
    
    // for the matrices/vectors
    int ndim=2;
	sa_size_t mydim[ndim];
	
	// matrices/vectors used in this method
	TArray<int> idx, label, u;
    TArray<double> data, dataDash, tmp1, tmp2, tmp3, R;
    // set their sizes
    mydim[0]=1; mydim[1]=n;
	data.SetSize(ndim, mydim);
	mydim[0]=nGaussian_; mydim[1]=1;
	idx.SetSize(ndim, mydim);
    mydim[0]=1; mydim[1]=nGaussian_;
	dataDash.SetSize(ndim, mydim);
	mydim[0]=nGaussian_; mydim[1]=n;
	tmp1.SetSize(ndim, mydim);
	mydim[0]=nGaussian_; mydim[1]=1;
	tmp2.SetSize(ndim, mydim);
	mydim[0]=nGaussian_; mydim[1]=n;
	tmp3.SetSize(ndim, mydim);
	mydim[0]=1; mydim[1]=n;
	label.SetSize(ndim, mydim);
	mydim[0]=1; mydim[1]=nGaussian_;
	u.SetSize(ndim, mydim);
    mydim[0]=n; mydim[1]=nGaussian_;
    R.SetSize(ndim, mydim);
    
    // fill data array
    for (int j=0; j<n; j++)
        data(0,j) = data_[j];
    
    // uniform sample nGaussian_ times between 1:n
    // pick out the randomly chosen nGaussian_ data points
    for (int i=0; i<nGaussian_; i++) {
        int rval = floor(1 + (n+1-1)*rg_.Flat01());
        idx(i,0)=(rval-1); //
        dataDash(0,i)=(data_[idx(i,0)]);
        //cout << dataDash(0,i) << "  ";
        }
    //cout << endl;

    // matrix multiply nGaussian_ data points with all data points
    TArray<double> dataDashTranspose=Transpose(dataDash); 
    TArray<double> tmp1a = Mult(dataDashTranspose,data);

    if (tmp1a.SizeX()!=tmp1.SizeX() || tmp1a.SizeY()!=tmp1.SizeY() )
        throw ParmError("Error tmp1 size wrong");
       
        
    // diag(dataDash'*dataDash)/2: calculate "magnitude" of data vector
    for (int i=0; i<nGaussian_; i++)
        tmp2(i,0) = dataDash(0,i)*dataDash(0,i)/2;

    // subtract magnitude from each dataDash'*data_
    for (int i=0; i<nGaussian_; i++)
        for (int j=0; j<n; j++)
            tmp3(i,j) = tmp1a(i,j)-tmp2(i,0);
    
    // find max value of each dataDash'*data-mag and store the index of that max
    for (int j=0; j<n; j++) {
        vector<double> vals;
        for (int i=0; i<nGaussian_; i++) {
            //cout << tmp3(i,j) << "  ";
            vals.push_back(tmp3(i,j));
            }
        
        int iMax;
        double maxval = findMaximum(vals, iMax);
        label(0,j) = iMax;
        }
    
    // find the unique values of label and all their indices
    // (label should remain unchanged)
    // (u can take any value between 1 and nGaussian_)
    vector<int> vec;
    for (int j=0; j<n; j++) {
        vec.push_back(label(0,j));
        }
    vector<int> ids;
    vector<int> uvec = uniqueVector(vec,ids);
    
    int nu = uvec.size();
    if (nu != nGaussian_)
        cout << "eep: nu != nGaussian_"<<endl;
    
    for (int i=0; i<nGaussian_; i++)
        u(0,i) = uvec[i];
   
    
    // size (n,nGaussian_) filled with 1's and 0's
    // if 1 is in column i (where 1<=i<=nGaussian_) then data point j (where 1<=j<=n)
    // is from Gaussian number 
    
    for (int j=0; j<n; j++) {
        int val = label(0,j);
        if (val == 0 ) {
            R(j,0)= 1.; R(j,1)= 0.; R(j,2)= 0.; }
        else if (val == 1 ) {
            R(j,0)= 0.; R(j,1)= 1.; R(j,2)= 0.; }
        else if (val == 2 ) {
            R(j,0)= 0.; R(j,1)= 0.; R(j,2)= 1.; }
        }    
        

    return R;
    
}
  
  


void EMAlgorithm::maximization(TArray<double> R)
{
    // number of data points
    int n = data_.size();
    
    // temporary model vectors
    vector<double> w(nGaussian_), mu(nGaussian_), variance(nGaussian_);
    
    // find number of each type of Gaussian
    vector<double> nk(nGaussian_);
    for (int i=0; i<nGaussian_; i++)
        nk[i]=0.; // initialize to zero
    for (int j=0; j<n; j++) {
        nk[0] += R(j,0);
        nk[1] += R(j,1);
        nk[2] += R(j,2);
        }
    
    /*cout << "     Number of each type of Gaussian "<<endl;
    for (int i=0; i<nGaussian_; i++)    
        cout << "     "<< nk[i] <<"  ";
    cout << endl;*/
    
    // the current weight, the proportion of each Gaussian component    
    for (int i=0; i<nGaussian_; i++){
        double val = (double)nk[i]/n;
        w[i] = val;
        }
        
    // mean: data_*R matrix multiplication, every element divided by proportion
    // of each Gaussian
    for (int i=0; i<nGaussian_; i++) {
        
        double sum = 0;
        for (int j=0; j<n; j++)
            sum += data_[j]*R(j,i);
            
        mu[i] = (sum/nk[i]);
        }
        
    // variance: sum[(data-mu)^2*R]/nGaussian
    for (int i=0; i<nGaussian_; i++) {
    
        double sum = 0;
        for (int j=0; j<n; j++)
            sum += (data_[j]-mu[i])*(data_[j]-mu[i])*R(j,i);
        sum /= nk[i];
        
        variance[i] = sum;
        }

    // update the model structure
    for (int i=0; i<nGaussian_; i++) {
        model_.mus[i] = mu[i];
        model_.sigmas[i] = sqrt(variance[i]);
        model_.weights[i] = w[i];
        }
};

TArray<double> EMAlgorithm::expectation()
{
    // number of data points
    int n = data_.size();

    // temporary model vectors
    vector<double> w(nGaussian_), mu(nGaussian_), sigma(nGaussian_);

    // retrieve current model
    for (int i=0; i<nGaussian_; i++) {
        mu[i] = model_.mus[i];
        sigma[i] = model_.sigmas[i];
        w[i] = model_.weights[i];
        }
    
    TArray<double> logRho,R,logR;
    TArray<int> Rinteg;
    int ndim=2;
	sa_size_t mydim[ndim];
	mydim[0]=n; mydim[1]=nGaussian_;
    logRho.SetSize(ndim,mydim);
    R.SetSize(ndim,mydim);
    logR.SetSize(ndim,mydim);
    Rinteg.SetSize(ndim,mydim);
    
    for (int j=0; j<n; j++) {
        for (int i=0; i<nGaussian_; i++) {
            // returns log Gaussian pdf for each Gaussian parameter trial
            logRho(j,i) = logGaussPDF(data_[j],mu[i],sigma[i]);
            
            logRho(j,i) += log(w[i]);
            }
        }
        
    // Sum over 2nd dimension (the Gaussians)
    TVector<double> T = computeLogSumExp(logRho,n); 
    
    double llh = 0; // log likelihood
    for (int j=0; j<n; j++) {
        llh += T(j)/n;
        
        for (int i=0; i<nGaussian_; i++){
            logR(j,i) = logRho(j,i) - T(j);
            R(j,i) = exp(logRho(j,i) - T(j));
            }
        }    
    
    // Add loglikelihood value to model
    model_.llhs.push_back(llh); 

    return R;
};


double EMAlgorithm::logGaussPDF(double datum, double mu, double sigma)
{
    double chisq = (datum-mu)/sigma;
    chisq *= chisq;
    double loggpdf = -0.5*log(2*PI*sigma*sigma) - 0.5*chisq;
    return loggpdf;
};

TVector<double> EMAlgorithm::computeLogSumExp(TArray<double> logRho, int n)
{
    // Make the new array
    TVector<double> logSumExp(n);

    // subtract the largest in each row
    for (int j=0; j<n; j++) {
    
        // find the maximum in the row
        double maxval = -1e10;
        for (int i=0; i<nGaussian_; i++) {
            if (logRho(j,i) > maxval)
                maxval = logRho(j,i);
            }
            
        // subtract the max in each row THEN sum over the exp of this row
        double sumexp = 0;
        for (int i=0; i<nGaussian_; i++) {
            sumexp += (exp(logRho(j,i)-maxval)); // subtract max & calc exp and sum up
            }
            
        logSumExp(j) = maxval + log(sumexp);
        
        if (my_isinf(logSumExp(j)))
            logSumExp(j) = maxval;
        }
        
    return logSumExp;
    
};

void EMAlgorithm::getLabels(TArray<double> R, int n)
{
    for (int j=0; j<n; j++) {
        vector<double> vals;
        for (int i=0; i<nGaussian_; i++)
            vals.push_back(R(j,i));
        int iMax;
        findMaximum(vals, iMax);
        model_.labels[j] = iMax;
        }
};
