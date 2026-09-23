/*
 *
 * Copyright (C) 2024 Juan Domingo, Esther Dura, Guillemo Ayala, Maite Leon ({Juan.Domingo,Esther.Dura,Guillermo.Ayala, Teresa.Leon}@uv.es)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*************************************************************************************
 *
 * Set of functions to implement the Standard GLM algortihm, as decribed in the book
 * "Negative Binomial Regression", Joseph M. Hilbe, 2nd ed., page 52
 *
 *************************************************************************************/
 
#include <blissowen.h>

//#include<fstream>

extern unsigned char DEB;

//' BlissOwen
//'
//' Gets a vector of p observations assuming they come from a negative binomial distribution with overdispersion
//' parameter Phi (var=mu+Phi*mu^2) and uses the Bliss-Owen method to get an estimator for Phi
//'
//' @param  y					A vector of observations, size p
//' @return VPhi["Phi","valid"] A list with a real value, the Bliss-Owen estimation for the overdispersion parameter, and a 
//' boolean mark to indicate if the estimation is considered valid. It will not be when statistics u, v or both are negative or null.
//' @examples
//' n = 200
//' mu = 30
//' size = .7
//' y = rnbinom(n = n, size = size, mu= mu)
//' phi_est = BlissOwen(y)
//' @export 
// [[Rcpp::export]]
Rcpp::List BlissOwen(Rcpp::NumericVector y)
{
 int p=y.length();
  
 double m=0.0,sum2=0.0;
 for (int i=0; i<p; i++)
 {
  m += y[i];
  sum2 += y[i]*y[i];
 }
 m /= double(p);
 double s2 = (sum2-double(p)*m*m)/double(p-1);
 double u=(m*m)-(s2/double(p));
 double v=s2-m;
 
 bool valid = (u>0.0) && (v>0.0);
 /*
 if (u<=0.0 && v>0.0)
  Rcpp::warning("Calculation of single Phi with BlissOwen: u-statistic is negative or null (u=%f)\n",u);
 if (u>0.0 && v<=0.0)
  Rcpp::warning("Calculation of single Phi with BlissOwen: v-statistic is negative or null (v=%f)\n",v);
 if (u<=0.0 && v<=0.0)
  Rcpp::warning("Calculation of single Phi with BlissOwen: u and v-statistics are negative or null (u=%f, v=%f)\n",u,v);
 */
 
 Rcpp::List ret=Rcpp::List::create(
                                   Rcpp::Named("Phi")=v/u,
 								   Rcpp::Named("valid")=valid);
 
 return(ret);
}

//' CommonBlissOwen
//'
//' Gets a matrix of n samples of p observations each given as a nxp matrix (one sample per row) and, assuming that all 
//' of them come from a negative binomial ditribution with a common overdispersion Phi value (var=mu+Phi*mu^2), uses
//' the Bliss-Owen method to get a common estimator for Phi using all the samples for which estimation is possible.
//' These samples are marked as true in an array of boolean values which is returned. Estimation is considered impossible
//' if the v-statistic (s^2-mu) is lower or equal than a threshold (0 by default). 
//' 
//' When estimation is possible, the method proceeds iteratively until either the difference of estimated values between
//' successive iterations is below the tolerance or until the fixed number of iterations has been reached.
//'
//' @param  y		 A matrix of sample, size n x p
//' @param  maxit	 Maximum number of allowed iterations. Default: 100
//' @param  tol      Tolerance to exit if change is below its value. Default: 1e-9
//' @param  posthres Threhold to consider a sample valid if its v-value is strictly above it. Default: 0.0
//' @return VPhi["Phi","valid"] A list of two element, a real number with the Phi estimated value and a  boolean array with length=nrows(y).
//' @examples
//' n = 100
//' N = 400
//' mu = 1:N
//' size = .7
//' y = t(sapply(mu,function(i) rnbinom(n = n, size = size, mu= i)))
//' VPhi=CommonBlissOwen(y)
//' cat("Estimated value: ",VPhi$Phi," using ",sum(VPhi$valid)," of all the samples.\n")
//' @export
// [[Rcpp::export]]
Rcpp::List CommonBlissOwen(Rcpp::NumericMatrix y,int maxit=100,double tol=1e-9,double posthres=0.0)
{
 int n=y.nrow();
 int p=y.ncol();
 
 if (n<2)
  Rcpp::stop("Incorrect size of input data matrix (less than 2 rows).\n");
 
 if (maxit<1)
  Rcpp::stop("Incorrect value for maximum iterations (less than 1).\n");
 
 if (tol<0.0)
  Rcpp::stop("Incorrect value for tolerance. It must be a possitive number.\n");
    
 if (DEB)
  Rcpp::Rcout << "Tolerance set to " << ((tol==1e-9) ? "its default value, " : "") << tol << "\n";
 
 // Vector with means of each row
 Rcpp::NumericVector m(n);
 for (int row=0; row<n; row++)
 {
  m[row]=0.0;
  for (int col=0; col<p; col++)
   m[row] += y(row,col);
  m[row] /= double(p);
 }
 
 // Vector with variances of each row
 Rcpp::NumericVector s2(n);
 // Vectors with statistics u and v
 Rcpp::NumericVector u(n);
 Rcpp::NumericVector v(n);
 // Vector with means of only good values
 Rcpp::NumericVector mgood(n);
 Rcpp::LogicalVector valid(n);
 
 // Same loop used to fill all vectors...
 double sum2;
 double ut,vt,s2t;
 int ng=0;
 for (int row=0; row<n; row++)
 {
  sum2=0.0;
  for (int col=0; col<p; col++) 
   sum2 += y(row,col)*y(row,col);
   
  // Var = (sum(y^2)-n*mean(y)^2)/(n-1) 
  s2t = (sum2 - double(p)*m[row]*m[row])/double(p-1);
  // u = mean(y)^2-(s2/n)
  ut = m[row]*m[row]-(s2t/double(p));
  // v = s2-mean(y)
  vt = s2t-m[row];
  
  valid[row] = (vt>posthres);
  if (valid[row])
  {
   s2[ng] = s2t;
   u[ng] = ut;
   v[ng] = vt;
   mgood[ng] = m[row];
   ng++;
  }
 }
 
 if (DEB)
 {
  if (ng<n)
   Rcpp::Rcout << n-ng << " samples have been discarded because v value was below " << posthres <<"\n";
  else
   Rcpp::Rcout << "All samples have been accepted (no negative values for v)\n";
 }
  
 // From now on, and to compare with the original notation of BlissOwen58:
 // BlissOwen	  --> This program
 // k             --> k
 // \overline{u}  --> mgood
 // \overline{x'} --> u (yes, our u for each row is really averaged...)
 // \overline{y'} --> v (as for u)
 // m             --> Bliss suggests using what he calls \overline{u} (which is our mgood).
 // m^2           --> Bliss suggests using \overline{x'}-s^2/N (which is our u)
 
 // Initial estimation of k: mean(u/v)
 double k=0.0;
 for (int row=0; row<ng; row++)
  k += (u[row]/v[row]);
  
 k /= double(ng);
 
 // Value of Phi
 double Phi=1/k;
 double initial_Phi=Phi;
 
 if (DEB)
  Rcpp::Rcout << "Starting value: Phi=" << initial_Phi << "\n";
  
 // Auxiliary variable to contain k^2
 double k2;

 // Value of Phi in previous iteration
 double Phi_old;
 
 int iter=0;
 
 // Variables to contain sum(w*u*u) and sum(w*u*v)
 double u1s,v1s;
 // The weight at each point
 double w;
 // A factor used as part of the weight that must be calculated for each iteration,
 // but it is common for all rows
 double fac;
 
 // Current update of Phi between consecutive iterations
 double upd;
 
 // Loop starts here. It makes always at least one iteration and end either when the limit of iterations has been
 // reached, or when the estimated Phi changes less than a small tolerance
 do
 {  
  Phi_old = Phi;
  
  // Just auxiliary to simplify formula
  k2=k*k;
  
  // The first factor of Eq. (8) of BlissOwen58 which does not depend on the row and therefore can be precalculated
  // outside of the loop: 0.5(N-1)k^4/(k(k+1)+(2k-1)/N-3/N^2)
  fac = 0.5*(double(ng)-1)*k2*k2/(k*(k+1)-((2*k-1)/double(ng))-(3/double(ng*ng)));
  
  // sum(x'^2) and sum(x'y') are initialized to zero and incremented in the loop which runs over all rows
  u1s=v1s=0.0;
  for (int row=0; row<ng; row++)
  {
   w = fac/(mgood[row]*mgood[row]*(mgood[row]+k)*(mgood[row]+k));
   u1s += w*u[row]*u[row];  // This is w*sum(x'*x') in Eq. (10) of BlissOwen58. Look at the former equivalences of symbols
   v1s += w*u[row]*v[row];  // This is w*sum(x'*y') in same eq.
  }
  // Eq. (10), inverted (we calculate Phi, not k)
  Phi = v1s/u1s;
  k = 1/Phi;
  
  // how much Phi has changed since former iteration
  upd=fabs(Phi-Phi_old);
  
  if (DEB)
   Rcpp::Rcout << "Iteration " << iter << ": Phi=" << Phi << "; update=" << upd << "\n";
   
  iter++;
 }
 while ((iter<maxit) && (upd>tol));

 Rcpp::List ret=Rcpp::List::create(
                                   Rcpp::Named("Phi")=Phi,
								   Rcpp::Named("valid")=true);
 
 if (DEB)
 {
  Rcpp::Rcout << "Common value of Phi estimated by Bliss-Owen method for all the cases. Value started at " << initial_Phi << ". ";
  Rcpp::Rcout << "Final value after " << iter << " iterations " << ((iter>=maxit) ? "(the limit)" : "") << " is " << Phi << "\n";
 }
 return(ret);
}

// This is an auxiliary function, not exported to R interface, that implements the same algorithm as in CommonBlissOwen
// but only for a group of numsel selected observations which are supposed to belong to the same 'group' (i.e.: to have
// been generated by a NB-distribution with a common or at least similar Phi value)
double AuxiliaryCommonBlissOwen(int n,Rcpp::NumericVector u,Rcpp::NumericVector v,Rcpp::NumericVector m,const std::vector<bool> &selec,int numsel,int maxit,double tol)
{
 // Initial estimation of k: mean(u/v)
 // As the complete function, but sum and mean is only calculated for the selected 'numsel' cases marked as true in the selec array
 double k=0.0;
 for (int row=0; row<n; row++)
  if (selec[row])
   k += (u[row]/v[row]);
  
 k /= double(numsel);
 // Value of Phi
 double Phi=1/k;
 
 // Auxiliary variable to contain k^2
 double k2;

 // Value of Phi in previous iteration
 double Phi_old;
 
 int iter=0;
 double u1s,v1s;
 // The weight at each point
 double w;
 // A factor used as part of the weight
 double fac;
 // Current update
 double upd;
 do
 {  
  Phi_old = Phi;
  
  k2=k*k;
  // The first factor of Eq. (8) of BlissOwen58 which does not depend on the row and therefore can be precalculated
  // outside of the loop: 0.5(N-1)k^4/(k(k+1)+(2k-1)/N-3/N^2)
  fac = 0.5*(double(numsel)-1)*k2*k2/(k*(k+1)-((2*k-1)/double(numsel))-(3/double(numsel*numsel)));
  
  u1s=v1s=0.0;
  for (int row=0; row<n; row++)
  {
   //Again, this is like the general function but only for the selected observations
   if (selec[row])
   {
    w = fac/(m[row]*m[row]*(m[row]+k)*(m[row]+k));
    u1s += w*u[row]*u[row];  // This is w*sum(x'*x') in Eq. (10) of BlissOwen58
    v1s += w*u[row]*v[row];  // This is w*sum(x'*y') in same eq.
   }
  }
  // Eq. (10), inverted (we calculate Phi, not k)
  Phi = v1s/u1s;
  k = 1/Phi;
  
  upd=fabs(Phi-Phi_old);
   
  iter++;
 }
 while ((iter<maxit) && (upd>tol));
 
 return(Phi);
}

int SortAndExpurge(const std::vector<double> &differences,std::vector<bool> &ingroup,double prop,int currently_sel,bool &expurged)
{
 int nsel=int(prop*differences.size()+0.5);
 if (currently_sel<=nsel)
 {
  expurged=false;
  return(currently_sel);
 } 
 
 std::vector<size_t> idx(differences.size());
 iota(idx.begin(), idx.end(), 0);
 stable_sort(idx.begin(), idx.end(), [&differences](size_t i1, size_t i2) {return differences[i1] < differences[i2];});
 for (size_t i=nsel; i<differences.size(); i++)
  ingroup[idx[i]]=false; 
 expurged=true;
 
 return(nsel);
}

//' GroupedCommonBlissOwen
//'
//' 
//' Gets a matrix of n samples of p observations each given as a nxp matrix (one sample per row) and, assuming that all 
//' of them come from several negative binomial ditributions, each with a different overdispersion Phi value
//' (var_i = mu_i+Phi_i*mu_i^2), uses the Bliss-Owen method to get estimations for the different Phi values when they
//' can be estimated (which is determined by the fact that its v-statistic is strictly above a threshold, 0 by default).
//' 
//' The estimation of each Phi proceeds by estimating first the value for each row with the simple BlissOwen method and
//' then, grouping for each sample all the others whose absolute value of difference of Phi is less than the tolerance and
//' using them to estimate the value through the same method used in the CommonBlissOwen function programmed in this package.\cr
//' The default value for tolphi, which is -1.0, does not mean -1 (which would make no sense) but a mark to indicate that\cr
//' it must be estimated using from formula (5) of BlissOwen[58], which depends on the sample. If you set any other value\cr
//' for tolphi that value will be used for all rows instead of individually estimating the tolerance for each sample.
//'
//' @param  y		A matrix of samples, size n x p
//' @param  maxit	Maximum number of allowed iterations for each Phi estimation. Default: 100 
//' @param  tol     Tolerance to exit if change in each Phi estimation is below its value. Default: 1e-9
//' @param  tolphi  Maximum absolute difference in Phi value to be part of the estimation around a given sample. Default: -1.0, which means 'Estimate from the sample'
//' @param  posthres Threshold to consider a sample valid if its v-value is strictly above it. Default: 0.0
//' @param  down_ndev Integer with the number of stdev to limit interval of grouping by the lower side. Default: 3
//' @param  up_ndev Integer with the number of stdev to limit interval of grouping by the upper side. Default: 2
//' @param  down_maxprop Maximum proportion of the total samples that, being in the interval [mean-down_ndev*stdev,mean], must be considered. Default: 1.0 (consider all the samples that fall in the interval)
//' @param  up_maxprop Maximum proportion of the total samples that, being in the interval [mean,mean+up_ndev*stdev], must be considered. Default: 1.0 (consider all the samples that fall in the interval)
//' @param  niter   Integer with the number of iteration smoothing steps. Default: 2
//' @return A list L of n lists, each one (i.e.: L[[row]]) with the following keys:
//' \itemize{
//'           \item Phi    - A numeric value with the estimated Phi value for L[[row]]
//'           \item valid  - A boolean indicating if row has been kept or discarded according to the v-threshold criterion
//'           \item nsel   - An integer indicating how many samples were selected to be averaged with the sample at this row
//'           \item tolphi - A numeric value with the estimated tolerance for sample grouping for this row
//' }
//' If L[[i]]$valid is false the values of L[[i]]$Phi$ and of L[[i]]$tolphi will be NaN and the value of L[[i]]$nsel will be 0.
//' @examples
//' n = 50
//' size = seq(.2,2,length.out=5)
//' mu = 30
//' y = t(sapply(size,function(i) rnbinom(n = n, size = i, mu= mu)))
//' VPhi=GroupedCommonBlissOwen(y)
//' @export
// [[Rcpp::export]]
Rcpp::List GroupedCommonBlissOwen(Rcpp::NumericMatrix y,int maxit=100,double tol=1e-9,double tolphi=-1.0,double posthres=0.0,int down_ndev=3,int up_ndev=2,double down_maxprop=1.0,double up_maxprop=1.0,int niter=2)
{
 int n=y.nrow();
 int p=y.ncol();
 
 if (n<2)
  Rcpp::stop("Incorrect size of input data matrix (less than 2 rows).\n");
 
 if (maxit<1)
  Rcpp::stop("Incorrect value for maximum iterations (less than 1).\n");
 
 if (tol<0.0)
  Rcpp::stop("Incorrect value for tolerance. It must be a possitive number.\n");
 
 if (DEB)
  Rcpp::Rcout << "Tolerance set to " << ((tol==1e-9) ? "its default value, " : "") << tol << "\n";
    
 if ((tolphi<0.0) && (tolphi!=-1.0))
  Rcpp::stop("Incorrect value for tolphi. It must be a possitive number, or -1 to indicate 'Estimate from the sample'.\n");
 
 if (DEB)
 {
  if (tolphi==-1.0)
   Rcpp::Rcout << "|Phi| difference will be set to an appropriate estimation, different for each row.\n";
  else
   Rcpp::Rcout << "|Phi| difference will be set to the fixed value " << tolphi << ".\n";
 }
 
 if (down_ndev<0)
  Rcpp::stop("Incorrect value for down_ndev. It must be a possitive integer.\n");
 if (up_ndev<0)
  Rcpp::stop("Incorrect value for down_ndev. It must be a possitive integer.\n"); 
 
 if (DEB)
  Rcpp::Rcout << "Interval for grouping: [value-" << down_ndev << "*stdev, value+" << up_ndev << "*stdev]\n"; 
  
 if ((down_maxprop<0.0) || (down_maxprop>1.0))
  Rcpp::stop("Incorrect value for parameter down_maxprop: must be a real number in [0.0 .. 1.0]\n");
 
 if ((up_maxprop<0.0) || (up_maxprop>1.0))
  Rcpp::stop("Incorrect value for parameter up_maxprop: must be a real number in [0.0 .. 1.0]\n");  
 
 // Vector with means of each row
 Rcpp::NumericVector m(n);
 for (int row=0; row<n; row++)
 {
  m[row]=0.0;
  for (int col=0; col<p; col++)
   m[row] += y(row,col);
  m[row] /= double(p);
 }
 
 // Vector with variances of each row
 Rcpp::NumericVector s2(n);
 // Vectors with statistics u and v
 Rcpp::NumericVector u(n);
 Rcpp::NumericVector v(n);
 // Vector with means of only good values
 Rcpp::NumericVector mgood(n);
 
 // Same loop used to fill all vectors...
 if (DEB)
  Rcpp::Rcout << "Calculating u,v, s2 and local Phi for up to " << n << " samples...";
  
 double sum2;
 double ut,vt,s2t;
 int ng=0;
 Rcpp::LogicalVector valid(n);
 Rcpp::NumericVector phi(n);
 
 for (int row=0; row<n; row++)
 {
  if (DEB && !(row%1000))
   Rcpp::Rcout << " " << row;
  if (DEB && row==n-1)
   Rcpp::Rcout << " " << row << "\n";
    
  sum2=0.0;
  for (int col=0; col<p; col++) 
   sum2 += y(row,col)*y(row,col);
   
  // Var = (sum(y^2)-n*mean(y)^2)/(n-1) 
  s2t = (sum2 - double(p)*m[row]*m[row])/double(p-1);
  // u = mean(y)^2-(s2/n)
  ut = m[row]*m[row]-(s2t/double(p));
  // v = s2-mean(y)
  vt = s2t-m[row];
 
  valid[row] = (vt>posthres);
  
  if (valid[row])
  {
   s2[ng] = s2t;
   u[ng] = ut;
   v[ng] = vt;
   mgood[ng] = m[row];
   phi[ng] = vt/ut;
   ng++;
  }
 }
 
 if (n!=ng)
  Rcpp::warning("%d samples have been discarded because statistic v was below %f\n",n-ng,posthres);
 if (DEB && ng==n)
  Rcpp::Rcout << "All samples have been accepted.\n";
 
 // Next loop is to find, for each sample, which other samples are close to it in terms of Phi estimation
 // and use them to calculate a common value for Phi
 
 // The marks to signal the group of samples which are close to the current one
 std::vector<bool> ingroup;
 ingroup.resize(ng);
 std::vector<bool> inlowgroup;
 inlowgroup.resize(ng);
 std::vector<bool> inupgroup;
 inupgroup.resize(ng);
 
 // The differences in absolute value between the samples in the current group and the current sample for
 // samples below the current one. It will be a huge number for samples not in the current group
 std::vector<double> difbelow;
 bool limit_down = (down_maxprop != 1.0);
 if (limit_down)
  difbelow.resize(ng);
 
 // The differences in absolute value between the samples in the current group and the current sample for
 // samples above the current one. It will be a huge number for samples not in the current group
 std::vector<double> difabove;
 bool limit_up = (up_maxprop != 1.0);
 if (limit_up)
  difabove.resize(ng);
 
 //std::ofstream f("res.m");
 //f << "M=[\n";
 
 // Also, the variables inf,sup,delta,nsel and dummyphi are reused in each loop
 double inf,sup;
 Rcpp::NumericVector dummyphi(ng);
 Rcpp::NumericVector delta(ng);
 Rcpp::IntegerVector nsel(ng);
 
 int nslow,nsup,nslowc,nsupc;
 int up_expurged,down_expurged;
 bool expurged;
 
 for (int l=0; l<niter; l++)
 { 
  up_expurged=0;
  down_expurged=0;
  for (int row=0; row<ng; row++)
  {
   // This is the estimator for tolphi, unless the user has chosen a fixed value.
   // It comes from Bliss-Owen, p. 38, eq. 5, but corrected.
   // It is printed as (2k(k+1)/(Nk^3))((m+k)/m)^2 which really should be
   //(2(k+1)/(Nk^3))((m+k)/m)^2 = (2/N)*(1+(1/k))*((1/k)+(1/m))^2
   // The stdev (taking square root) is sqrt((2/N)*(1+(1/k))) * ((1/k)+(1/m))
   // in which we approximately take k as 1/phi and m as mu which results in
   // sigma=(Phi+1/mu)*sqrt((2/N)(1+Phi))
   delta[row] = (tolphi==-1.0) ? ( (phi[row]+(1.0/mgood[row]))*sqrt(2.0*(1.0+phi[row])/double(p)) ) : tolphi;
   inf=phi[row]-down_ndev*delta[row];
   sup=phi[row]+up_ndev*delta[row];
   
   // No need to initialize inlowgroup and inupgroup at each loop, the loop through all values sets each to true or false
   // Also, ingroup is built later as logical or of both
 
   // On the contrary, we need to initialize difbelow and difabove, if they are to be used...
   if (limit_down)
    for (int row2=0; row2<ng; row2++)
     difbelow[row2] = DBL_MAX;
     
   if (limit_up)
    for (int row2=0; row2<ng; row2++) 
     difabove[row2] = DBL_MAX;
   
   nslow=nsup=0;
   
   for (int row2=0; row2<ng; row2++)
   {
    // If difference in Phi (which is v/u) is bigger than the lower threshold and lower than the upper threshold the sample is included in the group...
    if ( (phi[row2]>=inf) && (phi[row2]<=sup) )
    {
     // .. either to be in the lower part....
     if (phi[row2]<=phi[row])
     {
      inlowgroup[row2] = true;
      nslow++;
      // differences are stored only if we plan to use them later
      if (limit_down)
       difbelow[row2] = phi[row]-phi[row2];
     }
     // ... or in the upper part
     if (phi[row2] > phi[row])
     {
      inupgroup[row2] = true;
      nsup++;
      if (limit_up)
       difabove[row2] = phi[row2]-phi[row];
     }
    }
    else
     inlowgroup[row2] = inupgroup[row2] = false;
   }
   
   if (limit_down)
   {
    nslowc=SortAndExpurge(difbelow,inlowgroup,down_maxprop,nslow,expurged);
    if (DEB && nslowc != nslow)
     Rcpp::Rcout << "Row " << row << ", lower interval contains " << nslow << " samples, trimmed to " << nslowc << std::endl;
    nslow=nslowc;
    if (expurged)
     down_expurged++;
   } 
    
   if (limit_up)
   {
    nsupc=SortAndExpurge(difabove,inupgroup,up_maxprop,nsup,expurged);
    if (DEB && nsupc != nsup)
     Rcpp::Rcout << "Row " << row << ", upper interval contains " << nsup << " samples, trimmed to " << nsupc << std::endl;
    nsup=nsupc;
    if (expurged)
     up_expurged++;
   }
   
   for (int row2=0; row2<ng; row2++)
    ingroup[row2] = inlowgroup[row2] || inupgroup[row2];
    
   nsel[row] = nslow + nsup;
   
   // If there is only one observation close to our case, it is the case itself and its estimation is just v/u
   // Otherwise, estimate using the group of nsel selected observations with the function done for that
   //f << phi[row] << " " << delta[row] << " " << inf << " " << sup << "\n";
   dummyphi[row] = (nsel[row]==1) ? phi[row] : AuxiliaryCommonBlissOwen(ng,u,v,mgood,ingroup,nsel[row],maxit,tol);
  }
  for (int row=0; row<ng; row++)
   phi[row] = dummyphi[row];
  if (DEB && (limit_up || limit_down))
   Rcpp::Rcout << "Iteration " << l << ": ";
  if (DEB && limit_up) 
   Rcpp::Rcout << up_expurged << " rows with at least one sample trimmed in upper interval.\n";
  if (DEB && limit_down) 
   Rcpp::Rcout << down_expurged << " rows with at least one sample trimmed in lower interval.\n";  
 }
 //f << "];\n";
 //f.close();
  
 // A list of lists to be returned
 Rcpp::List ret(n);
 
 double minPhi=1e10;
 double maxPhi=0.0;
 for (int row=0,row2=0; row<n; row++)
 {
  if (valid[row])
  {
   Rcpp::List oneres=Rcpp::List::create(
                                        Rcpp::Named("Phi") = phi[row2],
                                        Rcpp::Named("valid") = true,
                                        Rcpp::Named("nsel") = nsel[row2],
                                        Rcpp::Named("tolphi") = delta[row2]);
   ret[row]=oneres;
                                        
   if (phi[row2]>maxPhi)
    maxPhi=phi[row2];
   if (phi[row2]<minPhi)
    minPhi=phi[row2]; 
   row2++;
  }
  else
  {
   Rcpp::List oneres=Rcpp::List::create(
                                        Rcpp::Named("Phi") = NAN,
                                        Rcpp::Named("valid") = false,
                                        Rcpp::Named("nsel") = 0,
                                        Rcpp::Named("tolphi") = NAN);
   ret[row]=oneres;
  }
 }
 
 if (DEB)
  Rcpp::Rcout << "Values of Phi estimated by grouped Bliss-Owen method. They belong to the interval [" << minPhi << "," << maxPhi << "]\n";

 return(ret);
}


