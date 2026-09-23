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
 * Set of functions to build a model matrix with interactions between neighbours and apply the GLM to it
 *
 *************************************************************************************/
 
#include <rowadjust.h>

extern unsigned char DEB;

// You should have added constants and names (all the name in lower case) for more distributions in
// future versions and increased the NPOSSDIST constant as needed in file rowadjust.h
// ALSO, MODIFY the array procestnames here.
// It must be such that procestnames[CONSTANT] be the corresponding name
std::string procestnames[N_PHI_EST]={"none","single","common","grouped","common_edger","tagwise_edger"};

int CheckEstprocName(std::string procname)
{ 
 // Convert to lower case for standarization
 std::transform(procname.begin(), procname.end(), procname.begin(),[](unsigned char c){ return std::tolower(c); });
 
 int i=0;
 while (i<N_PHI_EST && procname!=procestnames[i]) 
  i++;
 if (i>=N_PHI_EST)
  Rcpp::stop("%s is an incorrect Phi estimation procedure name.\n",procname);
 
 if (DEB)
  Rcpp::Rcout << "Procedure for Phi estimation set to " << procname << " (procedure number " << i << ")\n";
    
 return(i);
}

//' FitRowGLM
//' 
//' Fits a Poisson or General Linear Model to the requested row of a data matrix using as predictors
//' those included in the commonX model matrix
//'
//' Error management: This function manages the errors in the same way as StGLM, i.e.: it returns an emtpy vector for beta,
//' the value 0 for X2, sums and Dev and the empty matrix for CovM if the model is numerically unstable. This MUST be checked
//' by the caller function and it is to allow it to call a different model if the current one cannot be fitted
//'
//' @param    dat     Numeric matrix of counts with n individuals/genes as rows and c studies/cases as columns
//' @param    ng      Number of the individual/gene for which the model must be adjusted
//' @param    model   A string with the model to fit, either "poisson","binlogit","negbin" or "binlogitw"
//' @param    commonX Model matrix with c rows and with the columns common to all individuals/genes. Usually generated as model.matrix(...)
//' @param    offset  Numeric vector with offsets for the GLM. Its length must be the number of colums of the counts matrix, i.e. c. If the data are normalized this vector should be passed filled with 0.
//' @param    Phi     Estimated value for Phi in the case of negative binomial model. It can't be 0.
//' @param    retmu   Boolean to indicate whether or not the vector of estimated means must be returned. Default: FALSE
//' @return   A list with the following keys:
//' \itemize{
//' \item model - A string with the model used to adjust. A copy of the model parameter, useful to keep track of it when the result is included in a list 
//' \item beta  - A numeric vector with the value of the model parameters assuming the requested distribution
//' \item muest - A numeric vector with the estimated means or NULL if parameter retmu was FALSE
//' \item X2    - The numeric value of the X2 statistic (see help of StGLM)
//' \item sums  - The numeric value of the sums statistic (see help of StGLM)
//' \item Dev   - The numeric value of the deviance of the model
//' \item CovM  - Covariance matrix of the maximum likelihood estimator of the model coefficients
//' }
//' @examples
//' # To be done
//' @export
// [[Rcpp::export]]
Rcpp::List FitRowGLM(const Rcpp::NumericMatrix dat,int ng,std::string model,const Rcpp::NumericMatrix commonX,const Rcpp::NumericVector offset,double Phi=0.0,bool retmu=false)
{
 // Number of individuals/genes
 int n=dat.nrow();
 // Number of cases/studies
 int c=dat.ncol();
 
 if (commonX.nrow()!=c)
  Rcpp::stop("Number of columns of dat different to number of rows of commonX.\n");
 
 if (ng<1 || ng>n)
  Rcpp::stop("Number of row must be in the range of rows of dat matrix (i.e.: 1 to %d)\n",n);
   
 int nummodel=CheckDistName(model);
 if ((nummodel==DNEGBIN) && (Phi<=0.0))
  Rcpp::stop("Parameter Phi for negbinomial distribution must be passed as a strictly possitive number.\n");
 
 Rcpp::NumericVector y(c);
 // Notice that we get row ng-1 because in dat matrix we go back again to C++-index convention
 for (int q=0; q<c; q++)
  y[q]=dat(ng-1,q);
 
 Rcpp::List r=StGLM(y,commonX,model,offset,false,Phi,retmu);
 
 Rcpp::NumericVector beta=r["beta"];

 if (beta.length()==0)
  Rcpp::warning("Row %d: the base model could not be fitted.\n",ng+1);
  
 return(r);
}

//' FitExtendedRowGLM
//' 
//' Fits a Poisson or Negative Binomial General Linear Model to the requested row of a data matrix using as predictors
//' those included in the commonX model matrix plus the columns of those rows in the data matrix which are neighbours of the current row.\cr
//' The concept of neighbourhood is encoded in a nxp matrix of closeness which contains int its i-th row\cr
//' the index of the p cases which are closest to the i-th case. The dimension of the matrix is used, so p is not explicitly passed as parameter.
//'
//' @param    dat     Numeric matrix of counts with n individuals/genes as rows and c studies/cases as columns
//' @param    Cl      Numeric matrix of affinity with n individuals/genes as rows and the integer indexes of the p closest individuals in its p columns
//' @param    ng      Number of the individual/gene for which the complete model matrix has to be built
//' @param    model   A string with the model to fit, either "poisson", "binlogit", "negbin" or "binlogitw"
//' @param    commonX Partial model matrix with c rows and with the columns common to all individuals/genes. Usually generated as model.matrix(...)
//' @param    offset  Numeric vector with offsets for the GLM. Its length must be the number of colums of the counts matrix, i.e. c. If the data are normalized this vector should be passed filled with 0.
//' @param    Phi     Estimated value for Phi in the case of negative binomial model. It can't be 0.
//' @param    retmu   Boolean to indicate whether or not the vector of estimated means must be returned
//' @return   A list with the following keys:
//' \itemize{
//' \item model  - A string with the model used to adjust. A copy of the model parameter, useful to keep track of it when the result is included in a list 
//' \item beta   - A numeric vector with the value of the model parameters assuming the requested distribution
//' \item muest  - A numeric vector with the estimated means or NULL if parameter retmu was FALSE
//' \item X2     - The numeric value of the X2 statistic (see help of StGLM)
//' \item sums   - The numeric value of the sums statistic (see help of StGLM)
//' \item numnei - The number of neighbours that has been used to obtain a valid fit, either p or a lower value
//' \item Dev    - The numeric value of the deviance of the model
//' \item CovM   - Covariance matrix of the maximum likelihood estimator of the model coefficients
//' }
//' @examples
//' # To be done
//' @export
// [[Rcpp::export]]
Rcpp::List FitExtendedRowGLM(const Rcpp::NumericMatrix dat,const Rcpp::IntegerMatrix Cl,int ng,std::string model,const Rcpp::NumericMatrix commonX,const Rcpp::NumericVector offset,double Phi=0.0,bool retmu=false)
{
 // Number of individuals/genes
 int n=dat.nrow();
 // Number of cases/studies
 int c=dat.ncol();
 // Number of neighbours of each individual/gene
 int p=Cl.ncol();
 // Number of columns needed to accomodate the phenotipic variables and the offset term in a GLM (R knows about this)
 int prevc=commonX.ncol();
 
 if (Cl.nrow()!=n)
  Rcpp::stop("Number of rows of dat and Cl parameters are not the same.\n");
 
 if (commonX.nrow()!=c)
  Rcpp::stop("Number of columns of dat different to number of rows of commonX.\n");
 
 if (ng<1 || ng>n)
  Rcpp::stop("Number of row must be in the range of rows of dat matrix (i.e.: 1 to %d)\n",n);
 
 int nummodel=CheckDistName(model);
 if ((nummodel==DNEGBIN) && (Phi<=0.0))
  Rcpp::stop("Parameter Phi for negbinomial distribution must be passed as a strictly possitive number.\n");
 
 int num_used_variables=prevc+p;
 // The new model matrix has c rows (as many as columns in the data matrix) and as many columns as there are in
 // the common model, plus the variables we decide to add
 Rcpp::NumericMatrix X(c,num_used_variables);
 
 // The common part is copied (copied, not referenced)
 for (int col=0; col<prevc; col++)
  for (int row=0; row<c; row++)
   X(row,col)=commonX(row,col);
 
 // Now, the appropriate rows of the dat matrix are copied as columns of the model matrix
 // Notice the use of ng-1 instead of ng: ng has come from R as a R-index type (from 1). dat is indexed in C/Rcpp mode (from 0).
 // For the same reason, the VALUE of the Cl matrix, which is a R-like index, too, must be diminished by 1.
 for (int col=prevc; col<num_used_variables; col++)
  for (int row=0; row<c; row++)
   X(row,col)=dat(Cl(ng-1,col-prevc)-1,row);
 
 Rcpp::NumericVector y(c);
 // Notice that we get row ng-1 because in dat matrix we go back again to C++-index convention
 for (int q=0; q<c; q++)
  y[q]=dat(ng-1,q);
 
 bool fitted=false;
 Rcpp::NumericVector beta;
 double X2=0.0;
 double sums=0.0;
 double Dev=0.0;
 Rcpp::NumericMatrix CovM;
 Rcpp::NumericVector muest;
 int iters=0;
 do
 {   
  // This extracts from X all the rows, and some columns, one less each time, as long as the model cannot be fitted.
  Rcpp::NumericMatrix::Sub Xa=X(Rcpp::Range(0,c-1),Rcpp::Range(0,num_used_variables-1));

  Rcpp::List stret=StGLM(y,Xa,model,offset,false,Phi,retmu);
 
  beta=stret["beta"];
  // If the model is numerically unstable, we try to reduce the numer of added variables
  if (beta.length()==0)
   num_used_variables--;
  else
  {
   X2=stret["X2"];
   muest=stret["muest"];
   sums=stret["sums"];
   Dev=stret["Dev"];
   iters=stret["iters"];
   CovM=Rcpp::as<Rcpp::NumericMatrix>(stret["CovM"]);
   fitted=true;
   // Inside a thread DEB will be false, even if originally set to true
   if (DEB)
    Rcpp::Rcout << "adjusted with " << num_used_variables-prevc << " neighbours.\n";
  }
 }
 while (!fitted && num_used_variables>prevc);
 
 if (DEB && !fitted)
  Rcpp::Rcout << "could not be fitted even with one neighbour.\n";
 if (!fitted)
  Rcpp::warning("Row %d: the extended model could not be fitted even with one neighbour.\n",ng+1);
  
 // Notice that, if we have left the loop because num_used_variables==prevc, we would be in the base model
 // but we haven't adjusted it. beta would have the value returned by the call with one neighbour, i.e. a
 // vector of length 1, and X2 and sums will not have been changed, i.e. they will be both zero.
 // This will be corrected by the caller function when noticing that numnei=0. We could correct it here by
 // trying the fit with the base model but that was already done in the caller and it succedeed. Otherwise,
 // this function would have never been called.
 Rcpp::List ret=Rcpp::List::create(
 		Rcpp::Named("model")=model,
        Rcpp::Named("beta")=beta,
        Rcpp::Named("muest")=muest,
        Rcpp::Named("X2")=X2,
        Rcpp::Named("sums")=sums,
        Rcpp::Named("numnei")=num_used_variables-prevc,
        Rcpp::Named("Dev")=Dev,
        Rcpp::Named("CovM")=CovM,
        Rcpp::Named("iters")=iters);
 
 return(ret);
}

//' FitGLMWithNei
//'
//' Fits a General Linear Model to the requested rows of a data matrix by using as predictors those
//' included in the commonX model matrix and also those plus the columns of some rows in the data matrix which are neighbours of the current row.
//'
//' The concept of neighbourhood between rows is encoded in a nxp matrix of closeness which contains in its i-th row
//' the index of the p cases which are closest to the i-th case. The function tries to adjust two models, one
//' with and one without the added neighbours and returns its parameters and measures of goodnes-of-fit.\cr
//' It might be possible that no model can be adjusted if too many predictors are incorporated in the model matrix,
//' or even without incorporating any additional predictors. This manifests as the lack of convergence of the numerical
//' algorithms. The parameter nnei returns the maximum number of neighbours (from p and downward) that allows a stable
//' adjustment. Value 0 would indicate that only the common model could be adjusted and in this case beta1 is NULL and
//' the Wedderburn and related parameters are those of the base model.
//' Value -1 is used to indicate that not even the common model fitted. In this case values of beta0 and beta1
//' are the empty vectors and overdispersion estimators are 0.
//'
//' @param    dat     Numeric matrix of counts with n individuals/genes as rows and c studies/cases as columns
//' @param    normalized Boolean value to indicate if data are passed normalized (TRUE) or not (FALSE). If not, offset  vector is calculated inside.
//' @param    Cl      Numeric matrix of affinity with n individuals/genes as rows and the integer indexes of the p closest individuals in its p columns
//' @param    selected Boolean array of n positions indicating to which of the n rows the selected model should be fitted. Notice that calculation of
//'                    the offset vector, if requested by the normalized parameter, is done using all the individuals in the dat matrix, not only the selected ones.
//'                    Similarly, Cl must have the number of the closest neighbours, either if they are in the selected ones or not.
//' @param    model   A string with the model to fit, either "poisson","binlogit","negbin" or "binlogitw"
//' @param    commonX Partial model matrix with c rows and with the columns common to all individuals/genes. Usually generated as model.matrix(...)
//' @param    Phi_estproc Procedure to estimate the Phi parameter in the case of negbinomial model. Possible values are "none","single","common","grouped","common_edger" or "tagwise_edger".
//' @param    Phi_external Numeric vector with estimated Phi values if the estimation procedure is "common_edger" or "tagwise_edger". If a vector with correct length is passed, but estimation procedures has been set to any different to the former ones, the vector is used and the procedure is set to "common_edger". Default: empty numeric vector. 
//' @param    retmu  Boolean to indicat wheter or not the vectors of estimated means must be returned for both models. Default: FALSE
//' @return   A list of lists, each one with the following keys:
//' \itemize{
//' \item model - A string with the model used to adjust. A copy of the model parameter, useful to keep track of it when the result is included in a list, or none if selected was FALSE for that row
//' \item beta0 - A numeric vector with the value of the model parameters assuming a Poisson distribution and commonX as model matrix
//' \item Dev0  - A numeric value with the deviance of the base model
//' \item CovM0 - A numeric matrix, the covariance matrix of the maximum likelihood estimator of the model coefficients for the base model
//' \item muest0 - A numeric vector with the estimated means of the base model, or NULL if parameter retmu was FALSE
//' \item iters0 - An integer with the number of iterations until convergence of the base model
//' \item beta1 - A numeric vector with the value of the model parameters assuming a Poisson distribution and commonX plus as much as p added neighbours as model matrix.
//' \item Dev1  - A numeric value with the deviance of the extended model (commonX plus as much as p added neighbours)
//' \item CovM11 - A numeric matrix, the covariance matrix of the maximum likelihood estimator of the model coefficients for the extended model
//' \item muest1 - A numeric vector with the estimated means of the extended model, or NULL if parameter retmu was FALSE
//' \item iters1 - An integer with the number of iterations until convergence of the extended model
//' \item Bliss - The numeric value of Phi obtained using the requested Bliss-Owen method (single, common or grouped) for fits to the negative binomial model or 0.0 otherwise
//' \item Wedderburn - The numeric value of the cuasidispersion parameter of the adjusted model according to the Wedderburn estimator (1979)
//' \item Farrington - The numeric value of the cuasidispersion parameter of the adjusted model according to the Farrington estimator (1996)
//' \item Fletcher   - The numeric value of the cuasidispersion parameter of the adjusted model according to the Fletcher estimator (2012)
//' \item nnei  - An integer with the number of neigbours that still generate a valid adjusted model
//' }
//' 
//' @examples
//' # To be done
//' @export
// [[Rcpp::export]]
Rcpp::List FitGLMWithNei(const Rcpp::NumericMatrix dat,bool normalized,const Rcpp::IntegerMatrix Cl,std::string model,const Rcpp::LogicalVector selected,const Rcpp::NumericMatrix commonX,std::string Phi_estproc,Rcpp::NumericVector Phi_external=Rcpp::NumericVector::create(),bool retmu=false)
{
 // Number of columns and rows of the data
 int ncolsdat=dat.ncol();
 int nrowsdat=dat.nrow();
 // Number of columns of the initial model
 int ncols_basemod=commonX.ncol();
 
 int nummodel=CheckDistName(model);
 
 int numestproc=CheckEstprocName(Phi_estproc);
 
 if ((nummodel==DPOISS) && (numestproc!=PHI_EST_NOT_NEEDED))
 {
  Rcpp::warning("Poisson model does not admit a Phi estimation procedure, since it has no Phi parameter. Estimation procedure set to 'none'\n");
  numestproc=PHI_EST_NOT_NEEDED;
 }
 
 if (nrowsdat!=selected.length())
  Rcpp::stop("Boolean vector of selected cases has not the same length as number of rows of data matrix (%d instead of %d)\n",selected.length(),nrowsdat);
  
 int nsel=0;
 for (int row=0; row<nrowsdat; row++)
  if (selected[row])
   nsel++;
   
 if (Phi_external.length()==nrowsdat && numestproc!=PHI_EST_COMMON_EDGER && numestproc!=PHI_EST_TAGWISE_EDGER)
 {
  Rcpp::warning("You have passed a vector of externally estimated values for Phi but you have set the estimation procedure to %s.\nThe only estimation procedures where you must pass a vector of externally estimated Phi values are %s and %s.\nFrom now on, we will set the procedure to %s.\n",Phi_estproc,procestnames[PHI_EST_COMMON_EDGER],procestnames[PHI_EST_TAGWISE_EDGER],procestnames[PHI_EST_COMMON_EDGER]);
  numestproc=PHI_EST_COMMON_EDGER;
 }
 
 if (DEB)
 {
  Rcpp::Rcout << "Fitting model number " << nummodel << " (" << model << ") to " << nsel << " cases.\n";
  Rcpp::Rcout << "Filling vector of offsets...  ";
 }
 
 // offset vector is the sum of the selected data, by columns
 Rcpp::NumericVector offset(ncolsdat);
 for (int col=0; col<ncolsdat; col++)
 {
  offset[col]=0.0;
  if (normalized==false)
   for (int row=0; row<nrowsdat; row++)
    if (selected[row])
     offset[col]+=dat(row,col);
 }
 
 // The list to be returned has as many rows as dat, even if they have not been selected. Those not selected will return with empty values in ret[row]
 Rcpp::List ret(nrowsdat);
 
 DifftimeHelper Dt;
 
 Dt.StartClock("End of GLM application (serial version)."); 
 double Dev0,Dev1;
 int iters0,iters1;
 Rcpp::NumericMatrix CovM0;
 Rcpp::NumericMatrix CovM1;

 // Declared as different because internally they have a different structure. 
 // Phi_single_estimation and Phi_common_estimation are returned as a list with two fields: valid and Phi whereas
 // Phi_estimations is returned as list of lists of two fields (valid and Phi again).
 Rcpp::List Phi_single_estimation;
 Rcpp::List Phi_common_estimation;
 Rcpp::List Phi_estimations;
 
 // This is the vector of Phi estimations. Values will be:
 // -1 if estimation was not needed (row not selected or we are in Poisson model),
 //  0 if estimation has failed (BlissOwen does not work) 
 //  the estimated value otherwise
 Rcpp::NumericVector Phi(nrowsdat);
 
 switch (numestproc)
 {
  case PHI_EST_NOT_NEEDED:
   for (int row=0; row<nrowsdat; row++)
    Phi[row]=-1.0;
  break;
  case PHI_EST_SINGLE:
  {
   int nbad=0;
   for (int row=0; row<nrowsdat; row++)
   {
    if (selected[row])
    {
     Phi_single_estimation=BlissOwen(dat.row(row));
     if (Phi_single_estimation["valid"])
      Phi[row]=Phi_single_estimation["Phi"];
     else
     {
      Phi[row]=0.0;
      nbad++;
     }
    }
    else
     Phi[row]=-1.0;
   }
   if (((float)nbad/(float)nsel)>=0.5)
    Rcpp::stop("BlissOwen single estimation of Phi failed in more than half of the samples (%d out of %d).\n",nbad,nsel);
  }
  break;
  case PHI_EST_COMMON:
  {
   Rcpp::NumericMatrix datsel(nsel,ncolsdat);
   int s=0;
   for (int row=0; row<nrowsdat; row++)
    if (selected[row])
    {
     for (int col=0; col<ncolsdat; col++)
      datsel(s,col)=dat(row,col);
     s++;
    }
   Phi_common_estimation=CommonBlissOwen(datsel);
   if (!Phi_common_estimation["valid"])
   {
    Rcpp::stop("Phi cannot be estimated with the Bliss-Owen method for the selected sample (estimation procedure: common).\n");
    //for (int row=0; row<nrowsdat; row++)
    // Phi[row]=0.0;
   }
   else
    for (int row=0; row<nrowsdat; row++)
     Phi[row] = (selected[row]) ? Phi_common_estimation["Phi"] : -1.0;
  }
  break;
  case PHI_EST_GROUPED:
  {
   Rcpp::NumericMatrix datsel(nsel,ncolsdat);
   int s=0;
   for (int row=0; row<nrowsdat; row++)
    if (selected[row])
    {
     for (int col=0; col<ncolsdat; col++)
      datsel(s,col)=dat(row,col);
     s++;
    }
   Phi_estimations=GroupedCommonBlissOwen(datsel);
   int nbad=0;
   s=0;
   for (int row=0; row<nrowsdat; row++)
   {
    if (selected[row])
    {
     Rcpp::List r=Phi_estimations[s];
     if (r["valid"])
      Phi[row]=r["Phi"];
     else
     {
      nbad++;
      Phi[row]=0.0;
     }
     s++;
    }
    else
     Phi[row]=-1.0;
   }
   if (((float)nbad/(float)nsel)>=0.5)
    Rcpp::stop("BlissOwen grouped estimation of Phi failed in more than half of the samples (%d out of %d).\n",nbad,nsel);
  }
  break;
  case PHI_EST_COMMON_EDGER:
  case PHI_EST_TAGWISE_EDGER:
  {
   if (Phi_external.length()!=nrowsdat)
    Rcpp::stop("External estimation of Phi comes with a vector that has a wrong number of components (%d instead of %d).\n",Phi_external.length(),nrowsdat);
     
   for (int row=0; row<nrowsdat; row++)
   {
    if (selected[row])
     Phi[row]=Phi_external[row];
    else
     Phi[row]=-1.0;
   }
  }
  break;
  default: 
   Rcpp::stop("Wrong Phi estimation method passed to switch in function FitGLMWithNei. We should not be here...\n");
   break;
 }
 
 for (int row=0; row<nrowsdat; row++)
 {
  // Print the row in R-convention (from 1)
  if (DEB)
   Rcpp::Rcout << "Row " << row+1 << ": ";
  

  if (!selected[row] || Phi[row]==0.0 )
  {
   if (DEB)
   {
    if (!selected[row])
     Rcpp::Rcout << "not selected.\n";
    else
     Rcpp::Rcout << "Phi not estimated.\n";
   }
   
   Rcpp::List L2=Rcpp::List::create(
  		Rcpp::Named("model")="model",
        Rcpp::Named("beta0")=Rcpp::NumericVector(),
        Rcpp::Named("Dev0")=NULL,
        Rcpp::Named("CovM")=Rcpp::NumericMatrix(),
        Rcpp::Named("muest0")=NULL,
        Rcpp::Named("iters0")=NULL,
        Rcpp::Named("beta1")=Rcpp::NumericVector(),
        Rcpp::Named("Dev1")=NULL,
        Rcpp::Named("CovM")=Rcpp::NumericMatrix(),
        Rcpp::Named("muest1")=NULL,
        Rcpp::Named("iters1")=NULL,
        Rcpp::Named("nnei")=-1,
        Rcpp::Named("Bliss")=NULL,
        Rcpp::Named("Wedderburn")=NULL,
        Rcpp::Named("Farrington")=NULL,
        Rcpp::Named("Fletcher")=NULL);
   ret[row]=L2;
  }
  else
  {
   int numnei;
   
   // Here row is passed in R-convention (from 1).
   // This is because FitRowGLM is a function with R interface and could be called also from R
   // L0 contains the result of fitting the distribution to the base model, i.e. without additional neighbours but
   // with the independent term (commonX was returned by the R function model.matrix)
   Rcpp::List L0=FitRowGLM(dat,row+1,model,commonX,offset,Phi[row],retmu); 
   Rcpp::NumericVector beta0=L0["beta"];
   Dev0=L0["Dev"];
   CovM0=Rcpp::as<Rcpp::NumericMatrix>(L0["CovM"]);
   iters0=L0["iters"];
   Rcpp::NumericVector muest0=L0["muest"];
 
   if (beta0.length()==0)
   {
    numnei=-1;
    Rcpp::List L2=Rcpp::List::create(
  		Rcpp::Named("model")="model",
        Rcpp::Named("beta0")=Rcpp::NumericVector(),
        Rcpp::Named("Dev0")=NULL,
        Rcpp::Named("CovM0")=Rcpp::NumericMatrix(),
        Rcpp::Named("muest0")=NULL,
        Rcpp::Named("iters0")=NULL,
        Rcpp::Named("beta1")=Rcpp::NumericVector(),
        Rcpp::Named("Dev1")=NULL,
        Rcpp::Named("CovM1")=Rcpp::NumericMatrix(),
        Rcpp::Named("muest1")=NULL,
        Rcpp::Named("iters1")=NULL,
        Rcpp::Named("nnei")=numnei,
        Rcpp::Named("Bliss")=NULL,
        Rcpp::Named("Wedderburn")=NULL,
        Rcpp::Named("Farrington")=NULL,
        Rcpp::Named("Fletcher")=NULL);
    ret[row]=L2;
   }
   else
   { 
    // This function tries to fit the requested model with all neighbours in the Cl matrix, which are as many as columns in it
    // If the fit fails, it tries with one less, and so on, until 1. If this fails, too, numnei will return as 0
    // indicating that only the base model fitted (which is true, otherwise we would not be here...)
    Rcpp::List L1 = FitExtendedRowGLM(dat,Cl,row+1,model,commonX,offset,Phi[row],retmu);
     
    numnei=L1["numnei"];
    
    double X2;
    double sums;
    Rcpp::NumericVector beta1;
    Rcpp::NumericVector muest1;
  
    if (numnei != 0)
    {
     beta1 = L1["beta"];
     X2 = L1["X2"];
     sums = L1["sums"];
     Dev1 = L1["Dev"];
     CovM1 = Rcpp::as<Rcpp::NumericMatrix>(L1["CovM"]);
     iters1 = L1["iters"];
     muest1 = L1["muest"];
    }
    else
    {
     // Special case: the base model fitted but not even one neighbour can be used.
     Rcpp::warning("Row %d: the base model fitted, but not even the model with one neighbour could fit. Data for the extended model will be NULL.\n",row+1);
     
     beta1 = Rcpp::NumericVector();
     CovM1 = Rcpp::NumericMatrix();
     Dev1 = 0.0;
     iters1 = 0;
     muest1 = Rcpp::NumericVector();
     X2 = L0["X2"];
     sums = L0["sums"];
    }
       
    // The following formulae are taken from Fletcher (2012)
   
    // The number of degrees of freedom, denoted in Fletcher as p, is the number of columns in the base
    // model matrix plus the number of added columns due to addition of some neighbours.
      
    // Wedderburn: 
    // Denoted in Fletcher (2021), page 231 line 3 as \hat{Phi}_P where:
    // P is our X2
    // n is our ncolsdat
    // p is in our case would be ncols_in_base_model+numnei that, if only base model could be fitted,
    // would be base_model since numnei would be 0
    double Wed = X2/double(ncolsdat-(ncols_basemod+numnei));
   
    // Farrington: it is equal to the Wedderburn estimator, minus n*overline{s}/(n-num_deg_freedom)
    // Denoted in Fletcher (2021), page 231, eq. 1 as \hat{Phi}_F where:
    // n*overline{s} is our ncolsdat*average(sums)=ncolsdat*sums/ncolsdat=sums
    // n-p is in our case as in Wedderburn
    // WARNING: this could be negative, which makes no sense, so we clip it to zero
    double Farr = Wed-(sums/double(ncolsdat-(ncols_basemod+numnei)));
    if (Farr<0.0)
     Farr=0.0;
   
    // Fletcher: it is equal to the Wedderburn estimator divided by 1 plus overline{s}
    // denoted in Fletcher (2021), page 232, eq. 5 as \hat{Phi} where:
    // overline{s} is our sums/ncolsdat
    double Fle = Wed/(1.0+(sums/double(ncolsdat)));
    
    // The fields in the returned list are for the biggest model that fitted. They are
    // the coefficient vector beta, the two statistics X2 and sums, defined as
    // X2=\sum_{i=1}^{ncolsdat} (y_i-\hat{mu}_i)^2/\hat{V}_i
    // sums=\sum_{i_1}^{ncolsdat} (\frac{\hat{V}^'_i}{\hat{V}_i}(y_i-\hat{mu}_i)
    // the deviance and the matrix CovM
     
    Rcpp::List L2=Rcpp::List::create(
  		Rcpp::Named("model")=model,
        Rcpp::Named("beta0")=beta0,
        Rcpp::Named("Dev0")=Dev0,
        Rcpp::Named("CovM0")=CovM0,
        Rcpp::Named("muest0")=muest0,
        Rcpp::Named("iters0")=iters0,
        Rcpp::Named("beta1")=beta1,
        Rcpp::Named("Dev1")=Dev1,
        Rcpp::Named("CovM1")=CovM1,
        Rcpp::Named("muest1")=muest1,
        Rcpp::Named("iters1")=iters1,
        Rcpp::Named("nnei")=numnei,
        Rcpp::Named("Bliss")=Phi[row],
        Rcpp::Named("Wedderburn")=Wed,
        Rcpp::Named("Farrington")=Farr,
        Rcpp::Named("Fletcher")=Fle);
        
    ret[row]=L2;
   }
  }
 }
 
 Dt.EndClock(DEB);
 return(ret);
}

//' FitGLMInd
//'
//' Fits a General Linear Model to the requested rows of a data matrix by using as predictors those
//' included in the commonX model matrix. This is a particular case of FitGLMWithNei when no neighbours are to be considered.
//'
//' It might be possible that the model cannot be adjusted. This manifests as the lack of convergence of the numerical
//' algorithms. In such a case the vector of coefficients beta0 is returned as a 0-length vector.
//'
//' @param    dat     Numeric matrix of counts with n individuals/genes as rows and c studies/cases as columns
//' @param    normalized Boolean value to indicate if data are passed normalized (TRUE) or not (FALSE). If not, offset  vector is calculated inside.
//' @param    model   A string with the model to fit, either "poisson","binlogit","negbin" or "binlogitw"
//' @param    selected Boolean array of n positions indicating to which of the n rows the selected model should be fitted. Notice that calculation of
//'                    the offset vector, if requested by the normalized parameter, is done using all the individuals in the dat matrix, not only the selected ones.
//' @param    commonX Partial model matrix with c rows and with the columns common to all individuals/genes. Usually generated as model.matrix(...)
//' @param    Phi_estproc Procedure to estimate the Phi parameter in the case of negbinomial model. Possible values are "none","single","common","grouped","common_edger" or "tagwise_edger".
//' @param    Phi_external Numeric vector with estimated Phi values if the estimation procedure is "common_edger" or "tagwise_edger". If a vector with correct length is passed, but estimation procedures has been set to any different to the former ones, the vector is used and the procedure is set to "common_edger". Default: empty numeric vector. 
//' @param    retmu Boolean to indicate wheter or not the vector of estimated means must be returned. Default: FALSE
//' @return   A list of lists, each one with the following keys:
//' \itemize{
//' \item model - A string with the model used to adjust. A copy of the model parameter, useful to keep track of it when the result is included in a list, or none if selected was FALSE for that row
//' \item beta0 - A numeric vector with the value of the model parameters
//' \item Dev0  - A numeric value with the deviance of the base model
//' \item CovM0 - A numeric matrix, the covariance matrix of the maximum likelihood estimator of the model coefficients for the base model
//' \item muest0 - A numeric vector with the estimated means, or NULL if parameter retmu was FALSE
//' \item iters0 - An integer with the number of iterations until convergence of the base model
//' \item Bliss - The numeric value of Phi obtained using the requested Bliss-Owen method (single, common or grouped) for fits to the negative binomial model or 0.0 otherwise
//' \item Wedderburn - The numeric value of the cuasidispersion parameter of the adjusted model according to the Wedderburn estimator (1979)
//' \item Farrington - The numeric value of the cuasidispersion parameter of the adjusted model according to the Farrington estimator (1996)
//' \item Fletcher   - The numeric value of the cuasidispersion parameter of the adjusted model according to the Fletcher estimator (2012)
//' }
//' 
//' @examples
//' # To be done
//' @export
// [[Rcpp::export]]
Rcpp::List FitGLMInd(const Rcpp::NumericMatrix dat,bool normalized,std::string model,const Rcpp::LogicalVector selected,const Rcpp::NumericMatrix commonX,std::string Phi_estproc,Rcpp::NumericVector Phi_external=Rcpp::NumericVector::create(),bool retmu=false)
{
 // Number of columns and rows of the data
 int ncolsdat=dat.ncol();
 int nrowsdat=dat.nrow();
 // Number of columns of the initial model
 int ncols_basemod=commonX.ncol();
 
 int nummodel=CheckDistName(model);
 
 int numestproc=CheckEstprocName(Phi_estproc);
 
 if ((nummodel==DPOISS) && (numestproc!=PHI_EST_NOT_NEEDED))
 {
  Rcpp::warning("Poisson model does not admit a Phi estimation procedure, since it has no Phi parameter. Estimation procedure set to 'none'\n");
  numestproc=PHI_EST_NOT_NEEDED;
 }
 
 if (nrowsdat!=selected.length())
  Rcpp::stop("Boolean vector of selected cases has not the same length as number of rows of data matrix (%d instead of %d)\n",selected.length(),nrowsdat);
  
 int nsel=0;
 for (int row=0; row<nrowsdat; row++)
  if (selected[row])
   nsel++;
   
 if (Phi_external.length()==nrowsdat && numestproc!=PHI_EST_COMMON_EDGER && numestproc!=PHI_EST_TAGWISE_EDGER)
 {
  Rcpp::warning("You have passed a vector of externally estimated values for Phi but you have set the estimation procedure to %s.\nThe only estimation procedures where you must pass a vector of externally estimated Phi values are %s and %s.\nFrom now on, we will set the procedure to %s.\n",Phi_estproc,procestnames[PHI_EST_COMMON_EDGER],procestnames[PHI_EST_TAGWISE_EDGER],procestnames[PHI_EST_COMMON_EDGER]);
  numestproc=PHI_EST_COMMON_EDGER;
 }
 
 if (DEB)
 {
  Rcpp::Rcout << "Fitting model number " << nummodel << " (" << model << ") to " << nsel << " cases.\n";
  Rcpp::Rcout << "Filling vector of offsets...  ";
 }
 
 // offset vector is the sum of the selected data, by columns
 Rcpp::NumericVector offset(ncolsdat);
 for (int col=0; col<ncolsdat; col++)
 {
  offset[col]=0.0;
  if (normalized==false)
   for (int row=0; row<nrowsdat; row++)
    if (selected[row])
     offset[col]+=dat(row,col);
 }
 
 // The list to be returned has as many rows as dat, even if they have not been selected. Those not selected will return with empty values in ret[row]
 Rcpp::List ret(nrowsdat);
 
 DifftimeHelper Dt;
 
 Dt.StartClock("End of GLM application (serial version)."); 
 double Dev0;
 int iters0;
 Rcpp::NumericMatrix CovM0;
 Rcpp::NumericMatrix CovM1;

 // Declared as different because internally they have a different structure. 
 // Phi_single_estimation and Phi_common_estimation are returned as a list with two fields: valid and Phi whereas
 // Phi_estimations is returned as list of lists of two fields (valid and Phi again).
 Rcpp::List Phi_single_estimation;
 Rcpp::List Phi_common_estimation;
 Rcpp::List Phi_estimations;
 
 // This is the vector of Phi estimations. Values will be:
 // -1 if estimation was not needed (row not selected or we are in Poisson model),
 //  0 if estimation has failed (BlissOwen does not work) 
 //  the estimated value otherwise
 Rcpp::NumericVector Phi(nrowsdat);
 
 switch (numestproc)
 {
  case PHI_EST_NOT_NEEDED:
   for (int row=0; row<nrowsdat; row++)
    Phi[row]=-1.0;
  break;
  case PHI_EST_SINGLE:
  {
   int nbad=0;
   for (int row=0; row<nrowsdat; row++)
   {
    if (selected[row])
    {
     Phi_single_estimation=BlissOwen(dat.row(row));
     if (Phi_single_estimation["valid"])
      Phi[row]=Phi_single_estimation["Phi"];
     else
     {
      Phi[row]=0.0;
      nbad++;
     }
    }
    else
     Phi[row]=-1.0;
   }
   if (((float)nbad/(float)nsel)>=0.5)
    Rcpp::stop("BlissOwen single estimation of Phi failed in more than half of the samples (%d out of %d).\n",nbad,nsel);
  }
  break;
  case PHI_EST_COMMON:
  {
   Rcpp::NumericMatrix datsel(nsel,ncolsdat);
   int s=0;
   for (int row=0; row<nrowsdat; row++)
    if (selected[row])
    {
     for (int col=0; col<ncolsdat; col++)
      datsel(s,col)=dat(row,col);
     s++;
    }
   Phi_common_estimation=CommonBlissOwen(datsel);
   if (!Phi_common_estimation["valid"])
   {
    Rcpp::stop("Phi cannot be estimated with the Bliss-Owen method for the selected sample (estimation procedure: common).\n");
    //for (int row=0; row<nrowsdat; row++)
    // Phi[row]=0.0;
   }
   else
    for (int row=0; row<nrowsdat; row++)
     Phi[row] = (selected[row]) ? Phi_common_estimation["Phi"] : -1.0;
  }
  break;
  case PHI_EST_GROUPED:
  {
   Rcpp::NumericMatrix datsel(nsel,ncolsdat);
   int s=0;
   for (int row=0; row<nrowsdat; row++)
    if (selected[row])
    {
     for (int col=0; col<ncolsdat; col++)
      datsel(s,col)=dat(row,col);
     s++;
    }
   Phi_estimations=GroupedCommonBlissOwen(datsel);
   int nbad=0;
   s=0;
   for (int row=0; row<nrowsdat; row++)
   {
    if (selected[row])
    {
     Rcpp::List r=Phi_estimations[s];
     if (r["valid"])
      Phi[row]=r["Phi"];
     else
     {
      nbad++;
      Phi[row]=0.0;
     }
     s++;
    }
    else
     Phi[row]=-1.0;
   }
   if (((float)nbad/(float)nsel)>=0.5)
    Rcpp::stop("BlissOwen grouped estimation of Phi failed in more than half of the samples (%d out of %d).\n",nbad,nsel);
  }
  break;
  case PHI_EST_COMMON_EDGER:
  case PHI_EST_TAGWISE_EDGER:
  {
   if (Phi_external.length()!=nrowsdat)
    Rcpp::stop("External estimation of Phi comes with a vector that has a wrong number of components (%d instead of %d).\n",Phi_external.length(),nrowsdat);
     
   for (int row=0; row<nrowsdat; row++)
   {
    if (selected[row])
     Phi[row]=Phi_external[row];
    else
     Phi[row]=-1.0;
   }
  }
  break;
  default: break;
 }
 
 for (int row=0; row<nrowsdat; row++)
 {
  // Print the row in R-convention (from 1)
  if (DEB)
   Rcpp::Rcout << "Row " << row+1 << ": ";
  
  if (!selected[row] || Phi[row]==0.0)
  {
   if (DEB)
   {
    if (!selected[row])
     Rcpp::Rcout << "not selected.\n";
    else
     Rcpp::Rcout << "Phi not estimated.\n";
   }
   Rcpp::List L2=Rcpp::List::create(
  		Rcpp::Named("model")="model",
        Rcpp::Named("beta0")=NULL,
        Rcpp::Named("Dev0")=NULL,
        Rcpp::Named("CovM0")=NULL,
        Rcpp::Named("muest0")=NULL,
        Rcpp::Named("iters0")=NULL,
        Rcpp::Named("Bliss")=NULL,
        Rcpp::Named("Wedderburn")=NULL,
        Rcpp::Named("Farrington")=NULL,
        Rcpp::Named("Fletcher")=NULL);
   ret[row]=L2;
  }
  else
  {
   // Here row is passed in R-convention (from 1).
   // This is because FitRowGLM is a function with R interface and could be called also from R
   // L0 contains the result of fitting the distribution to the base model, i.e. without additional neighbours but
   // with the independent term (commonX was returned by the R function model.matrix)
   Rcpp::List L0=FitRowGLM(dat,row+1,model,commonX,offset,Phi[row],retmu);
   
   Rcpp::NumericVector beta0=L0["beta"];
   if (beta0.length()==0)
   {
    Rcpp::List L2=Rcpp::List::create(
  		Rcpp::Named("model")="model",
        Rcpp::Named("beta0")=NULL,
        Rcpp::Named("Dev0")=NULL,
        Rcpp::Named("CovM0")=NULL,
        Rcpp::Named("muest0")=NULL,
        Rcpp::Named("iters0")=NULL,
        Rcpp::Named("Bliss")=NULL,
        Rcpp::Named("Wedderburn")=NULL,
        Rcpp::Named("Farrington")=NULL,
        Rcpp::Named("Fletcher")=NULL);
    ret[row]=L2;
   }
   else
   {
    Dev0=L0["Dev"];
    CovM0=Rcpp::as<Rcpp::NumericMatrix>(L0["CovM"]);
    Rcpp::NumericVector muest0=L0["muest"];
    iters0=L0["iters"];
  
    double X2=L0["X2"];
    double sums=L0["sums"];
  
    // The following formulae are taken from Fletcher (2012)
   
    // The number of degrees of freedom, denoted in Fletcher as p, is the number of columns in the base
    // model matrix plus the number of added columns due to addition of some neighbours.
      
    // Wedderburn: 
    // Denoted in Fletcher (2021), page 231 line 3 as \hat{Phi}_P where:
    // P is our X2
    // n is our ncolsdat
    // p is in our case would be ncols_in_base_model+numnei
    double Wed = X2/double(ncolsdat-(ncols_basemod));
   
    // Farrington: it is equal to the Wedderburn estimator, minus n*overline{s}/(n-num_deg_freedom)
    // Denoted in Fletcher (2021), page 231, eq. 1 as \hat{Phi}_F where:
    // n*overline{s} is our ncolsdat*average(sums)=ncolsdat*sums/ncolsdat=sums
    // n-p is in our case as in Wedderburn
    // WARNING: this could be negative, which makes no sense, so we clip it to zero
    double Farr = Wed-(sums/double(ncolsdat-(ncols_basemod)));
    if (Farr<0.0)
     Farr=0.0;
   
    // Fletcher: it is equal to the Wedderburn estimator divided by 1 plus overline{s}
    // denoted in Fletcher (2021), page 232, eq. 5 as \hat{Phi} where:
    // overline{s} is our sums/ncolsdat
    double Fle = Wed/(1.0+(sums/double(ncolsdat)));
   
    Rcpp::List L2=Rcpp::List::create(
  		Rcpp::Named("model")=model,
        Rcpp::Named("beta0")=beta0,
        Rcpp::Named("Dev0")=Dev0,
        Rcpp::Named("CovM0")=CovM0,
        Rcpp::Named("muest0")=muest0,
        Rcpp::Named("iters0")=iters0,
        Rcpp::Named("Bliss")=Phi[row],
        Rcpp::Named("Wedderburn")=Wed,
        Rcpp::Named("Farrington")=Farr,
        Rcpp::Named("Fletcher")=Fle);
        
    ret[row]=L2;     
   }
  }
 }
 
 Dt.EndClock(DEB);
 return(ret);
}

