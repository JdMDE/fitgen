/*
 *
 * Copyright (C) 2024 Juan Domingo, Guillemo Ayala, Maite Leon ({Juan.Domingo,Guillermo.Ayala, Teresa.Leon}@uv.es)
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
 
#include "stglm.h"

// You should have added constants and names (all the name in lower case) for more distributions in
// future versions and increased the NPOSSDIST constant as needed in file stglm.h
// ALSO, MODIFY the array distnames here.
// It must be such that distnames[CONSTANT] be the corresponding name
std::string distnames[NPOSSDIST]={"poisson","binlogit","negbin","binlogitw"};


extern unsigned char DEB;

const double ToleranceForLSSolve=1e-12;

void PrintVector(Rcpp::NumericVector v,std::string name)
{
 Rcpp::Rcout << name << "=[ ";
 for (int i=0; i<v.length(); i++)
  Rcpp::Rcout << v(i) << " ";
 Rcpp::Rcout << "];\n";
}

void PrintMatrix(Rcpp::NumericMatrix M,std::string name)
{
 int n=M.nrow();
 int m=M.ncol();
 
 Rcpp::Rcout << name << "=[ ";
 for (int i=0; i<n; i++)
 {
  for (int j=0; j<m; j++)
   Rcpp::Rcout << M(i,j) << " ";
  if (i<(n-1))
   Rcpp::Rcout << "\n";
  else
   Rcpp::Rcout << "];\n";
 }
}

int CheckDistName(std::string model)
{ 
 // Convert to lower case for standarization
 std::transform(model.begin(), model.end(), model.begin(),[](unsigned char c){ return std::tolower(c); });
 
 int i=0;
 while (i<NPOSSDIST && model!=distnames[i]) 
  i++;
 if (i>=NPOSSDIST)
  Rcpp::stop("%s is an incorrect model name.\n",model);
 
 if (DEB)
  Rcpp::Rcout << "Model type set to " << model << " (model number " << i << ")\n";
    
 // TO BE REMOVED:
 if (i == DBINLOGITW)
  Rcpp::stop("Sorry, grouped binomial with logistic distributions has not yet been implemented.\n");
  
 return(i);
}

//' SolveLinearSystem
//' 
//' Solves the linear system Ax=b where A is a numeric matrix of size m x n (with n>1), b is a column vector of length m
//' and x (the solution) is a column vector of length n.
//'  
//' The solution is obtained in the minimal square error sense (i.e.: sum_i |Ax_i-b|^2 is a minimum) using
//' non-explicit QR factorization of A. Therefore, the number of equations must be bigger or equal than the number
//' of unknows so m must be >= n. A is supposed to be well conditioned (i.e.: rank n). Numerical problems may arise
//' for badly conditioned A (see Error management)\cr
//'
//' WARNING: matrix A IS CHANGED inside this function. By efficiency, we have opted for NOT copying the argument in
//' a new, temporary matrix inside the function. If we want to keep your original A make a deep copy with\cr
//'                   Acopy<-matrix(A,nrow=m,ncol=n)\cr
//' before calling this function.\cr
//' 
//' After the call A will contain the upper-triangular R matrix in its upper-diagonal portion and the second-to-final
//' components of each Housholder vector v_i in the lower-diagonal part of column i. The first component is not stored
//' (no place to keep it) but it can be recovered since the Housholder vectors are of modulus 1. The Householder
//' matrix H_i can be constructed as\cr
//'                  H_i = I(m)-2*v_i'*v_i\cr
//' being I(m) the identity matrix of size n x n and v_i' the transposed of v_i.\cr
//' Finally, if the Q matrix is needed, it is\cr
//'                  Q=H(n-1) H(n-2) ... H_0 (in that order).\cr
//' Indeed, system is solved writing A*x*b as\cr
//'                  Q*R*x=b ==> R*x = Q'*b = H(n-1)*H(n-2)* ... *H_0 b\cr
//' (notice that inv(Q)=Q' and that H'(i)=H(i)\cr
//' Multiplications are done accumulatively as\cr
//'                  b(1)=H(1)*b, b(2)=H(2)*b(1), b(n-1)=H(n-1)*b(n-2), ...\cr
//' and the final triangular system\cr
//'                  R*x=b(n-1)\cr
//' is solved by back-susbtitution.
//' 
//' Error management:\cr
//' This function calls Rcpp::stop to indicate argument passing errors but it returns an empty vector to indicate lack
//' of solution because of numerical unstability, i.e.: the norm of a Housholder vector is below the tolerance. This
//' is to allow other functions try with different models if one of them cannot be adjusted. This means that the caller
//' function SHOULD check the length of the returned solution vector. 
//'
//' @param  A      Numeric matrix (system matrix)
//' @param  b      Numeric vector (independent term)
//' @param  tol    Numeric tolerance for the norm of the Housholder vector. Possitive number. Default: 1e-9
//' @return x      A numeric vector of length n with the solution\cr
//' 
//' @examples
//' m=9
//' n=6
//' A=matrix(nrow=m,ncol=n)
//' for (i in c(1:n))
//' {
//'  A[,i]=c(rep(0,i-1),c(1:(9-i+1)))
//' }
//' b=as.vector(c(m:1))
//' x=SolveLinearSystem(A,b)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector SolveLinearSystem(Rcpp::NumericMatrix &A,const Rcpp::NumericVector &b,double tol=1e-9)
{
 int m=A.nrow();
 int n=A.ncol();
 
 if (n<2)
  Rcpp::stop("Number of columns of matrix A must be bigger or equal than 2.\n");
 
 if (m<n)
  Rcpp::stop("Number of rows of matrix A must be bigger or equal than number of columns.\n");

 if (b.length()!=m)
  Rcpp::stop("Length of vector b must be equal to the number of rows of matrix A.\n");
 
 // Set default value for tolerance
 if (tol<0.0)
  Rcpp::stop("Tolerance must be a possitive (small) number).\n");
 else
 {
  if (tol>1e-6)
   Rcpp::stop("Maximal tolerance is 1e-6.\n");
 } 
      
 // This variable will be used to store the first element of each Householder vector, the one that would be
 // stored in the diagonal of A if we could overwite the elements of the diagonal (which we CANNOT)
 Rcpp::NumericVector vdiag(m); 
 
 // These are local variables to store norms and dot products of the used vectors.
 double vnorm, vTa, vpartdot;
 for (int i = 0; i < n; i++)
 {
    // A temporary vector to hold a part of A we need and that we are going to change
    Rcpp::NumericVector vh(m);
    /* Set vh equal to subvector A[i : m][i] (the lower-diagonal part of column i) */
    for (int k=0; k<(m-i); k++)
     vh[k]=A(i+k,i);

    /* vpartdot = ||vh||^2 - vh[0] * vh[0]; since vpartdot 
       is unaffected by the change in vh[0], storing this value 
       prevents the need to recalculate the entire norm of vh 
       after updating vh[0] in the following step              */
    vpartdot=0.0;
    for (int k=1; k<(m-i); k++)
     vpartdot += vh[k]*vh[k];
        
    /* Set vh[0] = vh[0] + sign(vh[0]) * ||v[i]|| */
    if (vh[0] < 0.0)
     vh[0] -= sqrt(vh[0] * vh[0] + vpartdot);
    else
     vh[0] += sqrt(vh[0] * vh[0] + vpartdot);
        
    /* Normalize vh */
    vnorm = sqrt(vh[0] * vh[0] + vpartdot);
    if (vnorm<tol)
    {
     Rcpp::NumericVector empty(0);
     return(empty);
    }
    
    for (int k=0; k<(m-i); k++)
     vh[k] /= vnorm;
       
    /* Make the pivot-reduction step to get 0 in the lower diagonal of column i
       (without setting anything not zero in the formerly nullified lower-diagonal columns of A)
    */   
    for(int j = i; j < n; j++)
    {
        /* set A[i:m][j] = A[i:m][j] - 2 * (vh^T A[i:m][j]) * vh */
        vTa = 0.0;
        for (int k=0; k<(m-i); k++)
          vTa += A(k+i,j)*vh[k];
           
        vTa *= 2;
        for (int k=0; k<(m-i); k++)
          A(k+i,j) -= vTa*vh[k];
    }
    
    /* Now, the Householder vector is stored in the lower-diagonal of this row
       of A. It is OK to overwrite this space, since in the former loop it should
       has become 0's (that's precisely the objective of this loop: to transform
       gradually A into an upper-diagonal matrix). The first element of the Householder
       version cannot be stored in A, too, since it would ovewrite the main diagonal of
       A that we will need later. This is why we use another vector, vdiag
       
       Notice that it is possible to have a probably clearer version of this program
       by declaring a matrix v to store the Householder vectors but it would need twice
       the space and in this way we work "on the place". But even in this case, A would
       be being ovewritten, anyway.
    */
    vdiag[i]=vh[0];
    for (int k=1; k<m-i; k++)
     A(i+k,i)=vh[k];
 }

 Rcpp::NumericVector bmod(m);
 for (int i=0; i<m; i++)
  bmod[i]=b[i];
 
 double vhk,vhl; 
 for (int i=0; i < n; i++)
 {
  /*
    More inefficient (but clearer) version of the program:
    
    The lower-diagonal part of A is copied back to the Housholder vector (and the first component is restored, too).
    
    Rcpp::NumericVector vh(m);
    vh[0]=vdiag[i];
    for (int k=1; k<m-i; k++)
      vh[k]=A(k+i,i);
 
   
    Even the Householder matrix can be explictly constructed
     
    Rcpp::NumericMatrix H(m-i,m-i);
    for (int r=0; r<m-i; r++)
     for (int c=0; c<m-i; c++)
      if (r==c)
       H(r,c) = 1.0-2.0*vh[r]*vh[c];
      else
       H(r,c) = -2.0*vh[r]*vh[v]; 
       
  */
  
  // The modified part of the vector of independent terms
  Rcpp::NumericVector newb(m-i);
  
  for (int k=0; k<m-i; k++)
  {
   newb[k]=0.0;
   for (int l=0; l<m-i; l++)
   {
    /* The next line if we would have constructed the Householder matrix would be
     
       newb[k] += H(k,l)*bmod[l+i];
     
       Without constructing explicitly the H matrix this line would be replaced by:
       
       newb[k] -= 2.0*vh[k]*vh[l]*bmod[l+i];
       if (k==l)
        newb[k] += bmod[l+i];
        
       and these, getting rid even of the vh vector, becomes:        
    */
    vhk=(k==0) ? vdiag[i] : A(k+i,i);
    vhl=(l==0) ? vdiag[i] : A(l+i,i);
    
    newb[k] -= 2.0*vhk*vhl*bmod[l+i];
    if (k==l)
     newb[k] += bmod[l+i];
   }
  }
 
  /*
    The modified vector of independent terms is copied back to b for the next multiplication.
    If the original b vector is b_0, at each step,
    b_i = H_i * b_{i-1}
    being H_i the i´-th Housholder matrix.
    I think it could have been possible to make some tricks to avoid the temporary vector newb and
    save space and assignments but it is just a vector, not a matrix, and it doesn't pay for...
  */
  for (int k=0; k<m-i; k++)
   bmod[k+i]=newb[k]; 
 }
  
 // Finally, declare the solution vector...
 Rcpp::NumericVector x(n);

 // ... and fill it by back-substitution using the upper-diagonal part of A (i.e.: R)
 double sk;
 for (int k=n-1; k>=0; k--)
 {
  sk=0.0;
  for (int l=n-1; l>=k; l--)
   sk += A(k,l)*x[l];
  x[k] = (bmod[k]-sk)/A(k,k);
 }
 
 return(x);
}

// Auxiliary function to CalculateInverseFromModifiedA. Not to be exported.
Rcpp::NumericMatrix InverseFromUpperTriangular(Rcpp::NumericMatrix R)
{
 if (DEB==LL_DEB)
  PrintMatrix(R,"Rpassed");
  
 int n=R.nrow();
 
 Rcpp::NumericVector d1(n);
 for (int i=0; i<n; i++)
  d1(i)=1.0/R(i,i);
 
 if (DEB==LL_DEB)
  PrintVector(d1,"d");
 
 Rcpp::NumericMatrix N(n,n);
 Rcpp::NumericMatrix Np(n,n);
 Rcpp::NumericMatrix IN1(n,n);
 for (int i=0; i<n; i++)
 {
  for (int j=i; j<n; j++)
  {
   N(i,j) += d1(i)*R(i,j);
   if (i==j)
    N(i,j) -= 1.0;
  }
  for (int j=0; j<n; j++)
  {
   Np(i,j)=N(i,j);
   if (i==j)
    IN1(i,j)=1.0;
  }
 }
 
 if (DEB==LL_DEB)
 {
  PrintMatrix(Np,"Np");
  PrintMatrix(IN1,"IN1st");
 }
 
 for (int iter=0; iter<(n-1); iter++)
 {
  for (int i=0; i<n; i++)
   for (int j=0; j<n; j++)
    IN1(i,j) += (iter%2) ? Np(i,j) : -Np(i,j);
  Rcpp::NumericMatrix Ntr(n,n);
  for (int i=0; i<n; i++)
   for (int j=i; j<n; j++)
    for (int k=i; k<=j; k++)
     Ntr(i,j) += Np(i,k)*N(k,j);
  for (int i=0; i<n; i++)
   for (int j=0; j<n; j++)
    Np(i,j) = Ntr(i,j);
 }
 
 if (DEB==LL_DEB)
  PrintMatrix(IN1,"IN1");
 
 Rcpp::NumericMatrix Ri(n,n);
 for (int i=0; i<n; i++)
  for (int j=0; j<n; j++)
   Ri(i,j) += IN1(i,j)*d1(j);
  
 if (DEB==LL_DEB)
  PrintMatrix(Ri,"Ri");
 
 return(Ri);
}

// Auxiliary function to SolveLinearSystemWithInverse. Not to be exported.
Rcpp::NumericMatrix CalculateInverseFromModifiedA(Rcpp::NumericMatrix T,Rcpp::NumericVector vdiag)
{
 // T is the system matrix, normally called A, modified by the QR decomposition so it
 // implicitly contains the Q and R parts. Therefore, Q and R must be rebuilt. This is
 // the same piece of code as in QRDec 
 
 int m=T.nrow();
 
 if (DEB==LL_DEB)
  PrintMatrix(T,"Amod");
 
 // Temporary matrix Qt initialized to identity
 Rcpp::NumericMatrix Qt(m,m);
 for (int d=0; d<m; d++)
  Qt(d,d)=1.0;
 
 double t;
 Rcpp::NumericMatrix Q(m,m);
 // We have here always square matrices so n=m
 for (int i=0; i<m; i++)
 {
  // Householder matrix created and destroyed at every loop pass
  Rcpp::NumericMatrix Hi(m,m); 
  // The first i places of the diagonal are ones
  for (int k=0;k<i;k++)
   Hi(k,k)=1.0;
  
  // The others are built recovering the values in the lower part of the corresponding column of A
  // Notice that the diagonal of A does not contain the values we need know, we stored them in vdiag 
  for (int k=i; k<m; k++)
   for (int l=i; l<m; l++)
   {
    if (k==l)
     Hi(k,l) = 1.0-2.0*((k==i) ? vdiag[i] : T(k,i))*((l==i) ? vdiag[i] : T(l,i));
    else
     Hi(k,l) = -2.0*((k==i) ? vdiag[i] : T(k,i))*((l==i) ? vdiag[i] : T(l,i));
   }
   
  // This is the product Qt*H_i, which is directly stored in Q...
  for (int k=0; k<m; k++)
   for (int l=0; l<m; l++)
   {
    t=0.0;
    for (int r=0; r<m; r++)
     t += Qt(k,r)*Hi(r,l);
    Q(k,l)=t;
   }
   
  // ... and, for each loop pass except the last one, Q is copied into Qt to set it prepared for next pass
  if (i<m-1)
  {
   for (int k=0; k<m; k++)
    for (int l=0; l<m; l++)
     Qt(k,l)=Q(k,l);
  }
 
 }

 // Finally, R is filled with the upper-triangular part (including main diagonal) of T
 Rcpp::NumericMatrix R(m,m);
 for (int col=0; col<m; col++)
  for (int row=0; row<m; row++)
   R(row,col)=(row>col) ? 0.0 : T(row,col);
   
 if (DEB==LL_DEB)
 {
  PrintMatrix(Q,"Q");
  PrintMatrix(R,"R");
 }
 
 // and the inverse of R is calculated using appropriate decomposition (see code of the function)
 Rcpp::NumericMatrix R1=InverseFromUpperTriangular(R);
 
 // Finally, our desired inverse is R1*Q^t, since Q^-1=Q^t
 Rcpp::NumericMatrix Inverse(m,m);
 
 double s;
 for (int r=0; r<m; r++)
  for (int c=0; c<m; c++)
  {
   s=0.0;
   for (int k=0; k<m; k++)
    s += R1(r,k)*Q(c,k);    // Q(c,k) and not Q(k,c), as in usual matrix product, because we want Q^t
   Inverse(r,c)=s; 
  }
 
 if (DEB==LL_DEB)
  PrintMatrix(Inverse,"Inv");
 
 return(Inverse);
}

//' SolveLinearSystemWithInverse
//' 
//' Solves the linear system Ax=b where A is a square numeric matrix of size n x n (with n>1), b is a column vector of length n
//' and x (the solution) is a column vector of length n. It also returns the inverse of matrix A which, if used for estimation
//' of a linear model, is the variance/covariance matrix. Inversion is done using the QR decomposition of A used to solve the
//' linear system, so it should be faster than using simply SolveLinearSystem and invert the A matrix later.
//'  
//' The solution is obtained using QR factorization of A. A is supposed to be well conditioned (i.e.: rank n). Numerical problems
//' may arise for badly conditioned A (see Error management)\cr
//'
//' WARNING: matrix A IS CHANGED inside this function. By efficiency, we have opted for NOT copying the argument in
//' a new, temporary matrix inside the function. If you want to keep your original A make a deep copy with\cr
//'                   Acopy<-matrix(A,nrow=m,ncol=n)\cr
//' before calling this function.\cr
//' 
//' After the call A will contain the upper-triangular R matrix in its upper-diagonal portion and the second-to-final
//' components of each Housholder vector v_i in the lower-diagonal part of column i. The first component is not stored
//' (no place to keep it) but it can be recovered since the Housholder vectors are of modulus 1. The Householder
//' matrix H_i can be constructed as\cr
//'                  H_i = I(m)-2*v_i'*v_i\cr
//' being I(m) the identity matrix of size n x n and v_i' the transposed of v_i.\cr
//' Finally, if the Q matrix is needed, it is\cr
//'                  Q=H(n-1) H(n-2) ... H_0 (in that order).\cr
//' Indeed, system is solved writing A*x*b as\cr
//'                  Q*R*x=b ==> R*x = Q'*b = H(n-1)*H(n-2)* ... *H_0 b\cr
//' (notice that inv(Q)=Q' and that H'(i)=H(i)\cr
//' Multiplications are done accumulatively as\cr
//'                  b(1)=H(1)*b, b(2)=H(2)*b(1), b(n-1)=H(n-1)*b(n-2), ...\cr
//' and the final triangular system\cr
//'                  R*x=b(n-1)\cr
//' is solved by back-substitution.
//' 
//' Inverse of A is done using the QR-decomposition calculated to solve the system. Indeed, if\cr
//' A=Q*R then\cr
//' A^-1=R^-1*Q^-1.\cr
//' Being Q orthogonal, Q^-1 = Q^T, and for the upper-triangular matrix R, it can be written as\cr
//' R=D(I+N) where D is a diagonal matrix and N is nilpotent (i.e.: N^n=0). In such a case,\cr
//' R^-1=(I+N)^-1*D^-1.\cr
//' where, since D is diagonal, calculating D^-1 is trivial and a matrix lemma stablishes that, for nilpotent matrices N,\cr
//' (I+N)^-1=I+sum_p=1..(n-1) N^p\cr
//' 
//' Error management:\cr
//' This function calls Rcpp::stop to indicate argument passing errors but it returns an empty vector to indicate lack
//' of solution because of numerical unstability, i.e.: the norm of a Housholder vector is below the tolerance. This
//' is to allow other functions try with different models if one of them cannot be adjusted. This means that the caller
//' function SHOULD check the length of the returned solution vector. 
//'
//' @param  A      Numeric matrix (system matrix)
//' @param  b      Numeric vector (independent term)
//' @param  tol    Numeric tolerance for the norm of the Housholder vector. Possitive number. Default: 1e-9
//' @return A list with the following keys:
//' \itemize{
//' \item x      A numeric vector of length n with the solution\cr
//' \item Inv    A numeric matrix with the inverse of A\cr
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List SolveLinearSystemWithInverse(Rcpp::NumericMatrix &A,const Rcpp::NumericVector &b,double tol=1e-9)
{
 int m=A.nrow();
 int n=A.ncol();
 
 if (n<2)
  Rcpp::stop("Number of columns of matrix A must be bigger or equal than 2.\n");
 
 if (m<n)
  Rcpp::stop("Number of rows of matrix A must be bigger or equal than number of columns.\n");

 if (b.length()!=m)
  Rcpp::stop("Length of vector b must be equal to the number of rows of matrix A.\n");
 
 // Set default value for tolerance
 if (tol==0.0)
  tol=1e-9;
 else
 {
  if (tol<0.0)
   Rcpp::stop("Tolerance must be a possitive (small) number).\n");
  else
  {
   if (tol>1e-6)
    Rcpp::stop("Maximum tolerance is 1e-6.\n");
  } 
 }
     
 if (DEB==LL_DEB)
  PrintMatrix(A,"A");
  
 // This variable will be used to store the first element of each Householder vector, the one that would be
 // stored in the diagonal of A if we could overwite the elements of the diagonal (which we CANNOT)
 Rcpp::NumericVector vdiag(m); 
 
 // These are local variables to store norms and dot products of the used vectors.
 double vnorm, vTa, vpartdot;
 for (int i = 0; i < n; i++)
 {
    // A temporary vector to hold a part of A we need and that we are going to change
    Rcpp::NumericVector vh(m);
    // Set vh equal to subvector A[i : m][i] (the lower-diagonal part of column i) 
    for (int k=0; k<(m-i); k++)
     vh[k]=A(i+k,i);

    // vpartdot = ||vh||^2 - vh[0] * vh[0]; since vpartdot 
    // is unaffected by the change in vh[0], storing this value 
    // prevents the need to recalculate the entire norm of vh 
    // after updating vh[0] in the following step            
    vpartdot=0.0;
    for (int k=1; k<(m-i); k++)
     vpartdot += vh[k]*vh[k];
        
    // Set vh[0] = vh[0] + sign(vh[0]) * ||v[i]||
    if (vh[0] < 0.0)
     vh[0] -= sqrt(vh[0] * vh[0] + vpartdot);
    else
     vh[0] += sqrt(vh[0] * vh[0] + vpartdot);
        
    // Normalize vh
    vnorm = sqrt(vh[0] * vh[0] + vpartdot);
    if (vnorm<tol)
    {
     Rcpp::NumericVector emptyv(0);
     Rcpp::NumericMatrix emptym(0,0);
     Rcpp::List ret = Rcpp::List::create(
                                         Rcpp::Named("x") = emptyv,
                                         Rcpp::Named("Inv") =  emptym);
     return(ret);
    }
    
    for (int k=0; k<(m-i); k++)
     vh[k] /= vnorm;
       
    // Make the pivot-reduction step to get 0 in the lower diagonal of column i
    // (without setting anything not zero in the formerly nullified lower-diagonal columns of A)
      
    for(int j = i; j < n; j++)
    {
        // set A[i:m][j] = A[i:m][j] - 2 * (vh^T A[i:m][j]) * vh
        vTa = 0.0;
        for (int k=0; k<(m-i); k++)
          vTa += A(k+i,j)*vh[k];
           
        vTa *= 2;
        for (int k=0; k<(m-i); k++)
          A(k+i,j) -= vTa*vh[k];
    }
    
    // Now, the Householder vector is stored in the lower-diagonal of this row
    // of A. It is OK to overwrite this space, since in the former loop it should
    // has become 0's (that's precisely the objective of this loop: to transform
    // gradually A into an upper-diagonal matrix). The first element of the Householder
    // version cannot be stored in A, too, since it would ovewrite the main diagonal of
    // A that we will need later. This is why we use another vector, vdiag
       
    // Notice that it is possible to have a probably clearer version of this program
    // by declaring a matrix v to store the Householder vectors but it would need twice
    // the space and in this way we work "on the place". But even in this case, A would
    // be being ovewritten, anyway.
    
    vdiag[i]=vh[0];
    for (int k=1; k<m-i; k++)
     A(i+k,i)=vh[k];
 }

 Rcpp::NumericVector bmod(m);
 for (int i=0; i<m; i++)
  bmod[i]=b[i];
 
 double vhk,vhl; 
 for (int i=0; i < n; i++)
 {
  
  // More inefficient (but clearer) version of the program:
    
  // The lower-diagonal part of A is copied back to the Housholder vector (and the first component is restored, too).
    
  // Rcpp::NumericVector vh(m);
  // vh[0]=vdiag[i];
  // for (int k=1; k<m-i; k++)
  //  vh[k]=A(k+i,i);
 
  // Even the Householder matrix can be explictly constructed
     
  // Rcpp::NumericMatrix H(m-i,m-i);
  // for (int r=0; r<m-i; r++)
  //  for (int c=0; c<m-i; c++)
  //   if (r==c)
  //    H(r,c) = 1.0-2.0*vh[r]*vh[c];
  //   else
  //    H(r,c) = -2.0*vh[r]*vh[v]; 
       
 
  // The modified part of the vector of independent terms
  Rcpp::NumericVector newb(m-i);
  
  for (int k=0; k<m-i; k++)
  {
   newb[k]=0.0;
   for (int l=0; l<m-i; l++)
   {
    // The next line if we would have constructed the Householder matrix would be
     
    // newb[k] += H(k,l)*bmod[l+i];
     
    // Without constructing explicitly the H matrix this line would be replaced by:
       
    // newb[k] -= 2.0*vh[k]*vh[l]*bmod[l+i];
    // if (k==l)
    //  newb[k] += bmod[l+i];
        
    // and these, getting rid even of the vh vector, becomes:        
    
    vhk=(k==0) ? vdiag[i] : A(k+i,i);
    vhl=(l==0) ? vdiag[i] : A(l+i,i);
    
    newb[k] -= 2.0*vhk*vhl*bmod[l+i];
    if (k==l)
     newb[k] += bmod[l+i];
   }
  }
 
  // The modified vector of independent terms is copied back to b for the next multiplication.
  // If the original b vector is b_0, at each step,
  // b_i = H_i * b_{i-1}
  // being H_i the i´-th Housholder matrix.
  // I think it could have been possible to make some tricks to avoid the temporary vector newb and
  // save space and assignments but it is just a vector, not a matrix, and it doesn't pay for...
  
  for (int k=0; k<m-i; k++)
   bmod[k+i]=newb[k]; 
 }
  
 // Finally, declare the solution vector...
 Rcpp::NumericVector x(n);

 // ... and fill it by back-substitution using the upper-diagonal part of A (i.e.: R)
 double sk;
 for (int k=n-1; k>=0; k--)
 {
  sk=0.0;
  for (int l=n-1; l>=k; l--)
   sk += A(k,l)*x[l];
  x[k] = (bmod[k]-sk)/A(k,k);
 }
 
 Rcpp::NumericMatrix Ai=CalculateInverseFromModifiedA(A,vdiag);
 return Rcpp::List::create(
     Rcpp::Named("x") = x,
     Rcpp::Named("Inv") = Ai );
}


//' QRdec
//' 
//' Makes the QRdecomposition of a matrix A as A=Q*R using the Householder method.\cr
//'
//' WARNING: if you are trying to solve a linear system of the form Ax=B, DO NOT use this function\cr
//' Use instead SolveLinearSystem. Both operate with the same algorithm and are equally efficient in
//' the decomposition but this one explicitly constructs Q and R which uses more memory.\cr
//'
//' WARNING: the passed matrix IS CHANGED inside this function. By efficiency, we have opted for NOT copying
//' the argument in a new, temporary matrix inside the function. If we want to keep your original A make a 
//' deep copy with\cr
//'                   Acopy<-matrix(A,nrow=m,ncol=n)\cr
//' before calling this function.\cr
//'
//'
//' @param  A	  A (m x n) numeric matrix to be decomposed with m>=n
//' @param  tol   Numeric tolerance for the norm of the Housholder vector. Possitive number. Default: 1e-9
//' @return ret   A list with two elements
//'               ret["Q"] The orthogonal matrix Q
//'               ret["R"] The upper-triangular matrix R
//' @examples
//' m=10
//' n=8
//' A=matrix(runif(m*n),nrow=m,ncol=n)
//' Acop=matrix(A,nrow=m,ncol=n)
//' ret <- QRdec(A)
//' Arec <- ret$Q  %*% ret$R
//' error <- max(abs(Arec-Acop))
//' @export
// [[Rcpp::export]]
Rcpp::List QRdec(Rcpp::NumericMatrix &A,double tol=1e-9)
{
 int m=A.nrow();
 int n=A.ncol();
 
 if (m<n)
  Rcpp::stop("Number of rows of matrix A must be bigger or equal than number of columns.\n");
 
 // Set default value for tolerance
 if (tol==0.0)
  tol=1e-9;
 else
  if (tol<0.0)
   Rcpp::stop("Tolerance must be a possitive (small) number).\n");
  else
   if (tol>1e-6)
    Rcpp::stop("Maximal tolerance is 1e-6.\n");
      
 // This variable will be used to store the first element of each Householder vector, the one that would be
 // stored in the diagonal of A if we could overwite the elements of the diagonal (which we CANNOT)
 Rcpp::NumericVector vdiag(m); 
 
 // These are local variables to store norms and dot products of the used vectors.
 double vnorm, vTa, vpartdot;
 for (int i = 0; i < n; i++)
 {
    // A temporary vector to hold a part of A we need and that we are going to change
    Rcpp::NumericVector vh(m);
    /* Set vh equal to subvector A[i : m][i] (the lower-diagonal part of column i) */
    for (int k=0; k<(m-i); k++)
     vh[k]=A(i+k,i);

    /* vpartdot = ||vh||^2 - vh[0] * vh[0]; since vpartdot 
       is unaffected by the change in vh[0], storing this value 
       prevents the need to recalculate the entire norm of vh 
       after updating vh[0] in the following step              */
    vpartdot=0.0;
    for (int k=1; k<(m-i); k++)
     vpartdot += vh[k]*vh[k];
        
    /* Set vh[0] = vh[0] + sign(vh[0]) * ||v[i]|| */
    if (vh[0] < 0.0)
     vh[0] -= sqrt(vh[0] * vh[0] + vpartdot);
    else
     vh[0] += sqrt(vh[0] * vh[0] + vpartdot);
        
    /* Normalize vh */
    vnorm = sqrt(vh[0] * vh[0] + vpartdot);
    if (vnorm<tol)
     Rcpp::stop("System matrix is ill-conditioned.\n");
     
    for (int k=0; k<(m-i); k++)
     vh[k] /= vnorm;
       
    /* Make the pivot-reduction step to get 0 in the lower diagonal of column i
       (without setting anything not zero in the formerly nullified lower-diagonal columns of A)
    */   
    for(int j = i; j < n; j++)
    {
        /* set A[i:m][j] = A[i:m][j] - 2 * (vh^T A[i:m][j]) * vh */
        vTa = 0.0;
        for (int k=0; k<(m-i); k++)
          vTa += A(k+i,j)*vh[k];
           
        vTa *= 2;
        for (int k=0; k<(m-i); k++)
          A(k+i,j) -= vTa*vh[k];
    }
    
    /* Now, the Householder vector is stored in the lower-diagonal of this row
       of A. It is OK to overwrite this space, since in the former loop it should
       has become 0's (that's precisely the objective of this loop: to transform
       gradually A into an upper-diagonal matrix). The first element of the Householder
       version cannot be stored in A, too, since it would ovewrite the main diagonal of
       A that we will need later. This is why we use another vector, vdiag
       
       Notice that it is possible to have a probably clearer version of this program
       by declaring a matrix v to store the Householder vectors but it would need twice
       the space and in this way we work "on the place". But even in this case, A would
       have been overwritten, anyway.
    */
    vdiag[i]=vh[0];
    for (int k=1; k<m-i; k++)
     A(i+k,i)=vh[k];
 }

 /* The next lines make iteratively the multiplication of the Householder matrices and
    store the results in a temporaty matrix, Qt, which is at every loop pass copied in Q, i.e.
    Qt <- Identity
    for i=0..n-1
     Q <- Qt*H_i
    Qt <- Q
    so, for i==0, Q is H0; for i==1, Q is H0*H1, etc.
 */    
 
 // Temporary matrix Qt initialized to identity
 Rcpp::NumericMatrix Qt(m,m);
 for (int d=0; d<m; d++)
  Qt(d,d)=1.0;
 
 double t;
 Rcpp::NumericMatrix Q(m,m);
 for (int i=0; i<n; i++)
 {
  // Householder matrix created and destroyed at every loop pass
  Rcpp::NumericMatrix Hi(m,m); 
  // The first i places of the diagonal are ones
  for (int k=0;k<i;k++)
   Hi(k,k)=1.0;
  
  // The others are built recovering the values in the lower part of the corresponding column of A
  // Notice that the diagonal of A does not contain the values we need know, we stored them in vdiag 
  for (int k=i; k<m; k++)
   for (int l=i; l<m; l++)
   {
    if (k==l)
     Hi(k,l) = 1.0-2.0*((k==i) ? vdiag[i] : A(k,i))*((l==i) ? vdiag[i] : A(l,i));
    else
     Hi(k,l) = -2.0*((k==i) ? vdiag[i] : A(k,i))*((l==i) ? vdiag[i] : A(l,i));
   }
   
  // This is the product Qt*H_i, which is directly stored in Q...
  for (int k=0; k<m; k++)
   for (int l=0; l<m; l++)
   {
    t=0.0;
    for (int r=0; r<m; r++)
     t += Qt(k,r)*Hi(r,l);
    Q(k,l)=t;
   }
   
  // ... and, for each loop pass except the last one, Q is copied into Qt to set it prepared for next pass
  if (i<n-1)
  {
   for (int k=0; k<m; k++)
    for (int l=0; l<m; l++)
     Qt(k,l)=Q(k,l);
  }
 
 }

 // Finally, R is filled with the upper-triangular part (including main diagonal) of A
 Rcpp::NumericMatrix R(m,n);
 for (int col=0; col<n; col++)
  for (int row=0; row<=col; row++)
   R(row,col)=A(row,col);
   
 return Rcpp::List::create(
     Rcpp::Named("Q") = Q,
     Rcpp::Named("R") = R );

}

//' StGLM
//'
//' Gets a vector of observations, a type of link function, the phi parameter, the offset
//' vector and the Phi-star parameter and returns the vector of estimators beta, the value of X^2 and the final
//' values of the deviance and of the covariance matrix of the maximum likelihood estimator of the model coefficients
//'
//' Error management:\cr
//' This function calls Rcpp::stop only to signal error in pass of parameters. Errors due to numerical unstability
//' are not raised (at most, a warning is called, see parameter warn_me) but function returns an empty vector as beta value,
//' 0 for X2, sums and Dev and the empty matrix for CovM. This MUST be checked by the caller function and it is to allow it
//' to test a different, simpler model if the current one cannot be fitted.
//'
//' @param  y	  		A numeric vector of observations, length m
//' @param  X	  		A numeric matrix of predictors, size m x n
//' @param  model  		A string with the model to use. Possible values: 'Poisson', 'BinLogit', 'NegBin' or 'BinLogitW'\cr
//'                     for Possion, Binomial with logistic link, Negbinomial, and Grouped Binomial with logistic link, respectively.
//' @param  offset 		A numeric vector of offsets, length m
//' @param  warn_me     A logical value indicating if a warning must be raised in the case of numerical unstability 
//' @param  Phi    		A real value, the overdispersion parameter. Used only for Negbinomial. Default: 0.0
//' @param  retmu       A boolean value to indicate whether or not the vector of estimated means must be returned. Default: FALSE
//' @return A list with the following keys:
//' \itemize{
//' \item beta  		A numeric vector with the estimated model coefficients, or empty vector if they cannot be estimated
//' \item muest         A numeric vector with the estimated means, or NULL if parameter retmu was FALSE
//' \item X2            The statistic X^2, defined as sum((y_i-muest_i)^2/Vest_i). 0 if the model cannot be estimated
//' \item sums          The statistic sums, defined as sum((y_i-muest_i)*Vest'_i/Vest_i). 0 if the model cannot be estimated
//' \item Dev			The deviance of the model. 0 if the model cannot be estimated
//' \item CovM          Variance-covariance matrix of the maximum likelihood estimator of the model coefficients. Empty matrix if the model cannot be estimated
//' }
//' @examples
//' # To be done
//' @export
// [[Rcpp::export]]
Rcpp::List StGLM(Rcpp::NumericVector y,Rcpp::NumericMatrix X,std::string model,Rcpp::NumericVector offset,bool warn_me,double Phi=0.0,bool retmu=false)
{
 int m=y.size();
 int n=X.ncol();
 
 // Some cheks of validity of input parameters
 if (X.nrow()!=m || offset.length()!=m )
  Rcpp::stop("Incorrect dimension for either matrix X, vector y or offset vector in StGLM.\n");
 
 int nummodel=CheckDistName(model);
 
 if ((nummodel==DNEGBIN) && (Phi==0.0))
  Rcpp::stop("For negative binomial model, Phi cannot be 0. Please, provide an estimated value.\n");
 
 if ((nummodel==DBINLOGIT) || (nummodel==DBINLOGITW))
 {
  int i=0;
  while (i<m && ((y[i]==1.0) || (y[i]==0.0)) )
   i++;
  if (i<m)
   Rcpp::stop("For both types of logistic regressions the response vector must be composed only by 0 and 1. In this case y[%d] is %f.\n",i,y[i]);
 }
   
 // Calculation of the mean of the data  
 double meany=0.0;
 for (int i=0; i<m; i++)
  meany+=y[i];
 meany/=double(m);

 // This is only for low level debug (TRUE,TRUE)
 if (DEB==LL_DEB)
 {
  Rcpp::Rcout << "y=[";
  for (int i=0;i<m;i++)
   Rcpp::Rcout << " " << y[i];
  Rcpp::Rcout << "];\n";
  Rcpp::Rcout << "meany=" << meany << "\n";
 }
 
 Rcpp::NumericVector mu(m);
 Rcpp::NumericVector eta(m);
 // Vector of estimated means is initialized according to the family
 switch (nummodel)
 {
  case DPOISS:
  	   for (int i=0; i<m; i++)
       {
        mu[i]=(y[i]+meany)/2.0;
        eta[i]=log(mu[i]);
       }
       break;
  case DNEGBIN:                   // Initialization is the same for Poisson and negative binomial
       for (int i=0; i<m; i++)
       {
        mu[i]=(y[i]+meany)/2.0;
        eta[i]=log(mu[i]);
       }
       break;
  case DBINLOGIT:
       for (int i=0; i<m; i++)
       {
        mu[i]=(y[i]+0.5)/2.0;
        eta[i]=log(mu[i]/(1.0-mu[i]));
       }
       break;
  case DBINLOGITW:
       Rcpp::stop("Sorry: grouped binomial with logistic link function has not yet been implemented.\n");
       break;  
  default: 
       Rcpp::stop("Incorrect model number??? (We should not have arrived here...)\n");
       break;
 }
 
 if (DEB==LL_DEB)
 {
  Rcpp::Rcout << "mu=[";
  for (int i=0;i<m;i++)
   Rcpp::Rcout << " " << mu[i];
  Rcpp::Rcout << "];\n";
  Rcpp::Rcout << "eta=[";
  for (int i=0;i<m;i++)
   Rcpp::Rcout << " " << eta[i];
  Rcpp::Rcout << "];\n";
 }
 
 double Dev=0.0,oldDev=0.0;
 double tolerance=1e-10;
 double DeltaDev=tolerance+1.0;
 
 // W is indeed a diagonal matrix but we store only the main diagonal
 // That's why we use w instead of W
 Rcpp::NumericVector w(m);

 Rcpp::NumericVector z(m);
 
 // This is matrix X^T * W * X
 Rcpp::NumericMatrix XtWX(n,n);
 
 // This is vector X^T * W * z
 Rcpp::NumericVector Xtwz(n);
 
 // This is the vector of coefficients that will be filled
 
 Rcpp::NumericVector beta(n);
 int iteration=0;
 
 Rcpp::List retfromSolve;
 Rcpp::NumericMatrix emptymatrix(0,0);
 
 while (fabs(DeltaDev) > tolerance && iteration<101)
 {
  if (DEB)
   Rcpp::Rcout << "Starting iteration " << iteration << " with fabs(DeltaDev)=" << fabs(DeltaDev) << " and tolerance " << tolerance << " for model number " << nummodel << "\n";
   
  // First, fill vectors w and z and "matrix" W (even we keep only its diagonal in a vector)
  // This MUST be done in each iteration since mu and eta will change at each iteration
  
  // Check here that the update formulae corresponds to each family
  switch (nummodel)
  {
   case DPOISS:
      for (int i=0; i<m; i++)
      { 
       w[i] = mu[i];
       z[i] = eta[i]+(y[i]-mu[i])/mu[i];
      }
      break;
   case DNEGBIN:
      for (int i=0; i<m; i++)
      { 
       w[i] = mu[i]/(1.0+Phi*mu[i]);
       z[i] = eta[i]+((y[i]-mu[i])/mu[i])-offset[i];
      }
      break;
   case DBINLOGIT:
      for (int i=0; i<m; i++)
      {
       w[i] = mu[i]*(1.0-mu[i]);
       z[i] = eta[i]+((y[i]-mu[i])/(mu[i]*(1.0-mu[i])))-offset[i]; 
      }
      break;
   case DBINLOGITW:
      Rcpp::stop("Sorry: grouped binomial with logistic link function has not yet been implemented.\n");
      break;
   default:
    Rcpp::stop("Incorrect model number??? (We should not have arrived here...)\n");
    break;   
  }
  
  // Only for low level debug...
  if (DEB==LL_DEB)
  {
   Rcpp::Rcout << "X=[\n";
   for (int i=0; i<m; i++)
   {
    for (int j=0; j<n; j++)
     Rcpp::Rcout << X(i,j) << " ";
    Rcpp::Rcout << ((i<m-1) ? "; " : "];\n");
   }
   Rcpp::Rcout << "w=[";
   for (int i=0; i<m; i++)
    Rcpp::Rcout << w[i] << " ";
   Rcpp::Rcout << "];\n";
   Rcpp::Rcout << "z=[";
   for (int i=0; i<m; i++)
    Rcpp::Rcout << z[i] << " ";
   Rcpp::Rcout << "];\n";
  }
  
  
  // First, get the matrix X'*W*X
  // Since W is diagonal, the multiplication can be done in a just one matrix product, not two.
  // This must be correct, since it is the same for all families and it works with Poisson
  for (int i=0; i<n; i++)
   for (int j=0; j<n; j++)
   {
    XtWX(i,j) = 0.0;
    for (int k=0; k<m; k++)
     XtWX(i,j) += X(k,i)*X(k,j)*w(k);
   }
   
  // Now, get the vector X'*W*z. Again, being W diagonal allows just one matrix-by-vector product
  // Same comment on correctness as below
  for (int i=0; i<n; i++)
  {
   Xtwz[i] = 0.0;
   for (int j=0; j<m; j++)
    Xtwz[i] += X(j,i)*w[j]*z[j];
  }
 
  if (DEB==LL_DEB)
  {
   Rcpp::Rcout << "XtWX=[\n";
   for (int i=0; i<n; i++)
   {
    for (int j=0; j<n; j++)
     Rcpp::Rcout << XtWX(i,j) << " ";
    Rcpp::Rcout << ((i<n-1) ? "\n" : "]\n");
   }
   Rcpp::Rcout << "Xtwz=[ ";
   for (int i=0; i<n; i++)
    Rcpp::Rcout << Xtwz[i] << " ";
   Rcpp::Rcout << "]\n";
  } 
  
  // The system we have to solve is XtWX * beta = Xtwz
  // Interestingly, note that, in the same way as for any matrix X'*X is symmetric, if W is diagonal,
  // X'*W*X is symmetric, too. This means the system can be solved either by LU or by SVD decomposition,
  // but we'll not do that but a solution based on QR decomposition. This is the solution used
  // by R, so certainly there must be a good reason...
  // Notice that SolveLinearSystems operates "on the place", so matrix XtWX is overwritten with its QR
  // decomposition. But this is not important since what we really want is its inverse, which is calculated
  // inside the function using precisely the QR decompostion (which is an additional reason to use it).
  
  retfromSolve=SolveLinearSystemWithInverse(XtWX,Xtwz,ToleranceForLSSolve);
  beta=retfromSolve["x"];
    
  if (beta.length()==0)
  {
   if (warn_me)
    Rcpp::warning("Solution numerically unstable.\n");

   Rcpp::List ret=Rcpp::List::create(
                                   Rcpp::Named("beta")=beta,
                                   Rcpp::Named("muest")=NULL,
                                   Rcpp::Named("X2")=0.0,
                                   Rcpp::Named("sums")=0.0,
                                   Rcpp::Named("Dev")=0.0,
                                   Rcpp::Named("CovM")=emptymatrix);
 
   return(ret);
  } 
  if (DEB)
  {
   Rcpp::Rcout << "    beta=[ ";
   for (int i=0; i<n; i++)
    Rcpp::Rcout << beta(i) << " ";
   Rcpp::Rcout << " ]\n";
  }
  
  // Finally, once we have got beta we use it to update vector eta, since eta = X*beta
  // Again, this product is the same for all families
  for (int i=0; i<m; i++)
  {
   eta[i] = 0.0;
   for (int j=0; j<n; j++)
    eta[i] += beta[j]*X(i,j);
  }

  // Check now that the formulae are correct for each distribution.
  // Notice that currently eta contains X*beta
  switch (nummodel)
  {
   case DPOISS:  
     for (int i=0; i<m; i++)
      mu[i] = exp(eta[i]); 
     break;
   case DNEGBIN:
     for (int i=0; i<m; i++) 
     {
      // For the negative binomial case, eta = X*beta + offset and mu=e^eta so we just add offset
      eta[i] += offset[i];
      mu[i] = exp(eta[i]);
     }
     break;
   case DBINLOGIT:
     for (int i=0; i<m; i++) 
     {
      // For the logistic non-grouped case, eta = X*beta + offset and mu=1/(1+e^(-eta))
      eta[i] += offset[i];
      mu[i] = 1.0/(1.0+exp(-eta[i]));
     }
     break;
   case DBINLOGITW:
     Rcpp::stop("Sorry: grouped binomial with logistic link function has not yet been implemented.\n");
     break;
   default:
     Rcpp::stop("Incorrect model number??? (We should not have arrived here...)\n");
     break;   
  }
  
  if (DEB==LL_DEB)
  {
   Rcpp::Rcout << "mu_updated=[ \n";
   for (int i=0; i<m; i++)
    Rcpp::Rcout << mu(i) << " ";
   Rcpp::Rcout << " ]\n";
  }
  
  oldDev = Dev;
  Dev=0.0;
  switch (nummodel)
  {
   // In Poisson and negative binomial the undeterminacy y*log(y/mu) when y is 0 is taken to be zero (its limit as y --> 0)
   case DPOISS:
     for (int i=0; i<m; i++)
      Dev += ((y[i]!=0.0) ? y[i]*log(y[i]/mu[i]) : 0.0 ) - (y[i]-mu[i]);
     Dev = 2*Dev;
     break;
   case DNEGBIN:
     for (int i=0; i<m; i++)
      Dev += ((y[i]!=0.0) ? y[i]*log(y[i]/mu[i]) : 0.0 ) - (y[i]+(1.0/Phi))*log((1+Phi*y[i])/(1+Phi*mu[i]));
     Dev = 2*Dev;
     break;
   case DBINLOGIT:
     for (int i=0; i<m; i++)
     {
      if (y[i]==1.0)
       Dev += (mu[i]!=1.0) ? log(1.0/mu[i]) : 0.0;
      else  // This will be always 0.0. The check for 0/1 of y was done before
       Dev += (mu[i]!=0.0) ? log(1.0/(1.0-mu[i])) : 0.0;
     }
     Dev = 2*Dev;
     break;
   case DBINLOGITW:
     Rcpp::stop("Sorry: grouped binomial with logistic link function has not yet been implemented.\n");
     break;
   default:
     Rcpp::stop("Incorrect model number??? (We should not have arrived here...)\n");
     break;
  }
  
  DeltaDev = Dev - oldDev;
  
  if (DEB)
   Rcpp::Rcout << "    Deviance calculated. Dev=" << Dev <<". DeltaDev=" << DeltaDev << "\n";
     
  iteration++;
 }
 
 double X2=0.0;
 double sums=0.0;
 switch (nummodel)
 {
  case DPOISS:
    for (int i=0; i<m; i++)
    {
     X2 += (y[i]-mu[i])*(y[i]-mu[i])/mu[i];
     sums += (y[i]-mu[i])/mu[i];
    }
    break;
  case DNEGBIN:
    for (int i=0; i<m; i++)
    {
     X2 += (y[i]-mu[i])*(y[i]-mu[i])/(mu[i]+Phi*mu[i]*mu[i]);
     sums += (y[i]-mu[i])*(1.0+2.0*Phi*mu[i])/(mu[i]+Phi*mu[i]*mu[i]);
    }
    break; 
  case DBINLOGIT:
    for (int i=0; i<m; i++)
    {
     X2 += (y[i]-mu[i])*(y[i]-mu[i])/(mu[i]*(1.0-mu[i]));
     sums += (y[i]-mu[i])*(1.0-2.0*mu[i])/(mu[i]*(1.0-mu[i]));
    } 
    break;
  case DBINLOGITW:
     Rcpp::stop("Sorry: grouped binomial with logistic link function has not yet been implemented.\n");
     break; 
  default:
    Rcpp::stop("Incorrect model number??? (We should not have arrived here...)\n");
    break;
 }
 
 Rcpp::List ret=Rcpp::List::create(
                                    Rcpp::Named("beta")=beta,
                                    Rcpp::Named("muest")=(retmu ? mu : NULL),
                                    Rcpp::Named("X2")=X2,
                                    Rcpp::Named("sums")=sums,
                                    Rcpp::Named("Dev")=Dev,
                                    Rcpp::Named("CovM")=retfromSolve["Inv"],
                                    Rcpp::Named("iters")=iteration);
  
 return(ret);
}
 
