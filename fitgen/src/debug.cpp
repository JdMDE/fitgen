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

#include <debug.h>

// This is the only place where this variable is declared.
// It is global for the full package and should only be changed by FitGeneSetDebug
unsigned char DEB=0;

//' FitGenSetDebug
//' 
//' Sets debugging in fitgen package to ON (with TRUE) or OFF (with FALSE) or even sets\cr
//' a lower level of debug to see inside the numerical algorithms.\cr
//' On package load the default status of both debugging levels is OFF.\cr
//' Setting any debug level to ON shows a message. Setting to OFF does not show anything (since debugging is OFF...)
//'
//' @param deb   boolean, TRUE to generate normal debug messages and FALSE to turn them off.
//' @param lldeb boolean, TRUE to generate low level debug messages (very verbose) and FALSE to turn them off. Setting this parameter
//' to TRUE sets also deb to TRUE. Default: FALSE
//' @return      No return value, called for side effects (internal boolean flag changed)
//' @examples
//' FitGenSetDebug(TRUE)
//' FitGenSetDebug(FALSE)
//' FitGenSetDebug(TRUE,TRUE)
//' FitGenSetDebug(TRUE,FALSE)
//' @export
// [[Rcpp::export]]
void FitGenSetDebug(bool deb,bool lldeb = false)
{
 if (lldeb)
 {
  DEB=LL_DEB;
  Rcpp::Rcout << "Normal and low level debugging for fitgen package both set to ON.\n";
 }
 else
 {
  if (deb)
  {
   DEB=NOR_DEB;
   Rcpp::Rcout << "Debugging for fitgen package set to ON.\n";
  }
  else
   DEB=0;
 }
}

//' FitGenGetDebug
//' Returns the current value of the debug flag (TRUE or FALSE).\cr
//' Used by the R functions of the package to know if they should emit messages or not.
//'
//' @return Boolean value with the current value of the internal DEB flag.
//' @examples
//' if (FitGenGetDebug()) {
//'  cat("This is a message emitted because we are in debug mode.\n")
//' }
//' @export
// [[Rcpp::export]]
bool FitGenGetDebug()
{
 // We return the state of the "normal" debugging flag. R functions really do not need to know the state of the low level debug flag,
 // which is mainly for internal debugging.
 return(DEB!=0);
}

