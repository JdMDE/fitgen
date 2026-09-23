#' @importFrom stats rnbinom
NULL

#' rmatNB2
#' @param N Number of features (genes)
#' @param x Model matrix
#' @param common_beta Coefficients correspondings to \code{x} i.e.
#' the common predictors for all rows
#' @param neighbours_beta Coefficients of neighbours
#' @param phi The dispersion parameters (Nx1)
#' @param seed The seed for the random generator
#' @export
#' 
rmatNB2 = function(N,x,common_beta,neighbours_beta=NULL,phi,
                   seed=1234){
    set.seed(seed)
    n = nrow(x)
    mu0 = x %*% common_beta
    ## Independent rows
    y1 = t(sapply(1:N,function(i) stats::rnbinom(length(mu0),mu=mu0,
                                                 size=1/phi[i])))

    ## Neighbours
    y2 = NULL
    if(!is.null(neighbours_beta)){
        neighbours = t(replicate(N,sample(N,length(neighbours_beta))))
        
        ## Dependent using neighbours
        y2 = t(sapply(1:N,function(i)
            rnbinom(length(mu0),mu=mu0 +t(y1[neighbours[i,],])%*%
                                    neighbours_beta,size=1/phi[i])))
    }
    list(y1=y1,y2=y2,neighbours=neighbours)
}
