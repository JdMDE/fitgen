#' @importFrom  methods new
#' @importFrom stats coefficients fitted pnorm pchisq
NULL

#' Find neighbours for each row of a matrix
#' 
#' Identifies a fixed number of neighbouring rows for every row of a
#' numeric matrix. The computation is performed using the binary matrix
#' utilities provided by the \pkg{scellpam} package.
#' 
#' @param y A numeric \code{matrix}. Rows represent the objects for which
#' neighbours are to be identified and columns represent observations 
#' or features.
#' @param neighbours_n A non-negative integer giving the maximum number of
#' neighbours to return for each row.
#' @param neighbours_method Character string specifying the method used to
#' select neighbours. Available methods are \code{"trivial"}, \code{"absvalue"} 
#' and \code{"FDR"}.
#' Currently, only \code{"trivial"} is implemented.
#' @return
#' A \code{matrix} with one row for each row of \code{y}. Each row contains the
#' indices of the rows of \code{y} selected as neighbours. The number of columns
#' is determined by \code{neighbours_n}.
#' 
#' @seealso 
#' \code{\link{glmRows}}, \code{\link{NB2Rows}},\code{\link{FitRows}}.
#' 
#' @export
find_neighbours = function(y,neighbours_n=3,
                           neighbours_method=c("trivial","absvalue","FDR")){
    neighbours_method = match.arg(neighbours_method)

    ## y  
    if(!is.matrix(y)) stop("y is not a matrix \n")
    if(!is.numeric(y)) stop("y should be numeric \n")

    ## neighbours_n  
    if (!is.numeric(neighbours_n) ||
        length(neighbours_n) != 1L ||
        is.na(neighbours_n) ||
        !is.finite(neighbours_n) ||
        neighbours_n < 0 ||
        neighbours_n != floor(neighbours_n))
        stop("neighbours_n should be a single non-negative integer")

    if (neighbours_n >= nrow(y))
        stop("neighbours_n should be smaller than the number of rows of y")
    
    ## Evaluating the neighbours
    countbinfile=tempfile(pattern="Counts_full",fileext=".bin")   
    scellpam::JWriteBin(y,countbinfile,dtype="float",
                        dmtype="full",comment=" ")
    ## Evaluating the closest cases up to a maximum of neighbours_n
    if (fitgen::FitGenGetDebug())
        cat("Evaluating the closest cases up to a maximum of",neighbours_n,"\n")
    neighbours = scellpam::ClosestCases(countbinfile,q=neighbours_n,
                                        method=neighbours_method)
    on.exit(unlink(countbinfile), add = TRUE)
    neighbours
}


#' Fit a generalized linear model to a single response vector
#' 
#' @description
#' Fits a generalized linear model to a single response vector. Currently,
#' binomial models with a logit link and Poisson models with a log link are
#' supported.
#' 
#' @param x A numeric model \code{matrix} with one row per observation and 
#' one column per regression coefficient.
#' @param y A numeric response \code{vector}. Its length must be equal to
#' \code{nrow(x)}.
#' @param family A description of the random component and link function
#' to be used in the model. It can be specified as a character string, a
#' family function or a family object. For example,
#' the Poisson family can be specified as \code{poisson()},
#' \code{poisson} or \code{"poisson"}. Currently, only
#' \code{binomial(link = "logit")} and \code{poisson(link = "log")} are
#' supported.
#' @param offset An optional numeric \code{vector} of offsets with the same
#' length as \code{y}. If \code{NULL}, a vector of zeros is used. 
#' @param is_fitted_values Logical. If \code{TRUE}, fitted values are
#' computed and included in the returned object.
#' @param is_unstable Logical. If \code{TRUE}, warnings associated with
#' unstable model fitting are issued.
#'  
#' @return
#' An object of class \code{fitgen_glm} containing the estimated
#' coefficients, model-fitting information and, if requested,
#' the fitted values.
#' 
#' @details
#' The model is fitted by \code{\link{StGLM}}. The resulting object is
#' validated and assigned class \code{fitgen_glm} by \code{new_fitgen_glm()}.
#' 
#' This function fits a model to a single response vector. To fit the same
#' model separately to the rows of a response matrix, use
#' \code{\link{glmRows}}.
#' 
#' @seealso
#' \code{\link{glmRows}}, \code{\link{StGLM}},
#' \code{\link[stats]{glm}}
#' 
#' @export
#' 
glmRow = function(x,y,family,offset=NULL,is_fitted_values=FALSE,
                  is_unstable=TRUE){
    ## x, y 
    if (!is.numeric(y) || !is.null(dim(y)))
        stop("y should be a numeric vector")
    if(!is.matrix(x)) stop("x should be a matrix")
    if(!is.numeric(x)) stop("x should be numeric")
    if (nrow(x) != length(y))
        stop("The number of rows of x should be equal to the length of y")

    ## family 
    if (is.character(family)) 
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family = family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    ## offset
    if(is.null(offset)){
        offset = rep(0,length(y))
    } else {
        if(!is.numeric(offset) || !is.null(dim(offset)))
            stop("offset should be a numeric vector")
        if(!(length(y) == length(offset)))
            stop("offset and y should have the same length")
        if (anyNA(offset) || any(!is.finite(offset)))
            stop("offset should contain finite values")
    }
    
    ## is_fitted_values  
    if (!is.logical(is_fitted_values) ||
        length(is_fitted_values) != 1L ||
        is.na(is_fitted_values))
        stop("is_fitted_values should be TRUE or FALSE")
    
    ## is_unstable  
    if (!is.logical(is_unstable) ||
        length(is_unstable) != 1L ||
        is.na(is_unstable))
        stop("is_unstable should be TRUE or FALSE")

    ## Fitting the models  
    model1 = (family$family == "binomial" & family$link == "logit")
    if(model1)
        fit  = fitgen::StGLM(y=y,X=x,model="binlogit",offset=offset,
                             warn_me=is_unstable,retmu=is_fitted_values)

    
    model2 = (family$family == "poisson" & family$link == "log")
    if(model2)
        fit  = fitgen::StGLM(y=y,X=x,model="poisson",offset=offset,
                             warn_me=is_unstable,retmu=is_fitted_values)
    
    if(!any(c(model1,model2))) stop("Not implemented model \n")
    new_fitgen_glm(fit)
}

#' Fit generalized linear models by rows
#' @description
#' Fits a generalized linear model separately to each selected row of a
#' response matrix. Each row of \code{y} is treated as a response vector.
#' Currently, binomial models with a logit link and 
#' Poisson models with a log link are supported.
#' 
#' The models share a common model matrix and may additionally include
#' row-specific predictors obtained from neighbouring rows.
#' 
#' @param x A numeric model \code{matrix} common to all fitted models, with
#' one row per observation and one column per common regression coefficient.
#' Its number of rows must be equal to \code{ncol(y)}.
#' @param y A numeric response \code{matrix}  whose rows are response vectors
#' and whose columns correspond to observations.
#' @param family A description of the random component and link function
#' to be used in the models. It can be specified as a character string, a
#' family function or a family object. For example,
#' the Poisson family can be specified as \code{poisson()},
#' \code{poisson} or \code{"poisson"}. Currently, only
#' \code{binomial(link = "logit")} and \code{poisson(link = "log")} are
#' supported.
#' @param is_normalized Logical. Indicates whether the response data have
#' already been normalized.
#' @param selected A logical \code{vector} with one element per row of
#' \code{y}. Rows corresponding to \code{TRUE} are fitted.
#' @param is_fitted_values Logical. If \code{TRUE}, fitted values are
#' computed and included in the returned object.
#' @param neighbours An optional integer \code{matrix} with one row per row
#' of \code{y}. Each row contains the indices of the neighbouring rows to
#' be included as additional predictors.
#' @param neighbours_n A non-negative integer giving the number of
#' neighbours. It is used to compute the neighbourhood matrix when
#' \code{neighbours = NULL}.
#' @param neighbours_method Character string specifying the method used to
#' select neighbours. Currently, only \code{"trivial"} is implemented.
#' @return
#' An object of class \code{fitgen_Lglm} containing one fitted generalized
#' linear model for each selected row of \code{y}. The common predictors
#' are given by \code{x}. When neighbours are used, the responses of the
#' corresponding neighbouring rows are included as additional predictors.
#' @export
#' 
glmRows = function(x,y,family,
                       is_normalized=TRUE,
                       selected =rep(TRUE,nrow(y)),
                       is_fitted_values=FALSE,
                       neighbours = NULL,
                       neighbours_n = 0,
                       neighbours_method="trivial"){

    neighbours_method = match.arg(neighbours_method)

    ## x, y 
    if(!is.matrix(x)) stop("x should be a matrix")
    if(!is.numeric(x)) stop("x should be numeric")
    if(!is.matrix(y)) stop("y should be a matrix")
    if(!is.numeric(y)) stop("y should be numeric")
    if(!(nrow(x) == ncol(y))) stop("Number of rows of the model
    matrix should be equal to the number of columns of y")
    
    ## family  
    if (is.character(family)) 
        family = get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family = family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    ## is_normalized
    if (!is.logical(is_normalized) ||
        length(is_normalized) != 1L ||
        is.na(is_normalized))
        stop("is_normalized should be TRUE or FALSE")

    ## selected  
    if(!is.logical(selected)) stop("selected should be logical")
    if(!(length(selected) == nrow(y)))
        stop("The length of selected should be equal to the number of
               rows of y")
    if (anyNA(selected))
        stop("selected contains NAs")
    
    ## is_fitted_values  
    if (!is.logical(is_fitted_values) ||
        length(is_fitted_values) != 1L ||
        is.na(is_fitted_values))
        stop("is_fitted_values should be TRUE or FALSE")

    ## neighbours_n
    if(!is.numeric(neighbours_n) ||
       length(neighbours_n) != 1L ||
       is.na(neighbours_n) ||
       !is.finite(neighbours_n) ||
       neighbours_n < 0 ||
       neighbours_n != floor(neighbours_n))
        stop("neighbours_n should be a single non-negative integer")
    
    ## neighbours  
    if (!is.null(neighbours)) {
        if (!is.matrix(neighbours))
            stop("neighbours should be a matrix")
        if (!is.numeric(neighbours))
            stop("neighbours should be numeric")
        if (anyNA(neighbours) ||
            any(!is.finite(neighbours)) ||
            any(neighbours != floor(neighbours))) 
            stop("neighbours should contain finite integer indices")
        if (any(neighbours < 1L | neighbours > nrow(y)))
            stop("neighbours contains indices outside the rows of y")
        neighbours_n = ncol(neighbours)
        if(any(neighbours == row(neighbours))) 
            stop("a row cannot be included among its own neighbours")
    }
    if(is.null(neighbours) & (neighbours_n >0)){
        if (fitgen::FitGenGetDebug())
            cat("Evaluating the neighbours \n")
        neighbours = find_neighbours(y=y,neighbours_n=neighbours_n,
                                     neighbours_method=neighbours_method)
    }
    if (!is.null(neighbours) && nrow(neighbours) != nrow(y))
        stop("The number of rows of neighbours should be equal to the
              number of rows of y")

    model1 = (family$family == "binomial" & family$link == "logit")
    if(model1 & neighbours_n >0){
        if (fitgen::FitGenGetDebug())
        cat("Fitting logistic regression by rows with neighbours\n")
        fits = fitgen::FitGLMWithNei(dat=y,
                            normalized=is_normalized,        
                            Cl=neighbours,
                            model="binlogit",                     
                            selected=selected,
                            commonX = x,
                            Phi_estproc="none",
                            retmu=is_fitted_values)
    }
    if(model1 & neighbours_n == 0){
        if (fitgen::FitGenGetDebug())
            cat("Fitting logistic regression by rows without neighbours \n")
        fits = fitgen::FitGLMInd(dat=y,
                                 normalized=is_normalized,        
                                 model="binlogit",                     
                                 selected=selected,
                                 commonX = x,
                                 Phi_estproc="none",
                                 retmu=is_fitted_values)
    }
    
    model2 = (family$family == "poisson" & family$link == "log")
    
    if(!any(c(model1,model2)))
      stop("Not implemented model")
    
    if(model2 & neighbours_n >0){
        if (fitgen::FitGenGetDebug())
            cat("Fitting Poisson regression by rows with neighbours \n")
        fits = fitgen::FitGLMWithNei(dat=y,
                                     normalized=is_normalized,        
                                     Cl=neighbours,
                                     model="poisson",                     
                                     selected=selected,
                                     commonX = x,
                                     Phi_estproc="none",
                                     retmu=is_fitted_values)
    }

    if(model2 & neighbours_n == 0){
        if (fitgen::FitGenGetDebug())
            cat("Fitting Poisson regression by rows without neighbours \n")
        fits = fitgen::FitGLMInd(dat=y,
                                 normalized=is_normalized,        
                                 model="poisson",                     
                                 selected=selected,
                                 commonX = x,
                                 Phi_estproc="none",
                                 retmu=is_fitted_values)
    }
    class(fits) = "fitgen_Lglm"
    fits
}



#' Fit a negative binomial regression to a single response vector
#' 
#' @description
#' Fits a negative binomial regression model to a single response vector
#' using iteratively reweighted least squares (IRLS) and a supplied dispersion
#' parameter
#' 
#' @param x A numeric model \code{matrix} with one row per observation and
#' one column per regression coefficient. 
#' @param y A numeric response \code{vector}. Its length must be equal to
#' \code{nrow(x)}.
#' @param phi A positive numeric value giving the dispersion parameter of
#' the negative binomial model.
#' @param offset An optional numeric \code{vector} of offsets with the same
#' length as \code{y}. If \code{NULL}, a vector of zeros is used.
#' @param is_fitted_values Logical. If \code{TRUE}, fitted values are
#' computed and included in the returned object.
#' @param is_unstable Logical. If \code{TRUE}, warnings associated with
#' unstable model fitting are issued.
#' 
#' @return
#' An object of class \code{fitgen_glm} containing the estimated
#' coefficients, model-fitting information and, if requested, the fitted
#' values.
#' 
#' @details
#' The model is fitted by \code{\link{StGLM}} using IRLS,
#' \code{model = "negbin"} and the dispersion parameter supplied in
#' \code{phi}. The resulting object is validated and assigned class
#' \code{fitgen_glm} by \code{new_fitgen_glm()}.
#' This function fits a model to a single response vector. To fit negative
#' binomial models separately to the rows of a response matrix, use
#' \code{\link{NB2Rows}}.
#' 
#' @seealso 
#' \code{\link{NB2Rows}}, \code{\link{NB2_phi}},
#' \code{\link{StGLM}}
#' 
#' @export
#' 
NB2Row = function(x,y,phi,offset=NULL,is_fitted_values=FALSE,
                  is_unstable=TRUE){

    ## x, y
    if(!is.matrix(x)) stop("x should be a matrix")
    if(!is.numeric(x)) stop("x should be numeric")
    if (!is.numeric(y) || !is.null(dim(y)))
        stop("y should be a numeric vector")
    if (nrow(x) != length(y))
        stop("The number of rows of x should be equal to the length of y")  

    ## phi  
    if(!(is.numeric(phi) && length(phi) == 1L &&
         !is.na(phi) && is.finite(phi) && phi > 0)) 
        stop("phi should be a single, finite, strictly positive number.")

    ## offset  
    if(!is.null(offset) & !is.numeric(offset))
        stop("offset should be numeric")
    if(!is.null(offset) & !(length(y) == length(offset)))
        stop("offset and y should have the same length")
    
    ## is_fitted_values  
    if (!is.logical(is_fitted_values) ||
        length(is_fitted_values) != 1L ||
        is.na(is_fitted_values))
        stop("is_fitted_values should be TRUE or FALSE")
    
    ## is_unstable  
    if (!is.logical(is_unstable) ||
        length(is_unstable) != 1L ||
        is.na(is_unstable))
        stop("is_unstable should be TRUE or FALSE")
    

    ## offset  
    if(is.null(offset)){
        fit  = fitgen::StGLM(y=y,X=x,model="negbin",offset=rep(0,length(y)),
                             warn_me=is_unstable,retmu=is_fitted_values,
                             Phi=phi)
    } else {
        fit  = fitgen::StGLM(y=y,X=x,model="negbin",offset=offset,
                             warn_me=is_unstable,retmu=is_fitted_values,
                             Phi=phi)
    }
    new_fitgen_glm(fit)
    
}

#' Estimate negative binomial dispersion parameters
#'  
#' @description
#' Estimates negative binomial dispersion parameters for selected rows of
#' a numeric data matrix using procedures provided by the \pkg{edgeR} package.
#' The input matrix is assumed to have been previously preprocessed and
#' prepared for statistical modelling.
#' 
#' @param y A numeric count \code{matrix} whose rows represent features and
#' whose columns represent observations. 
#' @param method Character string specifying the dispersion estimation
#' method. Available methods are \code{"common_edger"}, which estimates a
#' common dispersion using \code{\link[edgeR]{estimateCommonDisp}}, and
#' \code{"tagwise_edger"}, which estimates row-specific dispersions using
#' \code{\link[edgeR]{estimateTagwiseDisp}}.
#' @param selected A logical \code{vector} with one element per row of
#' \code{y}. A value of \code{TRUE} indicates that the corresponding row
#' is to be used for dispersion estimation. At least one row must be
#' selected.
#' @return
#' A numeric \code{vector} with one element per row of \code{y}. Estimated
#' dispersion parameters are returned for the selected rows, whereas
#' unselected rows contain \code{NA}. When \code{method = "common_edger"}, 
#' the common dispersion parameter is assigned to every selected row.
#' 
#' @seealso 
#' \code{\link{NB2Row}}, \code{\link{NB2Rows}},
#' \code{\link[edgeR]{estimateCommonDisp}},
#' \code{\link[edgeR]{estimateTagwiseDisp}}
#' 
#' @export
#'
NB2_phi = function(y,method=c("common_edger","tagwise_edger"),selected){
    method = match.arg(method)
    ## y
    if(!is.matrix(y)) stop("y should be a matrix")
    if(!is.numeric(y)) stop("y should be numeric")
    
    ## selected  
    if (!is.logical(selected))
        stop("selected must be a logical vector.")
    if(length(selected) != nrow(y))
        stop("The length of selected should be equal to the number of rows of y")
    if (anyNA(selected))    
      stop("selected contains NAs")
    if (!any(selected))
        stop("At least one row in y must be selected.")
    selected = which(selected)
    
    y_selected = y[selected, , drop = FALSE]
    if (anyNA(y_selected) || any(!is.finite(y_selected)))
        stop("Selected rows of y should contain finite values")
    if (any(y_selected < 0))
        stop("Selected rows of y should contain non-negative counts")
    
    ret = rep(NA_real_,nrow(y))
    
    switch(method,
           "common_edger" = {
               dge = edgeR::DGEList(counts=y_selected)
               dge.c = edgeR::estimateCommonDisp(dge)
               phi = dge.c$common.dispersion
               if (fitgen::FitGenGetDebug())
                   cat("Common phi value estimated per row using ",
                       method,"method is",phi,"\n")
           },
           "tagwise_edger" = {
               dge = edgeR::DGEList(counts=y_selected)
               dge.c = edgeR::estimateCommonDisp(dge)
               dge.t = edgeR::estimateTagwiseDisp(dge.c)
               phi = dge.t$tagwise.dispersion
               if (fitgen::FitGenGetDebug())
                   cat("phi values estimated per rows using ",method,
                       "method. Their mean is",mean(phi),"\n")
           }
           )
    ret[selected] = phi
    ret 
}


#' Fit negative binomial models by rows 
#' 
#' @description
#' Fits a negative binomial regression model separately to each selected
#' row of a response matrix. Each row of \code{y} is treated as a response
#' vector.
#' 
#' All models share a common model matrix \code{x}. They may additionally
#' include row-specific predictors obtained from neighbouring rows of
#' \code{y}.
#' 
#' @param x A numeric model \code{matrix} common to all fitted models, with
#' one row per observation and one column per common regression coefficient.
#' Its number of rows must be equal to \code{ncol(y)}.
#' @param y A numeric count \code{matrix} whose rows are response vectors
#' and whose columns correspond to observations.
#' @param phi An optional numeric \code{vector} containing externally
#' supplied negative binomial overdispersion parameters. 
#'  If \code{NULL}, the overdispersion parameters are obtained according to
#' \code{phi_method}.
#' @param phi_method Character string specifying the procedure used to
#' obtain the negative binomial overdispersion parameters when
#' \code{phi = NULL}. See \code{Details}.
#' @param is_normalized Logical. Indicates whether the count data have
#' already been normalized.
#' @param selected A logical \code{vector} with one element per row of
#' \code{y}. Rows corresponding to \code{TRUE} are fitted.
#' @param is_fitted_values Logical. If \code{TRUE}, fitted values are
#' computed and included in the returned model objects.
#' @param neighbours_method Character string specifying the method used to
#' select neighbours. Currently, only \code{"trivial"} is implemented.
#' @param neighbours An optional integer \code{matrix} with one row per row
#' of \code{y}. Each row contains the indices of the neighbouring rows to
#' be included as additional predictors.
#' @param neighbours_n A non-negative integer giving the number of
#' neighbours. It is used to compute the neighbours when
#' \code{neighbours = NULL}.
#' 
#' @return
#' An object of class \code{fitgen_Lglm} containing one fitted negative
#' binomial model for each selected row of \code{y}. The common predictors
#' are given by \code{x}. When neighbours are used, the responses of the
#' corresponding neighbouring rows are included as additional predictors. 
#' 
#' @details
#' The negative binomial model uses the NB2 parametrization, in which the
#' variance increases quadratically with the mean. The parameter
#' \code{phi} controls the amount of overdispersion.
#' 
#' When \code{phi = NULL}, the available overdispersion estimation
#' procedures are:
#' \itemize{
#'   \item \code{"none"} No estimation is performed. This option is used
#'    when Poisson regression is used
#'   \item \code{"single"}: computes a separate moment-based estimate
#'   of the overdispersion parameter for each row of
#'   \code{y}.
#'   \item \code{"common"}: estimates a single common overdispersion
#'   parameter using the Bliss--Owen procedure and the rows for which
#'   estimation is valid.
#'   \item \code{"grouped"}: obtains row-specific overdispersion estimates
#'   by applying the common Bliss--Owen procedure locally to groups of
#'   rows with similar dispersion profiles.
#'   \item \code{"common_edger"}: estimates a common overdispersion
#'   parameter using \pkg{edgeR}.
#'   \item \code{"tagwise_edger"}: estimates row-specific overdispersion
#'   parameters using \pkg{edgeR}.
#' }
#' 
#' The procedures  \code{"single"}, \code{"common"} and
#' \code{"grouped"} are handled internally by the row-wise fitting
#' routines. For \code{"common_edger"} and \code{"tagwise_edger"}, the
#' overdispersion parameters are first estimated by
#' \code{\link{NB2_phi}} and then supplied to the fitting routines.
#' If \code{phi_method} is set to \code{"common_edger"} or
#' \code{"tagwise_edger"} then the provided \code{phi} are used. Otherwise
#' \code{phi} is ignored and \code{phi_method} is used internally.
#'  
#' Models without neighbours are fitted by
#' \code{\link[fitgen]{FitGLMInd}}, whereas models with neighbours are
#' fitted by \code{\link[fitgen]{FitGLMWithNei}}.
#' 
#' @seealso
#' \code{\link{NB2Row}}, \code{\link{NB2_phi}},
#' \code{\link{find_neighbours}},
#' \code{\link[fitgen]{BlissOwen}},
#' \code{\link[fitgen]{CommonBlissOwen}},
#' \code{\link[fitgen]{GroupedCommonBlissOwen}}
#' 
#' @export
#' 
NB2Rows  =  function(x,y,phi = NULL,
                    phi_method=c("single","common","grouped",
                                 "common_edger","tagwise_edger"),
                     is_normalized=TRUE,
                     selected,
                     is_fitted_values=FALSE,
                     neighbours_method = "trivial",
                     neighbours = NULL,
                    neighbours_n = 0){
    phi_method = match.arg(phi_method)
    neighbours_method = match.arg(neighbours_method)

    ## x, y 
    if(!is.matrix(x)) stop("x should be a matrix")
    if(!is.numeric(x)) stop("x should be numeric")
    if(!is.matrix(y)) stop("y should be a matrix")
    if(!is.numeric(y)) stop("y should be numeric")
    if(!(nrow(x) == ncol(y))) stop("Number of rows of the model
                      matrix should be equal to the number of columns of y")

    ## is_normalized
    if (!is.logical(is_normalized) ||
        length(is_normalized) != 1L ||
        is.na(is_normalized))
        stop("is_normalized should be TRUE or FALSE")
    
    ## selected  
    if(!is.logical(selected)) stop("selected should be logical")
    if(!(length(selected) == nrow(y)))
        stop("The length of selected should be equal to the number of
               rows of y")
    if (anyNA(selected))
        stop("selected contains NAs")
    if (!any(selected))
        stop("At least one row of y should be selected")

    ## is_fitted_values  
    if (!is.logical(is_fitted_values) ||
        length(is_fitted_values) != 1L ||
        is.na(is_fitted_values))
        stop("is_fitted_values should be TRUE or FALSE")
    
    ## phi  
    if (!is.null(phi)) {
        if (!is.numeric(phi) || length(phi) == 0L) {
            stop("phi should be a non-empty numeric vector.")
        }
        if (!(length(phi) %in% c(1L, nrow(y)))) {
            stop("phi should have length one or one element per row of y.")
        }
        phi_selected = if (length(phi) == 1L) {
                           phi
                       } else {
                           phi[selected]
                       }
        if (anyNA(phi_selected) ||
            any(!is.finite(phi_selected)) ||
            any(phi_selected <= 0)) {
            stop(
                paste(
                    "phi should contain finite, strictly positive",
                    "values for all selected rows."
                )
            )
        }
    }
    ## neighbours  
    if (!is.null(neighbours)) {
        if (!is.matrix(neighbours))
            stop("neighbours should be a matrix")
        if (!is.numeric(neighbours))
            stop("neighbours should be numeric")
        if (anyNA(neighbours) ||
            any(!is.finite(neighbours)) ||
            any(neighbours != floor(neighbours))) 
            stop("neighbours should contain finite integer indices")
        if (any(neighbours < 1L | neighbours > nrow(y)))
            stop("neighbours contains indices outside the rows of y")
        neighbours_n = ncol(neighbours)
        if(any(neighbours == row(neighbours))) 
            stop("a row cannot be included among its own neighbours")
    }
    ## neighbours_n 
    if (!is.numeric(neighbours_n) ||
        length(neighbours_n) != 1L ||
        is.na(neighbours_n) ||
        !is.finite(neighbours_n) ||
        neighbours_n < 0 ||
        neighbours_n != floor(neighbours_n))
      stop("neighbours_n should be a single non-negative integer")
    
    
    if(is.null(neighbours) & (neighbours_n >0)){
        if (fitgen::FitGenGetDebug())
            cat("Evaluating the neighbours \n")
        neighbours = find_neighbours(y=y,neighbours_n=neighbours_n,
                                     neighbours_method=neighbours_method)
    }
    if (!is.null(neighbours) && nrow(neighbours) != nrow(y))
        stop("The number of rows of neighbours should be equal to the
               number of rows of y")


    ## phi  
    if(is.null(phi) & is.element(phi_method,c("common_edger","tagwise_edger"))){
        if (fitgen::FitGenGetDebug())
            cat("Estimating phi \n")
        phi = NB2_phi(y=y,method=phi_method,selected=selected)
    }

    if(neighbours_n > 0L && is.null(neighbours))
        stop("neighbours does not exists and neighbours_n is positive")
    
    if(neighbours_n > 0L){
        if (fitgen::FitGenGetDebug())
            cat("Fitting negative binomial regression by rows with
                 neighbours. \n")
        fits = fitgen::FitGLMWithNei(dat=y,
                                     normalized=is_normalized,        
                                     Cl=neighbours,
                                     model="negbin",                     
                                     selected=selected,
                                     commonX = x,
                                     Phi_estproc=phi_method,
                                     Phi_external = phi,
                                     retmu=is_fitted_values)
    } else {
        if (fitgen::FitGenGetDebug())
            cat("Fitting negative binomial regression by rows without
                 neighbours. \n")
        fits = fitgen::FitGLMInd(dat=y,
                                 normalized=is_normalized,        
                                 model="negbin",                     
                                 selected=selected,
                                 commonX = x,
                                 Phi_estproc=phi_method,
                                 Phi_external = phi,
                                 retmu=is_fitted_values)
    }
    
    class(fits) = "fitgen_Lglm"
    fits
}

#' @importFrom matrixStats rowMins
#' @importFrom stats p.adjust
NULL

#' Fit Poisson and negative binomial models by rows
#' 
#' @description
#' Fits a Poisson regression model and a negative binomial regression
#' model separately to each row of a response matrix. Both models are
#' fitted to every row of \code{y}.
#' 
#' All row-wise models share the common model matrix \code{x}. They may
#' additionally include row-specific predictors obtained from neighbouring
#' rows of \code{y}.
#' 
#' @param x A numeric model \code{matrix} common to all fitted models, with
#' one row per observation and one column per common regression coefficient.
#' Its number of rows must be equal to \code{ncol(y)}.
#' @param y A numeric response \code{matrix} whose rows are response vectors
#' and whose columns correspond to observations.
#' @param phi A numeric \code{vector} containing the negative
#' binomial dispersion parameters. It may have length one or one element
#' per row of \code{y}.
#' @param is_normalized Logical. Indicates whether the response data have
#' already been normalized.
#' @param is_fitted_values Logical. If \code{TRUE}, fitted values are
#' computed and included in the returned model objects.
#' @param neighbours_method Character string specifying the method used to
#' select neighbours. Currently, only \code{"trivial"} is implemented.
#' @param neighbours An optional integer \code{matrix} with one row per row
#' of \code{y}. Each row contains the indices of the neighbouring rows to
#' be included as additional predictors.
#' @param neighbours_n A non-negative integer giving the number of
#' neighbours. It is used to compute the neighbourhood matrix when
#' \code{neighbours = NULL}.
#' @param verbose Logical. If \code{TRUE}, progress messages are printed.
#' @return
#' A \code{list} with three components:
#' \itemize{
#'   \item \code{Poi}: an object of class \code{fitgen_Lglm} containing
#'   the row-wise Poisson fits;
#'   \item \code{NB2}: an object of class \code{fitgen_Lglm} containing
#'   the row-wise negative binomial fits;
#'   \item \code{neighbours}: the neighbourhood matrix used in both fits,
#'   or \code{NULL} when no neighbours are used.
#' }
#' 
#' @export
#' 
FitRows = function(x,y,phi,
                   is_normalized=TRUE,
                   is_fitted_values=FALSE,
                   neighbours_method="trivial",
                   neighbours = NULL,
                   neighbours_n = 0,
                   verbose = FALSE){
    neighbours_method = match.arg(neighbours_method)
    ## x, y  
    if(!is.matrix(y)) stop("y should be a matrix")
    if(!is.numeric(y)) stop("y should be numeric")
    if(!is.matrix(x)) stop("x should be a matrix")
    if(!is.numeric(x)) stop("x should be numeric")
    if(!(nrow(x) == ncol(y))) stop("Number of rows of the model
         matrix should be equal to the number of columns of y")

    ## phi  
    if (!is.numeric(phi) ||
        length(phi) == 0L ||
        anyNA(phi) ||
        any(!is.finite(phi)) ||
        any(phi <= 0)) {
        stop("phi should contain finite, strictly positive numeric values.")
    }
    if (!(length(phi) %in% c(1L, nrow(y)))) 
        stop("phi should have length one or one element per row of y.")

    ## is_normalized
    if (!is.logical(is_normalized) ||
        length(is_normalized) != 1L ||
        is.na(is_normalized))
        stop("is_normalized should be TRUE or FALSE")

    ## is_fitted_values  
    if (!is.logical(is_fitted_values) ||
        length(is_fitted_values) != 1L ||
        is.na(is_fitted_values))
        stop("is_fitted_values should be TRUE or FALSE")
    
    ## neighbours  
    if (!is.null(neighbours)) {
        if (!is.matrix(neighbours))
            stop("neighbours should be a matrix")
        if (!is.numeric(neighbours))
            stop("neighbours should be numeric")
        if (anyNA(neighbours) ||
            any(!is.finite(neighbours)) ||
            any(neighbours != floor(neighbours))) 
            stop("neighbours should contain finite integer indices")
        if (any(neighbours < 1L | neighbours > nrow(y)))
            stop("neighbours contains indices outside the rows of y")
        neighbours_n = ncol(neighbours)
        if(any(neighbours == row(neighbours))) 
            stop("a row cannot be included among its own neighbours")
    }
    
    ## neigbours.n
    if (!is.numeric(neighbours_n) ||
        length(neighbours_n) != 1L ||
        is.na(neighbours_n) ||
        !is.finite(neighbours_n) ||
        neighbours_n < 0 ||
        neighbours_n != floor(neighbours_n))
      stop("neighbours_n should be a single non-negative integer")
    
    if(is.null(neighbours) & (neighbours_n >0)){
        if (fitgen::FitGenGetDebug())
            cat("Evaluating the neighbours \n")
        neighbours = find_neighbours(y=y,neighbours_n=neighbours_n,
                                     neighbours_method=neighbours_method)
    }
    if (!is.null(neighbours) && nrow(neighbours) != nrow(y))
        stop("The number of rows of neighbours should be equal to the
              number of rows of y")
    
    ## Fitting Poisson models
    if(verbose)  cat("Fitting the Poisson models.\n")

    ## Fit Poisson models to all rows
    LPoi = glmRows(x=x,y=y,family="poisson",
                   is_normalized=is_normalized,
                   selected = rep(TRUE,nrow(y)), 
                   is_fitted_values=is_fitted_values,
                   neighbours = neighbours,
                   neighbours_n = neighbours_n,
                   neighbours_method = neighbours_method)
    
    if(verbose) cat("Done  \n")

    ## Fit negative binomial models to all rows
    if (verbose)
        cat("Fitting the Negative binomial models. \n")
    
    LNB2 = NB2Rows(x=x,
                   y=y,
                   phi=phi,
                   is_normalized=is_normalized,
                   selected = rep(TRUE,nrow(y)),
                   is_fitted_values=is_fitted_values,
                   neighbours_method = neighbours_method,
                   neighbours = neighbours,
                   neighbours_n=neighbours_n)
    if(verbose) cat("Done \n")
    ## Return both sets of fitted models and the neighbourhood matrix used
    list(Poi=LPoi,NB2=LNB2,neighbours=neighbours)
}

#' Score tests for departures from the Poisson variance
#' @description
#' Computes score-type diagnostics for linear or quadratic departures from
#' the Poisson variance assumption. The linear statistic coincides with the
#' standardized score test proposed by Dean and Lawless (1989).
#' 
#' @param y A numeric \code{vector} containing the observed responses.
#' @param muest A numeric \code{vector} containing the fitted means under
#' the Poisson model. It must have the same length as \code{y}.
#' @param type Character string specifying the form of the variance
#' expansion. Available values are \code{"linear"} and
#' \code{"quadratic"}.
#' @param alternative Character string specifying the alternative
#' hypothesis. Available values are \code{"greater"}, \code{"less"}
#' and \code{"two.sided"}.
#' 
#' @return
#' A \code{list} with two components:
#' \itemize{
#'   \item \code{T}: the value of the score statistic;
#'   \item \code{p_value}: the corresponding p-value.
#' }
#' 
#' @details
#' Let \code{mu} denote the fitted Poisson mean. The diagnostics assess
#' whether the conditional variance departs from the variance implied by
#' the Poisson model.
#' 
#' For \code{type = "linear"}, the departure term is proportional to the 
#' mean. The resulting statistic coincides with the
#' standardized score statistic considered by Dean and Lawless (1989).
#' 
#' For \code{type = "quadratic"}, the departure term is proportional to
#' the square of the mean.
#' 
#' Under the null hypothesis of Poisson equidispersion, both statistics
#' are asymptotically standard normal. The alternative
#' \code{"greater"} is used to detect overdispersion, whereas
#' \code{"less"} is used to detect underdispersion. The option
#' \code{"two.sided"} detects departures from the Poisson variance in
#' either direction.
#' 
#' These diagnostics are intended to characterize departures from the
#' Poisson variance assumption and need not be used as an automatic
#' model-selection rule.
#' 
#' @references 
#' Dean, C. and Lawless, J. F. (1989). Tests for detecting overdispersion
#' in Poisson regression models. \emph{Journal of the American Statistical
#' Association}, \bold{84}(406), 467--472. 
#' 
#' Cameron, A. C. and Trivedi, P. K. (2013).
#' \emph{Regression Analysis of Count Data}. Second edition.
#' Cambridge University Press.
#' 
#' @export
PoissonVarianceScore = function(y,muest,type= c("linear", "quadratic"),
                                alternative=c("greater", "less", "two.sided")){
    type = match.arg(type)
    alternative = match.arg(alternative)

    ## y, muest  
    if(!is.numeric(y)|| !is.null(dim(y))) stop("y should be numeric vector")
    if(!is.numeric(muest) || !is.null(dim(muest)))
        stop("muest should be numeric")
    if(length(y) != length(muest))
        stop("y and muest should have the same length")
    
    remains = which(is.finite(y) & is.finite(muest))
    if(length(remains) == 0L) stop("There is no valid data")
    y = y[remains]
    if(any(y < 0)) stop("y should contain non-negative values")
    muest = muest[remains]
    if(min(muest) <= 0) stop("muest should contain strictly positive values")
    switch(type,
           "linear" = {
               T = sum(((y-muest)^2 - y)/muest)/sqrt(2*length(y))
           },
           "quadratic" = {
               T = sum((y-muest)^2 - y) / sqrt(2*sum(muest^2))
           })
    switch(alternative,
           "less" = {
               p_value = pnorm(T)
           },
           "greater" = {
               p_value = 1 - pnorm(T)
           },
           "two.sided" = {
               p_value = 2*(1 - pnorm(abs(T)))
           })
    list(T=T,p_value=p_value)
}


#' Wald tests for regression coefficients in row-wise models
#' 
#' @description
#' Performs a Wald test of one or more selected regression coefficients
#' for each row-wise fitted model. The null hypothesis is that all
#' coefficients selected by \code{coef} are equal to zero.
#' 
#' Tests are computed under the Poisson, quasi-Poisson, NB2 and quasi-NB2
#' model specifications.
#' 
#' @param fits A \code{list} containing the row-wise fitted models. It must
#' have components \code{Poi} and \code{NB2}, such as an object returned
#' by \code{\link{FitRows}}.
#' @param coef An integer \code{vector} giving the positions of the
#' regression coefficients to be tested.
#' @param quasi_method Character string specifying the estimator of the
#' quasi-likelihood scale parameter. Available methods are
#' \code{"Fletcher"}, \code{"Wedderburn"} and \code{"Farrington"}.
#' @param is_neighbours Logical. If \code{FALSE}, the coefficients and
#' covariance matrix of the model without neighbours are used. If
#' \code{TRUE}, those of the neighbourhood-adjusted model are used.
#' @param verbose Logical. If \code{TRUE}, progress messages are printed.
#' 
#' @return
#' A named \code{list} with components \code{Poi}, \code{quasiPoi},
#' \code{NB2} and \code{quasiNB2}. Each component is a numeric
#' \code{vector} containing one p-value per fitted response. An
#' \code{NA} value is returned when a test cannot be computed.
#' 
#' @details
#' For each fitted model, the function tests the null hypothesis that the
#' coefficients selected by \code{coef} are equal to zero. If
#' \code{coef} contains more than one position, the corresponding
#' coefficients are tested jointly.
#' 
#' For the quasi-Poisson and quasi-NB2 specifications, the Wald statistic
#' is divided by the quasi-likelihood scale estimate selected by
#' \code{quasi_method}. The test is not computed when this estimate is
#' not finite or is not strictly positive.
#' 
#' When \code{is_neighbours = FALSE}, the function uses the components
#' \code{beta0} and \code{CovM0} of each fitted model. When
#' \code{is_neighbours = TRUE}, it uses the components \code{beta1} and
#' \code{CovM1}.
#' 
#' @seealso
#' \code{\link{FitRows}}, \code{\link{testRows_contrast}}
#' 
#' @export
#' 

testRows_coef = function(fits,coef,
                         quasi_method = c("Fletcher","Wedderburn","Farrington"),
                         is_neighbours=FALSE,verbose=FALSE){
    quasi_method = match.arg(quasi_method)
    if (!is.numeric(coef))
        stop("coef should be numeric")
    if (length(coef) == 0L)
      stop("coef should contain at least one index")
    if (anyNA(coef) ||
        any(!is.finite(coef)) ||
        any(coef != floor(coef))) 
        stop("coef should contain finite integer indices")
    if (any(coef < 1L))
      stop("coef should contain strictly positive indices")
    if (anyDuplicated(coef))
        stop("coef should not contain duplicated indices")
    
    ## is_neighbours
    if (!is.logical(is_neighbours) ||
        length(is_neighbours) != 1L ||
        is.na(is_neighbours))
        stop("is_neighbours should be TRUE or FALSE")

    ## verbose
    if (!is.logical(verbose) ||
        length(verbose) != 1L ||
        is.na(verbose))
        stop("verbose should be TRUE or FALSE")

    p1 = p2 = p3 = p4 = NULL
    if(!is_neighbours){
        if(verbose) cat("Poisson \n")
        p1 = sapply(fits$Poi,function(x){
            ifelse(length(x$beta0)>0 && max(coef) <= length(x$beta0),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta0[coef],
                                             solve(x$CovM0[coef,coef],x$beta0[coef]))),
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})

        if(verbose) cat("Quasi-Poisson \n")
        p2 = sapply(fits$Poi,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            ifelse(!is.na(quasidisp) && quasidisp > 0 && is.finite(quasidisp)
                   && length(x$beta0)>0 && max(coef) <= length(x$beta0),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta0[coef],
                                             solve(x$CovM0[coef,coef],x$beta0[coef])))/quasidisp,
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})
        
        if(verbose) cat("NB2 \n")
        p3 = sapply(fits$NB2,function(x){
            ifelse(length(x$beta0) >0 && max(coef) <= length(x$beta0),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta0[coef],
                                             solve(x$CovM0[coef,coef],x$beta0[coef]))),
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})
        
        if(verbose) cat("Quasi-NB2 \n")
        p4 = sapply(fits$NB2,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            ifelse(!is.na(quasidisp) && quasidisp > 0 && length(x$beta0) >0 &&
                   is.finite(quasidisp) && max(coef) <= length(x$beta0),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta0[coef],
                                             solve(x$CovM0[coef,coef],x$beta0[coef])))/quasidisp,
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})
    }

    if(is_neighbours){
        if(verbose) cat("Poisson \n")
        p1 = sapply(fits$Poi,function(x){
            ifelse(length(x$beta1) >0 && max(coef) <= length(x$beta1),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta1[coef],
                                             solve(x$CovM1[coef,coef],x$beta1[coef]))),
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})
        
        if(verbose) cat("Quasi-Poisson \n")
        p2 = sapply(fits$Poi,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            ifelse(!is.na(quasidisp) && quasidisp > 0 && is.finite(quasidisp) &&
                   length(x$beta1)>0 && max(coef) <= length(x$beta1),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta1[coef],
                                             solve(x$CovM1[coef,coef],x$beta1[coef])))/quasidisp,
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})

        if(verbose) cat("NB2 \n")
        p3 = sapply(fits$NB2,function(x)
            ifelse(length(x$beta1) > 0  && max(coef) <= length(x$beta1),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta1[coef],
                                             solve(x$CovM1[coef,coef],x$beta1[coef]))),
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA))
        
        if(verbose) cat("Quasi-NB2 \n")
        p4 = sapply(fits$NB2,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            ifelse(length(x$beta1) > 0 && !is.na(quasidisp) && quasidisp > 0 &&
                   is.finite(quasidisp) && max(coef) <= length(x$beta1),
                   tryCatch(
                       pchisq(drop(crossprod(x$beta1[coef],
                                             solve(x$CovM1[coef,coef],x$beta1[coef])))/quasidisp,
                              df=length(coef),lower.tail = FALSE),
                       error=function(e) NA,
                       warning=function(w) NA
                   ),NA)})
    }
    
    list(Poi = p1,quasiPoi = p2,NB2 = p3,quasiNB2 = p4)
}


#' Wald tests for linear contrasts in row-wise models
#'
#' @description
#' Performs a separate Wald test for each linear contrast specified in
#' \code{A} and for each row-wise fitted model.
#'
#' Tests are computed under the Poisson, quasi-Poisson, NB2 and quasi-NB2
#' model specifications.
#'
#' @param fits A \code{list} containing the row-wise fitted models. It must
#' have components \code{Poi} and \code{NB2}, such as an object returned
#' by \code{\link{FitRows}}.
#' @param A A numeric contrast \code{matrix} for the common predictors.
#' Each column defines a separate linear contrast. Its number of rows must
#' be equal to the number of common regression coefficients.
#' @param quasi_method Character string specifying the estimator of the
#' quasi-likelihood scale parameter. Available methods are
#' \code{"Fletcher"}, \code{"Wedderburn"} and \code{"Farrington"}.
#' @param is_neighbours Logical. If \code{FALSE}, the coefficients and
#' covariance matrix of the model without neighbours are used. If
#' \code{TRUE}, those of the neighbourhood-adjusted model are used.
#' @param neighbours_n A non-negative integer giving the number of
#' neighbours included in the neighbourhood-adjusted models. It is used
#' to extend the contrast matrix with zero coefficients for the
#' neighbouring-row predictors.
#' @param verbose Logical. If \code{TRUE}, progress messages are printed.
#'
#' @return
#' A named \code{list} with components \code{Poi}, \code{quasiPoi},
#' \code{NB2} and \code{quasiNB2}. Each component contains the p-values
#' for the corresponding model specification, with one row per fitted
#' response and one column per contrast. An \code{NA} value is returned
#' when a contrast cannot be computed.
#'
#' @details
#' Each column of \code{A} defines a separate test of whether the
#' corresponding linear combination of the common regression
#' coefficients is equal to zero.
#'
#' For the quasi-Poisson and quasi-NB2 specifications, the standard errors
#' are adjusted using the quasi-likelihood scale estimate selected by
#' \code{quasi_method}.
#'
#' When \code{is_neighbours = FALSE}, the function uses the components
#' \code{beta0} and \code{CovM0} of each fitted model. When
#' \code{is_neighbours = TRUE}, it uses the components \code{beta1} and
#' \code{CovM1}. In the latter case, the contrast matrix is extended with
#' zero rows so that the neighbouring-row coefficients are not included
#' in the contrasts.
#'
#' @seealso
#' \code{\link{FitRows}}, \code{\link{testRows_coef}}
#'
#' @export
#' 

testRows_contrast = function(fits,A,quasi_method = c("Fletcher","Wedderburn",
                                                     "Farrington"),
                             is_neighbours = FALSE,neighbours_n = 0,
                             verbose = FALSE){

    quasi_method = match.arg(quasi_method)

    if (!is.matrix(A))
        stop("A should be a matrix")
    if (!is.numeric(A))
        stop("A should be numeric")
    if (nrow(A) == 0L || ncol(A) == 0L)
        stop("A should have at least one row and one column")
    if (anyNA(A) || any(!is.finite(A)))
        stop("A should contain finite values")

    if (!is.logical(is_neighbours) ||
        length(is_neighbours) != 1L ||
        is.na(is_neighbours))
        stop("is_neighbours should be TRUE or FALSE")
    
    if (!is.numeric(neighbours_n) ||
        length(neighbours_n) != 1L ||
        is.na(neighbours_n) ||
        !is.finite(neighbours_n) ||
        neighbours_n < 0 ||
        neighbours_n != floor(neighbours_n))
        stop("neighbours_n should be a single non-negative integer")
    
    if (!is.logical(verbose) ||
        length(verbose) != 1L ||
        is.na(verbose))
        stop("verbose should be TRUE or FALSE")
    
    p1 = p2 = p3 = p4 = NULL
    
    if(!is_neighbours){
        if(verbose) cat("Poisson \n")
        p1 = t(sapply(fits$Poi,
                      function(x){
            if(length(x$beta0)>0){
                tryCatch(
                2*pnorm(-abs(t(A)%*%x$beta0/
                               sqrt(diag(t(A)%*% x$CovM0 %*%A)))),
                error=function(e) rep(NA,ncol(A)),
                warning=function(w) rep(NA,ncol(A))
                )
            } else {
                rep(NA,ncol(A))
            }}))
        
        if(verbose) cat("Quasi-Poisson \n")
        p2 = t(sapply(fits$Poi,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            if(!is.na(quasidisp) && length(x$beta0)>0)
                tryCatch(
                    2*pnorm(-abs(t(A)%*% x$beta0/
                                   (sqrt(quasidisp * diag(t(A)%*% x$CovM0 %*%A))))),
                    error=function(e) rep(NA,ncol(A)),
                    warning=function(w) rep(NA,ncol(A))
                )
            else
                rep(NA,ncol(A))
        }))
        
        if(verbose) cat("NB2 \n")
        p3 = t(sapply(fits$NB2,function(x){
            if(length(x$beta0)>0){
                tryCatch(
                    2*pnorm(-abs(t(A)%*%x$beta0/
                                   sqrt(diag(t(A)%*%x$CovM0%*%A)))),
                    error=function(e) rep(NA,ncol(A)),
                    warning=function(w) rep(NA,ncol(A))
                )
            } else {
                rep(NA,ncol(A))
                }}))
        
        if(verbose) cat("Quasi-NB2 \n")
        p4 = t(sapply(fits$NB2,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            if(!is.na(quasidisp) && length(x$beta0)>0)
                tryCatch(
                2*pnorm(-abs(t(A)%*%x$beta0/
                               (sqrt(quasidisp* diag(t(A)%*%x$CovM0%*%A))))),
                error=function(e) rep(NA,ncol(A)),
                warning=function(w) rep(NA,ncol(A))
                )
            else
                rep(NA,ncol(A))
        }))
    }

    if(is_neighbours){
        A1 = rbind(A,matrix(0,nrow=neighbours_n,ncol=ncol(A)))
        if(verbose) cat("Poisson \n")
        p1 = t(sapply(fits$Poi,function(x){
            if(length(x$beta1) == nrow(A1))
                tryCatch(
                2*pnorm(-abs(t(A1[1:length(x$beta1),,drop=FALSE])%*%x$beta1/
                               sqrt(diag(t(A1[1:length(x$beta1),,drop=FALSE])
                          %*%x$CovM1%*%A1[1:length(x$beta1),,drop=FALSE])))),
                error=function(e) rep(NA,ncol(A)),
                warning=function(w) rep(NA,ncol(A))
                )
            else
                rep(NA,ncol(A))
        }))
        
        if(verbose) cat("Quasi-Poisson \n")
        p2 = t(sapply(fits$Poi,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            if(!is.na(quasidisp) && length(x$beta1)>0)
                tryCatch(
                2*pnorm(-abs(t(A1[1:length(x$beta1),,drop=FALSE])%*%x$beta1/
              (sqrt(quasidisp*diag(t(A1[1:length(x$beta1),,drop=FALSE])%*%
                   x$CovM1%*% A1[1:length(x$beta1),,drop=FALSE]))))),
                error=function(e) rep(NA,ncol(A)),
                warning=function(w) rep(NA,ncol(A))
                )                
            else
                rep(NA,ncol(A1))
        }))
        
        if(verbose) cat("NB2 \n")
        p3 = t(sapply(fits$NB2,function(x){
            if(length(x$beta1)>0)
                tryCatch(
                2*pnorm(-abs(t(A1[1:length(x$beta1),,drop=FALSE])%*%x$beta1/
              sqrt(diag(t(A1[1:length(x$beta1),,drop=FALSE])%*%x$CovM1%*%
                      A1[1:length(x$beta1),,drop=FALSE])))),
                    error=function(e) rep(NA,ncol(A)),
                warning=function(w) rep(NA,ncol(A))
                )
            else
                rep(NA,ncol(A1))}))
        
        if(verbose) cat("Quasi-NB2 \n")
        p4 = t(sapply(fits$NB2,function(x){
            switch(quasi_method,
                   "Wedderburn" = {
                       quasidisp = x$Wedderburn
                   },
                   "Farrington" = {
                       quasidisp = x$Farrington
                   },
                   "Fletcher" = {
                       quasidisp = x$Fletcher
                   }
                   )
            if(!is.na(quasidisp) && length(x$beta1)>0)
                tryCatch(
                    2*pnorm(-abs(t(A1[1:length(x$beta1),,drop=FALSE])%*%x$beta1/
                 (sqrt(quasidisp* diag(t(A1[1:length(x$beta1),,drop=FALSE]) %*%
                           x$CovM1 %*% A1[1:length(x$beta1),,drop=FALSE]))))),
                    error=function(e) rep(NA,ncol(A)),
                    warning=function(w) rep(NA,ncol(A))
                )
            else
                rep(NA,ncol(A1))
        }))
        colnames(p1) = colnames(p2) = colnames(p3)=colnames(p4) = colnames(A)
        
    }
    list(Poi = p1,quasiPoi = p2,NB2 = p3,quasiNB2 = p4)
}


