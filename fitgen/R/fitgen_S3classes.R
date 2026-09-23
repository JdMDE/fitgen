#' A function to create the S3 object \code{fitgen_glm}
#' @title new_fitgen_glm
#' @param object Object
#' @param ... Additional arguments
#' @rdname fitgen_glm
#' @export new_fitgen_glm
new_fitgen_glm = function(object, ...){
    stopifnot(is.numeric(object$beta))
    stopifnot(is.numeric(object$muest))
    stopifnot(is.numeric(object$X2))
    stopifnot(is.numeric(object$sums))
    stopifnot(is.numeric(object$Dev))
    stopifnot(is.matrix(object$CovM))
    stopifnot(is.integer(object$iters))
    stopifnot(length(object) == 7)
    structure(object,class="fitgen_glm")
}

#'  Accesor for the coefficients
#' @param object Object of class \code{fitgen_glm}
#' @param ... Additional arguments
#' @export
coefficients.fitgen_glm = function(object,...)
    object$beta

#'  Accesor for the coefficients
#' @param object Object of class \code{fitgen_glm}
#' @param ... Additional arguments
#' @export
coef.fitgen_glm = function(object, ...)
    object$beta


#'  Accesor for the fitted values 
#' @param object Object of class \code{fitgen_glm}
#' @param ... Additional arguments
#' @export
fitted.fitgen_glm = function(object,...)
    object$muest

#'  Accesor for the deviance
#' @param object Object of class \code{fitgen_glm}
#' @param ... Additional arguments
#' @export
deviance.fitgen_glm = function(object,...)
    object$Dev

#'  Class fitgen_Lglm
#' @param object Object of class \code{fitgen_Lglm}
#' @param ... Additional arguments
#' @export
new_fitgen_Lglm = function(object,...){
     stopifnot(is.list(object))
    structure(object,class="fitgen_Lglm")
}


#' Method to extract the coefficients of a list
#' of glm's 

#'  Accesor for coefficients
#' @param object Object of class \code{fitgen_Lglm}
#' @param type The coefficients corresponds to the full model (with
#' neighbours)
#' or to the common model (only with the common predictors)
#' @param ... Additional arguments
#' @export
coef.fitgen_Lglm = function(object,type=c("full","common"),...){
    type = match.arg(type)
    switch(type,
           "common" = lapply(object,function(u) u$beta0),
           "full" =  lapply(object,function(u) u$beta1))
}

#'  Accesor for coefficients
#' @param object Object of class \code{fitgen_Lglm}
#' @param type The coefficients corresponds to the full model (with neighbours)
#' or to the common model (only with the common predictors)
#' @param ... Additional arguments
#' @export
coefficients.fitgen_Lglm = function(object,type=c("full","common"),...){
    type = match.arg(type)
    switch(type,
        "common" = lapply(object,function(u) u$beta0),
        "full" =  lapply(object,function(u) u$beta1))
}

#'  Accesor for fitted values 
#' @param object Object of class \code{fitgen_Lglm}
#' @param type The coefficients corresponds to the full model (with neighbours)
#' or to the common model (only with the common predictors)
#' @param ... Additional arguments
#' @export
fitted.fitgen_Lglm = function(object,type=c("full","common"),...){
    type = match.arg(type)
    switch(type,
        "common" = lapply(object,function(u) u$muest0),
        "full" =  lapply(object,function(u) u$muest1))
}

#'  Accesor for deviances
#' @param object Object of class \code{fitgen_Lglm}
#' @param type The coefficients corresponds to the full model (with neighbours)
#' or to the common model (only with the common predictors)
#' @param ... Additional arguments
#' @export
deviance.fitgen_Lglm = function(object,type=c("full","common"),...){
    type = match.arg(type)
    switch(type,
        "common" = lapply(object,function(u) u$Dev0),
        "full" =  lapply(object,function(u) u$Dev1))
}

#' A function to create the S3 method quasidispersion
#' @title quasidispersion
#' @param object Object of class \code{fitgen_Lglm}
#' @param type Method of estimation of  the quasidispersions
#' @param ... Additional arguments
#' @export
quasidispersion <- function(object,type,...) {
    UseMethod("quasidispersion")
}

#'  Accesor for quasidispersions
#' @param object Object of class \code{fitgen_Lglm}
#' @param type Method of estimation of the quasidispersions
#' @param ... Additional arguments
#' @export
quasidispersion.fitgen_Lglm = function(object,type=c("Wedderburn",
                                                     "Farrington",
                                                     "Fletcher"),...){
    type = match.arg(type)
    switch(type,
        "Wedderburn" = lapply(object,function(u) u$Wedderburn),
        "Farrington" = lapply(object,function(u) u$Farrington),
        "Fletcher" = lapply(object,function(u) u$Fletcher))
}

#' A function to create the S3 method dispersion
#' @title dispersion
#' @param object Object of class \code{fitgen_Lglm}
#' @param type Method to estimate the dispersions
#' @param ... Additional arguments
#' @export
dispersion <- function(object,type,...) {
    UseMethod("dispersion")
}


#' Accesors for dispersions
#' @title Accesor for dispersion 
#' @param object Object of class \code{fitgen_Lglm}
#' @param type Method to estimate the dispersions
#' @param ... Additional arguments
#' @export
dispersion.fitgen_Lglm = function(object,type=c("Bliss"),...){
    type = match.arg(type)
    switch(type,
        "Bliss" = lapply(object,function(u) u$Bliss))
}


