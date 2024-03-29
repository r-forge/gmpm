gmpCreate <- function(formula, family, data=parent.frame(), ivars=c(),
                gmpControl=gmpCtrl(), ...)
{
  gmpFormals <- setdiff(names(formals(gmpCreate)), "...")
  
  if (missing(formula)) {
    stop("Argument 'formula' must be supplied.")
  }
  
  if (missing(family)) {
    stop("Argument 'family' must be supplied.")
  }
  
  myclass <- "Gmp.user" # default

  if (is.character(family)) {
      if (pmatch(c(family), c("multinomial", "user"),0) > 0) {
          if (pmatch(c(family), c("multinomial"), 0) > 0) {
              myclass = "Gmp.mul"
          } else {
            myclass <- "Gmp.user"
          }
      } else {
        myclass <- "Gmp.glm"
      }
  } else {
    if (is.function(family)) {
      family <- substitute(family)
    }
    myclass <- "Gmp.glm"
  }
  
  mc <- match.call(expand.dots = TRUE)
  c2 <- call("new", Class=myclass,
             formula=formula, family=family,
             data=data, ivars=ivars, gmpControl=gmpControl)
  mn <- setdiff(1:length(names(mc)), c(1,2,match(gmpFormals, names(mc))))  
  mc2 <- mc[c(1L, mn)]
  if (length(mn) > 0) {
    mcl <- as.call(c(as.list(c2), as.list(mc2[2L:length(mc2)])))
  } else {mcl <- c2}
  
  x <- eval(as.call(mcl))

  # now, assign variables common to all classes
  if (!missing(ivars))
    {
      x@ivars <- ivars
    }
  else {}
  
  x@gmpControl <- gmpControl
  
  .parseFormula(x, formula)

  # create model frame
  fakeform1 <- all.vars(x@mform)
  fakeform2 <- paste(fakeform1[1],
                     paste(fakeform1[2:length(fakeform1)],
                           collapse="+"),
                     sep="~")
  
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"),
             names(mf), 0)
  mf <- mf[c(1, m)]  
  mf$drop.unused.levels <- TRUE

  mfl1 <- as.list(mf)
  mfl1$formula <- fakeform2
  mfl1[[1]] <- as.name("model.frame")
  df1 <- eval(as.call(mfl1), parent.frame())

  attr(df1, "terms") <- NULL
  
  mfl2 <- as.list(mf)
  mfl2$formula <- x@dform
  mfl2[[1]] <- as.name("model.frame")
  df2 <- eval(as.call(mfl2), parent.frame())
  attr(df1, "terms") <- attr(df2, "terms")
  x@df1 <- df1
  x@DVname <- rownames(attr(attr(df1,"terms"),"factors"))[1]

  .checkMultilevel(x)
  .getDesign(x)
  x@fitcall <- .buildFitCall(x, ocall=match.call(expand.dots=TRUE))
  x@IVcoef <- .getFactorLabelsFromFit(x)
  
  .preparePermScheme(x)
  x@pspace <- permSpace(x)
  x@coef0 <- coef(fitOnce(x))
  x@coefTerms <- .getCoefTerms(x)  

  x@ivix <- .getIVix(x)
 
  return(x)
}
  

gmp <- function(formula, family, data=parent.frame(), ivars=c(),
                gmpControl=gmpCtrl(), ...)
{
  mc <- match.call()
  mc[[1]] <- gmpCreate

  # create the object
  x <- eval(mc)

  # fit the object
  return(gmpFit(x))
}
  
gmpCtrl <- function(maxruns = 999, report.interval=10,
                        outfile = NULL)
{
  return(list(maxruns=maxruns, report.interval=report.interval,
              outfile=outfile))
}

.getIVtypes <- function(x, ivars, munit)
  {
    nIVs <- length(ivars)
    nIVlevels <- rep(NA, nIVs)
    ivType <- rep("between", nIVs)
    if (is.null(munit) || munit=="") {
      return(ivType)
    }
    nSubj <- length(unique(x[,munit]))
    for (i in 1:nIVs) {
      nIVlevels[i] <- length(levels(x[,ivars[i]]))
      vt.form <- as.formula(paste("~",munit,"+",
                                  paste(ivars[[i]],collapse="+"),
                                  sep=""))
      vt.df <- as.data.frame(xtabs(vt.form,data=x))
      vt.df <- vt.df[vt.df$Freq > 0, ]
      vt.nobs <- length(vt.df[,1])
      if (vt.nobs == (nIVlevels[i])*nSubj) {
        ivType[i] <- "within"
      } else {
        if (vt.nobs != nSubj) {
          stop("error in design; all multilevel units must have all levels of within-subject factor")
        }
      }
    }
    
    return(ivType)
  }


#.getFactorLabels <- function(x, varname)
#  {
#    return(paste(varname, 1:(length(levels(x))-1), sep=""))
#  }

.sortByWithin <- function(x, ivWithin) {
  if (length(ivWithin) > 0) {
    for (i in length(ivWithin):1) {
      for (j in length(ivWithin[[i]]):1) {
        x <- x[order(x[,ivWithin[[i]][j]]),]
      }
    }
  } else {
    warning("no within IVs supplied; function had no effect")
  }
  return(x)
}

.getSig <- function(x) {
  sigvec <- rep("", length(x))
  x[is.na(x)] <- 1
  for (i in 1:length(x)) {
    if (x[i] <= .10) {
      sigvec[i] <- "."
    } else {}
    if (x[i] <= .05) {
      sigvec[i] <- "*"
    } else {}
    if (x[i] <= .01) {
      sigvec[i] <- "**"
    } else {}
    if (x[i] <= .001) {
      sigvec[i] <- "***"
    } else {}
  }
  return(sigvec)
}
