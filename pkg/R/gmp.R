setMethod("getModelFrame",
          signature(object="Gmp"),
          function(object) {
            return(object@df1)
          })

setMethod("getFactorCodes",
          signature(object="Gmp"),
          function(object) {
            ivs <- object@ivars
            nIVs <- length(object@ivars)
            fcodes <- list()
            for (i in 1:nIVs) {
              mx <- attr(object@df1[,ivs[i]],"contrasts")
              colnames(mx) <- object@IVcoef[[ivs[i]]]
              fcodes[[ivs[i]]] <- mx
            }
            return(fcodes)
          })

setMethod(".getFactorLabelsFromFit",
          signature(object="Gmp"),
          function(object, ivar) {
            fcall <- as.list(object@fitcall)
            dform <- object@dform
            lhs <- strsplit(deparse(dform), "~")[[1]][1]
            for (i in 1:length(object@ivars)) {
              nform <- paste(lhs, "~", object@ivars[i], sep="")
              fcall$formula <- nform
              if (object@famtype == "multinomial") {
                sink(".multinom.txt")
                nn1 <- colnames(coef(eval(as.call(fcall))))
                sink(NULL)
              } else {
                nn1 <- names(coef(eval(as.call(fcall))))
              }
              object@IVcoef[[object@ivars[i]]] <- nn1[2:length(nn1)]
            }
            return(object@IVcoef)
          })

setMethod(".getPredictorsFromFaclist",
          signature(x="Gmp"),
          function(x, faclist, allvars, nTests, j) {
            if (nTests > 1) {
              ivinc <- allvars[faclist[,j]==1]
            } else {
              ivinc <- allvars
            }

            if (length(intersect(ivinc, x@covars)) > 1) {
                                        # there is more than one co-variate in the term.
                                        # we need to perform a union rather than an intersection
              cvartmp <- intersect(ivinc, x@covars)
              cvars <- c()
              ivartmp <- intersect(ivinc, x@ivars)
              for (i in 1:length(cvartmp)) {
                cvars <- c(cvars, .getIXfromIV(x, c(ivartmp, cvartmp[i]), TRUE))
              }
            } else {
              cvars <- .getIXfromIV(x, ivinc, TRUE)
            }

            return(cvars)
          })

setMethod("testDV",
          signature(x="Gmp.mul"),
          function(x, excludeLevels, byCovar=FALSE) {
            DVlevels <- .getDVlevels(x)
            if (!missing(excludeLevels)) {
              if (is.character(excludeLevels)) {
                excludeLevels <- c(excludeLevels)
              } else {}
              if (DVlevels[1] %in% excludeLevels) {
                stop("Can't exclude baseline level '", DVlevels[1], "' from analysis.")
              } else {}
              ltest <- setdiff(DVlevels[2:length(DVlevels)], excludeLevels)
              if (length(ltest) < 2) {
                stop("Only ", length(ltest), " regions to test.\nMinimum of 2 required.")
              }
            } else {
              ltest <- DVlevels[2:length(DVlevels)]
            }
            ff <- .prepareMainSum(x, byCovar)
            nTests <- ff$nTests
            nDiffs <- length(ltest)-1
            pwid <- dim(x@pmx)[3]
            mx <- matrix(nrow=dim(x@pmx)[1], ncol=nDiffs*pwid)
            colnames(mx) <- rep(dimnames(x@pmx)$coef, nDiffs)
            newDVname <- paste("(",
                               paste(c(DVlevels[1], ltest), collapse=",", sep=""),
                               ")", sep="")
            for (k in 1:nDiffs) {
              ix0 <- (k-1)*(pwid)+1
              ix1 <- ix0+pwid-1
              mx[,ix0:ix1] <- x@pmx[,ltest[1],]-x@pmx[,ltest[k+1],]
            }
            mlist <- list()
            for (j in 1:nTests) {
              cvars <- .getPredictorsFromFaclist(x, ff$faclist, ff$allvars, ff$nTests, j)
              ctest <- cvars + rep(rep(0:(nDiffs-1))*pwid, each=length(cvars))
              tname <- paste(colnames(ff$faclist)[j], ":", newDVname, sep="")
              mlist[[tname]] <- ctest
            }
            return(mdTest(x, mx, mx[1,], mlist))
          })

setMethod(".getDVlevels",
          signature(x="Gmp.mul"),
          function(x) {
            if(length(grep("cbind", x@DVname))>0) {
              f1 <- strsplit(x@DVname, "\\(")[[1]][2]
              f2 <- gsub(")", "", f1)
              f3 <- gsub(" ", "", f2)
              return(strsplit(f3, ",")[[1]])
            }
            
            if (is.matrix(x@df1[,x@DVname]))
              return(colnames(x@df1[,x@DVname]))
            else if (is.factor(x@df1[,x@DVname]))
              return(levels(x@df1[,x@DVname]))
            else {}

            stop("unsure of how DV '",x@DVname,"' is represented\n",
                 "(was not factor, matrix, or cbind.)")
          })

setMethod(".writeFit",
          signature(x="Gmp"),
          function(x, y, outfile, append) {
            cat(y, "\n", file=outfile, append=append)
          })

setMethod(".writeFit",
          signature(x="Gmp.mul"),
          function(x, y, outfile, append) {
            nRows <- dim(y)[1]
            for (i in 1:nRows) {
              cat(y[i,], " ", file=outfile, append=append)
            }
            cat("\n", file=outfile, append=TRUE)
          })

setMethod("appendToPmx",
          signature(x="Gmp", y="Gmp"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            if (length(dim(x)) != length(dim(y))) {
              cat("Gmp source object permutation matrix has dimensions ", dim(y),
                  "\n")              
              cat("Gmp destination object permutation matrix has dimensions ", dim(x),
                  "\n")
              stop("These Gmp objects do not look the same.")
            }
            pmxSrc <- getPermMx(y)[-1,]
            pmxThis <- getPermMx(x)
            if (dim(pmxThis)[1]==0) {
              x@pmx <- getPermMx(y)
            } else {
              x@pmx <- rbind(pmxThis, pmxSrc)
            }
            cat("appended ", length(pmxSrc[,1]), " rows\n")            
            x@ncomp <- dim(x@pmx)[1]-1
            warning("Error checking not implemented yet; \nPlease ensure these two Gmp objects have the same underlying model / data.")
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod("appendToPmx",
          signature(x="Gmp", y="character"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            ff <- read.table(y)
            pmxSrc <- as.matrix(ff[-1,])
            rownames(pmxSrc) <- NULL
            pmxThis <- getPermMx(x)
            if (dim(pmxThis)[1]==0) {
              pmxSrc <- as.matrix(ff)
              colnames(pmxSrc) <- colnames(x@coef0)
              x@pmx <- pmxSrc
            } else {
              colnames(pmxSrc) <- colnames(pmxThis)
              x@pmx <- rbind(pmxThis, pmxSrc)
            }
            x@ncomp <- dim(x@pmx)[1]-1
            cat("appended ", length(pmxSrc[,1]), " rows\n")
            
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod("appendToPmx",
          signature(x="Gmp.mul", y="Gmp"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            if (length(dim(x)) != length(dim(y))) {
              cat("Gmp source object permutation matrix has dimensions ", dim(y),
                  "\n")              
              cat("Gmp destination object permutation matrix has dimensions ", dim(x),
                  "\n")
              stop("These Gmp objects do not look the same.")
            }
            pmxSrc <- getPermMx(y)[-1,,]
            srclen <- dim(pmxSrc)[1]
            pmxThis <- getPermMx(x)
            destlen <- dim(pmxThis)[1]
            nDVlevels <- dim(x@coef0)[1]
            nCoef <- dim(x@coef0)[2]
            dmn <- c("run","dv","coef")
            if (destlen > 0) {
              x@pmx <- array(dim=c(srclen+destlen,nDVlevels,nCoef),
                             dimnames=dmn)
              x@pmx[1:destlen,,] <- pmxThis
              x@pmx[(destlen+1):(srclen+destlen),,] <- pmxSrc
            } else {
              pmxSrc <- getPermMx(y)
              srclen <- dim(pmxSrc)[1]
              x@pmx <- getPermMx(y)
            }
            cat("appended ", srclen, " rows\n")            
            x@ncomp <- dim(x@pmx)[1]-1
            warning("Error checking not implemented yet; \nPlease ensure these two Gmp objects have the same underlying model / data.")
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod("appendToPmx",
          signature(x="Gmp.mul", y="character"),
          function(x, y) {
            nameObject <- deparse(substitute(x))            
            ff <- read.table(y)
            pmxSrc <- as.matrix(ff[-1,])
            rownames(pmxSrc) <- NULL
            pmxThis <- getPermMx(x)
            colnames(pmxSrc) <- colnames(pmxThis)
            x@pmx <- rbind(pmxThis, pmxSrc)
            x@ncomp <- dim(x@pmx)[1]-1
            cat("appended ", length(pmxSrc[,1]), " rows\n")
            
            assign(nameObject, x, envir=parent.frame())            
            return(invisible(x))
          })

setMethod(".prepareMainSum",
          signature(x="Gmp"),
          function(x, byCovar=FALSE) {

            # figure out tests we need to run from model frame
            faclist <-
              attr(attr(x@df1,"terms"),"factors")[-1,]
            if (is.vector(faclist)) {
              allvars <- x@ivars
              nTests <- 1
            } else {
              allvars <- rownames(faclist)
              nTests <- dim(faclist)[2]
            }

            # do not test main effects of covars
            if (length(x@covars) > 0) {
              m <- match(x@covars, colnames(faclist))
              m <- m[!is.na(m)]
              faclist <- faclist[,-m]
            }

            if ((!byCovar) && (length(x@covars) > 1)) {
              fl1 <- faclist[x@covars[-1],]
              if (is.vector(fl1)) {
                test1 <- fl1==0    
              } else {
                test1 <- colSums(fl1)==0    
              }
              newfl <- faclist[,faclist[x@covars[1],]==1 | test1]
              r1 <- newfl[x@covars[1],]
              c1 <- names(r1)[r1==1]
              newfl[x@covars[-1],c1] <- 1
              flnames <- strsplit(colnames(newfl), ":")
              srep <- paste("(",paste(x@covars, collapse=","),")",sep="")
              ncn <- rep("", length(flnames))
              for (i in 1:length(flnames)) {
                flnames[[i]][flnames[[i]]==x@covars[1]] <- srep
                ncn[i] <- paste(flnames[[i]],collapse=":",sep="")
              }
              colnames(newfl) <- ncn
              faclist <- newfl
            }
            
            if (is.vector(faclist)) {
              nTests <- 1
            } else {
              nTests <- dim(faclist)[2]
            }
            
            return(list(faclist=faclist,
                        allvars=allvars,
                        nTests=nTests))
          })

setMethod(".regSumProc",
          signature(x="Gmp"),
          function(x, pmx, index, includeCovars=TRUE) {
            #print("~~~ in .regSumProc ~~~")
            coef0 <- pmx[1,]
            ctest <- 1:length(coef0) %in% index

            c2 <- coef0
            se <- rep(NA, length(c2))
            nexceed <- rep(NA, length(c2))
            pval <- rep(NA, length(c2))              
            for (i in 1:length(c2)) {
              if (ctest[i]) {
                se[i] <- sd(pmx[,i])
                nexceed[i] <- getNExceeding(pmx, i)
                pval[i] <- getPValue(pmx, i)
              }
            }

            c2 <- round(c2,4)
            se <- round(se,4)
            pval <- round(pval,4)
            gmpRegSum <- data.frame(Coef=names(c2),
                                    Estimate=c2, se, nexceed,
                                    pval, .getSig(pval))
            rownames(gmpRegSum) <- 1:length(c2)
            colnames(gmpRegSum) <- c("Coefficient", "Estimate",
                                     "Std. Error", "N>=orig", "p-value", " ")
            return(gmpRegSum)
          })

setMethod("getRegSummary",
          signature(x="Gmp"),
          function(x, includeCovars=TRUE) {
            ff <- list(.regSumProc(x, x@pmx, x@ivix, includeCovars))
            names(ff) <- "Main Regression"
            return(ff)
          })

setMethod("getRegSummary",
          signature(x="Gmp.mul"),
          function(x, includeCovars=TRUE) {
            #print("~~~ in getRegSummary (Gmp.mul) ~~~")
            mlist <- list()
            DVlevels <- .getDVlevels(x)
            nDVlevels <- dim(x@coef0)[1]            
            mnames <- paste(DVlevels[2:length(DVlevels)], DVlevels[1],
                  sep=" versus ")
            for (i in 1:nDVlevels) {
              mlist[[i]] <- .regSumProc(x, x@pmx[,i,], x@ivix, includeCovars)
            }
            names(mlist) <- mnames
            #print("... exiting getMainSummary (Gmp.mul) ...")
            return(mlist)
            
          })

setMethod(".mainSumProc",
          signature(x="Gmp"),
          function(x, faclist, allvars, nTests, pmx) {
            #print("~~~ in .mainSumProc ~~~")

            coef0 <- pmx[1,]
            nge <- rep(NA, nTests)
            pval <- rep(NA, nTests)
            mcoef <- rep(NA, nTests)

            for (j in 1:nTests) {
              cvars <- .getPredictorsFromFaclist(x, faclist, allvars, nTests, j)

              if (length(cvars) > 1) {
                mdt <- mdTest(x, pmx, coef0, cvars)
                nge[j] <- getResults(mdt, 1, "nge")
                pval[j] <- getResults(mdt, 1, "pval")
                mcoef[j] <- paste(getResults(mdt, 1, "ix"), collapse=",", sep="")
              } else {
                nge[j] <- getNExceeding(pmx, cvars[1])
                pval[j] <- getPValue(pmx, cvars[1])
                mcoef[j] <- cvars[1]
              }
              ##############################################
            }

            gmpMainSum <- data.frame(mcoef, nge, pval, .getSig(pval))
            colnames(gmpMainSum) <- c("Coef","N>=Orig","p-value", " ")
            if (nTests > 1) {
              rownames(gmpMainSum) <- colnames(faclist)
            } else {
              rownames(gmpMainSum) <- x@ivars
            }

            if (length(x@covars) > 0) {                                       
              vv <- faclist[!(rownames(faclist) %in% x@ivars),]
              if (!is.vector(vv)) {
                vv <- as.vector(colSums(vv))
                vv[vv>1] <- 1
              }
              gmpMainSum <- gmpMainSum[order(vv),]
            }

            #print("... exiting .mainSumProc ...")               
            return(gmpMainSum)
          })

setMethod("getMainSummary",
          signature(x="Gmp"),
          function(x, byCovar=FALSE) {
            gg <- .prepareMainSum(x, byCovar)
            faclist <- gg[["faclist"]]
            allvars <- gg[["allvars"]]
            nTests <- gg[["nTests"]]
            #print("~~~ in getMainSummary (Gmp) ~~~")
            ff <- list(.mainSumProc(x, faclist, allvars, nTests, x@pmx))
            names(ff) <- c("Main Results")
            #print("... exiting getMainSummary Gmp) ...")
            return(ff)
          })

setMethod("getMainSummary",
          signature(x="Gmp.mul"),
          function(x, byCovar=FALSE) {
            #print("~~~ in getMainSummary (Gmp.mul) ~~~")
            gg <- .prepareMainSum(x, byCovar)
            faclist <- gg[["faclist"]]
            allvars <- gg[["allvars"]]
            nTests <- gg[["nTests"]]
            mlist <- list()
            DVlevels <- .getDVlevels(x)
            nDVlevels <- dim(x@coef0)[1]            
            mnames <- paste(DVlevels[2:length(DVlevels)], DVlevels[1],
                  sep=" versus ")
            for (i in 1:nDVlevels) {
              mlist[[i]] <-
                .mainSumProc(x, faclist, allvars, nTests, x@pmx[,i,])
            }
            names(mlist) <- mnames
            #print("... exiting getMainSummary (Gmp.mul) ...")
            return(mlist)
          })

setMethod(".getIXfromIV",
          signature(x="Gmp"),
          function(x, ivinc, includeCovars=TRUE) {
            flab <- x@IVcoef
            if (includeCovars) {
              if (length(x@covars) > 0) {
                for (i in 1:length(x@covars)) {
                  flab[[x@covars[i]]] <- x@covars[i]
                }
              } else {}
            } else {}
            
            cvars <- c()
            for (k in 1:length(x@coefTerms)) {
              ctmp <- c()
              for (m in 1:length(ivinc)) {
                ivForms <- flab[[ivinc[m]]]
                isIn <- ivForms %in% x@coefTerms[[k]]
                if (sum(isIn) > 0) {
                  if (sum(isIn) != 1) {
                    stop("something's wrong here; bailing out.")
                  } else {
                    ctmp <- c(ctmp, ivForms[isIn])
                  }
                }
              }
              if ((length(ctmp) == length(ivinc)) &&
                  (length(ctmp) == length(x@coefTerms[[k]]))) {
                cvars <- c(cvars, k)
              } else {}
            }
            return(cvars)
          })

           
setMethod(".getCoefTerms",
          signature(x="Gmp"),
          function(x, mycoef) {
            if (missing(mycoef)) {
              mycoef <- x@coef0
            }
            if (is.matrix(mycoef)) {
              n <- colnames(mycoef)
            } else {
              n <- names(mycoef)
            }
            nl <- list()
            for (i in 1:length(n)) {
              nl[[i]] <- strsplit(n[i], ":")[[1]]
            }
            x@coefTerms <- nl
            return(nl)
          }
          )

setMethod(".storeFitResult",
          signature(x="Gmp"),
          function(x, fit, index) {
            nameObject <- deparse(substitute(x))

            myFit <- coef(fit)
            x@pmx[index,] <- myFit

            assign(nameObject, x, envir=parent.frame())            
          }
          )

setMethod(".storeFitResult",
          signature(x="Gmp.mul"),
          function(x, fit, index) {
            
            nameObject <- deparse(substitute(x))

            myFit <- coef(fit)
            x@pmx[index,,] <- myFit
            if (!is.null(fit$convergence)) {
              x@convergence <- if(fit$convergence==1) FALSE else TRUE
            }

            assign(nameObject, x, envir=parent.frame())            
          }
          )

setMethod(".reportProgress",
          signature(x="Gmp"),
          function(x, myFit, ix, maxruns, elapsed) {
            if (is.matrix(myFit)) {
              myFit <- as.vector(myFit[1,])
            }
            tmp <- myFit[x@ivix]
            tmp1 <- tmp[if(length(tmp)>3) 1:3 else length(tmp)]
            fmtstr1 <- paste("%",nchar(as.character(maxruns)),"d",
                            sep="")
            fmtstr2 <- paste(fmtstr1, "/", fmtstr1, ":", sep="")
            stmp1 <- sprintf(fmtstr2, ix, maxruns)
            stmp2 <- sprintf("% 1.3f", tmp1)
            if (length(elapsed) > 1) {
              efac <- mean(elapsed)
            } else {
              efac <- elapsed[1]
            }
            cat(stmp1, stmp2, if (length(tmp)>3) "..." else " ",
                sprintf("%1.3fs/run, ", efac))

            totsecs <- efac*(maxruns-ix)
            cat(sprintf("%02dh:", floor(totsecs / 3600)))
            if (totsecs >= 3600) {
              remn <- totsecs %% 3600
            } else {
              remn <- totsecs
            }
            cat(sprintf("%02dm:", floor(remn/60)))
            cat(sprintf("%02ds", round(remn %% 60)), "left\n")
          }
          )

setMethod(".createPermMx",
          signature(x="Gmp"),
          function(x, maxruns) {
            x@pmx <- 
              matrix(nrow=maxruns+1, ncol=length(x@coef0))
            x@pmx[1,] <- x@coef0
            colnames(x@pmx) <- names(x@coef0)
            return(x@pmx)
          }
          )

setMethod(".createPermMx",
          signature(x="Gmp.mul"),
          function(x, maxruns) {
            nameObject <- deparse(substitute(x))
            
            x@pmx <- array(dim=c(maxruns+1, dim(x@coef0)[1], dim(x@coef0)[2]))
            sink(".multinom.txt")
            f0 <- origFit(x)
            sink(NULL)
            x@pmx[1,,] <-coef(f0)
              
            dimnames(x@pmx) <- list(run=1:(maxruns+1),
                                    dv=rownames(x@coef0),
                                    coef=colnames(x@coef0))

            if (!is.null(f0$convergence)) {
              x@convergence <- rep(FALSE, maxruns+1)
              x@convergence[1] <- if(f0$convergence == 1) FALSE else TRUE
            } else {
              warning("You are using an older version of package 'nnet' ",
                      "that doesn't provide\n", "information about ",
                      "convergence.  Please consider updating.\n",
                      "(see ?update.packages() for instructions).")
              x@convergence <- rep(TRUE, maxruns+1)              
            }
            assign(nameObject, x, envir=parent.frame())            

            return(x@pmx)
          }
          )

setMethod("getNExceeding",
          signature(x="Gmp"),
          function(x, y) {
            if (missing(y)) {
              ix <- .getIVix(x)
            } else {
              ix <- y
            }
            if (is.null(x@coef0) || is.null(x@pmx)) {
              stop("Model must be fit (gmpFit) before extracting results.")
            }
            if (length(setdiff(ix, x@ivix)) > 0) {
              print(names(x@coef0)[setdiff(ix, x@ivix)])
              stop("Error: p-values only meaningful for independent variables and IV*covariate interactions.")
            }

            return(getNExceeding(x@pmx, ix))
          }
          )

setMethod("getNExceeding",
          signature(x="matrix"),
          function(x, y) {
            pmx <- x
            coef0 <- x[1,]
            ctmp <- c()
            for (i in 1:length(y)) {
              ctmp <- c(ctmp,
                        sum(abs(pmx[,y[i]])>=abs(coef0[y[i]])))
            }
            
            return(ctmp)
          })


setMethod("getPValue",
          signature(x="Gmp"),
          function(x, y) {
            if (missing(y)) {
              y <- .getIVix(x)
            }
            return(getPValue(getPermMx(x), y))
            #nge <- getNExceeding(x, ix)
            #return(nge/(x@ncomp+1))
          }
          )

setMethod("getPValue",
          signature(x="matrix"),
          function(x, y) {
            nge <- getNExceeding(x, y)
            ntot <- dim(x)[1]
            return(nge/ntot)
          })

setMethod("gmpCoef",
          signature(x="Gmp"),
          function(x) {
            if (is.null(x@coef0)) {
              x@coef0 <- coef(fitOnce(x))
            }
            return(x@coef0)
          }
          )

#setMethod("coef",
#          signature(x="Gmp"),
#          function(x) {
#            coefficients(x)
#          })

setMethod("getPermMx",
          signature(x="Gmp"),
          function(x) {return(x@pmx)})

setMethod("permute",
    signature(x = "Gmp"),
    function (x) 
          {
                                        ##print("~~~ in Permute (Gmp) ~~~")
            # if there are any between subjects variables, then
            # shift around blocks of predictors by subject
            if (x@nBetween) {
              if (x@nunits > 1) {
                ordVec <- sample(1:x@nunits,x@nunits,replace=FALSE)
              } else {
                ordVec <- sample(1:length(x@psBetween[,1]),
                                 length(x@psBetween[,1]),
                                 replace=FALSE)
              }
              #x@lastPerm$between <- ordVec

              for (k in 1:x@nBetween) {
                  x@df1[,x@ivBetween[k]] <-
                    C(rep(x@psBetween[ordVec,x@ivBetween[k]],
                          x@psBetween[,'Freq']),
                      attr(x@df1[,x@ivBetween[k]], "contrasts"))
              }
            }

            # if there are within subject variables, shift around conditions
            if (x@nWithin) {
              for (k in 1:x@nWithin) {
                nCols <- dim(x@psWithin[[x@ivWithin[k]]])[2]
                newv1 <- rep(sample(1:nCols,x@nunits,replace=TRUE),
                             each=attr(x@psWithin[[x@ivWithin[k]]],
                               "mtimes"))
                newv2 <- c(x@psWithin[[x@ivWithin[k]]][,newv1])
                newv <- rep(newv2, each=attr(x@psWithin[[x@ivWithin[k]]],
                                     "meach"))
                newv <- factor(newv,
                               levels=levels(x@df1[,x@ivWithin[k]]))
                newv <- rep(newv, x@psWithin$ps[,'Freq'])
                x@df1[,x@ivWithin[k]] <-
                  C(newv, attr(x@df1[,x@ivWithin[k]], "contrasts"))
              }
            }
            
            return(invisible(x))
          }
          )

###############################

setMethod("origFit",
          signature(object = "Gmp"),
          function (object) {
            return(eval(as.call(object@fitcall)))
          }
          )

setMethod(".getIVix",
          signature(x = "Gmp"),
          function (x) {
            return(object@ivix)
          }
          )

          #####################

setMethod(".parseFormula",
          signature(object = "Gmp"),
          function(object, formula)
          {
            #print("~~~ in .parseFormula ~~~")
            nameObject <- deparse(substitute(object))

            fstr <- deparse(formula)
            object@mform <- formula

            fterms <- strsplit(fstr, "~")
            if (length(fterms[[1]]) != 2) {
              stop("Formula", paste("'", fstr, "'", sep=""),
                   "is malformed.  Must have exactly two sides, separated by '~'.")
            }
            lhs <- fterms[[1]][1]
            rhs <- fterms[[1]][2]
            spos <- gregexpr("\\|", rhs)[[1]][1]
            if (spos > 0) {
              rhs <- substr(rhs, 1, spos-1)
            }
            object@dform <- as.formula(paste(lhs, rhs, sep="~"))
            
            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))
          }
          )

setMethod(".checkMultilevel",
          signature(object="Gmp"),
          function(object)
          {
            #print("~~~ in .checkMultilevel ~~~")            
            nameObject <- deparse(substitute(object))

            ## MULTILEVEL DATA ? ###################################
            object@nunits <- 1
            object@munit <- ""
            ff <- object@mform
            if ("|" %in% as.character(ff[[3]]))
              {
                object@munit <- as.character(ff[[3]])[3]
                object@nunits <- dim(xtabs(as.formula(
                                                      paste("~",
                                                            object@munit,
                                                            sep="")),
                                           data=object@df1))

                if (object@nunits == dim(object@df1)[1])
                  {
                    #cat("Warning: Only one observation per sampling unit; data are single-level.\nGrouping factor will be ignored.", "\n")
                    stop("only one observation per sampling unit; your data are single-level.\ngmpm not yet configured for single-level data; please wait until next version.")
                    object@nunits <- 1
                    object@mform <- object@dform
                    object@munit <- ""
                  }
              } else {
                stop("No grouping factor supplied.\n gmpm not configured to handle single-level data; please wait for next version.")
                cat("No conditioning variable specified in model; assuming data are not multilevel.\n")
              }

            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))
          }
          )

setMethod(".getDesign",
          signature(object="Gmp"),
          function(object)
          {
            #print("~~~ in .getDesign ~~~")
            nameObject <- deparse(substitute(object))
            
            mfactors <- attr(attr(object@df1, "terms"), "factors")
            IVnames <- object@ivars
            covars <- setdiff(rownames(mfactors)[-1], IVnames)
            if (length(covars) > 0) {
              object@covars <- covars
            }
            
            if (length(IVnames)==0)
              {
                cat("Warning: No independent variables specified;\n",
                    "will perform one-sample test.\n")
              }
            else
              {
                allfacts <- rownames(mfactors)[2:(dim(mfactors))[1]]
                if (length(intersect(IVnames, allfacts)) < length(IVnames)) {
                  stop("Malformed formula: one or more IVs in list",
                       "missing from formula.")
                }

                ## check whether IVs are all coded as factors
                IVinfo <- data.frame(nLevels=rep(NA, length(IVnames)),
                                     Type=rep(NA, length(IVnames)),
                                     Levels=rep(NA, length(IVnames)))
                rownames(IVinfo) <- IVnames


                for (i in 1:length(IVnames)) {
                  f1 <- object@df1[,IVnames[i]]
                  if (!is.factor(f1)) {
                    errstr <- paste("Warning: Converting ", IVnames[i],
                                    " to a factor")
                    cat(errstr, "\n")
                    object@df1[,IVnames[i]] <- factor(f1)
                  }
                  nLevels <- length(levels(object@df1[,IVnames[i]]))
                  IVinfo[IVnames[i],"nLevels"] <- nLevels
                  if (nLevels == 2)
                    {
                      cmx <- c(-.5,.5)
                    }
                  else
                    {
                      cmx <- matrix(nrow=nLevels, ncol=nLevels)
                      cmx[] <- 0
                      diag(cmx) <- 1
                      cmx <- cmx[,2:nLevels]
                      for (j in 2:nLevels) {
                        cmx[,(j-1)] <- cmx[,(j-1)]-mean(cmx[,(j-1)])
                      }
                    }
                  object@df1[,IVnames[i]] <- C(object@df1[,IVnames[i]], cmx)

                  IVinfo[IVnames[i],"Levels"] <-
                    paste(levels(object@df1[,IVnames[i]]),
                          collapse=",")

                }


                IVinfo$Type <- .getIVtypes(object@df1, object@ivars,
                                          object@munit)
                object@nWithin <- sum(IVinfo$Type=="within")
                object@nBetween <- sum(IVinfo$Type=="between")                
                object@IVinfo <- IVinfo
              }

            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))
          }          
          )

setMethod(".buildFitCall",
          signature(object="Gmp"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (Gmp)~~~")            
            # build function call (common to all classes) ################
            # this function is only invoked by the '.buildFitCall'
            # functions for the child classes (Gmp.mul, Gmp.user).
            #
            # The current function handles situations common to all
            # classes, and passes the result back to the child class.
            #
            # Note that this function does not store the call in the object;
            # this is handled by the 'child' class.

            gmp.args <- c("gmpControl", "ivars")
            # note: update the above line when new arguments are
            # introduced to GMP to avoid them being passed along
            # to fitting function.
            
            # get rid of gmp-specific arguments
            ocall <- ocall[c(TRUE, !(names(ocall[2:length(ocall)]) %in%
                               c(arg.exclude, gmp.args)))]
            mcl <- c(list(fun="user"), as.list(ocall[2:length(ocall)]))
            mcl$formula <- object@dform
            mcl$data <- quote(object@df1)
            #print(names(mcl))

            return(invisible(mcl))
          }
          )

setMethod(".buildFitCall",
          signature(object="Gmp.glm"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (Gmp.glm)~~~")                        

            ncall <- callNextMethod(object=object, arg.exclude=c(),
                           ocall=ocall)
            ncall[[1]] <- glm
            return(invisible(ncall))
          }            
          )

setMethod(".buildFitCall",
          signature(object="Gmp.mul"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (Gmp.mul)~~~")
            # build function call (common to all classes) ################

            ncall <- callNextMethod(object=object, arg.exclude=c("family"),
                           ocall=ocall)
            library(nnet)
            ncall[[1]] <- multinom
            return(invisible(ncall))
          }            
          )

setMethod(".buildFitCall",
          signature(object="Gmp.user"),
          function(object, arg.exclude=c(), ocall)
          {
            #print("~~~ in .buildFitCall (Gmp.user)~~~")

            ncall <- callNextMethod(object=object, arg.exclude=c("family"),
                           ocall=ocall)
            
            stop("user-defined gmp objects not yet implemented.\n",
                 "Please wait for next version of GMP package.")
          }            
          )

setMethod("fitOnce",
          signature(object="Gmp"),
          function(object)
          {
                                        ##print("~~~ in fitOnce (Gmp) ~~~")
            ff <- eval(as.call(object@fitcall))
            return(ff)
          }
          )

setMethod("fitOnce",
          signature(object="Gmp.mul"),
          function(object)
          {
                                        ##print("~~~ in fitOnce (Gmp) ~~~")
            sink(".multinom.txt")
            fit1 <- eval(as.call(object@fitcall))
            sink(NULL)
            if (!is.null(fit1$convergence)) {
              if (fit1$convergence == 1) {
                cat("Warning: multinom function did not converge.\n",
                    "Try increasing number of multinom iterations",
                    " (use 'maxit' argument).\n")
              } else {}
            } else {}
                
            return(fit1)
          }
          )

setMethod(".preparePermScheme",
          signature(object="Gmp"),
          function(object)
          {
            #print("~~~ in .preparePermScheme ~~~")
            nameObject <- deparse(substitute(object))            
            
            # prepare permutation scheme ####################
            # and re-sort the data set.

            # save typing/make code more readable
            nWithin <- object@nWithin
            nBetween <- object@nBetween
            rnames <- rownames(object@IVinfo)
            ivWithin <- rnames[object@IVinfo$Type=="within"]
            ivBetween <- rnames[object@IVinfo$Type=="between"]
            object@ivWithin <- ivWithin
            object@ivBetween <- ivBetween
            id <- object@munit
            x <- object@df1
            psWithin <- object@psWithin
            
            # sort data by within-subjects variables
            # (this makes perm tests cleaner)            
            if (nWithin > 0) {
              x <- .sortByWithin(x, ivWithin)
            }
            if (object@nunits > 1) {
              x <- x[order(x[,id]),]
            }

            # create summary table (xt1) for use in rearranging labels
            tmpIVvec <- c(ivWithin[nWithin:1], (ivBetween))
            tmpIVvec <- tmpIVvec[!is.na(tmpIVvec)]
            if (object@nunits > 1) {
              xtabs.form <-
                as.formula(paste("~",paste(tmpIVvec,collapse="+"),
                                 "+",id,sep=""))
            } else {
              xtabs.form <-
                as.formula(paste("~",paste(tmpIVvec,collapse="+"),sep=""))
            }
            xt1 <- as.data.frame(xtabs(xtabs.form, data=x))
            xt1 <- xt1[xt1$Freq>0,]

            if (nWithin > 0) {
              xt1 <- .sortByWithin(xt1,ivWithin)
            }

            if (object@nunits > 1) {
              xt1 <- xt1[order(xt1[,id]),]

              xtform2 <-
                as.formula(paste("Freq~",paste(ivBetween,collapse="+"),
                                 "+",id,sep=""))
              xt2 <- as.data.frame(xtabs(xtform2,data=xt1))
              xt2 <- xt2[xt2$Freq>0,]
              xt2 <- xt2[order(xt2[,id]),]
              nvec <- xt2[,'Freq']
              xtform2 <-
                as.formula(paste("~",paste(ivBetween,collapse="+"),
                                 "+",id,sep=""))
              xt3 <- as.data.frame(xtabs(xtform2,data=xt1))
              xt3 <- xt3[xt3$Freq>0,]
              xt3 <- xt3[order(xt3[,id]),]
              nvec1 <- xt3[,'Freq']

              if (sum(nvec1==nvec1[1]) != length(nvec1)) {
                print(xt3)
                stop("design is unbalanced.")
              }
            } else {
              xt2 <- xt1
              nvec <- 1
            }

            # now handle psWithin
            if (nWithin > 0) {
              psWithin$ps <- xt1
              recPerms <- function(v) {
                if (length(v)==1) {
                  return(c(v))
                }
                lx <- length(v)
                nPerms <- factorial(lx)
                mx <- matrix(ncol=lx, nrow=nPerms)
                f1 <- nPerms / lx

                if (lx > 8) {
                  stop("Cannot currently handle IVs with more than 8 levels.")
                }
                
                for (i in 1:lx) {
                  start <- (i-1)*f1+1
                  end <- start+f1-1
                  mx[start:end,1] <- v[i]
                  mx[start:end,2:lx] <- recPerms(setdiff(v,v[i]))
                }

                return(mx)
              }
              for (i in 1:length(ivWithin)) {
                psWithin[[ivWithin[i]]] <- t(recPerms(levels(x[,ivWithin[i]])))
                meach <- 1
                mtimes <- 1
                if (i > 1) {
                  for (m in 1:(i-1)) {
                    mtimes <- mtimes*length(levels(x[,ivWithin[m]]))
                  }
                }
                if (i<length(ivWithin)) {
                  for (m in (i+1):length(ivWithin)) {
                    meach <- meach*length(levels(x[,ivWithin[m]]))
                  }
                } else {}                  
                attr(psWithin[[ivWithin[i]]], "meach") <- meach
                attr(psWithin[[ivWithin[i]]], "mtimes") <- mtimes
              }
            }
            object@nCellsPerUnit <- nvec[1]
            object@psBetween <- xt2
            object@df1 <- x
            object@psWithin <- psWithin
            
            assign(nameObject, object, envir=parent.frame())
            return(invisible(object))            
          }
          )

setMethod("gmpFit",
          signature(object="Gmp"),
          function(object,gmpControl)
          {
            #print("~~~ in gmpFit (Gmp) ~~~")
            if (missing(gmpControl)) {
              gmpControl <- object@gmpControl
            }
            if (is.null(object@coef0)) {
              object@coef0 <- coef(fitOnce(object))
            }
            maxruns <- gmpControl$maxruns
            report.interval <- gmpControl$report.interval
            outfile <- gmpControl$outfile
            elapsed <- rep(NA, maxruns)

            object@pmx <- .createPermMx(object, maxruns)
            
            if (!is.null(outfile)) {
              stop("writing to outfile not yet supported; bailing out.\n")
              .writeFit(object, coef(fitOnce(object)), outfile, FALSE)
            } else {}
            
            for (i in 1:maxruns) {
              t1 <- proc.time()["elapsed"]
              x <- permute(object)

              ff <- fitOnce(x)
              myFit <- coef(ff)
              .storeFitResult(object, ff, i+1)

              if (!is.null(outfile)) {
                .writeFit(object, myFit, outfile, TRUE)
              } else {}

              t2 <- proc.time()["elapsed"]
              elapsed[i] <- t2-t1

              if (report.interval > 0) {
                if ((i == maxruns) ||
                    ((i %% report.interval)==0)) {
                  .reportProgress(x, myFit, i, maxruns, elapsed[1:i])
                } else {}
              } else {}
            }

            object@ncomp <- i
            #print("~~~ leaving gmpFit ~~~")            
            return(object)
            
          }

          )

setMethod(".getIVix",
          signature(x="Gmp"),
          function(x) {
            ivix <- c()
            n <- names(x@coef0)
            nl <- x@coefTerms

            ivix <- c()
            for (i in 1:length(nl)) {
              varsIn <- nl[[i]]
              for (j in 1:length(x@IVcoef)) {
                for (m in 1:length(x@IVcoef[[j]])) {
                  if (x@IVcoef[[j]][m] %in% varsIn) {
                    ivix <- c(ivix, i)
                  }
                }
              }
            }
            ivix <- sort(unique(ivix))
            
            return(ivix)
          }
          )


setMethod("permSpace",
          signature(object="Gmp"),
          function(object)  {
            psBetween <- object@psBetween
            nCellsPerUnit <- object@nCellsPerUnit
            ivBetween <- object@ivBetween
            ivWithin <- object@ivWithin
            nunits <- object@nunits
            recChoose <- function(v) {
              tot <- sum(v)
              if (length(v)==1) {
                return(1)
              } else {
                n1 <- choose(tot, v[1])
                return(n1*recChoose(v[2:length(v)]))
              }
            }
            
            if (length(ivBetween)) {
              betform <- as.formula(paste("Freq ~", paste(ivBetween,collapse="+",sep="")))
              xt1 <- as.data.frame(xtabs(betform, data=psBetween))
              xt1 <- xt1[xt1$Freq > 0,] # get rid of unused levels
              bfac <- recChoose(xt1$Freq)
            } else {
              bfac <- 1
            }

            if (length(ivWithin)) {
              cc <- 1
              for (j in 1:length(ivWithin)) {
                cc <- cc*length(levels(object@df1[,ivWithin[j]]))
              }
              wfac <- cc^nunits
            } else {
              wfac <- 1
            }

            return(bfac*wfac)
          }
          )

setMethod("coefNames",
          signature(x="Gmp"),
          function(x)  {
            return(names(x@coef0))
          })

setMethod("coefNames",
          signature(x="Gmp.mul"),
          function(x)  {
            return(colnames(x@coef0))
          })
