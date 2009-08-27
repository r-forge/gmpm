setGeneric(name="getModelFrame",
           def=function(object) {
             standardGeneric("getModelFrame")})

setGeneric(name="getFactorCodes",
           def=function(object) {
             standardGeneric("getFactorCodes")})

setGeneric(name=".getFactorLabelsFromFit",
           def=function(object, ivar) {
             standardGeneric(".getFactorLabelsFromFit")})

setGeneric(name=".getPredictorsFromFaclist",
           def=function(x, faclist, allvars, nTests, j) {
             standardGeneric(".getPredictorsFromFaclist")})

setGeneric(name="testDV",
           def=function(x, excludeLevels, byCovar=FALSE) {standardGeneric("testDV")})

setGeneric(name=".getDVlevels",
           def=function(x) {standardGeneric(".getDVlevels")})

setGeneric(name=".writeFit",
           def=function(x,y,outfile,append) {standardGeneric(".writeFit")})

setGeneric(name="appendToPmx",
           def=function(x,y) {standardGeneric("appendToPmx")})

setGeneric(name=".prepareMainSum",
           def=function(x, byCovar=FALSE) {standardGeneric(".prepareMainSum")})

setGeneric(name=".mainSumProc",
           def=function(x, faclist, allvars, nTests, pmx) {
             standardGeneric(".mainSumProc")
           })

setGeneric(name=".regSumProc",
           def=function(x, pmx, index, includeCovars=TRUE) {
             standardGeneric(".regSumProc")
           })

setGeneric(name="getMainSummary",
           def=function(x, byCovar=FALSE) {
             standardGeneric("getMainSummary")
           })

setGeneric(name="getRegSummary",
           def=function(x, includeCovars=TRUE) {
             standardGeneric("getRegSummary")
           })

setGeneric(name=".getIXfromIV",
           def=function(x, ivinc, includeCovars=TRUE) {
             standardGeneric(".getIXfromIV")
           })

setGeneric(name=".getCoefTerms",
           def=function(x, mycoef) {
             standardGeneric(".getCoefTerms")}
           )

setGeneric(name=".storeFitResult",
           def=function(x, fit, index) {
             standardGeneric(".storeFitResult")}
           )

setGeneric(name=".reportProgress",
           def=function(x, myFit, ix, maxruns, elapsed) {
             standardGeneric(".reportProgress")}
           )


setGeneric(name=".createPermMx",
           def=function(x, maxruns) {
             standardGeneric(".createPermMx")}
           )

setGeneric(name="gmpCoef",
           def=function(x){standardGeneric("gmpCoef")})

#setGeneric(name="coef",
#           def=function(x){standardGeneric("coef")})

setGeneric(
           name="permute",
           def=function(x){standardGeneric("permute")}
           )

setGeneric(
           name="getNExceeding",
           def=function(x, y){standardGeneric("getNExceeding")}
           )

setGeneric(
           name="getPValue",
           def=function(x, y){standardGeneric("getPValue")}
           )

setGeneric(
           name=".getIVix",
           def=function(x){standardGeneric(".getIVix")}
           )

setGeneric(
           name=".parseFormula",
           def=function(object, formula){standardGeneric(".parseFormula")}
           )

setGeneric(
           name=".checkMultilevel",
           def=function(object)
           {
             standardGeneric(".checkMultilevel")
           }
           )

setGeneric(
           name=".getDesign",
           def=function(object)
           {
             standardGeneric(".getDesign")
           }
           )

setGeneric(
           name=".buildFitCall",
           def=function(object, arg.exclude=c(), ocall)
           {
             standardGeneric(".buildFitCall")
           }
           )

setGeneric(
           name="fitOnce",
           def=function(object)
           {
             standardGeneric("fitOnce")
           }
           )

setGeneric(
           name="origFit",
           def=function(object)
           {
             standardGeneric("origFit")
           }
           )

setGeneric(
           name=".preparePermScheme",
           def=function(object)
           {
             standardGeneric(".preparePermScheme")
           }
           )

setGeneric(
           name="permSpace",
           def=function(object) {standardGeneric("permSpace")}
           )

setGeneric(
           name="gmpFit",
           def=function(object,gmpControl){standardGeneric("gmpFit")})

setGeneric(name="getPermMx",
           def=function(x){standardGeneric("getPermMx")})

setGeneric(name="coefNames",
           def=function(x){standardGeneric("coefNames")})
