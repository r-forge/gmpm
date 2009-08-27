setClass("Gmp",
         representation(
                        df1="data.frame", # the model.frame
                        dform="formula", # design formula (including covars)
                        mform="formula", # full model formula
                        munit="character", # multilevel sampling unit
                        nunits="numeric", # number of units sampled
                        gmpControl="list", # control of fitting functions
                        fitcall="list", # call to fitting function
                        famtype="character", # type of data
                        DVname="character", # name of DV
                        ivix="numeric", # indices of IVs
                        IVinfo="data.frame", # info about IVs
                        nWithin="numeric", # nWithin unit variables
                        nBetween="numeric", # nBetween unit variables
                        ivWithin="character", # list of within IVs
                        ivBetween="character", # list of between IVs
                        ivars="vector", # list of names of IVs
                        IVcoef="list", # names of factor vars in glm output
                        covars="character", # list of names of covars
                        coefTerms="list", # names of variables from fit output
                                          # w/interactions separated.

                        psBetween="data.frame", # permn scheme (Betw unit IVs)
                        psWithin="list", # permutation scheme (Within unit IVs)
                        pspace="numeric", # size of permutation space
                        pmx="matrix", # matrix of permutation coefficients
                        nCellsPerUnit="numeric", # nCells per sampling unit

                        ncomp="numeric", # n of runs completed
                        "VIRTUAL"), # factor matrix for the model
         
         prototype=prototype(
           nunits=1, nWithin=0, nBetween=0, ncomp=0),
         )

setClass(Class="GmpSummary",
         representation(
                        gmpInfo="list", # misc. info about Gmp object
                        gmpMainSum="list", # list of data frames
                                           # with summary info
                        gmpRegSum="list", # main regression
                        showReg="logical" # whether to show reg coef?
                        ),
         prototype(showReg=FALSE)
         )

setClass(
         Class="Gmp.glm",
         representation(
                        coef0="numeric", # vector of original coefficients
                        family="list"
                        ),
         contains="Gmp"
         )

setClass(
         Class="Gmp.mul",
         representation(
                        coef0="matrix", # vector of original coefficients
                        family="character",
                        pmx="array",
                        convergence="vector" # did it converge?
                        ),
         prototype(family="multinomial",famtype="multinomial"),
         contains="Gmp"
         )

setClass(
         Class="Gmp.user",
         representation(
                        family="character"
                        ),
         prototype(family="user",famtype="user"),
         contains="Gmp"
         )
