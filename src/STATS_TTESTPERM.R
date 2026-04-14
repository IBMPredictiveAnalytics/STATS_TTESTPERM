#/***********************************************************************
# * (C) Copyright Jon K Peck, 2024
# ************************************************************************/

# version 1.0.0

# history
# Mar-2025    Initial version



# helpers
gtxt <- function(...) {
    return(gettext(...,domain="STATS_TTESTPERM"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_TTESTPERM"))
}

loadmsg = "The R %s package is required but could not be loaded."
# tryCatch(suppressWarnings(suppressPackageStartupMessages(library(MKinfer, warn.conflicts=FALSE))), error=function(e){
#     stop(gtxtf(loadmsg,"MKinfer"), call.=FALSE)
    
specialinstall = function() {
    # try to install gmp and arrangements on first use on Intel Mac and other OS's
    # may need update when SPSS goes native on ARM
    # if failure, all is lost
    
    ips = row.names(installed.packages())
    os = Sys.info()["sysname"]
    mtype = Sys.info()["machine"]  #??  mtype == "x86-64" on Intel Mac arm64 on ARM
    ptype = ifelse(os == "Darwin", "mac.binary.big-sur-x86_64", "binary")
    if (!"gmp" %in% ips) {
        print("attempting to install package gmp")
        install.packages("gmp", repos="https://cloud.r-project.org", type=ptype)
    }
    if (!"arrangements" %in% ips) {
        print("attempting to install package arrangements")
        install.packages("arrangements", repos="https://cloud.r-project.org", type=ptype)
    }
}


specialinstall()

tryCatch(suppressWarnings(suppressPackageStartupMessages(library(arrangements, warn.conflicts=FALSE))), error=function(e){
    stop(gtxtf(loadmsg,"arrangements"), call.=FALSE)
}
)
# tryCatch(suppressWarnings(suppressPackageStartupMessages(library(R.utils, warn.conflicts=FALSE))), error=function(e){
#     stop(gtxtf(loadmsg,"R.utils"), call.=FALSE)
# }
# )

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment
    
    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.
        
        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 
        
        if (is.null(msg) || dostop) {
            spssdata.CloseDataConnection()
            lcl$display(inproc)  # display messages and end procedure state

            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any
        

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
                procok = TRUE
            }
        } else {
            procok = inproc
            if (!inproc) {
                procok =tryCatch({
                    spsspkg.StartProcedure(lcl$procname, lcl$omsid)
                    procok = TRUE
                },
                error = function(e) {
                    prockok = FALSE
                }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings and Messages","Warnings", isSplit=FALSE) # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row,
                                               gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory,
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}


casecorrect = function(vlist, vardict, warns) {
    # correct the case of variable names
    # vlist is a list of names, possibly including TO and ALL
    # vardict is a variable dictionary
    # unrecognized names are returned as is as the GetDataFromSPSS api will handle them

    dictnames = vardict["varName",]
    names(dictnames) = tolower(dictnames)
    dictnames['all'] = "all"
    dictnames['to'] = "to"
    correctednames = list()
    for (item in vlist) {
        lcitem = tolower(item)
        itemc = dictnames[[lcitem]]
        if (is.null(itemc)) {
            warns$warn(gtxtf("Invalid variable name: %s", item), dostop=TRUE)
        }
        correctednames = append(correctednames, itemc)
    }
    return(correctednames)
}

procname=gtxt("Perm")
warningsprocname = gtxt("Perm t Test Notes and Warnings")
omsid="STATSTTESTPERM"

# main worker
doperm<-function(variable1, variable2=NULL, groupvar=NULL, groupvals=NULL,
    alternative="two.sided", conflevel=0.95, paired=FALSE, timelimit=Inf,
    equalvar=FALSE, plotx=FALSE, qqplotx=TRUE, nperm=9999
    ) {
    # plotx does not currently work due to ggplot format issues, but it is left
    # in here and disabled in the xml in case that changes in the future
    
    if (alternative == "twosided") {alternative = "two.sided"}
    #DEBUG
    # sink(file="c:/temp/permmout.log", type="output")
    # f = file("c:/temp/permmsgs.log", open="w")
    # sink(file=f, type="message")
    
    domain<-"STATS_TTESTPERM"
    setuplocalization(domain)
    warns = Warn(procname=warningsprocname,omsid=omsid)
    
    spsspkg.StartProcedure(gtxt("Permutation t Test"),"STATSTTESTPERM")
    vardict = spssdictionary.GetDictionaryFromSPSS()
    if (!spsspkg.IsUTF8mode()) {
        warns$warn(gtxt("This procedure requires SPSS to be in Unicode mode"), dostop=TRUE)
    }

    
    if (length(intersect(vardict['varName',], c(variable1, variable2, groupvar))) != 2) {
        warns$warn(gtxt("Variable specification is incorrect"), dostop=TRUE)
    }
    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(gtxt("Case weights are ignored by this procedure"), dostop=FALSE)
    }
    if (paired) {
        if (is.null(variable2)) {
            warns$warn(gtxt("Two variables are required for paired tests"), dostop=TRUE)
        }
        if (!is.null(groupvar)) {
            warns$warn(gtxt("The group variable is not used for paired tests and will be ignored"),
                dostop=FALSE)
        }
    } else {
        if (is.null(groupvar)) {
            warns$warn(gtxt("The group variable is required for this test"), dostop=TRUE)
        }
        if (length(groupvals) != 2) {
            warns$warn(gtxt("Exactly two values must be supplied for the group specification"), 
                dostop=TRUE)
        }
        if (vardict["varMeasurementLevel", which(vardict["varName",] == groupvar)] == "scale") {
            warns$warn(gtxt("The group variable must have a nominal or ordinal measurement level."))
        }
    }

    allvars = c(variable1, variable2, groupvar)
    # check for valid names in R
    tryCatch(
        {
        for (ch in allvars) {
            xxx = str2lang(ch)
        }
        },
        error = function(e) {warns$warn(gtxtf("%s is not a valid variable name in R.  Please rename it to use this procedure.", ch),
            dostop=TRUE)}
    )

    # the get data api requires case match
    variable1 = unlist(casecorrect(list(variable1), vardict, warns))
    if (!is.null(variable2)) {variable2 = unlist(casecorrect(list(variable2), vardict, warns))}
    if (!is.null(groupvar)) {groupvar = unlist(casecorrect(list(groupvar), vardict, warns))}

    first = TRUE
    scount = 0
    while (!spssdata.IsLastSplit()) {
        scount = scount + 1
    tryCatch(
        {
        dta = spssdata.GetSplitDataFromSPSS(allvars, missingValueToNA=TRUE, factorMode="levels",
            keepUserMissing=FALSE)
        },
        error=function(e) {warns$warn(paste(gtxt("error fetching data"), e, sep="\n"), dostop=TRUE)}
    )

    dta = dta[complete.cases(dta),]
    # ensure that data or split uses only cases defined by groupvals for nonpaired test
    if (!paired) {
        groupvals = as.character(groupvals)
        dta = subset(dta, as.character(dta[[2]]) %in% groupvals)
        tt = table(dta[2])
        ss = sum(tt > 0)
        if (sum(tt > 0) != 2) {
            warns$warn(gtxtf("There are no complete or eligible cases in the dataset or split number %s .  Skipping",
                scount), dostop=FALSE)
            next
        }
        if (any(tt == 1)) {
            warns$warn(gtxtf("There are too few cases in a group in the dataset or split number %s.  Skipping",
            scount), dostop=FALSE)
            next
        }
        }
    starttime = Sys.time()
    # tryCatch(withTimeout(
    #     {
    tryCatch(
        {
        # if the dataset or split has no rows, an error will be generated but processing
        # it will continue with tne next split if any
            
        res = NULL
        if (paired) {
            res = perm.t.test(x=dta[[1]], y=dta[[2]], alternative=alternative,
                paired = paired, var.equal=equalvar, conf.level=conflevel,
                R=nperm, permStat=TRUE)
        } else {
            ff = as.formula(sprintf("%s ~ %s", variable1, groupvar))
            res = perm.t.test(ff,  alternative=alternative,
                  paired = paired, var.equal=equalvar, conf.level=conflevel,
                  R=nperm, data=dta, permStat=TRUE)
        }
        },
        error = function(e) {
            warns$warn(paste(gtxt("error estimating equation"), e, sep="\n"), dostop=FALSE)
        },
        warning = function(w) {
            warns$warn(paste(gtxt("error estimating equation"), w, sep="\n"), dostop=FALSE)
        }
    )
    #     },
    #     timeout=timelimit), error = function(e) {
    #         warns$warn(sprintf("Procedure stopped.  Time limit of %s seconds has expired",
    #                            timelimit))}
    # )
    ###save(res, paired, variable1, variable2, groupvar, dta, file="c:/temp/afterest.rdata")    
    if (is.null(res)) {
        next
    }
    displaytables(res, nrow(dta), paired, variable1, variable2, groupvar, groupvals, equalvar, 
        conflevel, vardict, starttime, first, warns)
    first = FALSE
    # The plot function produces a ggplot object, but the SPSS Viewer
    # cannot handle that object.
    # if (plotx) {
    #     tryCatch(
    #         {
    #         h0plot(res, lwd=1.5, main="")
    #         },
    #     error = function(e) {warns$warn(e, dostop=FALSE)}
    #     )
    # }
    if (qqplotx) {
        tryCatch(
            {
            doqqplot(dta, groupvar, groupvals, vardict)
            }, 
            error = function(e) {warns$warn(e, dostop=FALSE)}
        )
    }
    }   
    
    spssdata.CloseDataConnection()
    warns$display(inproc=FALSE)
    
    # DEBUG
    # sink(file=NULL, type="output")
    # sink(file=NULL, type="message")
}

r4x = function(x) {ifelse(x < .001, "<.001", sprintf("%.4f", x))}
r4 = function(x) {ifelse(is.infinite(x), x, sprintf("%.4f", x))}

displaytables = function(res, nrows, paired, variable1, variable2, groupvar, groupvals, equalvar, 
    conflevel, vardict, starttime, first, warns) {
    # display all the tables
    
    # in some cases, the test result may be NULL without an error message
    # likely to be more of an issue with splits
    if (is.null(res)) {
        return()
    }
    
    if (!is.null(groupvar)) {
        catdict = spssdictionary.GetCategoricalDictionaryFromSPSS(groupvar)
        lvls = catdict$dictionary[[1]]$levels
        lbls = catdict$dictionary[[1]]$labels
        # NA if value not in dict
        index1 = match(groupvals[[1]], catdict$dictionary[[1]]$levels)
        label1 = catdict$dictionary[[1]]$labels[[index1]]
        index2 = match(groupvals[[2]], catdict$dictionary[[1]]$levels)
        label2 = catdict$dictionary[[1]]$labels[[index2]]
    }
        
    # specifications table
    if (first) {
        dd = data.frame()
        if (paired) {
            dd[gtxt("First Paired Variable"),1] = variable1
            dd[gtxt("Second Paired Variable"),1] = variable2
        } else {
            dd[gtxt("Test Variable"),1] = variable1
            dd[gtxt("Group Variable"),1] = groupvar
            dd[gtxt("First Group"),1] = label1   #groupvals[[1]]
            dd[gtxt("Second Group"),1] = label2  #groupvals[[2]]
        }
        dd[gtxt("Alternative Hypothesis"),1] = res$alternative   #TODO translation
        dd[gtxt("Equal Variance Assumed"),1] = ifelse(equalvar, gtxt("Yes"), gtxt("No"))
        dd[gtxt("Number of Permutations"),1] = res$R
        names(dd) = gtxt("Settings")
        caption=gtxtf("Results computed by module MKinfer")
        # which group is the first is unclear
        # if (!is.null(groupvar) && res$alternative != "two.sided") {
        #     caption = paste(caption, gtxt("The group with the smaller mean is considered to be the first."),
        #         sep = "\n")
        # }
        spsspivottable.Display(
            dd,
            title=gtxt("t Test Settings"),
            isSplit=FALSE,
            templateName="STATSPERMSETTINGS",
            caption=caption
        )
    }
    
    dd = data.frame()
    dd[gtxtf("Number of Cases"),1] = nrows
    dd[gtxt("t Statistic"),1] = r4(res$statistic)
    dd[gtxt("Degrees of Freedom"),1] = r4(res$parameter)
    dd[gtxt("P Value"),1] = r4x(res$p.value)
    if (paired) {
        dd[gtxt("Estimated Mean Difference"),1] = r4(res$estimate)
    } 
    else {
        dd[gtxt("Estimated Mean: Group 1"),1] = r4(res$estimate[[1]])
        dd[gtxt("Estimated Mean: Group 2"),1] = r4(res$estimate[[2]])
    }
    dd[gtxt("Lower Difference Confidence Interval"),1] = r4(res$conf.int[[1]])
    dd[gtxt("Upper Difference Confifdence Interval"),1] = r4(res$conf.int[[2]])
    dd[gtxt("Permutation P Value"),1] = r4x(res$perm.p.value)
    dd[gtxt("Permutation Estimated Mean Difference"),1] = r4(res$perm.estimate)
    dd[gtxt("Lower Permutation Confidence Interval"),1] = r4(res$perm.conf.int[[1]])
    dd[gtxt("Upper Permutation Confidence Interval"),1] = r4(res$perm.conf.int[[2]])
    dd[gtxt("Standardized Permutation Effect Size"),1] = r4(p2ses(res$perm.p.value, alternative=res$alternative))
    ###dd[gtxt("Elapsed Time (per split) in Seconds"),1] = as.integer(Sys.time() - starttime)

    names(dd) = gtxt("Statistics")
    spsspivottable.Display(
        dd,
        title=gtxt("Permutation t Tests"),
        templateName="STATSPERMTESTS"
    )
}

doqqplot = function(dta, groupvar, groupvals, vardict) {
    #do qqplot
    # dta[1:2] will be variables to plot if no groupvar and one variable plus group variable otherwise
    
    vnames = names(dta)
    var1col = which(vardict['varName',] == vnames[1])
    if (!is.null(groupvar)) {  # one-variable plus group plot
        varg = which(vardict['varName',] == vnames[2])
        vardict['varName',]
        
        ylabel = paste(ifelse(vardict[["varLabel", var1col]] == "", vnames[[1]], vardict[["varLabel", var1col]]), 
                       gtxtf(" (%s)", groupvals[1]))
        xlabel = paste(ifelse(vardict[["varLabel", var1col]] == "", vnames[[2]], vardict[["varLabel", var1col]]), 
                       gtxtf(" (%s)", groupvals[2]))
        dta = split(dta[[1]], as.character(dta[[2]]))
            title = gtxt("Two Group Q-Q Plot")
    }
    else {  # two-variable plot
        ylabel = ifelse(vardict[["varLabel", var1col]] == "", vnames[[1]], vardict[["varLabel", var1col]])
        var2col = which(vardict["varName", ] == vnames[2])
        xlabel = ifelse(vardict[["varLabel", var2col]] == "", vnames[[2]], vardict[["varLabel", var2col]]) 
        title = gtxt("Q-Q Plot")
    }
    
    rr = range(as.matrix(dta), na.rm=TRUE, finite=TRUE)
    qqplot(dta[[2]], dta[[1]], ylab=ylabel, xlab=xlabel, 
           col="blue", pch=16, xlim=rr, ylim=rr,
           main=title)
    abline(a=0, b=1, lwd=2)
}


setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    if (!is.null(fpath)) {
        bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
    }
} 


Run<-function(args){

    cmdname = args[[1]]
    args <- args[[2]]

    # variable keywords are typed as varname instead of existingvarlist in
    # order to allow for case correction of names later, since the data fetching apis are
    # case sensitive

    oobj <- spsspkg.Syntax(templ=list(
        spsspkg.Template("VARIABLE1", subc="", ktype="varname", var="variable1", islist=FALSE),
        spsspkg.Template("VARIABLE2", subc="", ktype="varname", var="variable2", islist=FALSE),
        spsspkg.Template("GROUPVAR", subc="", ktype="varname", var="groupvar", islist=FALSE),
        spsspkg.Template("GROUPVALS", subc="", ktype="literal", var="groupvals", islist=TRUE),
        spsspkg.Template("ALTERNATIVE", subc="", ktype="str", var="alternative",
            vallist=list("twosided", "greater", "less"), islist=FALSE),
        spsspkg.Template("CONFLEVEL", subc="", ktype="float", var="conflev", 
            vallist=list(.00001, .99999)),
        spsspkg.Template("PAIRED", subc="", ktype="bool", var="paired", islist=FALSE),
        spsspkg.Template("EQUALVAR", subc="", ktype="bool", var="equalvar", islist=FALSE),
        spsspkg.Template("DISTPLOT", subc="", ktype="bool", var="plotx", islist=FALSE),
        spsspkg.Template("QQPLOT", subc="", ktype="bool", var="qqplot", islist=FALSE),
        spsspkg.Template("NPERM", subc="", ktype="int", var="nperm", 
            vallist=list(2, 99999), islist=FALSE),
        spsspkg.Template("TIMELIMIT", subc="", ktype="int", var="timelimit", islist=FALSE)
        ))

    if ("HELP" %in% attr(args,"names"))
        #writeLines(helptext)
        helper(cmdname)
    else {
        spsspkg.processcmd(oobj, args, "doperm")
    }
}


helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }

    if (exists("spsspkg.helper")) {
        assign("helper", spsspkg.helper)
    }
}

### following code extracted from the MKinfer to avoid
### the many unused dependencies in the MKinfer module

p2ses <- function(p, alternative = c("two.sided", "less", "greater")){
    alternative <- match.arg(alternative)
    stopifnot(all(p <= 1))
    stopifnot(all(p >= 0))
    
    if(alternative == "less"){
        ses <- qnorm(p)  
    }else if (alternative == "greater"){
        ses <- qnorm(p, lower.tail = FALSE)
    }else{
        ses <- qnorm(p/2, lower.tail = FALSE) 
    }
    ses
}

sample.perm <- function(x, k = NULL, R = 1000, replace = FALSE){
    if(is.null(k)) k <- length(x)
    max.R <- try(npermutations(x, k = k, replace = replace), silent = TRUE)
    if(inherits(max.R, "try-error")) max.R <- Inf
    if(max.R < R){
        warning("The requested number of permutations (", R, ") is larger than ",
                "the total number of possible permutations (", max.R, ").\n",
                "Hence all possible permutations are computed.")
        res <- permutations(x, k = k, replace = replace)
    }else{
        res <- NULL
        iter <- R
        repeat{
            res <- rbind(res, 
                         t(replicate(iter, sample(x, size = k, 
                                                  replace = replace))))
            if(all(!duplicated(res))) break
            iter <- R - sum(!duplicated(res))
            res <- res[!duplicated(res),]
        }
    }
    res
}
perm.t.test <- function (x, ...){ 
    UseMethod("perm.t.test")
}
perm.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                                mu = 0, paired = FALSE, var.equal = FALSE, 
                                conf.level = 0.95, R = 9999, symmetric = TRUE, 
                                permStat = FALSE, ...){
    alternative <- match.arg(alternative)
    if(!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if(!is.null(y)){
        dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
        if (paired) 
            xok <- yok <- complete.cases(x, y)
        else{
            yok <- !is.na(y)
            xok <- !is.na(x)
        }
        y <- y[yok]
    }else{
        dname <- deparse1(substitute(x))
        if (paired) 
            stop("'y' is missing for paired test")
        xok <- !is.na(x)
        yok <- NULL
    }
    x <- x[xok]
    if(paired){
        x <- x - y
        y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    if (is.null(y)) {
        if (nx < 2) 
            stop("not enough 'x' observations")
        df <- nx - 1
        stderr <- sqrt(vx/nx)
        if (stderr < 10 * .Machine$double.eps * abs(mx)) 
            stop("data are essentially constant")
        tstat <- (mx - mu)/stderr
        method <- if (paired) "Permutation Paired t-test" else "Permutation One Sample t-test"
        estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
        x.cent <- x - mx
        X <- abs(x.cent)*sample.perm(c(-1,1), k = nx, R = R, replace = TRUE)
        R.true <- nrow(X)
        MX <- rowMeans(X)
        VX <- rowSums((X-MX)^2)/(nx-1)
        STDERR <- sqrt(VX/nx)
        perm.stderr <- mean(STDERR)
        TSTAT <- MX/STDERR
        EFF <- MX+mx
        perm.estimate <- mean(EFF)
        if(paired){
            names(perm.estimate) <- "permutation mean of the differences" 
        }else{
            names(perm.estimate) <- "permutation mean of x"
        } 
    }else{
        ny <- length(y)
        if(nx < 1 || (!var.equal && nx < 2)) 
            stop("not enough 'x' observations")
        if(ny < 1 || (!var.equal && ny < 2)) 
            stop("not enough 'y' observations")
        if(var.equal && nx + ny < 3) 
            stop("not enough observations")
        my <- mean(y)
        vy <- var(y)
        method <- paste("Permutation", paste(if (!var.equal) "Welch", "Two Sample t-test"))
        estimate <- c(mx, my)
        names(estimate) <- c("mean of x", "mean of y")
        z <- c(x, y)
        Z <- sample.perm(z, k = nx+ny, R = R)
        R.true <- nrow(Z)
        X <- Z[,1:nx]
        Y <- Z[,(nx+1):(nx+ny)]
        MX <- rowMeans(X)
        MY <- rowMeans(Y)
        EFF <- (MX+mx) - (MY+my)
        if(var.equal){
            df <- nx + ny - 2
            v <- 0
            if (nx > 1) 
                v <- v + (nx - 1) * vx
            if (ny > 1) 
                v <- v + (ny - 1) * vy
            v <- v/df
            stderr <- sqrt(v * (1/nx + 1/ny))
            V <- (rowSums((X-MX)^2) + rowSums((Y-MY)^2))/df
            STDERR <- sqrt(V*(1/nx + 1/ny))
        }else{
            stderrx <- sqrt(vx/nx)
            stderry <- sqrt(vy/ny)
            stderr <- sqrt(stderrx^2 + stderry^2)
            df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
            VX <- rowSums((X-MX)^2)/(nx-1)
            VY <- rowSums((Y-MY)^2)/(ny-1)
            STDERR <- sqrt(VX/nx + VY/ny)
        }
        perm.stderr <- mean(STDERR)
        perm.estimate <- mean(EFF) 
        names(perm.estimate) <- "permutation difference of means"
        if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
            stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
        TSTAT <- (MX - MY)/STDERR
    }
    if (alternative == "less") {
        pval <- pt(tstat, df)
        perm.pval <- max(mean(TSTAT < tstat), 1/R.true)
        cint <- c(-Inf, tstat + qt(conf.level, df))
        perm.cint <- c(-Inf, quantile(EFF, conf.level))
    }else if(alternative == "greater") {
        perm.pval <- max(mean(TSTAT > tstat), 1/R.true)
        pval <- pt(tstat, df, lower.tail = FALSE)
        cint <- c(tstat - qt(conf.level, df), Inf)
        perm.cint <- c(quantile(EFF, 1-conf.level), Inf)
    }else{
        pval <- 2 * pt(-abs(tstat), df)
        if(symmetric)
            perm.pval <- max(mean(abs(TSTAT) > abs(tstat)), 1/R.true)
        else
            perm.pval <- max(2*min(mean(TSTAT <= tstat), mean(TSTAT > tstat)), 1/R.true)
        alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
        cint <- tstat + c(-cint, cint)
        perm.cint <- quantile(EFF, c(alpha/2, 1-alpha/2))
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y)) "difference in means" else "mean"
    attr(cint, "conf.level") <- conf.level
    attr(perm.cint, "conf.level") <- conf.level
    if(permStat){ 
        perm.statistic <- TSTAT
    }else{
        perm.statistic <- NULL
    }
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
                 perm.p.value = perm.pval, R = R, R.true = R.true, 
                 p.min = perm.pval == 1/R.true,
                 conf.int = cint, perm.conf.int = perm.cint,
                 estimate = estimate, perm.estimate = perm.estimate, 
                 null.value = mu, stderr = stderr, perm.stderr = perm.stderr,
                 alternative = alternative, method = method, data.name = dname,
                 perm.statistic = perm.statistic)
    class(rval) <- c("perm.htest", "htest")
    rval
}
perm.t.test.formula <- function (formula, data, subset, na.action, ...){
    if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                    "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L) 
        stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("perm.t.test", c(DATA, list(...)))
    y$data.name <- DNAME
    if (length(y$estimate) == 2L) 
        names(y$estimate) <- paste("mean in group", levels(g))
    y
}
print.perm.htest <- function (x, digits = getOption("digits"), prefix = "\t", ...) {
    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    cat("number of permutations:  ", x$R.true, "\n", sep = "")
    out <- character()
    if (!is.null(x$perm.p.value)) {
        bfp <- format.pval(x$perm.p.value, digits = max(1L, digits - 3L))
        if(x$R.true < x$R){
            if(x$p.min){
                cat("(Exact) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("<", bfp), "\n")
            }else{
                cat("(Exact) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
            }
        }else{
            if(x$p.min){
                cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("<", bfp), "\n")
            }else{
                cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
            }
        }
    }
    if (!is.null(x$perm.estimate)) {
        cat(paste(names(x$perm.estimate), "(SE) =", 
                  format(x$perm.estimate, digits = digits),
                  paste("(", format(x$perm.stderr, digits = digits), ")", sep = "")), 
            "\n")
    }
    if (!is.null(x$perm.conf.int)) {
        if(x$R.true < x$R){
            cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
                " percent (exact) permutation percentile confidence interval:\n", 
                " ", paste(format(x$perm.conf.int[1:2], digits = digits), 
                           collapse = " "), "\n", sep = "")
        }else{
            cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
                " percent (Monte-Carlo) permutation percentile confidence interval:\n", 
                " ", paste(format(x$perm.conf.int[1:2], digits = digits), 
                           collapse = " "), "\n", sep = "")
        }
    }
    cat("\nResults without permutation:\n")
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(x$statistic, 
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$parameter)) 
        out <- c(out, paste(names(x$parameter), "=", format(x$parameter, 
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = max(1L, digits - 
                                                      3L))
        out <- c(out, paste("p-value", 
                            if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative, two.sided = "not equal to", 
                                   less = "less than", greater = "greater than")
                cat("true ", names(x$null.value), " is ", alt.char, 
                    " ", x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, digits = digits, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n", 
            " ", paste(format(x$conf.int[1:2], digits = digits), 
                       collapse = " "), "\n", sep = "")
    }
    if (!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, digits = digits, ...)
    }
    cat("\n")
    invisible(x)
}

sample.perm <- function(x, k = NULL, R = 1000, replace = FALSE){
    if(is.null(k)) k <- length(x)
    max.R <- try(npermutations(x, k = k, replace = replace), silent = TRUE)
    if(inherits(max.R, "try-error")) max.R <- Inf
    if(max.R < R){
        warning("The requested number of permutations (", R, ") is larger than ",
                "the total number of possible permutations (", max.R, ").\n",
                "Hence all possible permutations are computed.")
        res <- permutations(x, k = k, replace = replace)
    }else{
        res <- NULL
        iter <- R
        repeat{
            res <- rbind(res, 
                         t(replicate(iter, sample(x, size = k, 
                                                  replace = replace))))
            if(all(!duplicated(res))) break
            iter <- R - sum(!duplicated(res))
            res <- res[!duplicated(res),]
        }
    }
    res
}
perm.t.test <- function (x, ...){ 
    UseMethod("perm.t.test")
}
perm.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                                mu = 0, paired = FALSE, var.equal = FALSE, 
                                conf.level = 0.95, R = 9999, symmetric = TRUE, 
                                permStat = FALSE, ...){
    alternative <- match.arg(alternative)
    if(!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if(!is.null(y)){
        dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
        if (paired) 
            xok <- yok <- complete.cases(x, y)
        else{
            yok <- !is.na(y)
            xok <- !is.na(x)
        }
        y <- y[yok]
    }else{
        dname <- deparse1(substitute(x))
        if (paired) 
            stop("'y' is missing for paired test")
        xok <- !is.na(x)
        yok <- NULL
    }
    x <- x[xok]
    if(paired){
        x <- x - y
        y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    if (is.null(y)) {
        if (nx < 2) 
            stop("not enough 'x' observations")
        df <- nx - 1
        stderr <- sqrt(vx/nx)
        if (stderr < 10 * .Machine$double.eps * abs(mx)) 
            stop("data are essentially constant")
        tstat <- (mx - mu)/stderr
        method <- if (paired) "Permutation Paired t-test" else "Permutation One Sample t-test"
        estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
        x.cent <- x - mx
        X <- abs(x.cent)*sample.perm(c(-1,1), k = nx, R = R, replace = TRUE)
        R.true <- nrow(X)
        MX <- rowMeans(X)
        VX <- rowSums((X-MX)^2)/(nx-1)
        STDERR <- sqrt(VX/nx)
        perm.stderr <- mean(STDERR)
        TSTAT <- MX/STDERR
        EFF <- MX+mx
        perm.estimate <- mean(EFF)
        if(paired){
            names(perm.estimate) <- "permutation mean of the differences" 
        }else{
            names(perm.estimate) <- "permutation mean of x"
        } 
    }else{
        ny <- length(y)
        if(nx < 1 || (!var.equal && nx < 2)) 
            stop("not enough 'x' observations")
        if(ny < 1 || (!var.equal && ny < 2)) 
            stop("not enough 'y' observations")
        if(var.equal && nx + ny < 3) 
            stop("not enough observations")
        my <- mean(y)
        vy <- var(y)
        method <- paste("Permutation", paste(if (!var.equal) "Welch", "Two Sample t-test"))
        estimate <- c(mx, my)
        names(estimate) <- c("mean of x", "mean of y")
        z <- c(x, y)
        Z <- sample.perm(z, k = nx+ny, R = R)
        R.true <- nrow(Z)
        X <- Z[,1:nx]
        Y <- Z[,(nx+1):(nx+ny)]
        MX <- rowMeans(X)
        MY <- rowMeans(Y)
        EFF <- (MX+mx) - (MY+my)
        if(var.equal){
            df <- nx + ny - 2
            v <- 0
            if (nx > 1) 
                v <- v + (nx - 1) * vx
            if (ny > 1) 
                v <- v + (ny - 1) * vy
            v <- v/df
            stderr <- sqrt(v * (1/nx + 1/ny))
            V <- (rowSums((X-MX)^2) + rowSums((Y-MY)^2))/df
            STDERR <- sqrt(V*(1/nx + 1/ny))
        }else{
            stderrx <- sqrt(vx/nx)
            stderry <- sqrt(vy/ny)
            stderr <- sqrt(stderrx^2 + stderry^2)
            df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
            VX <- rowSums((X-MX)^2)/(nx-1)
            VY <- rowSums((Y-MY)^2)/(ny-1)
            STDERR <- sqrt(VX/nx + VY/ny)
        }
        perm.stderr <- mean(STDERR)
        perm.estimate <- mean(EFF) 
        names(perm.estimate) <- "permutation difference of means"
        if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
            stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
        TSTAT <- (MX - MY)/STDERR
    }
    if (alternative == "less") {
        pval <- pt(tstat, df)
        perm.pval <- max(mean(TSTAT < tstat), 1/R.true)
        cint <- c(-Inf, tstat + qt(conf.level, df))
        perm.cint <- c(-Inf, quantile(EFF, conf.level))
    }else if(alternative == "greater") {
        perm.pval <- max(mean(TSTAT > tstat), 1/R.true)
        pval <- pt(tstat, df, lower.tail = FALSE)
        cint <- c(tstat - qt(conf.level, df), Inf)
        perm.cint <- c(quantile(EFF, 1-conf.level), Inf)
    }else{
        pval <- 2 * pt(-abs(tstat), df)
        if(symmetric)
            perm.pval <- max(mean(abs(TSTAT) > abs(tstat)), 1/R.true)
        else
            perm.pval <- max(2*min(mean(TSTAT <= tstat), mean(TSTAT > tstat)), 1/R.true)
        alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
        cint <- tstat + c(-cint, cint)
        perm.cint <- quantile(EFF, c(alpha/2, 1-alpha/2))
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y)) "difference in means" else "mean"
    attr(cint, "conf.level") <- conf.level
    attr(perm.cint, "conf.level") <- conf.level
    if(permStat){ 
        perm.statistic <- TSTAT
    }else{
        perm.statistic <- NULL
    }
    rval <- list(statistic = tstat, parameter = df, p.value = pval, 
                 perm.p.value = perm.pval, R = R, R.true = R.true, 
                 p.min = perm.pval == 1/R.true,
                 conf.int = cint, perm.conf.int = perm.cint,
                 estimate = estimate, perm.estimate = perm.estimate, 
                 null.value = mu, stderr = stderr, perm.stderr = perm.stderr,
                 alternative = alternative, method = method, data.name = dname,
                 perm.statistic = perm.statistic)
    class(rval) <- c("perm.htest", "htest")
    rval
}
perm.t.test.formula <- function (formula, data, subset, na.action, ...){
    if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                    "term.labels")) != 1L)) 
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L) 
        stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("perm.t.test", c(DATA, list(...)))
    y$data.name <- DNAME
    if (length(y$estimate) == 2L) 
        names(y$estimate) <- paste("mean in group", levels(g))
    y
}
print.perm.htest <- function (x, digits = getOption("digits"), prefix = "\t", ...) {
    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    cat("number of permutations:  ", x$R.true, "\n", sep = "")
    out <- character()
    if (!is.null(x$perm.p.value)) {
        bfp <- format.pval(x$perm.p.value, digits = max(1L, digits - 3L))
        if(x$R.true < x$R){
            if(x$p.min){
                cat("(Exact) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("<", bfp), "\n")
            }else{
                cat("(Exact) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
            }
        }else{
            if(x$p.min){
                cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("<", bfp), "\n")
            }else{
                cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
            }
        }
    }
    if (!is.null(x$perm.estimate)) {
        cat(paste(names(x$perm.estimate), "(SE) =", 
                  format(x$perm.estimate, digits = digits),
                  paste("(", format(x$perm.stderr, digits = digits), ")", sep = "")), 
            "\n")
    }
    if (!is.null(x$perm.conf.int)) {
        if(x$R.true < x$R){
            cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
                " percent (exact) permutation percentile confidence interval:\n", 
                " ", paste(format(x$perm.conf.int[1:2], digits = digits), 
                           collapse = " "), "\n", sep = "")
        }else{
            cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
                " percent (Monte-Carlo) permutation percentile confidence interval:\n", 
                " ", paste(format(x$perm.conf.int[1:2], digits = digits), 
                           collapse = " "), "\n", sep = "")
        }
    }
    cat("\nResults without permutation:\n")
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(x$statistic, 
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$parameter)) 
        out <- c(out, paste(names(x$parameter), "=", format(x$parameter, 
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = max(1L, digits - 
                                                      3L))
        out <- c(out, paste("p-value", 
                            if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative, two.sided = "not equal to", 
                                   less = "less than", greater = "greater than")
                cat("true ", names(x$null.value), " is ", alt.char, 
                    " ", x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, digits = digits, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n", 
            " ", paste(format(x$conf.int[1:2], digits = digits), 
                       collapse = " "), "\n", sep = "")
    }
    if (!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, digits = digits, ...)
    }
    cat("\n")
    invisible(x)
}

