# Some Rcmdr dialogs for the TeachingDemos package

# last modified: 2015-12-08 by J. Fox

# Note: the following function (with contributions from Richard Heiberger and Milan Bouchet-Valat)
# can be included in any Rcmdr plug-in package to cause the package to load
# the Rcmdr if it is not already loaded

.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    putRcmdr("slider.env", new.env())
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if (!pkgname %in% plugins) {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
        if("package:Rcmdr" %in% search()) {
            if(!getRcmdr("autoRestart")) {
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
            }
        }
        else {
            Commander()
        }
    }
}


inspect <- function(data, x, ...)UseMethod("inspect")
inspect.tbl_df <- function(data, x){
  tmp <- data[[as.character(x)]]
  var.lab <- attr(tmp, "label")
  if(is.null(var.lab)){var.lab <- "No Label Found"}
  val.labs <- attr(tmp, "labels")
  if(is.null(val.labs)){val.labs <- sort(unique(tmp))}
  tab <- cbind(freq = table(tmp), prop = round(table(tmp)/sum(table(tmp), na.rm=T), 3))
  out <- list(variable_label = var.lab, value_labels=t(t(val.labs)), freq_dist = tab)
  return(out)
}
inspect.data.frame <- function(data, x){
  var.lab <- attr(data, "var.label")[which(names(data) == x)]
  val.labs <- if(!is.null(levels(data[[x]]))){levels(data[[x]])}
    else {sort(unique(data[[x]]))}
  tab <- cbind(freq = table(data[[x]]), prop = round(table(data[[x]])/sum(table(data[[x]]), na.rm=T), 3))
  out <- list(variable_label = var.lab, value_labels=t(t(val.labs)), freq_dist = tab)
  return(out)
}

## concordant, discordant, tau.b, tau.c, ord.somers.d, ord.gamma come from the ryouready package
## Phi and V come from the DescTools package
concordant <- function (x) {
    x <- matrix(as.numeric(x), dim(x))
    mat.lr <- function(r, c) {
        lr <- x[(r.x > r) & (c.x > c)]
        sum(lr)
    }
    r.x <- row(x)
    c.x <- col(x)
    sum(x * mapply(mat.lr, r = r.x, c = c.x))
}
discordant <- function(x){
    x <- matrix(as.numeric(x), dim(x))
    mat.ll <- function(r, c) {
        ll <- x[(r.x > r) & (c.x < c)]
        sum(ll)
    }
    r.x <- row(x)
    c.x <- col(x)
    sum(x * mapply(mat.ll, r = r.x, c = c.x))
}

tau.b <- function (x) {
    x <- matrix(as.numeric(x), dim(x))
    c <- concordant(x)
    d <- discordant(x)
    n <- sum(x)
    SumR <- rowSums(x)
    SumC <- colSums(x)
    tau.b <- (2 * (c - d))/sqrt(((n^2) - (sum(SumR^2))) * ((n^2) -
        (sum(SumC^2))))
    tau.b
}

ord.gamma <- function(x){
    x <- matrix(as.numeric(x), dim(x))
    c <- concordant(x)
    d <- discordant(x)
    gamma <- (c - d)/(c + d)
    class(gamma) <- "ord.gamma"
    gamma
}

ord.somers.d <- function(x){
    x <- matrix(as.numeric(x), dim(x))
    c <- concordant(x)
    d <- discordant(x)
    n <- sum(x)
    SumR <- rowSums(x)
    SumC <- colSums(x)
    sd.cr <- (2 * (c - d))/((n^2) - (sum(SumR^2)))
    sd.rc <- (2 * (c - d))/((n^2) - (sum(SumC^2)))
    sd.s <- (2 * (c - d))/((n^2) - (((sum(SumR^2)) + (sum(SumC^2)))/2))
    res <- list(sd.cr, sd.rc, sd.s)
    names(res) <- c("sd.cr", "sd.rc", "sd.symmetric")
    class(res) <- "ord.somersd"
    res
}

lambda <- function(x){
  wmax <- apply(x, 2, which.max)
  wgmax <- which.max(rowSums(x))
  nullcc <- rowSums(x)[wgmax]
  nullerr <- sum(rowSums(x)[-wgmax])
  corrpred <- x[cbind(wmax, 1:ncol(x))]
  errpred <- colSums(x) - corrpred
  E1 <- nullerr
  E2 <- sum(errpred)
  (E1-E2)/E1
}

phi <- function(x){
   as.numeric(sqrt(suppressWarnings(chisq.test(x, correct = FALSE)$statistic)/sum(x)))
}

V <- function(x){
  chi2 <- chisq.test(x)$statistic
  sqrt(chi2/(sum(c(x)) * (min(nrow(x), ncol(x)) -1)))
}

simtable <- function(x,y, n=1000, stat=NULL){
  out <- lapply(1:n, function(i)table(x, sample(y, length(y), replace=F)))
  if(is.null(stat)){
    return(out)
  }
  else{
    sapply(out, stat)
  }

}

simrho <- function(x,y, n=1000){
  rho0 <- cor(x,y, use="pair", method="spearman")
  simrho <- sapply(1:n, function(i)cor(x, sample(y, length(y), replace=F), use="pair", method="spearman"))
  return(list(rho0 = rho0, simrho = simrho, pv = mean(abs(simrho) > abs(rho0))))
}

makeStats <- function(x,y, chisq=FALSE, phi=FALSE, cramersV=FALSE, lambda=FALSE,
   gamma=FALSE, d=FALSE, taub=FALSE, rho=FALSE, n=1000){

  tabs <- simtable(x,y,n)
  tab <- table(x,y)
allStats <- NULL
if(chisq){
  stat0 <- do.call('chisq.test', list(x=tab))$statistic
  stats <- sapply(tabs, function(x)chisq.test(x)$statistic)
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Chi-squared"
}
if(phi){
  stat0 <- do.call('phi', list(x=tab))
  stats <- sapply(tabs, function(x)phi(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Phi"
}
if(cramersV){
  stat0 <- do.call('V', list(x=tab))
  stats <- sapply(tabs, function(x)V(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Cramers V"
}
if(lambda){
  stat0 <- do.call('lambda', list(x=tab))
  stats <- sapply(tabs, function(x)lambda(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Lambda"
}
if(gamma){
  stat0 <- do.call('ord.gamma', list(x=tab))
  stats <- sapply(tabs, function(x)lambda(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Kruskal-Goodman Gamma"
}
if(d){
  stat0 <- do.call('ord.somers.d', list(x=tab))$sd.symmetric
  stats <- sapply(tabs, function(x)ord.somers.d(x)$sd.symmetric)
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Somers D"
}
if(taub){
  stat0 <- do.call('tau.b', list(x=tab))
  stats <- sapply(tabs, function(x)tau.b(x))
  pv <- mean(stats > stat0)
  allStats <- rbind(allStats, c(stat0, pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Tau-b"
}
if(rho){
  x2 <-as.numeric(x)
  y2 <- as.numeric(y)
  r <- simrho(x2,y2,n)
  allStats <- rbind(allStats, c(r$rho0, r$pv))[,,drop=F]
  rownames(allStats)[nrow(allStats)] <- "Spearmans Rho"
}
if(!is.null(allStats)){
  colnames(allStats) <- c("statistic", "p-value")
  w <- which(allStats[,1] == 0 & allStats[,2] == 0)
  if(length(w) > 0){
    allStats[w,2] <- 1.000
  }
  allStats <- round(allStats, 4)
}
return(allStats)
}

plotStdRes <- function(x, col=RColorBrewer::brewer.pal(10, "RdBu")){
  x2 <- chisq.test(x)
  res <- x2$stdres
  minres <- ifelse(min(c(res)) > -3.5, -3.5, min(c(res)))
  maxres <- ifelse(max(c(res)) <3.5, 3.5, max(c(res)))
  if(maxres > abs(minres)){
    minres <- -maxres
  }
  if(maxres < abs(minres)){
    maxres <- -minres
  }
  lattice::levelplot(res, col.regions=col, cuts=10, at=c(minres, -3, -2, -1, 0, 1, 2, 3, maxres))
}

plotCIgroup <- function(form, data, horiz=FALSE,
    arrowlen = 0, includeOverall=TRUE, distr=c("normal", "t"), conflevel = .95, las=2, ...){
    mf <- model.frame(form, data)
    resp <- mf[,1]
    fac <- mf[,2]
    if(includeOverall){
        lfac <- levels(fac)
        fac <- factor(c(as.character(fac), rep("All Obs", length(resp))), levels=c(lfac, "All Obs"))
        resp <- c(c(resp), c(resp))
    }
    ag <- do.call(rbind, by(resp, list(fac), function(x)confidenceInterval(x, distr=distr,
       confidence=conflevel)))
    ngroup <- nrow(ag)
    if(!horiz){
        yl <- range(c(ag[,2:3]))
        xd <- (ngroup-1)*.25
        xl <- c(1-xd, ngroup + xd)
        plot(xl, yl, axes=F, type="n", xlab="", ...)
        points(1:ngroup, ag[,1], ...)
        axis(1, at=1:ngroup, labels=rownames(ag), las=las)
        axis(2)
        box()
        arrows(1:ngroup, ag[,2], 1:ngroup, ag[,3], code=3, length=arrowlen)
        if(includeOverall){
            abline(v=(max(as.numeric(fac))-.5), lty=2)
        }
    }
    if(horiz){
        yl <- range(c(ag[,2:3]))
        xd <- (ngroup-1)*.25
        xl <- c(1-xd, ngroup + xd)
        plot(yl, xl, axes=F, type="n", ylab="", ...)
        points(ag[,1], 1:ngroup, ...)
        axis(2, at=1:ngroup, labels=rownames(ag), las=las)
        axis(1)
        box()
        arrows(ag[,2], 1:ngroup, ag[,3], 1:ngroup, code=3, length=arrowlen)
        if(includeOverall){
            abline(h=(max(as.numeric(fac))-.5), lty=2)
        }
    }
}
searchVarLabels <- function(dat, str) UseMethod("searchVarLabels")
searchVarLabels.data.frame <-
function (dat, str)
{
    if ("var.labels" %in% names(attributes(dat))) {
        vlat <- "var.labels"
    }
    if ("variable.labels" %in% names(attributes(dat))) {
        vlat <- "variable.labels"
    }
    ind <- sort(union(grep(str, attr(dat, vlat), ignore.case = T), grep(str, names(dat), ignore.case = T)))
    vldf <- data.frame(ind = ind, label = attr(dat, vlat)[ind])
    rownames(vldf) <- names(dat)[ind]
    vldf
}
searchVarLabels.tbl_df <-
function (dat, str)
{
    vlat <- unlist(sapply(1:ncol(dat), function(i)attr(dat[[i]], "label")))
    ind <- sort(union(grep(str, vlat, ignore.case = T), grep(str, names(dat), ignore.case = T)))
    vldf <- data.frame(ind = ind, label = vlat[ind])
    rownames(vldf) <- names(dat)[ind]
    vldf
}

readDTA <- read_dta
freqDist <- function(x){
  tab <- table(x)
  ntab <- names(tab)
  pct <- tab/sum(tab)*100
  cpct <- cumsum(pct)
  tab <- c(tab, sum(tab))
  names(tab) <- c(ntab, "Total")
  pct <- c(pct, 100)
  cpct <- c(cpct, 100)
  maxnum <- max(nchar(tab))
  fr <- sprintf(paste0("%", maxnum, ".0f"), tab)
  pc <- sprintf("%6.2f", pct)
  cp <- sprintf("%6.2f", cpct)
  cp[length(cp)] <- ""
  out <- cbind(fr, pc, cp)
  rownames(out) <- names(tab)
  colnames(out) <- c("Freq", "  %   ", " Cu % ")
  noquote(out)
}

histDiscrete <- function(x, ...){
    l <- max(x, na.rm=T)
    b <- seq(.5, l+.5, by=1)
    hist(x, breaks=b, ...)
}

unalike <- function(x){
  o <- outer(x, x, "!=")
  mean(c(o[lower.tri(o)]), na.rm=T)
}

GKGamma <- function (x, y = NULL, conf.level = NA, ...){
## Function taken from DescTools v0.99.22
    if (!is.null(y))
        tab <- table(x, y, ...)
    else tab <- as.table(x)
    x <- ConDisPairs(tab)
    psi <- 2 * (x$D * x$pi.c - x$C * x$pi.d)/(x$C + x$D)^2
    sigma2 <- sum(tab * psi^2) - sum(tab * psi)^2
    gamma <- (x$C - x$D)/(x$C + x$D)
    if (is.na(conf.level)) {
        result <- gamma
    }
    else {
        pr2 <- 1 - (1 - conf.level)/2
        ci <- qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + gamma
        result <- c(gamma = gamma, lwr.ci = max(ci[1], -1), ups.ci = min(ci[2],
            1))
    }
    class(result) <- "gkg"
    return(result)
}

confidenceInterval <- function (x, confidence = 0.95,  na.rm = TRUE, distr=c("normal", "t")){
    distr <- match.arg(distr)
    nobs <- sum(!is.na(x))
    est <- mean(x, na.rm = na.rm)
    stderr <- sd(x, na.rm = na.rm)/sqrt(nobs)
    alpha <- 1-confidence
    if(distr == "t"){
      ci.low <- est + qt(alpha/2, nobs - 1) * stderr
      ci.high <- est - qt(alpha/2, nobs - 1) * stderr
    }
    else{
      ci.low <- est + qnorm(alpha/2) * stderr
      ci.high <- est - qnorm(alpha/2) * stderr
    }
    retval <- c(Estimate = est, `CI lower` = ci.low, `CI upper` = ci.high,
        `Std. Error` = stderr)
    retval
}




print.gkg <- function(x){
  if(class(x) != "gkg")stop("Object must be of class gkg\n")
  if(length(x) == 1){
    cat("Goodman-Kruskal's Gamma = ", round(x,3), "\n", sep="")
  }
  if(length(x) == 3){
    cat("Goodman-Kruskal's Gamma = ", round(x,3), ", 95% CI (", round(x[2], 3), ", ", round(x[3],3),  ")\n", sep="")
  }
}
print.ktb <- function(x){
  if(class(x) != "ktb")stop("Object must be of class ktb\n")
  cat("Kendall's Tau-b = ", round(x,3), "\n", sep="")
}


freqDist.iscss <- function () {
  defaults <- list (initial.x = NULL, initial.goodnessOfFit = "0")
  dialog.values <- getDialog ("freqDist.iscss", defaults)
  initializeDialog(title = gettextRcmdr("Frequency Distributions"))
  xBox <- variableListBox(top, selectmode = "single",
                          title = gettextRcmdr("Variables (pick one)"),
                          initialSelection = NULL)
  optionsFrame <- tkframe(top)
  goodnessOfFitVariable <- tclVar(dialog.values$initial.goodnessOfFit)
  goodnessOfFitCheckBox <- ttkcheckbutton(optionsFrame, variable = goodnessOfFitVariable)
  onOK <- function() {
    x <- getSelection(xBox)
    if (length(x) == 0) {
      errorCondition(recall = freqDist.iscss, message = gettextRcmdr("You must select a variable."))
      return()
    }
    goodnessOfFit <- tclvalue(goodnessOfFitVariable)
    putDialog ("freqDist.iscss", list (initial.x = x, initial.goodnessOfFit = goodnessOfFit))
    if (length(x) > 1 && goodnessOfFit == "1") {
      errorCondition(recall = freqDist.iscss, message = gettextRcmdr("Goodness-of-fit test not available when more than one variable is selected."))
      return()
    }
    closeDialog()
    .activeDataSet <- ActiveDataSet()
      command <- paste("with(", .activeDataSet, ", freqDist(", x, "))", sep = "")
      if (goodnessOfFit != 1) {
        doItAndPrint(command)
      }
    env <- environment()
    subwin <- NULL
    if (goodnessOfFit == 1) {
      initializeDialog(subwin, title = gettextRcmdr("Goodness-of-Fit Test"))
      hypothesisFrame <- tkframe(subwin)
      levs <- eval(parse(text = paste("levels(", .activeDataSet,
                                      "$", x, ")", sep = "")))
      n.levs <- length(levs)
      assign(".entry.1", tclVar(paste("1/", n.levs, sep = "")),
             envir = env)
      make.entries <- "labelRcmdr(hypothesisFrame, text='Hypothesized probabilities:   ')"
      make.lev.names <- "labelRcmdr(hypothesisFrame, text='Factor levels:')"
      for (i in 1:n.levs) {
        entry.varname <- paste(".entry.", i, sep = "")
        assign(entry.varname, tclVar(paste("1/", n.levs,
                                           sep = "")), envir = env)
        make.entries <- paste(make.entries, ", ", "ttkentry(hypothesisFrame, width='5', textvariable=",
                              entry.varname, ")", sep = "")
        make.lev.names <- paste(make.lev.names, ", labelRcmdr(hypothesisFrame, text='",
                                levs[i], "')", sep = "")
      }
      eval(parse(text = paste("tkgrid(", make.lev.names,
                              ", sticky='w')", sep = "")), envir = env)
      eval(parse(text = paste("tkgrid(", make.entries,
                              ", stick='w')", sep = "")), envir = env)
      tkgrid(hypothesisFrame, sticky = "w")
      onOKsub <- function() {
        probs <- rep(NA, n.levs)
        for (i in 1:n.levs) {
          entry.varname <- paste(".entry.", i, sep = "")
          res <- try(entry <- eval(parse(text = eval(parse(text = paste("tclvalue(",
                                                                        entry.varname, ")", sep = "")), envir = env))),
                     silent = TRUE)
          if (class(res) == "try-error") {
            errorCondition(subwin, message = gettextRcmdr("Invalid entry."))
            return()
          }
          if (length(entry) == 0) {
            errorCondition(subwin, message = gettextRcmdr("Missing entry."))
            return()
          }
          opts <- options(warn = -1)
          probs[i] <- as.numeric(entry)
          options(opts)
        }
        probs <- na.omit(probs)
        if (length(probs) != n.levs) {
          errorCondition(subwin, message = sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number levels (%d)."),
                                                   length(probs), n.levs))
          return()
        }
        if (any(probs < 0)) {
          errorCondition(subwin, message = gettextRcmdr("Negative probabilities not allowed."))
          return()
        }
        if (abs(sum(probs) - 1) > 0.001) {
          Message(message = gettextRcmdr("Probabilities rescaled to sum to 1."),
                  type = "warning")
          probs <- probs/sum(probs)
        }
        closeDialog(subwin)
        command <- paste(command, "\n  .Probs <- c(", paste(probs, collapse = ","), ")", sep = "")
        command <- paste(command, "\n  chisq.test(.Table, p=.Probs)\n})")
        doItAndPrint(command)
      }
      subOKCancelHelp(subwin)
      tkgrid(subButtonsFrame, sticky = "w")
      dialogSuffix(subwin, onOK = onOKsub, focus = subwin, force.wait=TRUE)
    }
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject = "table", reset = "freqDist.iscss", apply="freqDist.iscss")
  tkgrid(getFrame(xBox), sticky = "nw")
  tkgrid(goodnessOfFitCheckBox,
         labelRcmdr(optionsFrame, text = gettextRcmdr("Chi-square goodness-of-fit test (for one variable only)")),
         sticky = "w")
  tkgrid(optionsFrame, sticky = "w")
  tkgrid(buttonsFrame, sticky = "w")
  dialogSuffix()
}

histDiscrete.iscss <- function () {
    defaults <- list (initial.variable = NULL, initial.xlab=gettextRcmdr("<auto>"),
                      initial.ylab=gettextRcmdr("<auto>"), initial.main=gettextRcmdr("<auto>"),
                      initial.labelorient="horizontal", initial.scale="frequency", initial.colors="default", initial.tab=0)
    dialog.values <- getDialog ("histDiscrete.iscss", defaults)
    initializeDialog(title = gettextRcmdr("Discrete Histogram"), use.tabs=TRUE)
    optionsFrame <- tkframe(optionsTab)
    optionsFrame2 <- tkframe(optionsTab)
    variablesFrame <- tkframe(dataTab)
    variableBox <- variableListBox(variablesFrame, selectmode="single", title = gettextRcmdr("Variable (pick one)"),
                                   initialSelection = NULL)
    parFrame <- ttklabelframe(optionsFrame, labelwidget=tklabel(optionsFrame, text = gettextRcmdr("Plot Labels"),
                                                                font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    xlabVar <- tclVar(dialog.values$initial.xlab)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    mainVar <- tclVar(dialog.values$initial.main)
    xlabEntry <- ttkentry(parFrame, width = "25", textvariable = xlabVar)
    xlabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(xlabEntry, ...))
    tkconfigure(xlabEntry, xscrollcommand = function(...) tkset(xlabScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("x-axis label")), xlabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), xlabScroll, sticky = "ew", padx=6)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)
    mainEntry <- ttkentry(parFrame, width = "25", textvariable = mainVar)
    mainScroll <- ttkscrollbar(parFrame, orient = "horizontal",
                               command = function(...) tkxview(mainEntry, ...))
    tkconfigure(mainEntry, xscrollcommand = function(...) tkset(mainScroll,
                                                                ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("Graph title")), mainEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text=""), mainScroll, sticky = "ew", padx=6)
    onOK <- function() {
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        variable <- getSelection(variableBox)
        scale <- tclvalue(scaleVariable)
        xlab <- trim.blanks(tclvalue(xlabVar))
        xlab <- if (xlab == gettextRcmdr("<auto>"))
            paste(", xlab=\"", variable, "\"", sep = "")
        else paste(", xlab=\"", xlab, "\"", sep = "")
        ylab <- trim.blanks(tclvalue(ylabVar))
        ylab <- if (ylab == gettextRcmdr("<auto>")){
            if (scale == "frequency")
                paste(", ylab=\"Frequency\"", sep = "")
            else paste(", ylab=\"Percent\"", sep = "")
        }
        else paste(", ylab=\"", ylab, "\"", sep = "")
        main <- trim.blanks(tclvalue(mainVar))
        main <- if (main == gettextRcmdr("<auto>"))
            ", main = ''"
        else paste(", main=\"", main, "\"", sep = "")
        colors <- tclvalue(colorsVariable)
        putDialog ("histDiscrete.iscss", list(initial.variable = variable, initial.xlab=tclvalue(xlabVar),
                                    initial.ylab=tclvalue(ylabVar), initial.main=tclvalue(mainVar),
                                    initial.scale=scale, initial.colors=colors, initial.tab=tab))
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = histDiscrete.iscss, message = gettextRcmdr("You must select a variable"))
            return()
        }
        scale <- if (scale == "frequency") ", freq = TRUE" else ', freq=FALSE'
        col <- if (colors == "default") "" else paste0(", col=", colors)
        command <- paste0("with(", ActiveDataSet(),", histDiscrete(", variable, xlab, ylab, main, col, scale, "))")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    radioButtons(optionsFrame2, name = "scale", buttons = c("frequency", "proportions"),
                 labels = gettextRcmdr(c("Frequency counts", "Proportions")),
                 title = gettextRcmdr("Axis Scaling"),
                 initialValue = dialog.values$initial.scale)
    radioButtons(optionsFrame2, name = "colors", buttons = c("default", "palette"),
                 labels = gettextRcmdr(c("Default", "From color palette")),
                 title = gettextRcmdr("Color Selection"),
                 initialValue = dialog.values$initial.colors)
    OKCancelHelp(helpSubject = "histDiscrete", reset = "histDiscrete.iscss", apply = "histDiscrete.iscss")
    tkgrid(getFrame(variableBox), sticky="w")
    tkgrid(tklabel(variablesFrame, text=""))
    tkgrid(scaleFrame, sticky="w")
    tkgrid(colorsFrame, sticky="w")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(parFrame, sticky = "nw")
    tkgrid(optionsFrame2, optionsFrame, sticky = "nw")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE)
}

numSumAll.iscss <- function(){
    Library("abind")
    Library("e1071")
    defaults <- list(initial.x=NULL, initial.mean="1", initial.sd="1", initial.se.mean="0", initial.IQR="1", initial.cv="0",
                     initial.quantiles.variable="1",
                     initial.quantiles="0, .25, .5, .75, 1",
                     initial.skewness="0", initial.kurtosis="0", initial.type="2",
                     initial.counts="0",
                     initial.group=NULL, initial.tab=0)
    dialog.values <- getDialog("numSumAll.iscss", defaults)
    initial.group <- dialog.values$initial.group
    initializeDialog(title=gettextRcmdr("Numerical Summaries"), use.tabs=TRUE, tabs=c("dataTab", "statisticsTab"))
    xBox <- variableListBox(dataTab,  selectmode="multiple", title=gettextRcmdr("Variables (pick one or more)"),
                            initialSelection=NULL)
    checkBoxes(window = statisticsTab, frame="checkBoxFrame", boxes=c("mean", "sd", "se.mean", "IQR", "cv", "counts"),
               initialValues=c(dialog.values$initial.mean, dialog.values$initial.sd, dialog.values$initial.se.mean,
                               dialog.values$initial.IQR, dialog.values$initial.cv, dialog.values$initial.counts),
               labels=gettextRcmdr(c("Mean", "Standard Deviation", "Standard Error of Mean", "Interquartile Range",
                                     "Coefficient of Variation", "Binned Frequency Counts")))
    skFrame <- tkframe(statisticsTab)
    checkBoxes(window = skFrame, frame="skCheckBoxFrame", boxes=c("skewness", "kurtosis"),
               initialValues=c(dialog.values$initial.skewness, dialog.values$initial.kurtosis),
               labels=gettextRcmdr(c("Skewness", "Kurtosis")))
    radioButtons(window = skFrame, name="typeButtons", buttons=c("b1", "b2", "b3"), values=c("1", "2", "3"),
                 initialValue=dialog.values$initial.type,
                 labels=gettextRcmdr(c("Type 1", "Type 2", "Type 3")))
    quantilesVariable <- tclVar(dialog.values$initial.quantiles.variable)
    quantilesFrame <- tkframe(statisticsTab)
    quantilesCheckBox <- tkcheckbutton(quantilesFrame, variable=quantilesVariable,
                                       text=gettextRcmdr("Quantiles:"))
    quantiles <- tclVar(dialog.values$initial.quantiles)
    quantilesEntry <- ttkentry(quantilesFrame, width="20", textvariable=quantiles)
    groupsBox(recall=numSumAll.iscss, label=gettextRcmdr("Summarize by:"),
              initialLabel=if (is.null(initial.group)) gettextRcmdr("Summarize by groups")
              else paste(gettextRcmdr("Summarize by:"), initial.group),
              initialGroup=initial.group, window = dataTab)
    onOK <- function(){
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        x <- getSelection(xBox)
        quants <- tclvalue(quantiles)
        meanVar <- tclvalue(meanVariable)
        sdVar <- tclvalue(sdVariable)
        se.meanVar <- tclvalue(se.meanVariable)
        IQRVar <- tclvalue(IQRVariable)
        cvVar <- tclvalue(cvVariable)
        countsVar <- tclvalue(countsVariable)
        quantsVar <- tclvalue(quantilesVariable)
        skewnessVar <- tclvalue(skewnessVariable)
        kurtosisVar <- tclvalue(kurtosisVariable)
        typeVar <- tclvalue(typeButtonsVariable)
        putDialog("numSumAll.iscss", list(
            initial.x=x, initial.mean=meanVar, initial.sd=sdVar, initial.se.mean=se.meanVar, initial.IQR=IQRVar,
            initial.cv=cvVar, initial.counts=countsVar,
            initial.quantiles.variable=quantsVar, initial.quantiles=quants,
            initial.skewness=skewnessVar, initial.kurtosis=kurtosisVar, initial.type=typeVar,
            initial.group=if (.groups != FALSE) .groups else NULL, initial.tab=tab
        ))
        if (length(x) == 0){
            errorCondition(recall=numSumAll.iscss, message=gettextRcmdr("You must select a variable."))
            return()
        }
        closeDialog()
        quants <- paste("c(", gsub(",+", ",", gsub(" ", ",", quants)), ")", sep="")
        .activeDataSet <- ActiveDataSet()
        vars <- if (length(x) == 1) paste('"', x, '"', sep="")
        else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
        ds.vars <- paste("sapply(", vars, ", function(i)as.numeric(", .activeDataSet, "[[i]]))", sep="")
        stats <- paste("c(",
                       paste(c('"mean"', '"sd"', '"se(mean)"', '"IQR"', '"quantiles"', '"cv"', '"skewness"', '"kurtosis"')
                             [c(meanVar, sdVar, se.meanVar, IQRVar, quantsVar, cvVar, skewnessVar, kurtosisVar) == 1],
                             collapse=", "), ")", sep="")
        if (stats == "c()" && countsVar != 1){
            errorCondition(recall=numSumAll.iscss, message=gettextRcmdr("No statistics selected."))
            return()
        }
        type.text <- if (skewnessVar == 1 || kurtosisVar == 1) paste(', type="', typeVar, '"', sep="") else ""
        if (.groups != FALSE) grps <- paste(.activeDataSet, "$", .groups, sep="")
        if (stats != "c()"){
            command <- if (.groups != FALSE) {
                paste("numSummary(", ds.vars, ", groups=", grps, ", statistics=", stats,
                      ", quantiles=", quants, type.text, ")", sep="")
            }
            else  paste("numSummary(", ds.vars, ", statistics=", stats,
                        ", quantiles=", quants, type.text, ")", sep="")
            doItAndPrint(command)
        }
        if (countsVar == 1){
            if (.groups != FALSE){
                levels <- eval(parse(text=paste0("levels(", grps, ")")), envir=.GlobalEnv)
                for (level in levels){
                    command <- paste0("binnedCounts(", .activeDataSet, "[", grps, " == ", "'", level, "', ",
                                      vars, ", drop=FALSE])\n  # ", .groups, " = ", level)
                    doItAndPrint(command)
                }
            }
            else {
                command <- paste0("binnedCounts(", ds.vars, ")")
                doItAndPrint(command)
            }
        }
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="numSummary", reset="numSumAll.iscss", apply ="numSumAll.iscss")
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(checkBoxFrame, sticky="w")
    tkgrid(skCheckBoxFrame, typeButtonsFrame, sticky="nw")
    tkgrid(skFrame, sticky="w")
    tkgrid(quantilesCheckBox, quantilesEntry, sticky="w")
    tkgrid(quantilesFrame)
    tkgrid(groupsFrame, sticky = "w", padx=6)
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tabs=c("dataTab", "statisticsTab"),
                 tab.names=c("Data", "Statistics"))
}

subsetDataSet.iscss <- function(){
    dataSet <- activeDataSet()
    initializeDialog(title=gettextRcmdr("Subset Data Set"))
    allVariablesFrame <- tkframe(top)
    allVariables <- tclVar("1")
    allVariablesCheckBox <- ttkcheckbutton(allVariablesFrame, variable=allVariables)
    variablesBox <- variableListBox(top, Variables(), selectmode="multiple",
                                    initialSelection=NULL, title=gettextRcmdr("Variables (select one or more)"))
    subsetVariable <- tclVar(gettextRcmdr("<all cases>"))
    subsetFrame <- tkframe(top)
    subsetEntry <- ttkentry(subsetFrame, width="20", textvariable=subsetVariable)
    subsetScroll <- ttkscrollbar(subsetFrame, orient="horizontal",
                                 command=function(...) tkxview(subsetEntry, ...))
    tkconfigure(subsetEntry, xscrollcommand=function(...) tkset(subsetScroll, ...))
    newDataSetName <- tclVar(gettextRcmdr("<same as active data set>"))
    dataSetNameFrame <- tkframe(top)
    dataSetNameEntry <- ttkentry(dataSetNameFrame, width="25", textvariable=newDataSetName)
    onOK <- function(){
        newName <- trim.blanks(tclvalue(newDataSetName))
        if (newName == gettextRcmdr("<same as active data set>")) newName <- ActiveDataSet()
        if (!is.valid.name(newName)){
            errorCondition(recall=subsetDataSet.iscss,
                           message=paste('"', newName, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
        }
        if (is.element(newName, listDataSets())) {
            if ("no" == tclvalue(checkReplace(newName, type=gettextRcmdr("Data set")))){
                closeDialog()
                subsetDataSet.iscss()
                return()
            }
        }
        selectVars <- if (tclvalue(allVariables) == "1") ""
        else {
            x <- getSelection(variablesBox)
            if (0 == length(x)) {
                errorCondition(recall=subsetDataSet.iscss,
                               message=gettextRcmdr("No variables were selected."))
                return()
            }
            paste(", select=c(", paste(x, collapse=","), ")", sep="")
        }
        closeDialog()
        cases <- tclvalue(subsetVariable)
        selectCases <- if (cases == gettextRcmdr("<all cases>")) ""
        else if (length(grep(gettextRcmdr("sample"), cases))>0){
            gpct <- grep("%", cases, fixed=T)
            if(length(gpct) > 0){
                nsamp <- floor(eval(parse(text=paste0("nrow(", ActiveDataSet(), ")")))* (as.numeric(gsub(".*\\s(\\d+)%$", "\\1", cases))/100))
            }
            else nsamp <- as.integer(gsub(".*\\s(\\d+)$", "\\1", cases))
            insamp <- sample(1:eval(parse(text=paste0("nrow(", ActiveDataSet(), ")"))), nsamp, replace=FALSE)
            paste0(", subset=1:nrow(", ActiveDataSet(), ") %in% ", paste0("c(", paste(insamp, collapse=", "), ")"))
          }
          else paste(", subset=", cases, sep="")
        if (selectVars == "" && selectCases ==""){
            errorCondition(recall=subsetDataSet.iscss,
                           message=gettextRcmdr("New data set same as active data set."))
            return()
        }
        command <- paste(newName, " <- subset(", ActiveDataSet(), selectCases, selectVars, ")",
                         sep="")
        logger(command)
        result <- justDoIt(command)
        if (class(result)[1] !=  "try-error") activeDataSet(newName)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="subset")
    tkgrid(allVariablesCheckBox, labelRcmdr(allVariablesFrame, text=gettextRcmdr("Include all variables")),
           sticky="w")
    tkgrid(allVariablesFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("   OR"), fg="red"), sticky="w")
    tkgrid(getFrame(variablesBox), sticky="nw")
    tkgrid(labelRcmdr(subsetFrame, text=gettextRcmdr("Subset expression")), sticky="w")
    tkgrid(subsetEntry, sticky="w")
    tkgrid(subsetScroll, sticky="ew")
    tkgrid(subsetFrame, sticky="w")
    tkgrid(labelRcmdr(dataSetNameFrame, text=gettextRcmdr("Name for new data set")), sticky="w")
    tkgrid(dataSetNameEntry, sticky="w")
    tkgrid(dataSetNameFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}


plotCIgroup.iscss <- function(){
    defaults <- list(initial.row=NULL, initial.column=NULL,initial.conflevel=95,
      initial.ylab=gettextRcmdr("<auto>"), initial.distr="normal", initial.horizontal=0, initial.includeOverall = 1, initial.arrowlen=0, initial.las=0)
    dialog.values <- getDialog("plotCIgroup.iscss", defaults)
    initializeDialog(title=gettextRcmdr("Plot CIs by Group"), use.tabs=FALSE)
    variablesFrame <- tkframe(top)
    optionsFrame <- tkframe(top)
    optionsFrame2 <- tkframe(top)
    parFrame <- ttklabelframe(optionsFrame2, labelwidget=tklabel(optionsFrame2,
      text = gettextRcmdr("Plot Labels"),font="RcmdrTitleFont", foreground=getRcmdr("title.color")))
    rowBox <- variableListBox(variablesFrame, selectmode="single",
      title=gettextRcmdr("Quantitative Variable (pick one)"), initialSelection=NULL)
    columnBox <- variableListBox(variablesFrame, selectmode="single",
      title=gettextRcmdr("Grouping Variable (pick one)"), initialSelection=NULL)
    conflevelVariable <- tclVar(gettextRcmdr("95"))
    conflevelFrame <- tkframe(top)
    conflevelEntry <- ttkentry(conflevelFrame, width="20", textvariable=conflevelVariable)
    arrowlenVariable <- tclVar(gettextRcmdr("0"))
    arrowlenFrame <- tkframe(top)
    arrowlenEntry <- ttkentry(arrowlenFrame, width="20", textvariable=arrowlenVariable)
    ylabVar <- tclVar(dialog.values$initial.ylab)
    ylabEntry <- ttkentry(parFrame, width = "25", textvariable = ylabVar)
    ylabScroll <- ttkscrollbar(parFrame, orient = "horizontal",
      command = function(...) tkxview(ylabEntry, ...))
    tkconfigure(ylabEntry, xscrollcommand = function(...) tkset(ylabScroll, ...))
    tkgrid(labelRcmdr(parFrame, text = gettextRcmdr("y-axis label")), ylabEntry, sticky = "ew", padx=6)
    tkgrid(labelRcmdr(parFrame, text =""), ylabScroll, sticky = "ew", padx=6)

    onOK <- function() {
      row <- getSelection(rowBox)
      column <- getSelection(columnBox)
      conflevel <- tclvalue(conflevelVariable)
      arrowlen <- tclvalue(arrowlenVariable)
      distr <- tclvalue(distrVariable)
      las <- tclvalue(lasVariable)
      las2 <- as.numeric(las)
      horizontal <- tclvalue(horizontalVariable)
      horizontal2 <- ifelse(horizontal == "0", FALSE, TRUE)
      includeOverall <- tclvalue(includeOverallVariable)
      includeOverall2 <- ifelse(includeOverall == "0", FALSE, TRUE)
      ylab <- trim.blanks(tclvalue(ylabVar))
      ylab2 <- if (ylab == gettextRcmdr("<auto>"))row
        else ylab
      putDialog("plotCIgroup.iscss", list(initial.row = row, intial.column = column,
        initial.conflevel=conflevel, initial.ylab=ylab2, initial.distr=distr, initial.horizontal=horizontal, initial.includeOverall=includeOverall, initial.arrowlen=arrowlen, initial.las=las))
      closeDialog()
        if (length(row) == 0 | length(col) == 0) {
            errorCondition(recall = plotCIgroup.iscss, message = gettextRcmdr("You must select a quantitative variable and a grouping variable"))
            return()
        }
        command <- paste("plotCIgroup(", row, "~", column, ", ", ActiveDataSet(), ", conflevel = ", as.numeric(conflevel)/100, ", horiz = ", horizontal2, ", includeOverall=", includeOverall2, ", arrowlen = ", as.numeric(arrowlen), ", distr = '", distr, "', ylab = '", ylab2, "', las = ", las2, ")", sep="")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())

  }
  OKCancelHelp(helpSubject = "plotCIgroup", reset = "plotCIgroup.iscss", apply = "plotCIgroup.iscss")
 	rightFrame <- tkframe(top)
	radioButtons(top, name = "distr", buttons = c("normal", "t"), values = c("normal", "t"),
	             labels = gettextRcmdr(c("Normal", "Student's T")),
               title = gettextRcmdr("Distribution"),
	             initialValue = dialog.values$initial.distr)
  checkBoxes(top, frame="optionsFrame", boxes=c("horizontal", "includeOverall"),
          initialValues=c(dialog.values$initial.horizontal, dialog.values$initial.includeOverall),
          labels=gettextRcmdr(c("Horizontal Bars", "Include Overall Mean CI")))
	radioButtons(top, name = "las", buttons = c("parallel", "horizontal", "perpendicular", "vertical"), values = c(0,1,2,3),
	             labels = gettextRcmdr(c("Always Parallel", "Always Horizontal", "Always Perpendicular", "Always Vertical")),
               title = gettextRcmdr("Label Orientation"),
	             initialValue = dialog.values$initial.las)
  tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(optionsFrame, sticky="w")
  tkgrid(labelRcmdr(conflevelFrame, text=gettextRcmdr("Confidence Level")), sticky="w")
  tkgrid(conflevelEntry, sticky="w")
  tkgrid(conflevelFrame, sticky="w")
  tkgrid(labelRcmdr(arrowlenFrame, text=gettextRcmdr("Arrow Length")), sticky="w")
  tkgrid(arrowlenEntry, sticky="w")
  tkgrid(arrowlenFrame, sticky="w")
  tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
  tkgrid(distrFrame, rightFrame, sticky = "w")
  tkgrid(lasFrame, rightFrame, sticky = "w")
  tkgrid(optionsFrame2, sticky="w")
  tkgrid(parFrame, sticky = "nw")
  tkgrid(buttonsFrame, sticky="w")
  dialogSuffix()
}


ci.iscss <- function () {
    defaults <- list (initial.variable = NULL, initial.conflevel=95, initial.distr="normal")
    dialog.values <- getDialog ("ci.iscss", defaults)
    initializeDialog(title = gettextRcmdr("Confidence Interval"), use.tabs=FALSE)
    variableBox <- variableListBox(top, selectmode="single", title =
      gettextRcmdr("Variable (pick one)"), initialSelection = NULL)
    conflevelVariable <- tclVar(gettextRcmdr("95"))
    conflevelFrame <- tkframe(top)
    conflevelEntry <- ttkentry(conflevelFrame, width="20", textvariable=conflevelVariable)
    onOK <- function() {
        variable <- getSelection(variableBox)
        conflevel <- tclvalue(conflevelVariable)
        distr <- tclvalue(distrVariable)
      putDialog ("ci.iscss", list(initial.variable = variable, initial.conflevel=conflevel, initial.distr = distr))
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = ci.iscss, message = gettextRcmdr("You must select a variable"))
            return()
        }
        command <- paste0("with(", ActiveDataSet(),", confidenceInterval(", variable, ", confidence= ", as.numeric(conflevel)/100,", distr = '", distr, "'))")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "confidenceInterval", reset = "ci.iscss", apply = "ci.iscss")
   	rightFrame <- tkframe(top)
  	radioButtons(top, name = "distr", buttons = c("normal", "t"), values = c("normal", "t"),
  	             labels = gettextRcmdr(c("Normal", "Student's T")),
                 title = gettextRcmdr("Distribution"),
  	             initialValue = dialog.values$initial.distr)
    tkgrid(getFrame(variableBox), sticky="nw")
    tkgrid(labelRcmdr(conflevelFrame, text=gettextRcmdr("Confidence Level")), sticky="w")
    tkgrid(conflevelEntry, sticky="w")
    tkgrid(conflevelFrame, sticky="w")
	  tkgrid(labelRcmdr(rightFrame, text = ""), sticky = "w")
	  tkgrid(distrFrame, rightFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}

unalike.iscss <- function () {
    defaults <- list (initial.variable = NULL)
    dialog.values <- getDialog ("unalike.iscss", defaults)
    initializeDialog(title = gettextRcmdr("Unalikability Measure"), use.tabs=FALSE)
    variableBox <- variableListBox(top, selectmode="single", title =
      gettextRcmdr("Variable (pick one)"), initialSelection = NULL)
    onOK <- function() {
        variable <- getSelection(variableBox)
        putDialog ("unalike.iscss", list(initial.variable = variable))
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = unalike.iscss, message = gettextRcmdr("You must select a variable"))
            return()
        }
        command <- paste0("with(", ActiveDataSet(),", unalike(", variable, "))")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "unalike", reset = "unalike.iscss", apply = "unalike.iscss")
    tkgrid(getFrame(variableBox), sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}

inspect.iscss <- function () {
    defaults <- list (initial.variable = NULL)
    dialog.values <- getDialog ("inspect.iscss", defaults)
    initializeDialog(title = gettextRcmdr("Inspect Variable"), use.tabs=FALSE)
    variableBox <- variableListBox(top, selectmode="single", title =
      gettextRcmdr("Variable (pick one)"), initialSelection = NULL)
    onOK <- function() {
        variable <- getSelection(variableBox)
        putDialog ("inspect.iscss", list(initial.variable = variable))
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = inspect.iscss, message = gettextRcmdr("You must select a variable"))
            return()
        }
        command <- paste0("inspect(", ActiveDataSet(),  ", '", variable, "')")
        doItAndPrint(command)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "inspect", reset = "inspect.iscss", apply = "inspect.iscss")
    tkgrid(getFrame(variableBox), sticky="nw")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}

asFactor.iscss <- function () {
    defaults <- list (initial.variable = NULL, initial.name = "variable")
    dialog.values <- getDialog ("asFactor.iscss", defaults)
    initializeDialog(title = gettextRcmdr("Change Tibble Variable to Factor"), use.tabs=FALSE)
    variableBox <- variableListBox(top, selectmode="single", title =
      gettextRcmdr("Variable (pick one)"), initialSelection = NULL)
    variablesFrame <- tkframe(top)
    newVariableName <- tclVar(dialog.values$initial.name)
    newVariable <- ttkentry(variablesFrame, width = "20", textvariable = newVariableName)

    onOK <- function() {
        variable <- getSelection(variableBox)
        dataSet <- ActiveDataSet()
        closeDialog()
        if (length(variable) == 0) {
            errorCondition(recall = asFactor.iscss, message = gettextRcmdr("You must select a variable"))
            return()
        }
        name <- trim.blanks(tclvalue(newVariableName))
        putDialog ("asFactor.iscss", list(initial.variable = variable, initial.name=name))
        command <- paste0(dataSet, "[['", name, "']] <- haven:::as_factor.labelled(", dataSet, "[['", variable, "']])")
        result <- doItAndPrint(command)
        if (class(result)[1] != "try-error") activeDataSet(ActiveDataSet(), flushModel = FALSE,
            flushDialogMemory = FALSE)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "as_factor", reset = "asFactor.iscss", apply = "asFactor.iscss")
    tkgrid(getFrame(variableBox), sticky="nw")
    tkgrid(labelRcmdr(variablesFrame, text = ""))
    tkgrid(labelRcmdr(variablesFrame, text = gettextRcmdr("New variable name: ")), newVariable, sticky = "w")
    tkgrid(variablesFrame, sticky = "w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}



importSTATA.iscss <- function() {
    initializeDialog(title=gettextRcmdr("Import STATA Data Set"))
    dsname <- tclVar("Dataset")
    dsnameFrame <- tkframe(top)
    entryDsname <- ttkentry(dsnameFrame, width="20", textvariable=dsname)
    onOK <- function(){
        closeDialog()
        setBusyCursor()
        on.exit(setIdleCursor())
        dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == ""){
            errorCondition(recall=importSTATA,
                           message=gettextRcmdr("You must enter the name of a data set."))
            return()
        }
        if (!is.valid.name(dsnameValue)){
            errorCondition(recall=importSTATA,
                           message=paste('"', dsnameValue, '" ', gettextRcmdr("is not a valid name."), sep=""))
            return()
        }
        if (is.element(dsnameValue, listDataSets())) {
            if ("no" == tclvalue(checkReplace(dsnameValue, gettextRcmdr("Data set")))){
                importSTATA.iscss()
                return()
            }
        }
        file <- tclvalue(tkgetOpenFile(
            filetypes=gettextRcmdr('{"All Files" {"*"}} {"STATA datasets" {".dta" ".DTA"}}')))
        if (file == "") {
            tkfocus(CommanderWindow())
            return()
        }
        command <- paste('readDTA("', file,'")', sep="")
        logger(paste(dsnameValue, " <- ", command, sep=""))
        result <- justDoIt(command)
        if (class(result)[1] !=  "try-error"){
            gassign(dsnameValue, result)
            activeDataSet(dsnameValue)
        }
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="read_dta")
    tkgrid(labelRcmdr(dsnameFrame, text=gettextRcmdr("Enter name for data set:  ")), entryDsname, sticky="w")
    tkgrid(dsnameFrame, columnspan=2, sticky="w")
    tkgrid(buttonsFrame, columnspan="2", sticky="ew")
    dialogSuffix(focus=entryDsname)
}



searchVarLabels.iscss <- function(){
    dataSet <- activeDataSet()
    initializeDialog(title=gettextRcmdr("Search Variable Labels"))
    searchVariable <- tclVar(gettextRcmdr("<search string>"))
    searchFrame <- tkframe(top)
    searchEntry <- ttkentry(searchFrame, width="20", textvariable=searchVariable)
    searchScroll <- ttkscrollbar(searchFrame, orient="horizontal",
                                 command=function(...) tkxview(searchEntry, ...))
    tkconfigure(searchEntry, xscrollcommand=function(...) tkset(searchScroll, ...))
    onOK <- function(){
        closeDialog()
        searchstr1 <- tclvalue(searchVariable)
        searchstr <- if (searchstr1 == gettextRcmdr("<search string>")) ""
          else searchstr1
        command <- paste("searchVarLabels(", ActiveDataSet(), ", '", searchstr, "')",
                         sep="")
#        logger(command)
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="searchVarLabels")
    tkgrid(labelRcmdr(searchFrame, text=gettextRcmdr("Search String")), sticky="w")
    tkgrid(searchEntry, sticky="w")
    tkgrid(searchScroll, sticky="ew")
    tkgrid(searchFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
}




twoWayTable.iscss <- function(){
    Library("abind")
    defaults <- list(initial.row=NULL, initial.column=NULL,
        initial.chisq=1, initial.chisqComp=0, initial.expected=0,
        initial.tab=0, initial.colpct = 1, initial.phi=0,
        initial.cramersV = 0, initial.lambda = 0, initial.gamma = 0, initial.d = 0,
        initial.taub = 0, initial.rho = 0, initial.plotStdRes = 0)
    dialog.values <- getDialog("twoWayTable.iscss", defaults)
    initializeDialog(title=gettextRcmdr("Two-Way Table"), use.tabs=TRUE)
    variablesFrame <- tkframe(dataTab)
    .factors <- Factors()
    rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
        initialSelection=varPosn(dialog.values$initial.row, "factor"))
    columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
        initialSelection=varPosn(dialog.values$initial.column, "factor"))
    onOK <- function(){
        tab <- if (as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        row <- getSelection(rowBox)
        column <- getSelection(columnBox)
        chisqComp <- tclvalue(chisqCompVariable)
        expected <- tclvalue(expFreqVariable)
        colpct <- tclvalue(colpctVariable)
        chisq <- tclvalue(chisqTestVariable)
        phi <- tclvalue(phiTestVariable)
        cramersV <- tclvalue(cramersVTestVariable)
        lambda <- tclvalue(lambdaTestVariable)
        gamma <- tclvalue(gammaTestVariable)
        d <- tclvalue(dTestVariable)
        taub <- tclvalue(taubTestVariable)
        rho <- tclvalue(rhoTestVariable)
        plotStdRes <- tclvalue(plotStdResVariable)
        putDialog("twoWayTable.iscss", list(
            initial.row=row,
            initial.column=column, initial.colpct=colpct,
            initial.chisq=chisq, initial.chisqComp=chisqComp,
            initial.expected=expected, initial.tab=tab,
            initial.phi = phi, initial.cramersV = cramersV, initial.lambda=lambda,
            initial.gamma = gamma, initial.d = d, initial.taub=taub,
            initial.rho = rho, initial.plotStdRes = plotStdRes))
        if (length(row) == 0 || length(column) == 0){
            errorCondition(recall=twoWayTable.iscss, message=gettextRcmdr("You must select two variables."))
            return()
        }
        if (row == column) {
            errorCondition(recall=twoWayTable.iscss, message=gettextRcmdr("Row and column variables are the same."))
            return()
        }
        closeDialog()
        command <- paste0("local({.out <- with(", ActiveDataSet(), ", gmodels::CrossTable(", row, ", ", column, ", prop.r = FALSE, prop.t=FALSE, sresid=TRUE, prop.c=", as.logical(as.numeric(colpct)), ", expected = ", as.logical(as.numeric(expected)), ", prop.chisq = ", as.logical(as.numeric(chisqComp)), ", chisq=F))\n")
        if(any(c(chisq, phi, cramersV, lambda, gamma, d, taub, rho) == 1)){
        command.2 <- paste0(".allStats <- with(", ActiveDataSet(), ", makeStats(",row, ", ", column, ", chisq=", as.logical(as.numeric(chisq)),
        ", phi = ", as.logical(as.numeric(phi)), ", cramersV = ", as.logical(as.numeric(cramersV)), ", lambda = ", as.logical(as.numeric(lambda)), ", gamma = ", as.logical(as.numeric(gamma)), ", d = ", as.logical(as.numeric(d)), ", taub = ", as.logical(as.numeric(taub)), ", rho = ", as.logical(as.numeric(rho)),
         ", 2500))\nprint(.allStats)")
         command <- paste0(command, "\n", command.2)
        }
        if(plotStdRes == 1){
            command <- paste(command, "\nplotStdRes(.out$t)\n")
        }
        command <- paste(command, "\n})", sep="")
        doItAndPrint(command)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="xtabs", reset="twoWayTable.iscss", apply="twoWayTable.iscss")
    checkBoxes(optionsTab, frame="percentsFrame",
        boxes=c("colpct", "expFreq", "chisqComp"),
        initialValue=c(dialog.values$initial.colpct, dialog.values$initial.expected, dialog.values$initial.chisqComp),
        labels=gettextRcmdr(c("Column percentages", "Expected Counts", "Chi-Square Contribution")),
        title=gettextRcmdr("Cell Entries"))
    checkBoxes(optionsTab, frame="testsFrame",
        boxes=c("chisqTest", "phiTest", "cramersVTest", "lambdaTest", "gammaTest", "dTest", "taubTest", "rhoTest"),
        initialValues=c(dialog.values$initial.chisq,  dialog.values$initial.phi,
            dialog.values$initial.cramersV, dialog.values$initial.lambda, dialog.values$initial.gamma, dialog.values$initial.d, dialog.values$initial.taub, dialog.values$initial.rho),
        labels=gettextRcmdr(c("Chi-square test of independence", "Phi", "Cramer's V", "Lambda", "Kruskal's Gamma", "Somer's D", "Kendall's Tau-b", "Spearman's Rho")))
    checkBoxes(optionsTab, frame="sdresFrame", boxes="plotStdRes", initialValues=dialog.values$initial.plotStdRes,
        labels=gettextRcmdr("Plot Standardized Residuals"))
    tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox), sticky="nw")
    tkgrid(variablesFrame, sticky="w")
    tkgrid(percentsFrame, sticky="w")
    tkgrid(labelRcmdr(optionsTab, text=gettextRcmdr("Hypothesis Tests"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(labelRcmdr(optionsTab, text=gettextRcmdr("Plot Residuals"), fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(sdresFrame, sticky="w")
    dialogSuffix(use.tabs=TRUE, grid.buttons=TRUE, tab.names=c("Data", "Statistics"))
}
