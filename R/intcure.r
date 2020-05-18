
intcure <- function(formula = formula(data), cureform = NULL,
	data = sys.parent(), subset, na.action, bt = NULL, gm = NULL,
	basepara = NULL, sigma = c(0, 0, 0),
	optimcfg = list(ndeps = 0.001, maxit = 1000, reltol = 0.00001,
		method = "Nelder-Mead", hessian = FALSE),
	intcfg = list(eps = 0.0001, lower = c(-5, -5), upper = c(5, 5)),
	basedist = c("exponential", "weibull", "piecewise"),
	npiece = 3, piececut = NULL, piececuttype = c("quantile", "even"),
	model = FALSE, y = TRUE, x = FALSE, z = FALSE, funval = FALSE,
	debug = c("intcure", "integration", "optim"))
{
	call <- match.call()
	formal.args <- formals(sys.function(sys.parent()))
	m <- match.call(expand.dots = FALSE)
	temp <- c("", "formula", "data", "subset", "na.action")
	m1 <- m[match(temp, names(m), nomatch = 0)]
	temp <- c("", "cureform", "data", "subset", "na.action")
	m2 <- m[match(temp, names(m), nomatch = 0)]

	if (missing(debug)) debug <- NULL
	else debug <- match.arg(debug, several.ok = TRUE)

	special <- "cluster"
	Terms1 <- if(missing(data)) terms(formula, special)
		else terms(formula, special, data = data)
	m1$formula <- Terms1

	m1[[1]] <- as.name("model.frame")
	if(m1[[2]][[1]] == as.name("Surv")) { 			
		m1[[2]][[2]] <- m1[[2]]
		m1[[2]][[1]] <- as.name("~")
		m1[[2]][[3]] <- 1
	}
	m1 <- eval(m1, sys.parent())
	rawY <- model.extract(m1, "response")
	if(!inherits(rawY, "Surv"))
		stop("Response must be a survival object")
	sorted <- order(rawY[, 1], 1-rawY[,2])				
	Y <- rawY[sorted, ]
	datalen <- nrow(Y)

	offsetx <- attr(Terms1, "offset")
	tt <- length(offsetx)
	offsetx <- if(tt == 0) rep(0, datalen) else if(tt == 1)
		m1[[offsetx]]
	else {
		ff <- m1[[offsetx[1]]]
		for(i in 2:tt)
			ff <- ff + m1[[offsetx[i]]]
		ff
	}
	offsetx <- offsetx[sorted]

	clusterx <- attr(Terms1, "specials")$cluster
	dropx <- NULL
	nclusterx <- 0
	if(length(clusterx)) {
		tempc <- untangle.specials(Terms1, "cluster", 1:10)
		ord <- attr(Terms1, "order")[tempc$terms]
		if(any(ord > 1))
			stop("Cluster can not be used in an interaction")
		clusterx <- strata(m1[, tempc$vars], shortlabel = TRUE)
		dropx <- tempc$terms
		clusterx <- as.numeric(clusterx)[sorted]
		nclusterx <- length(unique(clusterx))
	}
	else {
		clusterx <- 1:datalen
		nclusterx <- datalen
	}

	if(length(dropx))
		X <- model.matrix(Terms1[ - dropx], m1)[, -1, drop = FALSE]
	else X <- model.matrix(Terms1, m1)[, -1, drop = FALSE]
	X <- X[sorted, , drop = FALSE]

	nbt <- ncol(X)
	if (!length(X)) X <- NULL
	if (is.null(nbt)) nbt <- 0
	cofnames <- dimnames(X)[[2]]

	ngm <- 0
	Z <- NULL
	clusterz <- NULL
	nclusterz <- 0
	offsetz <- NULL
	if (!missing(cureform)) {
		Terms2 <- if(missing(data)) terms(cureform, special)
			else terms(cureform, special, data = data)
		m2$formula <- Terms2
		m2$cureform <- NULL
    
		m2[[1]] <- as.name("model.frame")
		curedata <- eval(m2, sys.parent())
		Terms2 <- attr(curedata, "terms")
    
		offsetz <- attr(Terms2, "offset")
		tt <- length(offsetz)
		offsetz <- if(tt == 0) rep(0, datalen) else if(tt == 1)
			curedata[[offsetz]]
		else {
			ff <- curedata[[offsetz[1]]]
			for(i in 2:tt)
				ff <- ff + curedata[[offsetz[i]]]
			ff
		}
		offsetz <- offsetz[sorted]
    
		clusterz <- attr(Terms2, "specials")$cluster
		dropz <- NULL
		if(length(clusterz)) {
			tempc <- untangle.specials(Terms2, "cluster", 1:10)
			ord <- attr(Terms2, "order")[tempc$terms]
			if(any(ord > 1))
				stop("Cluster can not be used in an interaction")
			clusterz <- strata(curedata[, tempc$vars], shortlabel = TRUE)
			dropz <- tempc$terms
			clusterz <- as.numeric(clusterz)[sorted]
		}
		else {
			clusterz <- 1:datalen
			nclusterz <- datalen
		}
    
		if(length(dropz))
			Z <- model.matrix(Terms2[ - dropz], curedata)[, , drop = FALSE]
		else Z <- model.matrix(Terms2, curedata)[, , drop = FALSE]
		Z <- Z[sorted, , drop = FALSE]
		ngm <- ncol(Z)
    
		cofnames <- c(cofnames,dimnames(Z)[[2]])
		nclusterz <- length(unique(clusterz))
	}

	basedist <- match.arg(basedist)
	piececuttype <- match.arg(piececuttype)

	names(sigma) = c("lsigmau", "lsigmav", "ftrho")
	if (any(is.na(sigma[1:2])) && !is.na(sigma[3]))
		stop("If the covariance is present, the two random effects must both be present.")

	intcfg.def = eval(formal.args[["intcfg"]], envir = sys.frame(sys.parent()))
	missed = names(intcfg.def)[!(names(intcfg.def) %in% names(intcfg))]
	intcfg[missed] = intcfg.def[missed]

	if (!is.null(debug) && any(debug=="integration")) intcfg$trace = TRUE
	else intcfg$trace = FALSE

	optimcfg.def = eval(formal.args[["optimcfg"]], envir = sys.frame(sys.parent()))
	missed = names(optimcfg.def)[!(names(optimcfg.def) %in% names(optimcfg))]
	optimcfg[missed] = optimcfg.def[missed]

	if (!is.null(debug) && any(debug=="optim")) optimcfg$trace = 2
	else optimcfg$trace = 0

	ptime = NULL
	if (basedist=="piecewise") {
		pdist = piecetime(Y, piececuttype, piececut, npiece)
		nbasepara = pdist$nbasepara
		ptime = pdist$ptime
	}
	else if (basedist=="exponential") nbasepara = 1
	else nbasepara = 2

	if(model) {
		list(Y = Y, Z = Z, X = X, ngm = ngm, nbt = nbt, nbasepara = nbasepara,
			 ptime = ptime)
	}
	else {
		para = c(initpara(formula, cureform, basedist, bt, gm, basepara, nbt,
						  ngm, nbasepara, data, rawY), list(sigma = sigma))

		data2 <- list(y = Y, x = X, z = Z, id = clusterx, ncluster = nclusterx,
			offsetx = offsetx, offsetz = offsetz)

		if (!is.null(debug) && any(debug=="intcure")) {
			fcat("Initial values:\n")
			fcat("sigma =", para$sigma, "\n")
			fcat("bt =", para$bt, "\n")
			fcat("gm =", para$gm, "\n")
			fcat("basepara =", para$basepara, "\n")
		}

		fit <- intcure.fit(para, data2, intcfg, optimcfg, basedist, funval, ptime)

		if (is.null(names(fit$bt))) names(fit$bt) <- dimnames(X)[[2]]
		if (is.null(names(fit$gm))) names(fit$gm) <- dimnames(Z)[[2]]

		fit$call <- call
		fit$basedist <- basedist
		fit$ptime <- ptime
		fit$n <- datalen
		class(fit) <- "intcure"
		fit
	}
}

parv2l <- function(vpara, lpara)
{
	nbt = length(lpara$bt)
	ngm = length(lpara$gm)
	nbasepara = length(lpara$basepara)
	np = length(vpara)

	if (nbt > 0) bt = vpara[1:nbt]
	else bt = NULL

	gm = vpara[nbt+(1:ngm)]
	basepara = vpara[nbt + ngm + (1:nbasepara)]

	sigma = lpara$sigma
	if (np > nbt + ngm + nbasepara)
		sigma[!is.na(sigma)] = vpara[(nbt + ngm + nbasepara + 1):np]

	list(bt = bt, gm = gm, basepara = basepara, sigma = sigma)
}

piecetime = function(Y, piececuttype, piececut, npiece)
{
	maxtime <- max(Y[Y[, 2]==1, 1])
	mintime <- min(Y[Y[, 2]==1, 1])

	if (is.null(piececut)) {
		nbasepara = npiece
		if (piececuttype == "even") 
			ptime <- mintime + (1:(npiece-1))*(maxtime - mintime)/npiece
		else ptime <- quantile(Y[Y[, 2]==1, 1],
						probs = seq(0, 1, 1/npiece))[-c(1, npiece+1)]
	}
	else {
		if (max(piececut) >= maxtime || min(piececut) <= mintime)
			stop("piececut values are beyond survival time range.")
		nbasepara = length(piececut) + 1
		ptime <- piececut
	}
	list(ptime = ptime, nbasepara = nbasepara)
}
	

baseline <- function(basepara, time, basedist, ptime)
{
	nbasepara = length(basepara)
	if (basedist == "exponential") {
		rate <- exp(basepara)
		h <- rep(rate, length(time))
		S <- exp(-rate*time)
		list(h = h, S = S)
	}
	else if (basedist == "weibull") {
		p <- exp(basepara[1])
		rate <- exp(basepara[2])
		h <- rate*p*time^(p-1.0)
		S <- exp(-rate*time^p)
		list(h = h, S = S)
	}
	else if (basedist == "piecewise") {
		hazard <- exp(basepara)
		num <- as.numeric(cut(time, breaks = c(0, ptime, Inf), right = FALSE,
			include.lowest = TRUE))
		h <- hazard[num]
		cumhaz <- c(0, cumsum(hazard[-nbasepara]*diff(c(0, ptime))))
		chaz <- h*(time - c(0, ptime[-nbasepara])[num])
		S <- exp(-(cumhaz[num] + chaz))
		list(h = h, S = S)
	}
	else stop("baseline distribution is not defined.")
}

basesurvfun = function(time, intcureobj)
{
	basedist = intcureobj$basedist
	basepara = intcureobj$basepara
	nbasepara = length(basepara)

	if (basedist == "exponential") {
		rate <- exp(basepara)
		exp(-rate*time)
	}
	else if (basedist == "weibull") {
		p <- exp(basepara[1])
		rate <- exp(basepara[2])
		exp(-rate*time^p)
	}
	else if (basedist == "piecewise") {
		ptime = intcureobj$ptime
		hazard <- exp(basepara)
		num <- as.numeric(cut(time, breaks = c(0, ptime), right = FALSE,
			include.lowest = TRUE))
		cumhaz <- c(0, cumsum(hazard*diff(c(0,
			ptime)))[-nbasepara])
		chaz <- hazard[num]*(time - c(0, ptime[-nbasepara])[num])
		exp(-(cumhaz[num] + chaz))
	}
	else stop("baseline distribution is not defined.")
}

clusterlik <- function(bx, gz, basesurv, basehaz, cens)
{
	uncure <- exp(gz)/(1+exp(gz))
	surv <- basesurv^(exp(bx))
	val <- sum(ifelse(cens==1, gz - log(1.0 + exp(gz))+log(basehaz) + bx +
		exp(bx)*log(basesurv),
		-log(1.0 + exp(gz)) + log(1.0 + exp(gz)*surv)))
	val
}

clusterliku <- function(u, bx, gz, basesurv, basehaz, cens, sigma)
{
	sapply(u, function(x, bx, gz, basesurv, basehaz, cens, sigma){
		uncure <- exp(gz)/(1+exp(gz))
		surv <- basesurv^(exp(bx+x))
		exp(sum(ifelse(cens==1, gz - log(1+exp(gz)) +
			log(basehaz) + bx + x + exp(bx+x)*log(basesurv),
			log(1-uncure+uncure*surv))) + dnorm(x, sd = exp(sigma), log = TRUE))
	}, bx, gz, basesurv, basehaz, cens, sigma)
}

clusterlikv <- function(v, bx, gz, basesurv, basehaz, cens, sigma)
{
	sapply(v, function(x, bx, gz, basesurv, basehaz, cens, sigma){
		uncure <- exp(gz + x)/(1+exp(gz + x))
		surv <- basesurv^(exp(bx))
		exp(sum(ifelse(cens==1, gz + x - log(1+exp(gz + x)) +
			log(basehaz) + bx + exp(bx)*log(basesurv),
			log(1-uncure+uncure*surv))) + dnorm(x, sd = exp(sigma), log = TRUE))
	}, bx, gz, basesurv, basehaz, cens, sigma)
}

clusterlikuv <- function(uv, bx, gz, basesurv, basehaz, cens, sigma)
{
	uncure <- exp(gz + uv[2])/(1+exp(gz + uv[2]))
	surv <- basesurv^(exp(bx+uv[1]))

	val <- sum(ifelse(cens==1, gz + uv[2] - log(1+exp(gz + uv[2])) +
			log(basehaz) + bx+uv[1] + exp(bx+uv[1])*log(basesurv),
			log(1-uncure+uncure*surv)))
	tmp <- dnorm(uv[1], sd = exp(sigma[1]), log = TRUE) +
				dnorm(uv[2], sd = exp(sigma[2]), log = TRUE)
	exp(val + tmp)
}

clusterlikuv2 <- function(uv, bx, gz, basesurv, basehaz, cens, sigma)
{
	uncure <- exp(gz + uv[2])/(1+exp(gz + uv[2]))
	surv <- basesurv^(exp(bx+uv[1]))
	rho = (exp(2.0*sigma[3])-1.0)/(exp(2.0*sigma[3])+1.0)
	cvar <- rho*exp(sigma[1]+sigma[2])
	val <- sum(ifelse(cens==1, gz + uv[2] - log(1+exp(gz + uv[2])) +
			log(basehaz) + bx+uv[1] + exp(bx+uv[1])*log(basesurv),
			log(1-uncure+uncure*surv)))
	tmp <- dmvnorm(uv, sigma = matrix(c(exp(2.0*sigma[1]), cvar, cvar,
				exp(2.0*sigma[2])), nrow = 2), log = TRUE)
	exp(val + tmp)
}

loglikfun <- function(vpara, data, basedist, intcfg, lpara, failinfo)
{
	pars <- parv2l(vpara, lpara)
	time <- data$y[, 1]

	if (is.null(data$x) || length(pars$bt)==0) bx = rep(0, length(time))
	else bx <- data$x %*% pars$bt

	gz <- data$z %*% pars$gm
	sigma <- pars$sigma
	
	hS <- baseline(pars$basepara, time, basedist, failinfo)

	if (intcfg$trace) fcat("Parameter =", vpara, "\n")

	val <- 0.0
	for (i in 1:data$ncluster) {
		member <- data$id==i
		cens <- data$y[member, 2]
		basesurv <- hS$S[member]
		basehaz <- hS$h[member]
		cbx <- bx[member]
		cgz <- gz[member]

		if (all(is.na(lpara$sigma))) {
			intv = list()
			intv$value = clusterlik(cbx, cgz, basesurv, basehaz, cens)
			intv$error <- intv$code <- NULL
		}
		else if (all(is.na(lpara$sigma)[2:3])) {
			intcfg$fun = "integrate"
			intv <- intfun(clusterliku, intcfg, bx = cbx, gz = cgz,
				basesurv = basesurv, basehaz = basehaz, cens = cens,
				sigma = sigma[1])
		}
		else if (all(is.na(lpara$sigma)[c(1, 3)])) {
			intcfg$fun = "integrate"
			intv <- intfun(clusterlikv, intcfg, bx = cbx, gz = cgz,
				basesurv = basesurv, basehaz = basehaz, cens = cens,
				sigma = sigma[2])
		}
		else if (is.na(lpara$sigma[3])) {
			if(is.null(intcfg$fun)) intcfg$fun = "adaptIntegrate"
			intv <- intfun(clusterlikuv, intcfg, bx = cbx, gz = cgz,
				basesurv = basesurv, basehaz = basehaz, cens = cens,
				sigma = sigma[1:2])
		}
		else {
			if(is.null(intcfg$fun)) intcfg$fun = "adaptIntegrate"
			intv <- intfun(clusterlikuv2, intcfg, bx = cbx, gz = cgz,
				basesurv = basesurv, basehaz = basehaz, cens = cens,
				sigma = sigma)
		}

		if (intcfg$trace)
			fcat("Cluster", i, "likelihood =", intv$value, "with code =",
				 intv$code, "and error =", intv$error, "\n")

		val <- val + intv$value
	}
	if (intcfg$trace) fcat("total likelihood =", val, "\n")
	val
}

intfun = function(fn, intcfg, bx, gz, basesurv, basehaz, cens, sigma)
{
	if (intcfg$fun == "integrate") {
		intv = integrate(fn, lower = intcfg$lower[1], upper = intcfg$upper[1],
				bx = bx, gz = gz, basesurv = basesurv, basehaz = basehaz,
				cens = cens, sigma = sigma, rel.tol = intcfg$eps)
		list(value = log(intv$value), error = intv$abs.error, code = intv$message)
	}
	else if (intcfg$fun == "adaptIntegrate") {
		if (!exists("adaptIntegrate"))
			stop("adaptIntegrate cannot be found. It requires cubature package\n")

		intv <- adaptIntegrate(fn, lowerLimit = intcfg$lower,
			upperLimit = intcfg$upper, tol = intcfg$eps, bx = bx, gz = gz,
			basesurv = basesurv, basehaz = basehaz, cens = cens, sigma = sigma)

		list(value = log(intv$integral), error = intv$error, code = intv$returnCode)
	}
	else stop("Unknown integration function in intcfg$fun\n")
}

intcure.fit <- function(lpara, data, intcfg, optimcfg, basedist, funval, failinfo)
{
	vpara <- c(lpara$bt, lpara$gm, lpara$basepara, lpara$sigma[!is.na(lpara$sigma)])

	if (funval) {
		val <- lpara
		val$loglik <- loglikfun(vpara, data = data, basedist = basedist,
			intcfg = intcfg, lpara = lpara, failinfo = failinfo)
		val
	}
	else {
		method = optimcfg$method
		hessian = optimcfg$hessian
		optimcfg$method <- optimcfg$hessian <- NULL
		if (length(optimcfg$ndeps)==1)
			optimcfg$ndeps = rep(optimcfg$ndeps, length(vpara))
		else {
			if (length(optimcfg$ndeps)!=length(vpara))
				stop("ndeps is of the wrong length")
		}
		fit <- optim(vpara, loglikfun, method = method,
			control = c(list(fnscale = -1), optimcfg), data = data,
			basedist = basedist, intcfg = intcfg, lpara = lpara,
			failinfo = failinfo, hessian = hessian)

		val <- parv2l(fit$par, lpara)

		if (all(is.na(lpara$sigma))) val$sigma = NULL
		else val$sigma = val$sigma[(1:3)[!is.na(lpara$sigma)]]

		val$loglik <- fit$value
		if (hessian) val$hessian <- fit$hessian
		val
	}
}


print.intcure <- function(x, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}

	cat("\nMixture cure model with", x$basedist, "baseline\n")
	cat("\nFixed effects in survival model:\n")
	print(x$bt, ...)
	cat("\nBaseline parameters in survival model:\n")
	print(x$basepara, ...)
	cat("\nFixed effects in logistic model:\n")
	print(x$gm, ...)
	cat("\nRandom effects:\n")
	print(x$sigma, ...)

	cat("\nMaximum log-likelihood:", x$loglik, "\n")

	invisible(x)
}

summary.intcure = function(object, ...)
{
	if (is.null(object$hessian)) stop("Hessian matrix is needed")

	varpar = solve(-object$hessian)
	stderr = parv2l(sqrt(diag(varpar)), object[c("bt", "gm", "basepara", "sigma")])
	object$statistics = list(
			bt = cbind(Estimate = object$bt, Stderr = stderr$bt,
				zvalue = object$bt/stderr$bt,
				pvalue = (1-pnorm(abs(object$bt/stderr$bt)))*2),
			gm = cbind(Estimate = object$gm, Stderr = stderr$gm,
				zvalue = object$gm/stderr$gm,
				pvalue = (1-pnorm(abs(object$gm/stderr$gm)))*2),
			basepara = cbind(Estimate = object$basepara, Stderr = stderr$basepara,
				zvalue = object$basepara/stderr$basepara,
				pvalue = (1-pnorm(abs(object$basepara/stderr$basepara)))*2),
			sigma = cbind(Estimate = object$sigma, Stderr = stderr$sigma))
	class(object) = "summary.intcure"
	object
}

print.summary.intcure = function(x, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}

	cat("\nMixture cure model with", x$basedist, "baseline\n")
	cat("\nFixed effects in survival model:\n")
	print(x$statistics$bt, ...)
	cat("\nBaseline parameters in survival model:\n")
	print(x$statistics$basepara, ...)
	cat("\nFixed effects in logistic model:\n")
	print(x$statistics$gm, ...)
	cat("\nRandom effects:\n")
	print(x$statistics$sigma, ...)

	cat("\nMaximum log-likelihood:", x$loglik, "\n")

	invisible(x)
}


coef.intcure <- function(object, ...)
{
	c(object$bt, object$gm, object$basepara, object$sigma)
}

loglik.intcure <- function(x)x$loglik

initpara = function(formula, cureform, basedist, bt, gm, basepara, nbt, ngm,
					nbasepara, data, Y)
{
	phcall <- gsub(" *[+] *cluster[(][^)]*[)] *[+]*", "",
		paste(deparse(formula, width.cutoff = 500), collapse = ""))
	logitcall <- gsub(" *[+] *cluster[(][^)]*[)] *[+]*", "",
		paste(deparse(cureform, width.cutoff = 500), collapse = ""))
	
	lfit = survreg(eval(parse(text = phcall)), data = data,
		dist = "exponential", subset = Y[, 1] <= max(Y[Y[, 2]==1, 1]))
	if (is.null(bt)) bt <- coef(lfit)[-1]
	else names(bt) = names(coef(lfit)[-1])

	maxtime = max(Y[Y[, 2]==1, 1])
	mintime = min(Y[Y[, 2]==1, 1])
	uncure = ifelse(Y[, 2] == 1 | Y[, 1] < mintime, 1,
				ifelse(Y[, 1] > maxtime, 0,
					   1-punif(Y[, 1], min = mintime, max = maxtime)))
	ifit = glm(eval(parse(text = paste("uncure", logitcall))),
			   data = cbind(uncure, data), family=binomial(link=logit))
	if (is.null(gm)) gm <- coef(ifit)
	else names(gm) = names(coef(ifit))
	
	if (length(bt)!=nbt || length(gm)!=ngm)
		stop("The length of initial values doesn't match needs.")
	
	if (is.null(basepara)) {
		if (basedist=="exponential")
			basepara <- c(lograte = unname(-coef(lfit)[1]))
		else if (basedist=="weibull")
			basepara <- c(logshape = 0, lograte = unname(-coef(lfit)[1]))
		else if (basedist=="piecewise")
			basepara <- structure(rep(-coef(lfit)[1], nbasepara),
							names = paste("loghazard", 1:nbasepara, sep = ""))
		else stop("basedist is not valid in baseinit.")
	}
	else {
		if (basedist=="exponential") names(basepara) <- "lograte"
		else if (basedist=="weibull") names(basepara) <- c("logshape",
														   "lograte")
		else if (basedist=="piecewise")
			names(basepara) <- paste("loghazard", 1:nbasepara, sep = "")
		else stop("basedist is not valid in baseinit.")
	}

	if (length(bt) != nbt)
		stop("The number of regression parameters in latency is not correct.")
	if (length(gm) != ngm)
		stop("The number of regression parameters in incidence is not correct.")
	if (length(basepara) != nbasepara)
		stop("The number of baseline parameters is not correct.")

	list(bt = bt, gm = gm, basepara = basepara)
}

fcat = function (...) 
{
    cat(...)
    flush.console()
}
