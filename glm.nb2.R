## This file contains a modification of the function from MASS to
## fit a negative binomial GLM function.  It is modified to default
## to fitting a poisson regression model when the estimated
## dispersion is very low.  This is a practical modification for the
## purpose of running simulations.  The modification is not guaranteed
## to be generally useful, and may have detrimental effects under
## some circumstancers.  Its general use is not recommended.


theta.ml2 <- function (y, mu, n = sum(weights), weights, limit = 10, eps = .Machine$double.eps^0.25, 
    trace = FALSE) 
{
    score <- function(n, th, mu, y, w) sum(w * (digamma(th + 
        y) - digamma(th) + log(th) + 1 - log(th + mu) - (y + 
        th)/(mu + th)))
    info <- function(n, th, mu, y, w) sum(w * (-trigamma(th + 
        y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu + 
        th)^2))
    if (inherits(y, "lm")) {
        mu <- y$fitted.values
        y <- if (is.null(y$y)) 
            mu + residuals(y)
        else y$y
    }
    if (missing(weights)) 
        weights <- rep(1, length(y))
    t0 <- n/sum(weights * (y/mu - 1)^2)
    it <- 0
    del <- 1
    if (trace) 
        message(sprintf("theta.ml: iter %d 'theta = %f'", it, 
            signif(t0)), domain = NA)
    while ((it <- it + 1) < limit && abs(del) > eps & t0<1000) {
        t0 <- abs(t0)
        del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, 
            mu, y, weights))
        t0 <- t0 + del
        if (trace) 
            message("theta.ml: iter", it, " theta =", signif(t0))
    }
    if (t0 < 0) {
        t0 <- 0
        warning("estimate truncated at zero")
        attr(t0, "warn") <- gettext("estimate truncated at zero")
    }
    if (it == limit) {
        warning("iteration limit reached")
        attr(t0, "warn") <- gettext("iteration limit reached")
    }
    attr(t0, "SE") <- sqrt(1/i)
    t0
}




glm.nb2<-function (formula, data, weights, subset, na.action, start = NULL, 
    etastart, mustart, control = glm.control(...), method = "glm.fit", 
    model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ..., 
    init.theta, link = log) 
{
    loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th + 
        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
        log(mu + (y == 0)) - (th + y) * log(th + mu)))
    link <- substitute(link)
    fam0 <- if (missing(init.theta)) 
        do.call("poisson", list(link = link))
    else do.call("negative.binomial", list(theta = init.theta, 
        link = link))
    mf <- Call <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    if (method == "model.frame") 
        return(mf)
    Y <- model.response(mf, "numeric")
    X <- if (!is.empty.model(Terms)) 
        model.matrix(Terms, mf, contrasts)
    else matrix(, NROW(Y), 0)
    w <- model.weights(mf)
    if (!length(w)) 
        w <- rep(1, nrow(mf))
    else if (any(w < 0)) 
        stop("negative weights not allowed")
    offset <- model.offset(mf)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    n <- length(Y)
    if (!missing(method)) {
        if (!exists(method, mode = "function")) 
            stop(gettextf("unimplemented method: %s", sQuote(method)), 
                domain = NA)
        glm.fitter <- get(method)
    }
    else {
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
    }
    if (control$trace > 1) 
        message("Initial fit:")
    fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart, 
        mustart = mustart, offset = offset, family = fam0, control = list(maxit = control$maxit, 
            epsilon = control$epsilon, trace = control$trace > 
                1), intercept = attr(Terms, "intercept") > 0)
    class(fit) <- c("glm", "lm")
    mu <- fit$fitted.values
  th <- as.vector(theta.ml2(Y, mu, sum(w), w, limit = control$maxit, 
        trace = control$trace > 2))
  if (control$trace > 1) 
        message(gettextf("Initial value for 'theta': %f", signif(th)), 
            domain = NA)
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, Y, w)
    Lm0 <- Lm + 2 * d1

    while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - 
        Lm)/d1 + abs(del)/d2) > control$epsilon ) {
      eta <- g(mu)
        fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta, 
            offset = offset, family = fam, control = list(maxit = control$maxit, 
                epsilon = control$epsilon, trace = control$trace > 
                  1), intercept = attr(Terms, "intercept") > 
                0)
        t0 <- th
        th <- theta.ml2(Y, mu, sum(w), w, limit = control$maxit, 
            trace = control$trace > 2)
        fam <- do.call("negative.binomial", list(theta = th, 
            link = link))
        mu <- fit$fitted.values
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, Y, w)
        if (control$trace) {
            Ls <- loglik(n, th, Y, Y, w)
            Dev <- 2 * (Ls - Lm)
            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
                iter, signif(th), signif(Dev)), domain = NA)
        }
    }
    if (!is.null(attr(th, "warn"))) 
        fit$th.warn <- attr(th, "warn")
    if (iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    if (length(offset) && attr(Terms, "intercept")) {
        null.deviance <- if (length(Terms)) 
            glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, 
                offset = offset, family = fam, control = list(maxit = control$maxit, 
                  epsilon = control$epsilon, trace = control$trace > 
                    1), intercept = TRUE)$deviance
        else fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit$contrasts <- attr(X, "contrasts")
    fit$xlevels <- .getXlevels(Terms, mf)
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    if(th>=10000){
      fit <- glm.fit(X,Y,family=poisson())
      attr(fit,"class") <- "glm"
    }
    fit
}
