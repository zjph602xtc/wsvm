tune.control <- function(random = FALSE,
                         nrepeat = 1,
                         repeat.aggregate = mean,
                         sampling = c("cross", "fix", "bootstrap"),
                         sampling.aggregate = mean,
                         sampling.dispersion = sd,
                         cross = 10,
                         fix = 2 / 3,
                         nboot = 10,
                         boot.size = 9 / 10,
                         best.model = TRUE,
                         performances = TRUE,
                         error.fun = NULL) {
    structure(list(random = random,
                   nrepeat = nrepeat,
                   repeat.aggregate = repeat.aggregate,
                   sampling = match.arg(sampling),
                   sampling.aggregate = sampling.aggregate,
                   sampling.dispersion = sampling.dispersion,
                   cross = cross,
                   fix = fix,
                   nboot = nboot,
                   boot.size = boot.size,
                   best.model = best.model,
                   performances = performances,
                   error.fun = error.fun
    ),
    class = "tune.control"
    )
}

tune_wsvm <- function(train.x, train.y = NULL, weight, use_zero_weight = FALSE, pre.check = TRUE, data = list(), validation.x = NULL, validation.y = NULL, validation.weight = NULL, weigthed.error = TRUE, ranges = NULL, predict.func = predict, tunecontrol = tune.control(), ...) {
    method <- 'wsvm'
    call <- match.call()
    weight <- as.vector(weight)
    if (!is.null(validation.weight))
        validation.weight <- as.vector(validation.weight)

    ## internal helper functions
    resp <- function(formula, data) {

        model.response(model.frame(formula, data))
    }

    # classAgreement <- function (tab) {
    #     n <- sum(tab)
    #     if (!is.null(dimnames(tab))) {
    #         lev <- intersect(colnames(tab), rownames(tab))
    #         p0 <- sum(diag(tab[lev, lev])) / n
    #     } else {
    #         m <- min(dim(tab))
    #         p0 <- sum(diag(tab[1:m, 1:m])) / n
    #     }
    #     p0
    # }

    ## parameter handling
    if (tunecontrol$sampling == "cross")
        validation.x <- validation.y <- NULL
    useFormula <- is.null(train.y)
    if (useFormula && (is.null(data) || length(data) == 0))
        data <- model.frame(train.x)
    if (is.vector(train.x)) train.x <- t(t(train.x))
    if (is.data.frame(train.y))
        train.y <- as.matrix(train.y)

    if (useFormula && length(weight)!=NROW(data))
        stop("The length of weight is not equal to the number of samples.")
    if (!useFormula && length(weight)!=NROW(train.x))
        stop("The length of weight is not equal to the number of samples.")
    if (!is.null(validation.x) && !is.null(validation.weight) && length(validation.weight)!=NROW(validation.x))
        stop("The length of validation.weight is not equal to the number of samples in validation.x.")

    ## use_zero_weight = F
    if (!use_zero_weight){
        if (useFormula)
            data <- data[weight>0,,drop=FALSE]
        else{
            train.x <- train.x[weight>0,,drop=FALSE]
            train.y <- train.y[weight>0]
        }
        weight <- weight[weight>0]
        if (!is.null(validation.x) & !is.null(validation.weight)){
            validation.x <- validation.x[validation.weight>0,,drop=FALSE]
            validation.y <- validation.y[validation.weight>0]
            validation.weight <- validation.weight[validation.weight>0]
        }
    }

    ## construct parameters
    parameters <- if (is.null(ranges))
        data.frame(dummyparameter = 0)
    else
        expand.grid(ranges)
    p <- nrow(parameters)
    if (!is.logical(tunecontrol$random)) {
        if (tunecontrol$random < 1)
            stop("random must be a strictly positive integer")
        if (tunecontrol$random > p) tunecontrol$random <- p
        parameters <- parameters[sample(1:p, tunecontrol$random),]
        p <- nrow(parameters)
    }
    model.variances <- model.errors <- c()

    ## par for the try model
    pars_tmp <- if (is.null(ranges))
        NULL
    else
        lapply(parameters[1,,drop = FALSE], unlist)


    ## prepare training indices
    if (!is.null(validation.x)) tunecontrol$fix <- 1
    n <- NROW(if (useFormula) data else train.x)
    if (tunecontrol$sampling == "cross") {
        if (tunecontrol$cross > n)
            stop(sQuote("cross"), " must not exceed sampling size!")
        if (tunecontrol$cross == 1)
            stop(sQuote("cross"), " must be greater than 1!")
    }
    if (pre.check){
        try_num <- 0
        repeat{
            perm.ind <- sample(n)
            train.ind <- if (tunecontrol$sampling == "cross")
                tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
            else if (tunecontrol$sampling == "fix")
                list(perm.ind[1:trunc(n * tunecontrol$fix)])
            else ## bootstrap
                lapply(1:tunecontrol$nboot,
                       function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))

            model_tmp <- if (useFormula)
                lapply(1:length(train.ind),function(sample)suppressWarnings(do.call(method, c(list(formula = train.x, weight = weight, data = data, subset = train.ind[[sample]]), pars_tmp, fitted = FALSE, list(...)))))
            else
                lapply(1:length(train.ind),function(sample)suppressWarnings(do.call(method, c(list(x = train.x[train.ind[[sample]],,drop=FALSE], y = train.y[train.ind[[sample]]], weight = weight[train.ind[[sample]]]), pars_tmp, fitted = FALSE, list(...)))))

            if (all(sapply(model_tmp, function(model){model$tot.nSV >= 1})))
                break

            try_num <- try_num + 1
            warning("Model fitting fails on the training data. Try to re-partition the data.")
            if (try_num > 10)
                stop("Model fitting fails on the training data. This may be due to misspecification of the parameters, too many zero weights in the training data, .....")
        }
    }else{
        perm.ind <- sample(n)
        train.ind <- if (tunecontrol$sampling == "cross")
            tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
        else if (tunecontrol$sampling == "fix")
            list(perm.ind[1:trunc(n * tunecontrol$fix)])
        else ## bootstrap
            lapply(1:tunecontrol$nboot,
                   function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))
    }

    ## - loop over all models
    for (para.set in 1:p) {
        sampling.errors <- c()

        ## - loop over all training samples
        for (sample in 1:length(train.ind)) {
            repeat.errors <- c()

            ## - repeat training `nrepeat' times
            for (reps in 1:tunecontrol$nrepeat) {

                ## train one model
                pars <- if (is.null(ranges))
                    NULL
                else
                    lapply(parameters[para.set,,drop = FALSE], unlist)

                model <- if (useFormula)
                    suppressWarnings(do.call(method, c(list(formula = train.x, weight = weight, data = data, subset = train.ind[[sample]]), pars, list(...))))
                else
                    suppressWarnings(do.call(method, c(list(x = train.x[train.ind[[sample]],,drop=FALSE], y = train.y[train.ind[[sample]]], weight = weight[train.ind[[sample]]]), pars, list(...))))

                ## predict validation set
                pred <- predict.func(model,
                                     if (!is.null(validation.x))
                                         validation.x
                                     else if (useFormula)
                                         data[-train.ind[[sample]],,drop = FALSE]
                                     else if (inherits(train.x, "matrix.csr"))
                                         train.x[-train.ind[[sample]],]
                                     else
                                         train.x[-train.ind[[sample]],,drop = FALSE]
                )

                ## compute performance measure
                true.y <- if (!is.null(validation.y))
                    validation.y
                else if (useFormula) {
                    if (!is.null(validation.x))
                        resp(train.x, validation.x)
                    else
                        resp(train.x, data[-train.ind[[sample]],])
                } else
                    train.y[-train.ind[[sample]]]
                if (is.null(true.y)) true.y <- rep(TRUE, length(pred))

                if (weigthed.error){
                    true.w <- if (!is.null(validation.weight))
                        validation.weight
                    else
                        weight[-train.ind[[sample]]]
                }
                if (is.null(true.w)) true.w <- rep(1, length(pred))

                repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun))
                    tunecontrol$error.fun(true.y, pred, true.w)
                else if ((is.logical(true.y) || is.factor(true.y)) && (is.logical(pred) || is.factor(pred) || is.character(pred))) ## classification error
                    1 - weighted.mean(pred == true.y, true.w, na.rm = T)
                else if (is.numeric(true.y) && is.numeric(pred)) ## mean squared error
                    weighted.mean((pred - true.y)^2, true.w, na.rm = T)
                else
                    stop("Dependent variable has wrong type!")
            }
            repeat.errors <- repeat.errors[!is.na(repeat.errors)]
            sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
        }
        model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
        model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors)
    }

    ## return results
    best <- which.min(model.errors)
    pars <- if (is.null(ranges))
        NULL
    else
        lapply(parameters[best,,drop = FALSE], unlist)
    structure(list(best.parameters  = parameters[best,,drop = FALSE],
                   best.performance = model.errors[best],
                   method           = if (!is.character(method))
                       deparse(substitute(method)) else method,
                   nparcomb         = nrow(parameters),
                   train.ind        = train.ind,
                   sampling         = switch(tunecontrol$sampling,
                                             fix = "fixed training/validation set",
                                             bootstrap = "bootstrapping",
                                             cross = if (tunecontrol$cross == n) "leave-one-out" else
                                                 paste(tunecontrol$cross,"-fold cross validation", sep="")
                   ),
                   performances     = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances),
                   best.model       = if (tunecontrol$best.model) {
                       modeltmp <- if (useFormula)
                           suppressWarnings(do.call(method, c(list(formula = train.x,
                                                                   weight = weight,
                                                                   data = data),
                                                              pars, list(...))))
                       else
                           suppressWarnings(do.call(method, c(list(x = train.x,
                                                                   y = train.y,
                                                                   weight = weight),
                                                              pars, list(...))))
                       call[[1]] <- as.symbol("best.tune_wsvm")
                       modeltmp$call <- call
                       modeltmp
                   }
    ),
    class = c("tune_wsvm","tune")
    )
}

best.tune_wsvm <- function(...) {
    call <- match.call()
    modeltmp <- tune_wsvm(...)$best.model
    modeltmp$call <- call
    modeltmp
}

print.tune_wsvm <- function(x, ...) {
    if (x$nparcomb > 1) {
        cat("\nParameter tuning of ", sQuote(x$method), ":\n\n", sep="")
        cat("- sampling method:", x$sampling,"\n\n")
        cat("- best parameters:\n")
        tmp <- x$best.parameters
        rownames(tmp) <- ""
        print(tmp)
        cat("\n- best performance:", x$best.performance, "\n")
        cat("\n")
    } else {
        cat("\nError estimation of ", sQuote(x$method), " using ", x$sampling, ": ",
            x$best.performance, "\n\n", sep="")
        cat("There is maybe only 1 pair of parameters, or some parameters are invalid. Check the 'range' parameter and your dataset. ")
    }
}

summary.tune_wsvm <- function(object, ...)
    structure(object, class = c("summary.tune_wsvm","summary_tune"))

print.summary.tune_wsvm <- function(x, ...) {
    print.tune_wsvm(x)
    if (!is.null(x$performances) && (x$nparcomb > 1)) {
        cat("- Detailed performance results:\n")
        print(x$performances)
        cat("\n")
    }
}

hsv_palette <- function(h = 2/3, from = 0.7, to = 0.2, v = 1)
    function(n) hsv(h = h, s = seq(from, to, length.out = n), v = v)

plot.tune_wsvm <- function(x,
                           type=c("contour","perspective"),
                           theta=60,
                           col="lightblue",
                           main = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           swapxy = FALSE,
                           transform.x = NULL,
                           transform.y = NULL,
                           transform.z = NULL,
                           color.palette = hsv_palette(),
                           nlevels = 20,
                           ...)
{
    if (is.null(x$performances))
        stop("Object does not contain detailed performance measures!")
    k <- ncol(x$performances)
    if (k > 4) stop("Cannot visualize more than 2 parameters")
    type = match.arg(type)

    if (is.null(main))
        main <- paste("Performance of `", x$method, "'", sep="")

    if (k == 3)
        plot(x$performances[,1:2], type = "b", main = main)
    else  {
        if (!is.null(transform.x))
            x$performances[,1] <- transform.x(x$performances[,1])
        if (!is.null(transform.y))
            x$performances[,2] <- transform.y(x$performances[,2])
        if (!is.null(transform.z))
            x$performances[,3] <- transform.z(x$performances[,3])
        if (swapxy)
            x$performances[,1:2] <- x$performances[,2:1]
        x <- xtabs(error~., data = x$performances[,-k])
        if (is.null(xlab)) xlab <- names(dimnames(x))[1 + swapxy]
        if (is.null(ylab)) ylab <- names(dimnames(x))[2 - swapxy]
        if (type == "perspective")
            persp(x=as.double(rownames(x)),
                  y=as.double(colnames(x)),
                  z=x,
                  xlab=xlab,
                  ylab=ylab,
                  zlab="accuracy",
                  theta=theta,
                  col=col,
                  ticktype="detailed",
                  main = main,
                  ...
            )
        else
            filled.contour(x=as.double(rownames(x)),
                           y=as.double(colnames(x)),
                           xlab=xlab,
                           ylab=ylab,
                           nlevels=nlevels,
                           color.palette = color.palette,
                           main = main,
                           x, ...)
    }
}


# tune.wsvm <- function(x, y = NULL, weight, data = NULL, degree = NULL, gamma = NULL,
#                      coef0 = NULL, cost = NULL, nu = NULL, class.weights = NULL,
#                      epsilon = NULL, ...) {
#     call <- match.call()
#     call[[1]] <- as.symbol("best.wsvm")
#     ranges <- list(degree = degree, gamma = gamma,
#                    coef0 = coef0, cost = cost, nu = nu,
#                    class.weights = class.weights, epsilon = epsilon)
#     ranges[sapply(ranges, is.null)] <- NULL
#     if (length(ranges) < 1)
#         ranges = NULL
#     modeltmp <- if (inherits(x, "formula"))
#         tune_wsvm("wsvm", train.x = x, weight = weight, data = data, ranges = ranges, ...)
#     else
#         tune_wsvm("wsvm", train.x = x, train.y = y, weight = weight, ranges = ranges, ...)
#     if (!is.null(modeltmp$best.model))
#         modeltmp$best.model$call <- call
#     modeltmp
# }
#
# best.wsvm <- function(x, tunecontrol = tune.control(), ...) {
#     call <- match.call()
#     tunecontrol$best.model = TRUE
#     modeltmp <- tune.wsvm(x, ..., tunecontrol = tunecontrol)$best.model
#     modeltmp$call <- call
#     modeltmp
# }
