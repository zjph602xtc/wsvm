wsvm <-
  function (x, ...)
    UseMethod ("wsvm")

wsvm.formula <-
  function (formula, weight = NULL, data = NULL, ..., subset, na.action = na.omit, scale = TRUE)
  {
    call <- match.call()
    if(is.null(weight)) stop(sQuote("weight"), " must not be NULL! If you don't have subject weights, use svm() in package e1071.")
    if (!inherits(formula, "formula"))
      stop("method is only for formula objects")
    m <- match.call(expand.dots = FALSE)
    if (inherits(eval.parent(m$data), "matrix"))
      m$data <- as.data.frame(eval.parent(m$data))
    m$... <- NULL
    m$scale <- NULL
    m[[1L]] <- quote(stats::model.frame)
    m$na.action <- na.action
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    x <- model.matrix(Terms, m)
    y <- model.extract(m, "response")
    weight <- model.extract(m, "weight")
    attr(weight,"na.action") <- attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")
    if (length(scale) == 1)
      scale <- rep(scale, ncol(x))
    if (any(scale)) {
      remove <- unique(c(which(labels(Terms) %in%
                                 names(attr(x, "contrasts"))),
                         which(!scale)
      )
      )
      scale <- !attr(x, "assign") %in% remove
    }
    ret <- wsvm.default (x = x, y = y, weight = weight, scale = scale, ..., na.action = na.action)
    ret$call <- call
    ret$call[[1]] <- as.name("wsvm")
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action")))
      ret$na.action <- attr(m, "na.action")
    class(ret) <- c("wsvm.formula", class(ret))
    return (ret)
  }

wsvm.default <-
  function (x,
            y           = NULL,
            weight      = NULL,
            scale       = TRUE,
            type        = NULL,
            kernel      = "radial",
            degree      = 3,
            gamma       = if (is.vector(x)) 1 else 1 / ncol(x),
            coef0       = 0,
            cost        = 1,
            nu          = 0.5,
            class.weights = NULL,
            cachesize   = 100,
            tolerance   = 0.001,
            epsilon     = 0.1,
            shrinking   = TRUE,
            cross       = 0,
            probability = FALSE,
            fitted      = TRUE,
            ...,
            subset,
            na.action = na.omit)
  {
    yorig <- y
    if(inherits(x, "Matrix")) {
      loadNamespace("SparseM")
      loadNamespace("Matrix")
      x <- as(x, "matrix.csr")
    }
    if(inherits(x, "simple_triplet_matrix")) {
      loadNamespace("SparseM")
      ind <- order(x$i, x$j)
      x <- new("matrix.csr",
               ra = x$v[ind],
               ja = x$j[ind],
               ia = as.integer(cumsum(c(1, tabulate(x$i[ind])))),
               dimension = c(x$nrow, x$ncol))
    }
    if (sparse <- inherits(x, "matrix.csr"))
      loadNamespace("SparseM")

    ## NULL parameters?
    if(is.null(degree)) stop(sQuote("degree"), " must not be NULL!")
    if(is.null(gamma)) stop(sQuote("gamma"), " must not be NULL!")
    if(is.null(coef0)) stop(sQuote("coef0"), " must not be NULL!")
    if(is.null(cost)) stop(sQuote("cost"), " must not be NULL!")
    if(is.null(nu)) stop(sQuote("nu"), " must not be NULL!")
    if(is.null(epsilon)) stop(sQuote("epsilon"), " must not be NULL!")
    if(is.null(tolerance)) stop(sQuote("tolerance"), " must not be NULL!")
    if(is.null(weight)) stop(sQuote("weight"), " must not be NULL! If you don't have subject weights, use svm() in package e1071.")

    weightorig <- weight
    xhold   <- if (fitted) x else NULL
    x.scale <- y.scale <- NULL
    formula <- inherits(x, "wsvm.formula")

    ## determine model type
    if (is.null(type)) type <-
      if (is.null(y)) "one-classification"
    else if (is.factor(y)) "C-classification"
    else "eps-regression"

    type <- pmatch(type, c("C-classification",
                           "nu-classification",
                           "one-classification",
                           "eps-regression",
                           "nu-regression"), 99) - 1

    if (type > 10) stop("wrong type specification!")

    kernel <- pmatch(kernel, c("linear",
                               "polynomial",
                               "radial",
                               "sigmoid",
                               "precomputed"), 99) - 1

    if (kernel > 10) stop("wrong kernel specification!")
    if (kernel==4){
      if (inherits(x, "kernelMatrix")) x <- x@.Data
      if (dim(x)[1]!=dim(x)[2]) stop("the kernel matrix must be square!")
      if (!missing(subset)) stop("cannot use the 'subset' argument when the kernel is precalculated")
      if (any(is.na(x))) stop("cannot have missing values in x when the kernel is precalculated")
      if (any(is.na(y))) stop("cannot have missing values in y when the kernel is precalculated")
      if (sparse)  stop("cannot use sparse x when the kernel is precalculated")
      scale <- FALSE
    }

    nac <- attr(x, "na.action")

    ## scaling, subsetting, and NA handling
    if (sparse) {
      scale <- rep(FALSE, ncol(x))
      if(!is.null(y)) na.fail(y)
      na.fail(weight)
      x <- SparseM::t(SparseM::t(x)) ## make shure that col-indices are sorted
    } else {
      x <- as.matrix(x)
      weight <- as.numeric(weight)


      ## subsetting and na-handling for matrices
      if (!formula) {
        if (length(weight)!=NROW(x)){
          stop("The length of weight must be the same as the length of the number of rows in x.")
        }
        if (!missing(subset)) {
          x <- x[subset,]
          y <- y[subset]
          weight <- weight[subset]
          if (!is.null(xhold))
            xhold <- as.matrix(xhold)[subset,]
        }
        if (is.null(y)){
          df <- na.action(data.frame(weight, x))
          x <- as.matrix(df[,-1], rownames.force = TRUE)
          weight <- df[,1]
          nac <-
            attr(x, "na.action") <-
            attr(weight, "na.action") <-
            attr(df, "na.action")
        }else {
          df <- na.action(data.frame(1:NROW(y),y, x))
          yorig <- df[,2]
          weightorig <- weight[df[,1]]
          if (!is.null(xhold))
            xhold <- as.matrix(df[,-(1:2)], rownames.force = TRUE)
          df <- na.action(data.frame(y, weight, x))
          y <- df[,1]
          x <- as.matrix(df[,-(1:2)], rownames.force = TRUE)
          weight <- df[,2]
          nac <-
            attr(x, "na.action") <-
            attr(y, "na.action") <-
            attr(weight, "na.action") <-
            attr(df, "na.action")
        }
      }

      ## scaling
      if (length(scale) == 1)
        scale <- rep(scale, ncol(x))
      if (any(scale)) {
        co <- !apply(x[,scale, drop = FALSE], 2, var)
        if (any(co, na.rm = T)) {
          warning(paste("Variable(s)",
                        paste(sQuote(colnames(x[,scale,
                                                drop = FALSE])[co]),
                              sep="", collapse=" and "),
                        "constant. Cannot scale data.")
          )
          scale <- rep(FALSE, ncol(x))
        } else if (any(is.na(co))){
          stop("x is not a design matrix which only contains numeric values. Transform x into a design matrix and set 'scale' to FALSE for categorical variables. (OR use a formula). See example XXXXX.")
        } else {
          xtmp <- scale(x[,scale])
          x[,scale] <- xtmp
          x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
          if (is.numeric(y) && (type > 2)) {
            # yorig <- y
            y <- scale(y)
            y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
            y <- as.vector(y)
          }
        }
      }else{
        co <- !apply(x[,scale, drop = FALSE], 2, var)
        if (any(is.na(co))){
          stop("x is not a design matrix which only contains numeric values. Transform x into a design matrix and set 'scale' to FALSE for categorical variables. (OR use a formula). See example XXXXX.")
        }
      }
    }

    ## further parameter checks
    nr <- nrow(x)
    if (cross > nr)
      stop(sQuote("cross"), " cannot exceed the number of observations!")

    ytmp <- y
    attributes(ytmp) <- NULL
    if (!is.vector(ytmp) && !is.factor(y) && type != 2)
      stop("y must be a vector or a factor.")
    if (type != 2 && length(y) != nr)
      stop("x and y don't match.")

    if (any(is.infinite(weight)))
      stop("weight contains infinity.")

    if (all(weight==0))
      warning("all weights are zeros.")

    if (any(weight<=0))
      warning("drop subjects with zero or negative weights.")

    if (cachesize < 0.1)
      cachesize <- 0.1

    if (type > 2 && !is.numeric(y))
      stop("Need numeric dependent variable for regression.")

    lev <- NULL
    weightlabels <- NULL

    ## in case of classification: transform factors into integers
    if (type == 2) # one class classification --> set dummy
      y <- rep(1, nr)
    else
      if (is.factor(y)) {
        lev <- levels(y)
        y <- as.integer(y)
      } else {
        if (type < 3) {
          if(any(as.integer(y) != y))
            stop("dependent variable has to be of factor or integer type for classification mode.")
          y <- as.factor(y)
          lev <- levels(y)
          y <- as.integer(y)
        } else lev <- unique(y)
      }

    if (type < 3 && !is.null(class.weights)) {
      if (is.character(class.weights) && class.weights == "inverse")
        class.weights <- 1 / table(y)
      if (is.null(names(class.weights)))
        stop("Weights have to be specified along with their according level names !")
      weightlabels <- match (names(class.weights), lev)
      if (any(is.na(weightlabels)))
        stop("At least one level name is missing or misspelled.")
    }

    nclass <- 2
    if (type < 2) nclass <- length(lev)

    if (type > 1 && length(class.weights) > 0) {
      class.weights <- NULL
      warning(sQuote("class.weights"), " are set to NULL for regression mode. For classification, use a _factor_ for ", sQuote("y"),
              ", or specify the correct ", sQuote("type"), " argument.")
    }

    err <- empty_string <- paste(rep(" ", 255), collapse = "")

    if (is.null(type)) stop("type argument must not be NULL!")
    if (is.null(kernel)) stop("kernel argument must not be NULL!")
    if (is.null(degree)) stop("degree argument must not be NULL!")
    if (is.null(gamma)) stop("gamma argument must not be NULL!")
    if (is.null(coef0)) stop("coef0 seed argument must not be NULL!")
    if (is.null(cost)) stop("cost argument must not be NULL!")
    if (is.null(nu)) stop("nu argument must not be NULL!")
    if (is.null(cachesize)) stop("cachesize argument must not be NULL!")
    if (is.null(tolerance)) stop("tolerance argument must not be NULL!")
    if (is.null(epsilon)) stop("epsilon argument must not be NULL!")
    if (is.null(shrinking)) stop("shrinking argument must not be NULL!")
    if (is.null(cross)) stop("cross argument must not be NULL!")
    if (is.null(sparse)) stop("sparse argument must not be NULL!")
    if (is.null(probability)) stop("probability argument must not be NULL!")

    if (kernel==4){
      x <- cbind(1:NROW(x), x)
    }
    cret <- .C (C_wsvmtrain,
                ## data
                as.double  (if (sparse) x@ra else t(x)),
                as.integer (nr), as.integer(ncol(x)),
                as.double  (y),
                as.double  (weight),
                ## sparse index info
                as.integer (if (sparse) x@ia else 0),
                as.integer (if (sparse) x@ja else 0),

                ## parameters
                as.integer (type),
                as.integer (kernel),
                as.integer (degree),
                as.double  (gamma),
                as.double  (coef0),
                as.double  (cost),
                as.double  (nu),
                as.integer (weightlabels),
                as.double  (class.weights),
                as.integer (length (class.weights)),
                as.double  (cachesize),
                as.double  (tolerance),
                as.double  (epsilon),
                as.integer (shrinking),
                as.integer (cross),
                as.integer (sparse),
                as.integer (probability),

                ## results
                nclasses = integer  (1),
                nr       = integer  (1), # nr of support vectors
                index    = integer  (nr),
                labels   = integer  (nclass),
                nSV      = integer  (nclass),
                rho      = double   (nclass * (nclass - 1) / 2),
                coefs    = double   (nr * (nclass - 1)),
                sigma    = double   (1),
                probA    = double   (nclass * (nclass - 1) / 2),
                probB    = double   (nclass * (nclass - 1) / 2),

                cresults = double   (cross),
                ctotal1  = double   (1),
                ctotal2  = double   (1),
                error    = err

    )

    if (cret$error != empty_string)
      stop(paste(cret$error, "!", sep=""))

    cret$index <- cret$index[1:cret$nr]

    ret <- list (
      call     = match.call(),
      type     = type,
      kernel   = kernel,
      cost     = cost,
      degree   = degree,
      gamma    = gamma,
      coef0    = coef0,
      nu       = nu,
      epsilon  = epsilon,
      sparse   = sparse,
      scaled   = scale,
      x.scale  = x.scale,
      y.scale  = y.scale,
      weight   = weightorig,

      nclasses = cret$nclasses, #number of classes
      levels   = lev,
      tot.nSV  = cret$nr, #total number of sv
      nSV      = cret$nSV[1:cret$nclasses], #number of SV in diff. classes
      labels   = cret$labels[1:cret$nclasses], #labels of the SVs.
      SV       = if (sparse) SparseM::t(SparseM::t(x[cret$index]))
      else t(t(x[cret$index,,drop = FALSE])), #copy of SV
      index    = cret$index,  #indexes of sv in x
      ##constants in decision functions
      rho      = cret$rho[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
      ##probabilites
      compprob = probability,
      probA    = if (!probability) NULL else
        cret$probA[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
      probB    = if (!probability) NULL else
        cret$probB[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
      sigma    = if (probability) cret$sigma else NULL,
      ##coefficiants of sv
      coefs    = if (cret$nr == 0) NULL else
        t(matrix(cret$coefs[1:((cret$nclasses - 1) * cret$nr)],
                 nrow = cret$nclasses - 1,
                 byrow = TRUE)),
      na.action = nac
    )

    ## cross-validation-results
    if (cross > 0)
      if (type > 2) {
        scale.factor     <- if (any(scale)) crossprod(y.scale$"scaled:scale") else 1;
        ret$MSE          <- cret$cresults * c(scale.factor);
        ret$tot.MSE      <- cret$ctotal1  * c(scale.factor);
        ret$scorrcoeff   <- cret$ctotal2;
      } else {
        ret$accuracies   <- cret$cresults;
        ret$tot.accuracy <- cret$ctotal1;
      }

    class (ret) <- c("wsvm", "svm")

    if (fitted) {
      ret$fitted <- na.action(predict(ret, xhold,
                                      decision.values = TRUE))
      ret$decision.values <- attr(ret$fitted, "decision.values")
      attr(ret$fitted, "decision.values") <- NULL
      if (type > 1) ret$residuals <- yorig - ret$fitted
    }

    ret
  }

predict.wsvm <-
  function (object, newdata,
            decision.values = FALSE,
            probability = FALSE,
            ...,
            na.action = na.omit)
  {
    if (missing(newdata))
      return(fitted(object))

    if (object$tot.nSV < 1)
      stop("Model is empty!")


    if(inherits(newdata, "Matrix")) {
      loadNamespace("SparseM")
      loadNamespace("Matrix")
      newdata <- as(newdata, "matrix.csr")
    }
    if(inherits(newdata, "simple_triplet_matrix")) {
      loadNamespace("SparseM")
      ind <- order(newdata$i, newdata$j)
      newdata <- new("matrix.csr",
                     ra = newdata$v[ind],
                     ja = newdata$j[ind],
                     ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))),
                     dimension = c(newdata$nrow, newdata$ncol))
    }

    sparse <- inherits(newdata, "matrix.csr")
    if (object$sparse || sparse)
      loadNamespace("SparseM")

    act <- NULL
    if ((is.vector(newdata) && is.atomic(newdata)))
      newdata <- t(t(newdata))
    if (sparse)
      newdata <- SparseM::t(SparseM::t(newdata))
    preprocessed <- !is.null(attr(newdata, "na.action"))
    rowns <- if (!is.null(rownames(newdata)))
      rownames(newdata)
    else
      1:nrow(newdata)
    if (!object$sparse) {
      if (inherits(object, "wsvm.formula")) {
        if(is.null(colnames(newdata)))
          colnames(newdata) <- colnames(object$SV)
        newdata <- model.matrix(delete.response(terms(object)),
                                as.data.frame(newdata))
        newdata <- na.action(newdata)
        act <- attr(newdata, "na.action")
      } else {
        newdata <- na.action(as.matrix(newdata))
        act <- attr(newdata, "na.action")
      }
    }

    if (!is.null(act) && !preprocessed)
      rowns <- rowns[-act]

    if (any(object$scaled))
      newdata[,object$scaled] <-
      scale(newdata[,object$scaled, drop = FALSE],
            center = object$x.scale$"scaled:center",
            scale  = object$x.scale$"scaled:scale"
      )

    if (object$kernel==4){
      newdata <- cbind(1:NROW(newdata), newdata)
    }

    if (ncol(object$SV) != ncol(newdata))
      stop ("test data does not match model !")



    ret <- .C(C_wsvmpredict,
              as.integer (decision.values),
              as.integer (probability),

              ## model
              as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
              as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
              as.integer (if (object$sparse) object$SV@ia else 0),
              as.integer (if (object$sparse) object$SV@ja else 0),
              as.double  (as.vector(object$coefs)),
              as.double  (object$rho),
              as.integer (object$compprob),
              as.double  (if (object$compprob) object$probA else 0),
              as.double  (if (object$compprob) object$probB else 0),
              as.integer (object$nclasses),
              as.integer (object$tot.nSV),
              as.integer (object$labels),
              as.integer (object$nSV),
              as.integer (object$sparse),

              ## parameter
              as.integer (object$type),
              as.integer (object$kernel),
              as.integer (object$degree),
              as.double  (object$gamma),
              as.double  (object$coef0),

              ## test matrix
              as.double  (if (sparse) newdata@ra else t(newdata)),
              as.integer (nrow(newdata)),
              as.integer (if (sparse) newdata@ia else 0),
              as.integer (if (sparse) newdata@ja else 0),
              as.integer (sparse),

              ## decision-values
              ret = double(nrow(newdata)),
              dec = double(nrow(newdata) * object$nclasses * (object$nclasses - 1) / 2),
              prob = double(nrow(newdata) * object$nclasses)


    )

    ret2 <- if (is.character(object$levels)) # classification: return factors
      factor (object$levels[ret$ret], levels = object$levels)
    else if (object$type == 2) # one-class-classification: return TRUE/FALSE
      ret$ret == 1
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
      ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
    else
      ret$ret

    names(ret2) <- rowns
    ret2 <- napredict(act, ret2)

    if (decision.values) {
      colns = c()
      for (i in 1:(object$nclasses - 1))
        for (j in (i + 1):object$nclasses)
          colns <- c(colns,
                     paste(object$levels[object$labels[i]],
                           "/", object$levels[object$labels[j]],
                           sep = ""))
      attr(ret2, "decision.values") <-
        napredict(act,
                  matrix(ret$dec, nrow = nrow(newdata), byrow = TRUE,
                         dimnames = list(rowns, colns)
                  )
        )
    }

    if (probability && object$type < 2) {
      if (!object$compprob)
        warning("Weighted SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
      else
        attr(ret2, "probabilities") <-
          napredict(act,
                    matrix(ret$prob, nrow = nrow(newdata), byrow = TRUE,
                           dimnames = list(rowns, object$levels[object$labels])
                    )
          )
    }

    ret2
  }

print.wsvm <-
  function (x, ...)
  {
    cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
    cat("Parameters:\n")
    cat("   SVM-Type: ", c("C-classification",
                           "nu-classification",
                           "one-classification",
                           "eps-regression",
                           "nu-regression")[x$type+1], "\n")
    cat(" SVM-Kernel: ", c("linear",
                           "polynomial",
                           "radial",
                           "sigmoid",
                           "precomputed")[x$kernel+1], "\n")
    if (x$type==0 || x$type==3 || x$type==4)
      cat("       cost: ", x$cost, "\n")
    if (x$kernel==1)
      cat("     degree: ", x$degree, "\n")
    if (x$type==1 || x$type==2 || x$type==3)
      cat("      gamma: ", x$gamma, "\n")
    if (x$kernel==1 || x$kernel==3)
      cat("     coef.0: ", x$coef0, "\n")
    if (x$type==1 || x$type==2 || x$type==4)
      cat("         nu: ", x$nu, "\n")
    if (x$type==3) {
      cat("    epsilon: ", x$epsilon, "\n\n")
      if (x$compprob)
        cat("Sigma: ", x$sigma, "\n\n")
    }

    cat("\nNumber of Support Vectors: ", x$tot.nSV)
    cat("\n\n")

  }

summary.wsvm <-
  function(object, ...)
    structure(object, class="summary.wsvm")

print.summary.wsvm <-
  function (x, ...)
  {
    print.wsvm(x)
    if (x$type<2) {
      cat(" (", x$nSV, ")\n\n")
      cat("\nNumber of Classes: ", x$nclasses, "\n\n")
      cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels)
    }
    cat("\n\n")
    if (x$type==2) cat("\nNumber of Classes: 1\n\n\n")

    if ("MSE" %in% names(x)) {
      cat(length (x$MSE), "-fold cross-validation on training data:\n\n", sep="")
      cat("Total Mean Squared Error:", x$tot.MSE, "\n")
      cat("Squared Correlation Coefficient:", x$scorrcoef, "\n")
      cat("Mean Squared Errors:\n", x$MSE, "\n\n")
    }
    if ("accuracies" %in% names(x)) {
      cat(length (x$accuracies), "-fold cross-validation on training data:\n\n", sep="")
      cat("Total Accuracy:", x$tot.accuracy, "\n")
      cat("Single Accuracies:\n", x$accuracies, "\n\n")
    }
    cat("\n\n")
  }

scale.data.frame <-
  function(x, center = TRUE, scale = TRUE)
  {
    i <- sapply(x, is.numeric)
    if (ncol(x[, i, drop = FALSE])) {
      x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
      if(center || !is.logical(center))
        attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
      if(scale || !is.logical(scale))
        attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
  }

plot.wsvm <-
  function(x, data, formula = NULL, fill = TRUE,
           grid = 50, slice = list(), symbolPalette = palette(),
           svSymbol = "x", dataSymbol = "o", ...)
  {
    if (x$type < 3) {
      if (is.null(formula) && ncol(data) == 3) {
        formula <- formula(delete.response(terms(x)))
        formula[2:3] <- formula[[2]][2:3]
      }
      if (is.null(formula))
        stop("missing formula.")
      if (fill) {
        sub <- model.frame(formula, data)
        xr <- seq(min(sub[, 2]), max(sub[, 2]), length = grid)
        yr <- seq(min(sub[, 1]), max(sub[, 1]), length = grid)
        l <- length(slice)
        if (l < ncol(data) - 3) {
          slnames <- names(slice)
          slice <- c(slice, rep(list(0), ncol(data) - 3 -
                                  l))
          names <- labels(delete.response(terms(x)))
          names(slice) <- c(slnames, names[!names %in%
                                             c(colnames(sub), slnames)])
        }
        for (i in names(which(sapply(data, is.factor))))
          if (!is.factor(slice[[i]])) {
            levs <- levels(data[[i]])
            lev <- if (is.character(slice[[i]])) slice[[i]] else levs[1]
            fac <- factor(lev, levels = levs)
            if (is.na(fac))
              stop(paste("Level", dQuote(lev), "could not be found in factor", sQuote(i)))
            slice[[i]] <- fac
          }

        lis <- c(list(yr), list(xr), slice)
        names(lis)[1:2] <- colnames(sub)
        new <- expand.grid(lis)[, labels(terms(x))]
        preds <- predict(x, new)
        filled.contour(xr, yr,
                       matrix(as.numeric(preds),
                              nrow = length(xr), byrow = TRUE),
                       plot.axes = {
                         axis(1)
                         axis(2)
                         colind <- as.numeric(model.response(model.frame(x, data)))
                         dat1 <- data[-x$index,]
                         dat2 <- data[x$index,]
                         coltmp1 <- symbolPalette[colind[-x$index]]
                         coltmp2 <- symbolPalette[colind[x$index]]
                         points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
                         points(formula, data = dat2, pch = svSymbol, col = coltmp2)
                       },
                       levels = 1:(length(levels(preds)) + 1),
                       key.axes = axis(4, 1:(length(levels(preds))) + 0.5,
                                       labels = levels(preds),
                                       las = 3),
                       plot.title = title(main = "Weighted SVM classification plot",
                                          xlab = names(lis)[2], ylab = names(lis)[1]),
                       ...)
      }
      else {
        plot(formula, data = data, type = "n", ...)
        colind <- as.numeric(model.response(model.frame(x,
                                                        data)))
        dat1 <- data[-x$index,]
        dat2 <- data[x$index,]
        coltmp1 <- symbolPalette[colind[-x$index]]
        coltmp2 <- symbolPalette[colind[x$index]]
        points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
        points(formula, data = dat2, pch = svSymbol, col = coltmp2)
        invisible()
      }
    }
  }

# write.wsvm <-
#     function (object, svm.file = "Rdata.svm", scale.file = "Rdata.scale",
#               yscale.file = "Rdata.yscale", weight.file = "Weight")
#     {
#
#         ret <- .C (C_wsvmwrite,
#                    ## model
#                    as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
#                    as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
#                    as.integer (if (object$sparse) object$SV@ia else 0),
#                    as.integer (if (object$sparse) object$SV@ja else 0),
#                    as.double  (as.vector(object$coefs)),
#                    as.double  (object$rho),
#                    as.integer (object$compprob),
#                    as.double  (if (object$compprob) object$probA else 0),
#                    as.double  (if (object$compprob) object$probB else 0),
#                    as.integer (object$nclasses),
#                    as.integer (object$tot.nSV),
#                    as.integer (object$labels),
#                    as.integer (object$nSV),
#                    as.integer (object$sparse),
#
#                    ## parameter
#                    as.integer (object$type),
#                    as.integer (object$kernel),
#                    as.integer (object$degree),
#                    as.double  (object$gamma),
#                    as.double  (object$coef0),
#
#                    ## filename
#                    as.character(svm.file)
#
#
#         )$ret
#
#         write.table(data.frame(center = object$x.scale$"scaled:center",
#                                scale  = object$x.scale$"scaled:scale"),
#                     file=scale.file, col.names=FALSE, row.names=FALSE)
#
#         write.table(data.frame(weight = object$weight),
#                     file=weight.file, col.names=FALSE, row.names=FALSE)
#
#         if (!is.null(object$y.scale))
#             write.table(data.frame(center = object$y.scale$"scaled:center",
#                                    scale  = object$y.scale$"scaled:scale"),
#                         file=yscale.file, col.names=FALSE, row.names=FALSE)
#     }

coef.wsvm <- function(object, ...)
{
  if (object$kernel != 0 || object$nclasses > 2)
    stop("Only implemented for regression or binary classification with linear kernel.")
  ret <- if (inherits(object$SV, "matrix.csr")) {
    loadNamespace("SparseM")
    drop(crossprod(object$coefs, SparseM::as.matrix(object$SV)))
  } else
    drop(crossprod(object$coefs, object$SV))
  trm <- object$terms
  if(!is.null(trm))
    names(ret) <- labels(trm)
  c(`(Intercept)` = -object$rho, ret)
}
