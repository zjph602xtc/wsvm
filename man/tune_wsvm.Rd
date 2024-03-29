\name{tune_wsvm}
\alias{tune_wsvm}
\alias{best.tune_wsvm}
\alias{print.tune_wsvm}
\alias{summary.tune_wsvm}
\alias{print.summary.tune_wsvm}
\title{Parameter Tuning of Functions Using Grid Search}
\description{
  This generic function tunes hyperparameters of statistical methods
  using a grid search over supplied parameter ranges.
}
\usage{
tune_wsvm(train.x, train.y = NULL, weight, use_zero_weight = FALSE,
     pre.check = TRUE, data = list(), validation.x = NULL,
     validation.y = NULL, validation.weight = NULL,
     weigthed.error = TRUE, ranges = NULL, predict.func = predict,
     tunecontrol = tune.control(), ...)
best.tune_wsvm(...)
}
\arguments{
  \item{train.x}{either a formula or a '\emph{design} matrix' of predictors.}
  \item{train.y}{the response variable if \code{train.x} is a predictor
    matrix. Ignored if \code{train.x} is a formula.}
  \item{weight}{the weight of each subject. It should be in the same length of \code{train.y}.}
  \item{use_zero_weight}{if \code{FALSE}, any subjects in the training data and the validation data (if exist) with zero (or negative) weights will be removed.}
  \item{pre.check}{if \code{TRUE}, we prefit the model with partitioned training data using the first set of parameters in \code{range}. If fails (i.e., too many zero weight subjects in the partitioned training data), we re-partition the data and re-try the model for up to 10 times. This is useful when \code{use_zero_weight=TRUE} and there many zero weights subjects in the data. }
  \item{data}{data, if a formula interface is used. Ignored, if
    predictor matrix and response are supplied directly.}
  \item{validation.x}{an optional validation set. Depending on whether a
    formula interface is used or not, the response can be
    included in \code{validation.x} or separately specified using
    \code{validation.y}. Only used for bootstrap and fixed validation
    set (see \code{\link{tune.control}})}
  \item{validation.y}{if no formula interface is used, the response of
    the (optional) validation set. Only used for bootstrap and fixed validation
    set (see \code{\link{tune.control}})}
  \item{validation.weight}{the weight of each subject in the validation set. Will be set to 1, if the user does not provide.}
  \item{weigthed.error}{if \code{TRUE}, the preformance measure will be weighted.}
  \item{ranges}{a named list of parameter vectors spanning the sampling
    space. See \code{\link{wsvm}}. The vectors will usually be created by \code{seq}.}
  \item{predict.func}{optional predict function, if the standard \code{predict}
    behavior is inadequate.}
  \item{tunecontrol}{object of class \code{"tune.control"}, as created by the
    function \code{tune.control()}. In additon, \code{tune.control$error.fun} should be a
    function that takes three arguments: (true y, predicted y, weight). If omitted, \code{tune.control()}
    gives the defaults.}
  \item{\dots}{Further parameters passed to the training functions.}
}
\value{
  For \code{tune_wsvm}, an object of class \code{tune_wsvm}, including the components:
  \item{best.parameters}{a 1 x k data frame, k number of parameters.}
  \item{best.performance}{best achieved performance.}
  \item{performances}{if requested, a data frame of all parameter
    combinations along with the corresponding performance results.}
  \item{train.ind}{list of index vectors used for splits into
    training and validation sets.}
  \item{best.model}{if requested, the model trained on the complete training data
    using the best parameter combination.}

  \code{best.tune_wsvm()} returns the best model detected by \code{tune_wsvm}.
}
\details{
  As performance measure, the classification error is used
  for classification, and the mean squared error for regression. It is
  possible to specify only one parameter combination (i.e., vectors of
  length 1) to obtain an error estimation of the specified type
  (bootstrap, cross-classification, etc.) on the given data set.

  Cross-validation randomizes the data set before building the splits
  which---once created---remain constant during the training
  process. The splits can be recovered through the \code{train.ind}
  component of the returned object.
}
\author{
  David Meyer \cr
  Modified by Tianchen Xu \email{tx2155@columbia.edu}
}
\seealso{\code{\link{tune.control}}, \code{\link{plot.tune_wsvm}}}
\examples{
data(iris)

obj <- tune_wsvm(Species~., weight = c(rep(0.8, 50),rep(1,100)),
            data = iris, ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
            tunecontrol = tune.control(sampling = "fix"))

\donttest{
set.seed(11)
obj <- tune_wsvm(Species~., weight = c(rep(0, 52),rep(1,98)),
            data = iris, use_zero_weight = TRUE,
            ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
            tunecontrol = tune.control(sampling = "bootstrap"))
}


summary(obj)
plot(obj, transform.x = log2, transform.y = log2)
plot(obj, type = "perspective", theta = 120, phi = 45)

best.tune_wsvm(Species~.,weight = c(rep(0.08, 50),rep(1,100)),
            data = iris, ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
            tunecontrol = tune.control(sampling = "fix"))
}
\keyword{models}

