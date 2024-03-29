\name{predict.wsvm}
\alias{predict.wsvm}
\title{Predict Method for Subject Weighted Support Vector Machines}
\description{
  This function predicts values based upon a model trained by \code{wsvm}.
}
\usage{
\method{predict}{wsvm}(object, newdata, decision.values = FALSE,
probability = FALSE, ..., na.action = na.omit)
}

\arguments{
  \item{object}{Object of class \code{"wsvm"}, created by \code{wsvm}.}
  \item{newdata}{An object containing the new input data: either a
  matrix or a sparse matrix (object of class
    \code{\link[Matrix]{Matrix}} provided by the \pkg{Matrix} package,
    or of class \code{\link[SparseM]{matrix.csr}}
    provided by the \pkg{SparseM} package, or of class
    \code{\link[slam]{simple_triplet_matrix}} provided by the \pkg{slam}
    package). A vector will
    be transformed to a n x 1 matrix.}
  \item{decision.values}{Logical controlling whether the decision values
    of all binary classifiers computed in multiclass classification
    shall be computed and returned.}
  \item{probability}{Logical indicating whether class probabilities
    should be computed and returned. Only possible if the model was
    fitted with the \code{probability} option enabled.}
  \item{na.action}{A function to specify the action to be taken if \sQuote{NA}s are
    found. The default action is \code{na.omit}, which leads to rejection of cases
    with missing values on any required variable. An alternative
    is \code{na.fail}, which causes an error if \code{NA} cases
    are found. (NOTE: If given, this argument must be named.)}
  \item{\dots}{Currently not used.}
}
\value{
  A vector of predicted values (for classification: a vector of labels, for density
  estimation: a logical vector). If \code{decision.value} is
  \code{TRUE}, the vector gets a \code{"decision.values"} attribute
  containing a n x c matrix (n number of predicted values, c number of
  classifiers) of all c binary classifiers' decision values. There are k
  * (k - 1) / 2 classifiers (k number of classes). The colnames of
  the matrix indicate the labels of the two classes. If \code{probability} is
  \code{TRUE}, the vector gets a \code{"probabilities"} attribute
  containing a n x k matrix (n number of predicted values, k number of
  classes) of the class probabilities.
}
\note{
  If the training set was scaled by \code{wsvm} (done by default), the
  new data is scaled accordingly using scale and center of
  the training data.
}
\author{
  David Meyer (based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin) \cr
  Modified by Tianchen Xu \email{tx2155@columbia.edu}
}
\seealso{
  \code{\link{wsvm}}
}
\examples{
## load dataset
data(iris)
attach(iris)

## classification mode
# default with factor response:
model1 <- wsvm(Species ~ ., weight = rep(1,150), data = iris) # same weights
model2 <- wsvm(x = iris[,1:4], y = iris[,5],
              weight = c(rep(0.08, 50),rep(1,100))) # less weights to setosa
x <- subset(iris, select = -Species)
y <- iris$Species
model3 <- wsvm(x, y, weight = rep(10,150)) # similar to model 1, but larger weights for all subjects

# test with train data
pred <- predict(model1, iris[,1:4])
# (same as:)
pred <- fitted(model1)

# Check accuracy:
table(pred, y) # model 1, equal weights

# compute decision values and probabilities:
pred <- predict(model1, x, decision.values = TRUE)
attr(pred, "decision.values")[1:4,]


## try regression mode on two dimensions
# create data
x <- seq(0.1, 5, by = 0.05)
y <- log(x) + rnorm(x, sd = 0.2)

# estimate model and predict input values
model1 <- wsvm(x, y, weight = rep(1,99))
model2 <- wsvm(x, y,
        weight = seq(99,1,length.out = 99)) # decreasing weights

# visualize
plot(x, y)
points(x, log(x), col = 2)
points(x, fitted(model1), col = 4)
points(x, fitted(model2), col = 3) # better fit for the first few points
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}

