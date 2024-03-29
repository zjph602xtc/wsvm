\name{plot.wsvm}
\alias{plot.wsvm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plot WSVM Objects}
\description{
  Generates a scatter plot of the input data of a \code{wsvm} fit for
  classification models by highlighting the classes and support
  vectors. Optionally, draws a filled contour plot of the class regions.
}
\usage{
\method{plot}{wsvm}(x, data, formula, fill = TRUE, grid = 50, slice = list(),
symbolPalette = palette(), svSymbol = "x", dataSymbol = "o", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{wsvm}}
  \item{data}{data to visualize. Should be the same used for fitting.}
  \item{formula}{formula selecting the visualized two dimensions. Only
    needed if more than two input variables are used.}
  \item{fill}{switch indicating whether a contour plot for the class
    regions should be added.}
  \item{grid}{granularity for the contour plot.}
  \item{slice}{a list of named values for the dimensions held
    constant (only needed if more than two variables are
    used). The defaults for unspecified dimensions are 0 (for numeric
    variables) and the first level (for factors). Factor levels can
    either be specified as factors or character vectors of length 1.}
  \item{symbolPalette}{Color palette used for the class the data points and support
    vectors belong to.}
  \item{svSymbol}{Symbol used for support vectors.}
  \item{dataSymbol}{Symbol used for data points (other than support vectors).}
  \item{\dots}{additional graphics parameters passed to
    \code{filled.contour} and \code{plot}.}
}

\value{None}

\author{David Meyer \cr
  Modified by Tianchen Xu \email{tx2155@columbia.edu}}

\seealso{\code{\link{wsvm}}}

\examples{
## a simple example
data(cats, package = "MASS")
m <- wsvm(Sex~., data = cats, weight = rep(1,144))
plot(m, cats)

## more than two variables: fix 2 dimensions
data(iris)
m2 <- wsvm(Species~., data = iris, weight = rep(1,150))
plot(m2, iris, Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

## plot with custom symbols and colors
plot(m, cats, svSymbol = 1, dataSymbol = 2, symbolPalette = rainbow(4),
color.palette = terrain.colors)
}
\keyword{neural}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line
\keyword{nonlinear}% __ONLY ONE__ keyword per line
