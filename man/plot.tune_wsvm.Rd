\name{plot.tune_wsvm}
\alias{plot.tune_wsvm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plot Tuning Object}
\description{
  Visualizes the results of parameter tuning.
}
\usage{
\method{plot}{tune_wsvm}(x, type = c("contour", "perspective"), theta = 60,
          col = "lightblue", main = NULL, xlab = NULL, ylab = NULL,
          swapxy = FALSE, transform.x = NULL, transform.y = NULL,
          transform.z = NULL, color.palette = hsv_palette(),
          nlevels = 20, ...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{tune_wsvm}}
  \item{type}{choose whether a contour plot or a perspective plot is
    used if two parameters are to be visualized. Ignored if only one
    parameter has been tuned.}
  \item{theta}{angle of azimuthal direction.}
  \item{col}{the color(s) of the surface facets.  Transparent colors are
          ignored.}
  \item{main}{main title.}
  \item{xlab, ylab}{titles for the axes.  N.B. These must be character
          strings; expressions are not accepted.  Numbers will be
          coerced to character strings.}
  \item{swapxy}{if \code{TRUE}, the parameter axes are swaped (only used
    in case of two parameters).}
  \item{transform.x, transform.y, transform.z}{functions to transform
    the parameters (\code{x} and \code{y}) and the error measures
    (\code{z}). Ignored if \code{NULL}.}
  \item{color.palette}{color palette used in contour plot.}
  \item{nlevels}{number of levels used in contour plot.}
  \item{\dots}{Further graphics parameters.}
}
\value{None}
\author{
  David Meyer (based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin)\cr
  Modified by Tianchen Xu \email{tx2155@columbia.edu}
}
\seealso{\code{\link{tune_wsvm}}}

\examples{
data(iris)

obj <- tune_wsvm(Species~., weight = c(rep(0.8, 50),rep(1,100)),
            data = iris, ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)),
            tunecontrol = tune.control(sampling = "fix"))

summary(obj)
plot(obj, transform.x = log2, transform.y = log2)
plot(obj, type = "perspective", theta = 120, phi = 45)
}
\keyword{models}% __ONLY ONE__ keyword per line
