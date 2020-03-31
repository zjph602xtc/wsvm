\name{write.wsvm}
\alias{write.wsvm}
\title{Write WSVM Object to File}
\description{
  This function exports an WSVM object (trained by \code{wsvm}) to three
  specified files. One is in the format that the
  function 'svm_load_model' of \code{libsvm} can read, this dataset does NOT contain weight since the original \code{libsvm} does NOT support subject weighted svm. The second is for
  scaling data, containing a data with centers and scales for all variables. The last is the weights for all subjects.
}
\usage{
write.wsvm(object, svm.file = "Rdata.svm", scale.file = "Rdata.scale",
              yscale.file = "Rdata.yscale", weight.file = "Weight")
}
\arguments{
  \item{object}{Object of class \code{"wsvm"}, created by \code{wsvm}.}
  \item{svm.file}{filename to export the wsvm object to.}
  \item{scale.file}{filename to export the scaling data of the
    explanatory variables to.}
  \item{yscale.file}{filename to export the scaling data of the dependent
    variable to, if any.}
  \item{weight.file}{filename to export the weights info.}
}
\details{
  This function is NOT useful since the \code{libsvm} do not support subject weighted SVM. But maybe you can use if for the modified version \url{https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/#weights_for_data_instances}. The SVM model is saved in the standard format of
  libsvm. The scaling data are written to separate files because scaling
  data are not included in the standard format of libsvm. The format
  of the scaling data file is a n times 2 matrix: the n-th row
  corresponds to the n-th dimension of the data, the columns being formed
  of the corresponding mean and scale. If scaling information for the
  dependent variable exists (in case of regression models), it is stored
  in yet another file (1 times 2 matrix). The weight data just contains all weights.
}
\author{
  Tomomi TAKASHINA (based on 'predict.svm' by David Meyer)
  \email{t.takashina@computer.org}\cr
    Modified by Tianchen Xu \email{tx2155@columbia.edu}
}
\seealso{
  \code{\link{wsvm}}
}
\examples{
data(iris)
attach(iris)

## classification mode
# default with factor response:
model <- wsvm(Species ~ ., weight = rep(1,150), data = iris)

# export WSVM object to (temporary) files
svm_file <- tempfile()
scale_file <- tempfile()
weight_file <- tempfile()

write.wsvm(model, svm.file = svm_file, scale.file = scale_file, weight = weight_file)

# read scale file
# the n-th row is corresponding to n-th dimension. The 1st column contains the
# center value, the 2nd column is the scale value.
read.table(scale_file)

# clean up
unlink(svm_file)
unlink(scale_file)
unlink(weight_file)
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}


