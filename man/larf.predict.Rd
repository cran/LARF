\name{predict.larf}
\alias{predict.larf}

\title{Local Average Response Functions}
\description{
Predict new outcomes based on a regression fitted by \code{\link{larf}}.
}
\usage{
\method{predict}{larf}(object, newCov, newTreatment, \dots)
}

\arguments{
  \item{object}{an object of class \code{larf} as fitted by \code{\link{larf}}.}
  \item{newCov}{A matrix containing the new covariates.}
  \item{newTreatment}{A vector containing the new binary treatment variable.}
  \item{\dots}{currently not used.}
}

\details{
Predicted outcomes are based on the estimated coefficients and new covariates and/or new treatment variable. The predicted outcomes  are probabilities when the outcome is binary.
}

\value{
\item{predicted.values}{The function returns a vector of the predicted outcomes given the new data.}
}

\references{
Abadie, Alberto. 2003. "Semiparametric Instrumental Variable Estimation of Treatment Response Models." \emph{Journal of Econometric} 113: 231-263.
}

\author{
Weihua An and Xuefu Wang, Indiana University Bloomington
}


\seealso{
\code{\link{larf}}, \code{\link{larf.fit}}
}

\examples{
data(c401k)

# Choose the first 500 units for a small example
c401k <- c401k[1:500,]

Y <- c401k$pira
X <- as.matrix(c401k[c("inc", "incsq", "marr",  "male", "age", "agesq",  "fsize"  )])
D <- c401k$p401k
Z <- c401k$e401k

# high-level interface
est <- larf(Y~X, D, Z)

# make prediction based on new treatment status
D <- sample(D, 500, replace = TRUE)
prediction <- predict(est, cbind(1,X), D)
}
