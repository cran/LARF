\name{larf.fit}
\alias{larf.fit}

\title{Local Average Response Functions}
\description{
Provides instrumental variable estimation of treatment effects for compliers when both the treatment variable and instrumental variable are binary.
}
\usage{
larf.fit(Y, X, D, Z, outcome, method, discrete)
}

\arguments{
  \item{Y}{A vector containing the outcome variable. It can be either binary or continuous.}
  \item{X}{A matrix containing the covariates.}
  \item{D}{A vector containing the binary treatment variable.}
  \item{Z}{A vector containing the binary instrumental variable, for example, the treatment assignment.}
  \item{outcome}{A character string for the type of outcomes. The default is "binary". The other option is "continuous".}
  \item{method}{A character string for the estimation method to be used for binary outcome models. The default is "NLS", standing for nonlinear least squares. Maximum likelihood method (MLE) may be added in future versions.}
  \item{discrete}{A logical flag for whether marginal effects should be reported. The default is TRUE. This option is valid only for binary outcomes.}
}
\details{
\code{\link{larf.fit}} is the work-horse function for its high-level interface \code{larf}. A set of standard methods (including \code{print}, \code{summary}, \code{coef}, \code{vcov}, \code{fitted}, \code{resid}, \code{predict}) can be used to extract the corresponding information from a \code{\link{larf}} object.

The function provides instrumental variable estimation of treatment effects when both the treatment variable and instrumental variable are binary. First, pseudo weights are constructed from a probit regression of the instrumental variable on covariates. Then, the weights are used in the regression of the outcome on the treatment and covariates. When the outcome is binary, weighted nonlinear least squares method with a probit link is used. When the outcome is continuous, weighted least squares is used. More details of the methods can be found in Abadie (2003).
}

\value{
\item{coefficients}{Estimated coefficients for treatment and covariates.}
\item{StdErr}{Standard errors of the estimated coefficients.}
\item{MargEff}{Estimated marginal effects of treatment and covariates, available only for binary outcomes.}
\item{MargStdErr}{Standard errors of the estimated marginal effects, available only for binary outcomes.}
\item{vcov}{Variance covariance matrix between the parameters with the first entry being the treatment variable.}
\item{fitted.values}{The predicted outcomes based on the estimated coefficients and the observed covariates. They are probabilities when the outcome is binary.}
\item{residuals}{The differnece between the observed outcomes and the fitted values.}
}

\references{
Abadie, Alberto. 2003. "Semiparametric Instrumental Variable Estimation of Treatment Response Models." \emph{Journal of Econometric} 113: 231-263.
}
\author{
Weihua An and Xuefu Wang, Indiana University Bloomington
}

\note{
We derived part of the code from the Matlab code written by Professor Alberto Abadie, available at \url{http://www.hks.harvard.edu/fs/aabadie/larf.html}.
}

\seealso{
\code{\link{larf}}, \code{\link{c401k}}
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

# low-level interface
est2 <- larf.fit(Y, X, D, Z)
}