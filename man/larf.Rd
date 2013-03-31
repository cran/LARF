\name{larf}
\alias{larf}

\title{Local Average Response Functions}
\description{
Provides instrumental variable estimation of treatment effects for compliers when both the treatment variable and instrumental variable are binary.
}
\usage{
larf(Y, D, X, Z, outcome = "binary", intercept = TRUE, method = "NLS", discrete = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A vector containing the outcome variable. It can be either binary or continuous.}
  \item{D}{A vector containing the binary treatment variable.}
  \item{X}{A matrix containing the covariates.}
  \item{Z}{A vector containing the binary instrumental variable, for example, the treatment assignment.}
  \item{outcome}{A character string for the type of outcomes. The default is "binary". The other is "continuous".}
  \item{intercept}{A logical flag for whether an intercept should be included in the model. The default is TRUE.}
  \item{method}{A character string for the estimation method to be used for binary outcome models. The default is "NLS", standing for nonlinear least squares. Maximum likelihood method (MLE) may be added in future versions.}
  \item{discrete}{A logical flag for whether marginal effects should be reported. The default is TRUE. This option is valid only for binary outcomes.}
}
\details{
The function provides instrumental variable estimation of treatment effects when both the treatment variable and instrumental variable are binary. First, pseudo weights are constructed from a probit regression of the instrumental variable on covariates. Then, the weights are used in the regression of Y on the treatment and covariates. When the outcome is binary, weighted nonlinear least squares method with a probit link is used. When the outcome is continuous, weighted least squares is used. More details of the methods can be found in Abadie (2003).
}

\value{
\item{b}{Estimated coefficients for treatment and covariates}
\item{se}{Standard errors of the estimated coefficients}
\item{db}{Estimated marginal effects of treatment and covariates, available only for binary outcomes}
\item{dse}{Standard errors of the estimated marginal effects, available only for binary outcomes}
}
\references{
Abadie, Alberto. 2003. "Semiparametric Instrumental Variable Estimation of Treatment Response Models. \emph{Journal of Econometric} 113: 231-263.
}
\author{
Weihua An and Xuefu Wang, Indiana University Bloomington
}

\note{
We derived part of our codes from the Matlab codes shared by Professor Alberto Abadie, available at \url{http://www.hks.harvard.edu/fs/aabadie/larf.html}.
}

\seealso{
\code{\link{c401k}}
}

\examples{
data(c401k)

# Choose the first 100 units for a small example
c401k <- c401k[1:100,]

D <- c401k$p401k
Z <- c401k$e401k
X <- c401k[c("inc", "incsq", "marr",  "male", "age", "agesq",  "fsize"	)]

Y <- c401k$nettfa
estimate <- larf(Y, D, X, Z, outcome = "continuous")
}
