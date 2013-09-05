\name{summary.larf}
\alias{summary.larf}
\alias{print.summary.larf}


\title{Local Average Response Functions}
\description{
  Methods to standard generics for instrumental-variable regressions fitted by \code{\link{larf}}.
}

\usage{
\method{summary}{larf}(object, \dots)
}

\arguments{
  \item{object}{an object of class \code{"larf"} as fitted by \code{\link{larf}}.}
  \item{\dots}{currently not used.}
}

\details{
\code{\link{larf}} is the high-level interface to the work-horse function \code{\link{larf.fit}}. A set of standard methods (including \code{summary}, \code{vcov}, \code{predict}) is available.
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

summary(est)
}
