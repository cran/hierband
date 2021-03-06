% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/hierband.R
\name{hierband}
\alias{hierband}
\title{Solves main optimization problem for fixed lambda value}
\usage{
hierband(Sighat, lam, w = NULL, delta = NULL, maxiter = 100,
  tol = 1e-07)
}
\arguments{
\item{Sighat}{The sample covariance matrix}

\item{lam}{Non-negative penalty parameter.  Controls sparsity level.}

\item{w}{\code{(p-1)}-by-\code{(p-1)} lower-triangular matrix (above diagonal ignored).
\code{w[l,]} gives the \code{l} weights for g_l.
Defaults to \code{w[l,m]=sqrt(2 * l)/(l - m + 1)} for \code{m <= l}}

\item{delta}{Lower bound on eigenvalues.  If this is NULL (which is default), then no eigenvalue
constraint is included.}

\item{maxiter}{Number of iterations of blockwise coordinate descent to perform.}

\item{tol}{Only used when \code{delta} is non-NULL.  When no eigenvalue changes by more than
\code{tol} in BCD, convergence is assumed.}
}
\value{
Returns the convex banded estimate of covariance.
}
\description{
Solves the main optimization problem appearing in Bien, Bunea, & Xiao (2015):
\deqn{min_P || Sighat - P ||_F^2 + lam * sum_l || (W_l * P)_{g_l} ||_2}
where g_l are the outermost l(l+1) elements of a square matrix.
and || (W_l * P)_{g_l} ||^2 = sum_{m<=l} w_{lm}^2 ||P_{s_m}||^2.  If a non-NULL \code{delta} is provided,
then a constraint of the form $P >= delta I_p$ is included. Problem is solved by
performing blockwise coordinate descent on the dual problem.  See paper
for more explanation.
}
\examples{
set.seed(123)
p <- 100
n <- 50
K <- 10
true <- ma(p, K)
x <- matrix(rnorm(n*p), n, p) \%*\% true$A
Sighat <- cov(x)
fit <- hierband(Sighat, lam=0.4)
min(eigen(fit)$values)
fit2 <- hierband(Sighat, lam=0.4, delta=0.2)
min(eigen(fit2)$values)
# Use cross validation to select lambda:
path <- hierband.path(Sighat)
cv <- hierband.cv(path, x)
fit <- hierband(Sighat, lam=cv$lam.best)
}
\seealso{
\code{\link{hierband.path}} \code{\link{hierband.cv}}
}

