context("Test NB regression model")

# Set-Up
set.seed(42)
n = 10
p = 3
X = matrix(rnorm(n*p), nrow=n, ncol=p)
colnames(X) <- sprintf("par:%d", 1:p)
row.names(X) <- sprintf("r:%d", 1:n)
beta = rnorm(p)
mu = exp(X %*% beta)
y = rnbinom(n, mu=mu, size=1)

test_that("Simple NB regression fit", {
  model = nbreg.fit(X=X, y=y, phi=1, verbose=TRUE)
  expect_equal(model$converged, TRUE)
})

test_that("Simple NB regression fit - no convergence", {
  expect_warning(nbreg.fit(X=X, y=y, phi=1, tol=0))
})

test_that("NB regression - invalid data", {
  y = rnorm(n, mu, 1)
  expect_error(nbreg.fit(X=X, y=y, phi=1))
})

test_that("NB regression - low rank design matrix", {
  X.2 = cbind(X, X)
  beta.2 = rnorm(2 * p)
  mu.2 = exp(X.2 %*% beta.2)
  y.2 = rnbinom(n, mu, 1)
  expect_error(nbreg.fit(X=X.2, y=y.2, phi=1))
})