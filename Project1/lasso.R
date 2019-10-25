

lasso <- function (x, y)
{
  dimensions <- dim(x)
  rows <- dimensions[1]
  cols <- dimensions[2]
  im <- inactive <- seq(cols)
  one <- rep(1, rows)
  attrs <- dimnames(x)[[2]]

  alpha <- 2.220446e-16 # choose a small value for alpha

  meanx <- drop(one %*% x)/rows

  x <- scale(x, meanx, FALSE)
  mu <- mean(y)
  y <- drop(y - mu)

  normx <- sqrt(drop(one %*% (x^2)))
  nosignal <- normx/sqrt(rows) < alpha # remove the features which have lower value than alpha

  # add the features too low to the ignore list
  if (any(nosignal)) {
    ignores <- im[nosignal]
    inactive <- im[-ignores]
    normx[nosignal] <- alpha * sqrt(rows)
  }
  else
    ignores <- NULL

  names(normx) <- NULL

  # standardize the values in X(feature) matrix
  x <- scale(x, FALSE, normx)

  XtX <- t(x) %*% x #calculate product of transpose of X with X as it will help in further calculations
  Cvec <- drop(t(y) %*% x)
  ssy <- sum(y^2)
  residuals <- y

  max.steps <- 8 * min(cols, rows - 1)

  # initialize the coefficient matrix to all 0s
  beta <- matrix(0, max.steps + 1, cols)

  lambda = double(max.steps)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(cols)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0

  # loop through a maximum of max.steps to find the complete lasso path
  while ((k < max.steps) & (length(active) < min(cols - length(ignores),
                                                 rows - 1))) {
    action <- NULL
    C <- Cvec[inactive]

    # find the maximum correlation
    Cmax <- max(abs(C))

    # if the max correlation is less than alpha then we are done
    if (Cmax < alpha * 100) {
      break
    }

    k <- k + 1
    lambda[k] = Cmax
    if (!any(drops)) {
      new <- abs(C) >= Cmax - alpha
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {

        R <- updateR(XtX[inew, inew], R, drop(XtX[inew,
                                                      active]), Gram = TRUE, eps = alpha)

        if (attr(R, "rank") == length(active)) {
          n_rows <- seq(length(active))
          R <- R[n_rows, n_rows, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action, -inew)
        }
        else {
          if (first.in[inew] == 0)
            first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)
        }
      }
    }
    else action <- -drop.id
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL

    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A * Gi1

    if ((length(active) >= min(rows - 1, cols - length(ignores)))) {
      gamhat <- Cmax/A
    }
    else {

      # drop the inactive variables
      a <- drop(w %*% XtX[active, -c(active, ignores),
                             drop = FALSE])
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > alpha], Cmax/A)
    }
    drop.id <- NULL
    b1 <- beta[k, active]
    z1 <- -b1/w
    zmin <- min(z1[z1 > alpha], gamhat)
    if (zmin < gamhat) {
      gamhat <- zmin
      drops <- z1 == zmin
    }
    else drops <- FALSE
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w

    Cvec <- Cvec - gamhat * XtX[, active, drop = FALSE] %*% w

    Gamrat <- c(Gamrat, gamhat/(Cmax/A))
    arc.length <- c(arc.length, gamhat)
    if (any(drops)) {
      drop.id <- seq(drops)[drops]
      for (id in rev(drop.id)) {
        R <- downdateR(R, id)
      }
      drop.id <- active[drops]
      beta[k + 1, drop.id] <- 0 # set the coefficients of dropped variables to 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
    if (!is.null(attrs))
      names(action) <- attrs[abs(action)]
    actions[[k]] <- action
    inactive <- im[-c(active, ignores)]
  }
  beta <- beta[seq(k + 1), , drop = FALSE]
  lambda = lambda[seq(k)]
  dimnames(beta) <- list(paste(0:k), attrs)
  residuals <- y - x %*% t(beta)
  beta <- scale(beta, FALSE, normx)
  RSS <- apply(residuals^2, 2, sum)
  R2 <- 1 - RSS/RSS[1]
  actions = actions[seq(k)]
  netdf = sapply(actions, function(x) sum(sign(x)))
  df = cumsum(netdf)

  df = c(Intercept = 1, df + 1)
  rss.big = rev(RSS)[1]
  df.big = rows - rev(df)[1]
  if (rss.big < alpha | df.big < alpha)
    sigma2 = NaN
  else sigma2 = rss.big/df.big
  Cp <- RSS/sigma2 - rows + 2 * df
  attr(Cp, "sigma2") = sigma2
  attr(Cp, "rows") = rows
  object <- list(type = "LASSO",
                 actions = attrs[seq(k)],entry = first.in, beta = beta,
                 normx = normx, meanx = meanx)
  class(object) <- "lars"
  object
}

library(lars)
data(prostate, package = "faraway")

X <- as.matrix(prostate[1:8])
Y <- as.matrix(prostate[9])

res <- lasso(X,Y)

png("fig1.png")
plot(res)
dev.off()


