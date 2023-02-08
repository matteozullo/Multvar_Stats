#' Multivariate Statistics Codebank
#'
#' This script contains a collection of functions for performing common
#' multivariate statistics calculations, such as Hotelling's T-squared test,
#' Fisher's linear discriminant, and MANOVA.
#'
#' @author Matteo Zullo
#' @date December 7, 2022
#' 


#' Vector length
#'
#' This function calculates the length of a given vector. The length of a vector is the square root of the sum of the squares of its components.
#'
#' @param v A numeric vector.
#'
#' @return A numeric value representing the length of the vector.
#'
#' @export
#'
#' @examples
#' v.length(c(1,1))
v.length <- function(v) {
  length = sqrt(sum(v^2))
  return(length)
}


#' Radians to degrees
#'
#' This function converts an angle from radians to degrees using the formula
#' `degrees = 180 * radians / pi`.
#'
#' @param rad An angle measure in radians.
#'
#' @return An angle measure in degrees
#'
#' @export
#'
#' @examples
#' rad_to_deg(1)
#'
rad_to_deg <- function(rad) {
  # Calculate the angle in degrees using the formula
  # degrees = 180 * radians / pi
  deg <- 180 * rad / pi
  
  # Return the result
  return(deg)
}


#' Calculate the angle between two vectors
#'
#' This function calculates the angle between two vectors using the dot product and the
#' inverse cosine function. The angle is returned in degrees.
#'
#' @param v1 The first vector
#' @param v2 The second vector
#' @return The angle between the two vectors, in degrees
#' @export
#' @examples
#' v1 <- c(1, 0, 0)
#' v2 <- c(0, 1, 0)
#' vecs.angle(v1, v2)
#'
vecs.angle <- function(v1, v2) {
  # Check that the vectors are the same length
  if (length(v1) != length(v2)) {
    stop("Vectors must be the same length")
  }
  
  # Calculate the dot product of the vectors
  dot_prod <- v1 %*% v2
  
  # Calculate the lengths of the vectors
  length1 <- sqrt(v1 %*% v1)
  length2 <- sqrt(v2 %*% v2)
  
  # Calculate the angle in radians
  a_rad <- acos(dot_prod / (length1 * length2))
  
  # Convert the angle to degrees and return the result
  return(180 * a_rad / pi)
}


#' Calculate the projection of one vector onto another
#'
#' This function calculates the projection of the second vector onto the first vector using
#' the formula `proj_v1(v2) = ((v1 . v2) / (v1 . v1)) * v1`.
#'
#' @param v1 The first vector, onto which the second vector will be projected
#' @param v2 The second vector, which will be projected onto the first vector
#' @return A vector representing the projection of the second vector onto the first vector
#' @export
#' @examples
#' v1 <- c(1, 0, 0)
#' v2 <- c(0, 1, 0)
#' vecs.proj(v1, v2)
#'
vecs.proj <- function(v1, v2) {
  # Check that the first vector is not the zero vector
  if (all(v1 == 0)) {
    stop("First vector cannot be the zero vector")
  }
  
  # Calculate the projection using the formula
  # proj_v1(v2) = ((v1 . v2) / (v1 . v1)) * v1
  proj <- as.vector((v1 %*% v2) / (v1 %*% v1)) * v1
  
  # Return the result
  return(proj)
}


#' Generate random samples from a bivariate normal distribution
#'
#' This function generates random samples from a bivariate normal distribution with
#' given mean and covariance matrix.
#'
#' @param mu A vector of length 2 specifying the mean of the distribution
#' @param S A 2x2 covariance matrix for the distribution
#' @param n The number of samples to generate
#' @return A matrix of size n x 2 containing the generated samples
#' @export
#' @examples
#' mu <- c(0, 0)
#' S <- matrix(c(1, 0, 0, 1), nrow = 2)
#' bivnorm(mu, S, 10)
#'
bivnorm <- function(mu, S, n) {
  # Check that the mean vector has length 2
  if (length(mu) != 2) {
    stop("Mean vector must have length 2")
  }
  
  # Check that the covariance matrix is 2x2
  if (nrow(S) != 2 || ncol(S) != 2) {
    stop("Covariance matrix must be 2x2")
  }
  
  # Generate random samples using the given parameters
  X <- rep(mu, each = n) + matrix(rnorm(n*2), ncol = 2) %*% chol(S)
  
  # Return the result
  return(X)
}


#' Calculate the probability density for a bivariate normal distribution
#'
#' This function calculates the probability density at a given point for a bivariate
#' normal distribution with given mean and covariance matrix.
#'
#' @param mu A vector of length 2 specifying the mean of the distribution
#' @param S A 2x2 covariance matrix for the distribution
#' @param X A matrix of size n x 2 containing the points at which to calculate the density
#' @return A vector of length n containing the density values at the given points
#' @export
#' @examples
#' mu <- c(0, 0)
#' S <- matrix(c(1, 0, 0, 1), nrow = 2)
#' X <- matrix(c(1, 1, 2, 2), nrow = 2)
#' bivnorm.prob(mu, S, X)
#'
bivnorm.prob <- function(mu, S, X) {
  # Check that the input dimensions are valid
  if (length(mu) != 2) {
    stop("Mean vector must have length 2")
  }
  
  if (nrow(S) != 2 || ncol(S) != 2) {
    stop("Covariance matrix must be 2x2")
  }
  
  if (ncol(X) != 2) {
    stop("Points matrix must have 2 columns")
  }
  
  # Compute the probability
  mu1 = mu[1]
  mu2 = mu[2]
  x1 = X[,1]
  x2 = X[,2]
  s12 = S[1,2]
  s11 = S[1,1]
  s22 = S[2,2]
  num = exp(-(s22*(x1-mu1)^2+s11*(x2-mu2)^2-2*s12*(x1-mu1)*(x2-mu2))/(2*(s11*s22-s12^2)))
  den = (2*pi*sqrt(s11*s22-s12^2))
  prob = num / den
  return(prob)
}


#' Calculate the square root of a matrix
#'
#' This function calculates the square root of a matrix using its eigendecomposition.
#'
#' @param M A numeric matrix
#' @param inv A logical value indicating whether to calculate the inverse square root (TRUE) or the regular square root (FALSE)
#' @return A matrix containing the square root of M
#' @export
#' @examples
#' M = matrix(c(4, 1, 1, 4), nrow = 2)
#' msqrt(M)
#' msqrt(M, inv = TRUE)
msqrt <- function(M, inv = FALSE){
  # Calculate the eigendecomposition of M
  eigen = eigen(M)
  
  # If inv is FALSE, calculate M^(1/2).
  # If inv is TRUE, calculate M^(-1/2)
  if (inv == FALSE){
    # Calculate the square root of the eigenvalues of M
    sqrt_eigen_values = sqrt(eigen$values)
    
    # Reconstruct M^(1/2) from the eigenvectors,
    # the matrix of square root eigenvalues,
    # and the transpose of the eigenvectors
    M_out = eigen$vectors %*% diag(sqrt_eigen_values) %*% t(eigen$vectors)
  } else {
    # Calculate the inverse square root of the eigenvalues of M
    inv_sqrt_eigen_values = 1 / sqrt(eigen$values)
    
    # Reconstruct M^(-1/2) from the eigenvectors,
    # the matrix of inverse square root eigenvalues,
    # and the transpose of the eigenvectors
    M_out = eigen$vectors %*% diag(inv_sqrt_eigen_values) %*% t(eigen$vectors)
  }
  
  # Return the square root of M
  return(M_out)
}


#' Perform a multivariate Hotelling's T-squared test
#'
#' This function performs a multivariate Hotelling's T-squared test to test
#' whether the means of several variables are equal to a specified value.
#'
#' @param data A data frame containing the variables to be tested.
#' @param mu The vector of means to test against.
#' @param S The covariance matrix of the variables.
#' @param n The number of observations in the data.
#' @param method A string indicating which method to use for computing the confidence interval for the means. Valid values are "small", "large", and "bonf".
#' @param alpha The significance level for the test.
#' @param vars A logical vector indicating which columns in the data to use for the test.
#' @param N The number of points to use for plotting the ellipsoid.
#' @return A list containing the results of the test, including the mean, covariance, T-squared statistic, p-value, confidence interval, ellipsoid radius and axes, and eigenvalues and eigenvectors of the covariance matrix.
#' @export
#' @examples
#' data = mtcars[, c("mpg", "cyl", "disp")]
#' mu = c(20, 6, 200)
#' S = cov(data)
#' n = nrow(data)
#' multi.hotelling(data = data, mu = mu, S = S, n = n)
multi.hotelling <- function(
    data = NULL,
    mu = NULL,
    S = NULL,
    n = NULL,
    method = "small",
    alpha = 0.05,
    vars = TRUE,
    N = 1000){
  
  # Calculate the mean and covariance if not provided
  if(is.null(S)){
    data = data[,vars]
    n = nrow(data)
    mu = as.matrix(apply(data, 2, mean))
    S = cov(data)
  }
  
  # Ensure that mu is a column vector
  if(dim(mu)[2] != 1){
    mu = t(mu)
  }
  
  # Calculate the degrees of freedom
  k = length(mu)
  n = n
  
  # Calculate the inverse of the covariance matrix
  S.inv = solve(S)
  
  # Calculate the eigendecomposition of the covariance matrix
  eigen = eigen(S)
  
  # Calculate the Hotelling T-squared statistic
  T.sq = t(mu) %*% S.inv %*% mu
  
  # Transform T-squared into F-statistic for significance test
  F.stat = (n-k)/(k*(n-1)) * T.sq
  
  # Calculate the p-value for the test
  pval = pf(F.stat, k, n-k, lower.tail = FALSE)
  
  # Store variances
  SE = sqrt(diag(S)/n)
  
  # Calculate critical statistic with relevant method
  if (method == "small") { 
    sqrt = TRUE
    correction = k*(n-1)/(n-k)
    crit.stat = qf(p=1-alpha, df1=k, df2=n-k)
    
  } else if (method == "large") {
    sqrt = TRUE
    correction = 1
    crit.stat = qchisq(1-alpha, df = k)
    
  } else if (method == "bonf") {
    sqrt = FALSE
    correction = 1
    crit.stat = qt(p = 1 - alpha/(2*k), df = n-1)
  }
  
  # Confidence interval
  SE.stat = (correction * crit.stat)^(1/(sqrt+1)) * SE
  upr = mu + SE.stat
  lwr = mu - SE.stat
  CI = cbind(lwr, upr)
  colnames(CI) = c("lwr", "upr")
  
  # Calculate ellipsoid
  radius = sqrt(correction * crit.stat / n)
  axes = sqrt(eigen$values) * radius
  
  out = list(
    "mu" = mu,
    "S" = S,
    "T_sq" = T.sq,
    "pval" = pval,
    "CI" = CI,
    "radius" = radius,
    "axes" = axes, 
    "eigenvalues" = eigen$values,
    "eigenvectors" = eigen$vectors
  )
  
  return(out)
}


#' Calculate the pooled standard deviation
#'
#' This function calculates the pooled standard deviation for two groups
#' of observations.
#'
#' @param data A data frame containing the observations for the two groups.
#' @param grpvar A string or integer indicating the column in the data frame
#'   that contains the group labels.
#' @param grplabs A numeric vector of length 2 indicating the labels for the
#'   two groups.
#' @param vars A logical vector indicating which columns in the data frame
#'   to use for the calculation.
#' @return The pooled standard deviation for the two groups.
#' @export
#' @examples
#' data = mtcars
#' grpvar = "am"
#' grplabs = c(0, 1)
#' vars = c("mpg", "cyl")
#' sd.pooled(data, grpvar, grplabs, vars)
sd.pooled <- function(
    data,
    grpvar,
    grplabs = c(0,1),
    vars = NULL){
  
  # Extract the observations for each group
  data0 = as.matrix(data[data[,grpvar] == grplabs[1],vars])
  data1 = as.matrix(data[data[,grpvar] == grplabs[2],vars])
  
  # Calculate the mean, standard deviation, and sample size for each group
  mu0 <- mean(data0)
  mu1 <- mean(data1)
  s0 <- sd(data0)
  s1 <- sd(data1)
  n0 <- length(data0)
  n1 <- length(data1)
  
  # Calculate the pooled standard deviation
  sd.pooled <- sqrt(((n0-1)*s0^2 + (n1-1)*s1^2) / (n1+n0-2))
  
  # Return the pooled standard deviation
  return(sd.pooled)
}


#' Calculate the pooled covariance matrix
#'
#' This function calculates the pooled covariance matrix for two groups
#' of observations.
#'
#' @param data A data frame containing the observations for the two groups.
#' @param grpvar A string or integer indicating the column in the data frame
#'   that contains the group labels.
#' @param grplabs A numeric vector of length 2 indicating the labels for the
#'   two groups.
#' @param vars A character vector indicating which columns in the data frame
#'   to use for the calculation.
#' @return The pooled covariance matrix for the two groups.
#' @export
#' @examples
#' data = mtcars
#' grpvar = "am"
#' grplabs = c(0, 1)
#' vars = c("mpg", "cyl")
#' sd.pooled.matrix(data, grpvar, grplabs, vars)
sd.pooled.matrix <- function(
    data,
    grpvar,
    grplabs = c(0,1),
    vars = NULL){
  
  # Convert grpvar and vars to numeric indices if provided as strings
  if (is.numeric(grpvar) == FALSE){
    grpvar = which(names(data) %in% grpvar)
  }
  
  if (is.null(vars)){
    vars = colnames(data)[-grpvar]
  }
  
  if (is.numeric(vars) == FALSE){
    vars = which(colnames(data) %in% vars)
  }
  
  # Subset the data based on the group labels
  grp0 = data[,grpvar] == grplabs[1]
  grp1 = data[,grpvar] == grplabs[2]
  data0 = as.matrix(data[grp0,vars])
  data1 = as.matrix(data[grp1,vars])
  
  # Calculate the sample size, covariance matrix, and pooled covariance matrix for each group
  n0 = nrow(data0)
  n1 = nrow(data1)
  S0 = cov(data0)
  S1 = cov(data1)
  S.pooled = ((n0-1)*S0 + (n1-1)*S1)/((n0-1)+(n1-1))
  
  return(S.pooled)
}


#' Multivariate analysis of variance
#'
#' This function implements a two-way multivariate analysis of variance (MANOVA), which is a statistical technique used to test whether there are significant differences between two or more groups on multiple continuous response variables.
#'
#' @param data A data frame containing the data to be used for the analysis.
#' @param grpvar A variable indicating the group the data belongs to. Can be a numeric vector or a character vector containing the column names of the matrix data.
#' @param grplabs A list of labels for the groups.
#' @param vars A list of variables to be used for the analysis. Can be a numeric vector or a character vector containing the column names of the matrix data.
#'
#' @return A list containing the sum of squares for each population and the cross-products.
#'
#' @export
#'
#' @examples
#' mav(data, grpvar = 1, grplabs = c("A","B","C), vars = c(2,3,4))
mav <- function(data, grpvar, grplabs, vars){
  
  # no. of variables
  k = length(vars)
  
  # reshape input matrices
  n0 = nrow(data[data[,grpvar] == grplabs[1],])
  n1 = nrow(data[data[,grpvar] == grplabs[2],])
  n2 = nrow(data[data[,grpvar] == grplabs[3],])
  n.max = max(n0,n1,n2)
  
  pop0 = as.matrix(cbind(t(data[data[,grpvar] == grplabs[1],vars]),
                         matrix(, nrow = 2, ncol = n.max - n0)))
  pop1 = as.matrix(cbind(t(data[data[,grpvar] == grplabs[2],vars]),
                         matrix(, nrow = 2, ncol = n.max - n1)))
  pop2 = as.matrix(cbind(t(data[data[,grpvar] == grplabs[3],vars]),
                         matrix(, nrow = 2, ncol = n.max - n2)))
  
  # population (1)
  OBS1 = rbind(
    pop0[1,],
    pop1[1,],
    pop2[1,]
  )
  rownames(OBS1) = NULL
  colnames(OBS1) = NULL
  
  grand_avg1 = mean(OBS1, na.rm = TRUE)
  treat1 = rowMeans(OBS1, na.rm = TRUE) - grand_avg1
  
  # manova (1)
  RES1 = OBS1 - (grand_avg1 + treat1)
  AVG1 = OBS1
  TREAT1 = OBS1
  AVG1[,] <- ifelse(is.na(OBS1),NA,grand_avg1)
  TREAT1[,] <- ifelse(is.na(OBS1),NA,treat1)
  
  # sum of squares (1)
  SS_obs1 = sum(OBS1^2, na.rm = TRUE)
  SS_avg1 = sum(AVG1^2, na.rm = TRUE)
  SS_treat1 = sum(TREAT1^2, na.rm = TRUE)
  SS_res1 = sum(RES1^2, na.rm = TRUE)
  SS_tot1 = SS_obs1 - SS_avg1
  
  # population (2)
  OBS2 = rbind(
    pop0[2,],
    pop1[2,],
    pop2[2,]
  )
  rownames(OBS2) = NULL
  colnames(OBS2) = NULL
  
  grand_avg2 = mean(OBS2, na.rm = TRUE)
  treat2 = rowMeans(OBS2, na.rm = TRUE) - grand_avg2
  RES2 = OBS2 - (grand_avg2 + treat2)
  
  # manova (2)
  AVG2 = OBS2
  TREAT2 = OBS2
  AVG2[,] <- ifelse(is.na(OBS2),NA,grand_avg2)
  TREAT2[,] <- ifelse(is.na(OBS2),NA,treat2)
  
  # sum of squares (2)
  SS_obs2 = sum(OBS2^2, na.rm = TRUE)
  SS_avg2 = sum(AVG2^2, na.rm = TRUE)
  SS_treat2 = sum(TREAT2^2, na.rm = TRUE)
  SS_res2 = sum(RES2^2, na.rm = TRUE)
  SS_tot2 = SS_obs2 - SS_avg2
  
  # cross-products
  SS_obs = sum(OBS1 * OBS2, na.rm = TRUE)
  SS_avg = sum(AVG1 * AVG2, na.rm = TRUE)
  SS_treat = sum(TREAT1 * TREAT2, na.rm = TRUE)
  SS_res = sum(RES1 * RES2, na.rm = TRUE)
  SS_tot = SS_obs - SS_avg
  
  # MANOVA decomposition
  TREAT = cbind(c(SS_treat1,SS_treat),c(SS_treat,SS_treat2))
  RES = cbind(c(SS_res1,SS_res),c(SS_res,SS_res2))
  TOTAL = cbind(c(SS_tot1,SS_tot),c(SS_tot,SS_tot2))
  
  n = sum(!is.na(OBS1)) # no. of samples for each population
  df_treat = k - 1
  df_res = n - k
  df_total = n -1
  
  # Wilks lambda
  lambda = det(RES)/det(TOTAL)
  F_stat = (1- sqrt(lambda))/sqrt(lambda) * (n-k-1)/(k-1)
  pval = pf(F_stat, df1 = 2*(k-1) , df2 = 2*(n-k-1), lower.tail = FALSE)
  
  out = list(
    "manova" = list("treat" = TREAT, "residual" = RES, "total" = TOTAL),
    "wilks_stat" = lambda,
    "pval" = pval
  )
  return(out)
}


#' Multivariate t-test
#'
#' This function implements a multivariate t-test.
#'
#' @param data A data frame containing the data to be used for the analysis.
#' @param grpvar A variable indicating the group the data belongs to. Can be a numeric vector or a character vector containing the column names of the matrix data.
#' @param grplabs A list of labels for the groups.
#' @param vars A list of variables to be used for the analysis. Can be a numeric vector or a character vector containing the column names of the matrix data.
#' @param variance A string holding the variance matrix to be implemented, "pooled" or "unpooled".
#' @param method A string holding the method to be implemented, simultaneous or Bonferroni.
#' @param alpha The significance level for the test.
#'
#' @return A list containing the sum of squares for each population and the cross-products.
#'
#' @export
#'
#' @examples
#' multi.ttest(data = matrix(1:6, ncol = 2), grpvar = 1, grplabs = c("A","B"), vars = c(2,3))

multi.ttest <- function(
    data,
    grpvar,
    grplabs,
    vars = NULL,
    variance = "pooled",
    method = "simul",
    alpha = 0.05){
  
  # if group variable is not numeric, convert
  if (is.numeric(grpvar) == FALSE){
    grpvar = which(names(data) %in% grpvar)
  }
  
  # if variables are not defined, pull all other than group indicator
  if (is.null(vars) == TRUE){
    vars = names(data[,-grpvar])
  }
  
  # if variables are not numeric, convert
  if (is.numeric(vars) == FALSE){
    vars = which(colnames(data) %in% vars)
  }
  
  k = length(vars)  # number of variables
  
  # extract data
  grp0 = data[,grpvar] == grplabs[1]
  grp1 = data[,grpvar] == grplabs[2]
  data.grp = data[, -grpvar]
  data0 = as.matrix(data.grp[grp0,vars])
  data1 = as.matrix(data.grp[grp1,vars])
  
  # mean, sd, and sample size 
  mu0 <- apply(data0, 2, mean)
  mu1 <- apply(data1, 2, mean)
  mu_diff = colSums(rbind(mu0,(-1)*mu1))
  n0 <- nrow(data0)
  n1 <- nrow(data1)
  n = n1 + n0
  
  # Use either pooled or unpooled variances
  if (variance == "pooled") {
    
    S = sd.pooled.matrix(data, grpvar, grplabs, vars)
    SE = sqrt((1/n0 + 1/n1)*diag(S))
    
    # Hotelling's T
    T.sq = t(mu0-mu1) %*% solve(S*(1/n0 + 1/n1)) %*% (mu0-mu1)
    
    # Transform Hotelling into F stat
    F.stat = (n0+n1-k-1)/(k*(n0+n1-2)) * T.sq
    pval = pf(F.stat, k, n0+n1-k-1, lower.tail = FALSE)
    
    sqrt = TRUE
    correction = k*(n0+n1-2)/(n0+n1-k-1) 
    crit.stat = qf(1-alpha, k, n0+n1-k-1)
    
  } else if (variance == "unpooled") {
    
    S0 = cov(data0)
    S1 = cov(data1)
    SE = sqrt(diag(S0)*(1/n0) + diag(S1)*(1/n1))
    
    # Hotelling's T
    T.sq = t(mu0-mu1) %*% solve(S0*(1/n0) + S1*(1/n1)) %*% (mu0-mu1)
    pval = 1 - pchisq(T.sq, df = k)
    
    # Transform into chi-2
    sqrt = TRUE
    correction = 1
    crit.stat = qchisq(p=1-alpha, df=k)
    
  }
  
  # Calculate critical statistics using relevant method
  if (method == "simul"){
    sqrt = TRUE
    
  } else if (method == "bonf") {
    sqrt = FALSE
    correction = 1
    crit.stat = qt(p = 1-alpha/(2*k), df = n0+n1-2)
  }
  
  # Calculate confidence intervals
  SE.stat = (correction * crit.stat)^(1/(sqrt+1)) * SE
  upr = mu_diff + SE.stat
  lwr = mu_diff - SE.stat
  CI = rbind(upr, lwr)
  rownames(CI) = c("upr", "lwr")
  
  out = list(
    "mu0" = mu0,
    "mu1" = mu1,
    "mu_diff" = mu_diff,
    "T_sq" = T.sq,
    "pval" = pval,
    "CI" = CI
  )
  
  return(out)
}


#' Fisher's linear discriminant
#'
#' This function calculates Fisher's linear discriminant for the given data.
#'
#' @param data A matrix containing the data to be used for calculating the discriminant.
#' @param grpvar The variable indicating the group the data belongs to. Can be a numeric vector or a character vector containing the column names of the matrix data.
#' @param vars The variables to be used for calculating the discriminant. Can be a numeric vector or a character vector containing the column names of the matrix data. If no variables are provided, all columns in data except the group variable will be used.
#' @param p0 The prior probability of the first group.
#' @param p1 The prior probability of the second group.
#' @param plot A vector containing the indices of the variables to be used for plotting the discriminant line. If no variables are provided, no plot will be created.
#'
#' @return A list containing the following elements:
#' \itemize{
#' \item mu: A matrix containing the means of the two groups for each variable.
#' \item classification: A data frame containing the original data, the group the data belongs to, and the predicted group.
#' \item beta: A vector containing the un-scaled discriminant coefficients.
#' \item beta.scaled: A vector containing the scaled discriminant coefficients.
#' \item APER: The average percentage error rate of the classification.
#' \item threshold: The discriminant threshold.
#' \item lineplot: A list containing the intercept and slope of the line plotted using the specified variables. If no variables were provided, this element will not be included in the output.
#' }
#'
#' @export
#'
#' @examples
#' ld(data = matrix(1:6, ncol = 2), grpvar = c(1,2), p0 = 0.5, p1 = 0.5, plot = c("x","y))
ld <- function(
    data,
    grpvar,
    vars = NULL,
    p0 = 0.5,
    p1 = 0.5,
    plot_xy = NULL){
  
  # if variables are not numeric, convert
  if (is.numeric(grpvar) == FALSE){
    grpvar = which(names(data) %in% grpvar)
  }
  
  if (is.null(vars)){
    vars = colnames(data)[-grpvar]
  }
  
  if (is.numeric(vars) == FALSE){
    vars = which(colnames(data) %in% vars)
  }
  
  # send group to binary values
  data[,grpvar] = as.factor(data[,grpvar])
  flevels = levels(data[,grpvar])
  data[,grpvar] = as.numeric(data[,grpvar]) - 1
  
  # means
  data0 = data[data[,grpvar] == 0, vars]
  data1 = data[data[,grpvar] == 1, vars]
  mu0 = apply(data0,2,mean)
  mu1 = apply(data1,2,mean)
  mu = cbind(mu0, mu1)
  
  # regression matrices
  X = as.matrix(data[, vars])
  Y = as.matrix(data[, grpvar])
  varnames = colnames(data)[vars]
  grpname = colnames(data)[grpvar]
  
  # covariance matrix
  S.pooled = sd.pooled.matrix(data, grpvar, grplabs = c(0,1), vars)
  
  # discriminant coefficients
  B = solve(S.pooled) %*% (mu1 - mu0)
  B.scaled = B / drop(sqrt(t(B) %*% S.pooled %*% B))
  
  # threshold
  D = (as.numeric(t(B) %*% (mu1 + mu0)/2) - log(p1/p0))
  
  # classifications
  # higher prior of population 1 (p0), lower the threshold (right-hand side)
  # therefore, it increases the probability of a positive classification (1)
  classifications = as.integer(X %*% B - D >= 0)
  
  # error rate
  APER = 1 - sum(classifications == Y) / length(Y)
  
  # prepare output data
  outdata               = cbind(X, Y, classifications)
  outdata               = as.data.frame(outdata)
  colnames(outdata)     = c(varnames, grpname, "class")
  
  outdata[,c(grpname)] <- ifelse(outdata[,grpname] == 0, flevels[1], flevels[2])
  outdata[,c("class")] <- ifelse(outdata[,"class"] == 0, flevels[1], flevels[2])
  
  # plot
  lineplot = NULL
  if (!is.null(plot_xy)){
    
    # define line
    x = sum(mu[plot_xy[1],]) / 2
    y = sum(mu[plot_xy[2],]) / 2
    b = - B.scaled[plot_xy[1],] /B.scaled[plot_xy[2],]
    m = y - b*x
    lineplot = c(m, b)
    names(lineplot) = c("intercept", "slope")
    
  }
  out <- list(
    "mu" = mu,
    "classification" = outdata,
    "beta" = B,
    "beta.scaled" = B.scaled,
    "APER" = APER,
    "threshold" = D,
    "lineplot" = lineplot
  )
  return(out)
}


#' Function for performing Leave-One-Out Cross Validation (LOOCV) on a Linear Discriminant (LD) model
#'
#' @param data: a data frame containing the predictor and dependent variables
#' @param grpvar: a character string specifying the name of the dependent variable in the data
#' @param p0: a numeric value specifying the prior probability for group 0 (default is 0.5)
#' @param p1: a numeric value specifying the prior probability for group 1 (default is 0.5)
#'
#' @return a data frame containing the original data, the predicted classifications, and the true class labels for each sample
#'
#' @example ld.lach(data, grpvar = "class")
ld.lach <- function(
    data,
    grpvar,
    vars = NULL,
    p0 = 0.5,
    p1 = 0.5){
  
  # if variables are not numeric, convert
  if (is.numeric(grpvar) == FALSE){
    grpvar = which(names(data) %in% grpvar)
  }
  
  if (is.null(vars)){
    vars = colnames(data)[-grpvar]
  }
  
  if (is.numeric(vars) == FALSE){
    vars = which(colnames(data) %in% vars)
  }
  
  # send group to binary values
  data[,grpvar] = as.factor(data[,grpvar])
  flevels = levels(data[,grpvar])
  data[,grpvar] = as.numeric(data[,grpvar]) - 1
  
  # Extract the column names of the independent variables in the data set
  varnames = colnames(data)[vars]
  grpname = colnames(data)[grpvar]
  
  # Initialize an empty data frame to store the results of the cross validation
  outdata <- NULL
  
  # Loop over each row in the data set
  for (i in 1:nrow(data)){
    
    # Split the data into a training set and a classification set
    data.train = data[-i, ]
    data.pred.X = as.matrix(data[i, vars])
    data.pred.Y = as.matrix(data[i, grpvar])
    
    # Convert the columns in the data frame to numeric type
    data.pred.X <- as.numeric(data.pred.X)
    data.pred.Y <- as.numeric(data.pred.Y)
    
    # Train a LD model on the training set
    LD <- ld(data.train, grpvar, vars, p0 = 0.5, p1 = 0.5)
    
    # Use the trained model to classify the samples in the classification set
    data.class = as.integer(data.pred.X %*% LD$beta - LD$threshold >= log(p1/p0))
    data.pred = cbind(t(data.pred.X), data.pred.Y, data.class)
    
    # Append the classification results to the output data frame
    outdata <- rbind(outdata, data.pred)
  }
  
  # Return the output data frame, which contains the original data, the predicted classifications, and the true class labels for each sample
  outdata = as.data.frame(outdata)
  colnames(outdata) = c(varnames,grpname,"class")
  
  outdata[,c(grpname)] <- ifelse(outdata[,grpname] == 0, flevels[1], flevels[2])
  outdata[,c("class")] <- ifelse(outdata[,"class"] == 0, flevels[1], flevels[2])
  
  return(outdata)
}


#' Calculate canonical correlation
#'
#' This function calculates the canonical correlation of the given data.
#'
#' @param data A matrix containing the data to be used for calculating the canonical correlation.
#' @param p The number of columns in the matrix data.
#' @param n The number of rows in the matrix data. This parameter is optional and is only used if the correlation matrix is not provided.
#' @param corr A boolean value indicating whether the correlation matrix is provided. If corr is FALSE (the default value), the correlation matrix will be calculated from the data. If corr is TRUE, the provided matrix will be used as the correlation matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#' \item loadingsA: A matrix containing the loadings of the first canonical variate.
#' \item loadingsB: A matrix containing the loadings of the second canonical variate.
#' \item correlation: A matrix containing the canonical correlations.
#' \item chi2: A matrix containing the chi-squared values and corresponding p-values for each canonical variate.
#' \item communalities: A matrix containing the variance explained by each canonical variate.
#' }
#'
#' @export
#'
#' @examples
#' cca(data = matrix(1:6, ncol = 2), p = 2)
cca <- function(
    data,
    p,
    n = NULL,
    corr = FALSE){
  
  # if correlation is not supplied, obtain correlation matrix
  if (corr == FALSE){
    R = cor(data)
    n = nrow(data)
  } else {
    R = data
  }
  
  q = nrow(R) - p
  d = min(p,q)
  
  # partial matrices
  R11 = R[1:p,1:p]
  R12 = R[1:p,(p+1):nrow(R)]
  R21 = R[(p+1):nrow(R),1:p]
  R22 = R[(p+1):nrow(R),(p+1):nrow(R)]
  R11_inv_sqrt = msqrt(R11, inv = TRUE)
  R22_inv_sqrt = msqrt(R22, inv = TRUE)
  R22_inv = solve(R22)
  
  # rho matrix
  rho = R11_inv_sqrt %*% R12 %*% R22_inv %*% R21 %*% R11_inv_sqrt
  
  # calculate scores
  A = NULL
  B = NULL
  CORR = NULL
  CHI2 = NULL
  communalities = NULL
  
  for (i in 1:p){
    
    e = eigen(rho)$vectors[,i]
    
    # canonical variate A
    a = R11_inv_sqrt %*% e
    A = cbind(A,a)
    colnames(A)[i] = paste0("V", i)
    
    # canonical variate B
    b = R22_inv_sqrt %*% R22_inv_sqrt %*% R21 %*% R11_inv_sqrt %*% e
    b_std = b %*% (1/sqrt(t(b) %*% R22 %*% b))
    B = cbind(B, b_std)
    colnames(B)[i] = paste0("V", i)
    
    # canonical correlation
    cov_ab = t(a) %*% R12 %*% b
    var_a = t(a) %*% R11 %*% a
    var_b = t(b) %*% R22 %*% b
    corr_AB = cov_ab/(sqrt(var_a)*sqrt(var_b))
    CORR = rbind(CORR, corr_AB)
    rownames(CORR)[i] = paste0("V", i)
    
    # significance
    chi2 = -(n-1-0.5*(p+q+1))*log(prod(1-eigen(rho)$values[i:d]))
    df = (p-(i-1))*(q-(i-1))
    pval = 1-pchisq(chi2, df)
    CHI2 = cbind(CHI2,rbind(chi2, pval))
    colnames(CHI2)[i] = paste0("V", i)
    
    # variance explained
    communality_a = sum((t(a) %*% R11)^2)/p
    communality_b = sum((t(b_std) %*% R22)^2)/q
    communalities = cbind(communalities, rbind(communality_a, communality_b))
    colnames(communalities)[i] = paste0("V", i)
  }
  
  out = list(
    "loadingsA" = A,
    "loadingsB" = B,
    "correlation" = t(CORR),
    "chi2" = CHI2,
    "communalities" = communalities
  )
  
  return(out)
}


#' Function calculating euclidean or Jaccard distances
#'
#' @param data: a data frame containing the predictor and dependent variables
#' @param type: a character string specifying the distance calculation, "euclidean" or "jaccard"
#' @param std: a boolean indicating whether the data has to be standardized
#' @param lab.col: a number or string specifying the column holding row names (default is that the data has legacy row names)
#' @param wgt_true: an integer indicating the weight on true comparisons (for Jaccard, default is 1)
#' @param wgt_false: an integer indicating the weight on false comparisons (for Jaccard, default is 1)
#'
#' @return a data frame containing the original data, the predicted classifications, and the true class labels for each sample
#'
#' @example distx(data, type = "euclidean", std = TRUE, lab.col = 0)
distx <- function(
    data,
    type = "euclidean",
    std = F,
    lab.col = 0,
    wgt_true = 1,
    wgt_false = 1){
  
  # dimension of distance matrix
  n = nrow(data) 
  
  # save column and row labels
  # column: 1,..,n
  # row: 2,...,n+1
  if (type == "jaccard"){
    if(lab.col == 0){
      collabels = rownames(data)[1:(n-1)]
      rowlabels = rownames(data)[2:n]
    } else{
      collabels = data[1:(n-1),lab.col]
      rowlabels = data[2:n,lab.col]
      data = data[,-c(lab.col)]
    }
  }
  
  # Save labels for euclidean distance matrix
  if (type == "euclidean"){
    if(lab.col == 0){
      collabels = rownames(data)
      rowlabels = rownames(data)
    } else{
      collabels = data[1:n,lab.col]
      rowlabels = data[1:n,lab.col]
      data = data[,-c(lab.col)]
    }
  }
  
  # Scale if necessary
  if (std == T){
    D = scale(data, scale = T, center = T)
  } else{
    D = as.matrix(data)
  }
  
  # Calculate Jaccard distances
  if (type == "jaccard"){
    out <- matrix(0,n-1,n-1)
    for(r in 1:(n-1)) {
      for(c in 1:(n-1)) {
        out[r,c] = (wgt_true*sum(D[r + 1,] == D[c,]))/(wgt_true*sum(D[r + 1,] == D[c,]) + wgt_false*sum(D[r + 1,] != D[c,]))
      }
    }
  }
  
  # Calculate Euclidean distances
  if (type == "euclidean"){
    out <- matrix(0,n,n)
    for(r in 1:n) {
      for(c in 1:n) {
        out[r,c] = sqrt(sum((D[r,] - D[c,])^2))
      }
    }
  }
  
  # label
  colnames(out) = collabels
  rownames(out) = rowlabels
  
  return(out)
}