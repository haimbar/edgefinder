
#' Detect edges in co-expression datasets.
#'
#' Fit the L2N model to normalized correlation coefficients between pairs of
#' genes. The mixture model has three component - the null component follows
#' a normal distribution, and the two non-null components follow lognormal
#' distributions. An edge is in the graph if the correlation between the two
#' end-point genes is large enough and determined to be in one of the non-null
#' components.
#' @param Exprs A numeric matrix with normalized gene expression data. Rows
#' correspond to genes, and columns correspond to samples.
#' @param BHthr the Benjamini-Hochberg fasle discovery rate threshold to be
#' used to determine which pairs are strongly correlated. Default=0.01.
#' @param rndseed The random seed used to select a subset of the pairs.
#' @param maxLen The maximum number of pairs that will be randomly selected
#' to fit the L2N model. Default=20000.
#' @param LOvals the maximum log-odds ratio to be used to determine the
#' cut-off points to declare which correlations are significant.
#' The program will check which log-odds ratio (1,2,...,LOvals) results in
#' FDR less than or equal to the user-specified BHthr. Default=30.
#' @param ttl Title for the fitted-model plot. Default=""
#' @param trim Fraction of extreme values to exclude from the fitted-model
#' plot. Default=0 (show all the data).
#' @param verbose Whether to show progress message to the user. Default=FALSE.
#' @param plot.it Whether to show the fitted mixture plot to the user. Default=FALSE.
#' @return A list with the following elements
#' \itemize{
#' \item{G} {The total number of genes.}
#' \item{p1} {The proportion of genes in the right mixture component (positively correlated.)}
#' \item{p2} {The proportion of genes in the left mixture component (negtively correlated.)}
#' \item{p0} {The proportion of genes in the null component (un-correlated.)}
#' \item{m0, m1, m2, s0, s1, s2} {The location and scale parameters of the three mixture components.}
#' \item {thrtable} {A table with 6 columns: posterior probability ratio (ppr) between the non-null components and the null component), the right component cutoff corresponding to the ppr, the left component cutoff, the estimated probability of Type-I errors, the estimated power, the estimated FDR.}
#' \item {LogOddsRatio} {The log-odds ratio that yields FDR closest to the desired level.}
#' \item {fitted} {The fitted model (as returned by the EM function).}
#' \item {rmse} {The root mean-squared error of the fitted model.}
#' \item {rt, lt} {The significant edges (from the right and left mixture components.)}
#' \item {AdjMat} {The (sparse) adjacency matrix with edges corresponding to rt, lt.}
#' }
#' @export
#' @import stats
#' @importFrom Matrix Matrix
#' @importFrom grDevices rgb
#' @importFrom grDevices col2rgb
#' @importFrom grDevices colours
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#' }
edgefinder <- function(Exprs, BHthr = 0.01, rndseed=112211,
                       maxLen=20000, LOvals=30, ttl="",trim=0,
                       verbose=FALSE, plot.it=FALSE) {
  corM <- cor(t(Exprs), use = "pairwise.complete.obs")
  N <- ncol(Exprs)
  y <- atanh(corM[upper.tri(corM)])
  fix <- which(is.infinite(y))
  if (length(fix) > 0)
    y[fix] <- max(abs(y[-fix]))*(1 + runif(length(fix)))
  set.seed(rndseed)
  sset <- sample(1:length(y),size = min(maxLen,length(y)))
  y0 <- y[sset]
  fittedL2N <- EM(y0*sqrt(N-3))
  rmseL2N <- GoodnessOfFit(fittedL2N)
  if (plot.it)
    plotMixture(fittedL2N,gof=rmseL2N,trim=trim, ttl=ttl)
  if (verbose)
    cat("Calculating the posterior density...\n")
  B <- posteriorDensityL2N(fittedL2N, y*sqrt(N-3))
  p1L2N <- mean(fittedL2N$b1)
  p2L2N <- mean(fittedL2N$b2)
  p0L2N <- 1-(p1L2N+p2L2N)
  m0L2N <- fittedL2N$theta
  m1L2N <- fittedL2N$mu1
  m2L2N <- fittedL2N$mu2
  s0L2N <- fittedL2N$tau
  s1L2N <- fittedL2N$s1
  s2L2N <- fittedL2N$s2

  if (verbose)
    cat("Calculating the log-odds...\n")
  ret <- logoddsValues(fittedL2N$x,m0L2N,s0L2N,m1L2N,s1L2N,
                       m2L2N,s2L2N,p1L2N,p2L2N,vals=1:LOvals)
  if (length(which(ret[,6] < BHthr) > 0)) {
    LogOddsRatio <- max(min(which(ret[,6] < BHthr)),2)
  } else {
    LogOddsRatio <- LOvals
  }
  RtBFL2N <- which(B[[2]]/B[[1]] > LogOddsRatio)
  LtBFL2N <- which(B[[3]]/B[[1]] > LogOddsRatio)

  if (verbose)
    cat("Calculating the adjacency matrix...\n")
  G <- nrow(Exprs)
  sigW <- sort(union(RtBFL2N,LtBFL2N))
  tmpmat <- Matrix::Matrix(0,G, G)
  vec <- rep(0, choose(G,2))
  vec[sigW] <- 1
  tmpmat[upper.tri(tmpmat)] <- vec
  AdjMat <- tmpmat+Matrix::t(tmpmat)

  list(G=G, p1=p1L2N, p2=p2L2N, p0=p0L2N, m0=m0L2N, m1=m1L2N, m2=m2L2N,
       s0=s0L2N, s1=s1L2N, s2=s2L2N, thrtable=ret, LogOddsRatio=LogOddsRatio,
       fitted=fittedL2N, rmse=rmseL2N, rt=RtBFL2N, lt=LtBFL2N, AdjMat=AdjMat)
}


#' Print a short summary of the fitted mixture model.
#'
#' Show the number of nodes, the number of possible and detected edges, the estimated proportion of positively/negatively correlated pairs, and the estimated false discovery rate.
#' @param edgefinderobj The object (list) returned from the edgefinder function.
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    shortSummary(WTres)
#' }
shortSummary <- function(edgefinderobj) {
  with(edgefinderobj,{
    cat("No. nodes =", prettyNum(G,big.mark = ","),"\n")
    cat("Max no. edges =", prettyNum(choose(G, 2),big.mark = ","),"\n")
    cat("No. edges detected =", prettyNum(length(union(lt,rt)),big.mark = ","),"\n")
    cat("p1 =",format(p1,digits=3),"\n")
    cat("p2 =",format(p2,digits=3),"\n")
    cat("Est. FDR <=", format(thrtable[LogOddsRatio,6],digits=3),"\n")
  })
}

# The EM algorithm to fit the L2N model.
#
# Fit the L2N model to normalized correlation coefficients between pairs of genes. The mixture model has three component - the null component follows a normal distribution, and the two non-null components follow lognormal distributions. An edge is in the graph if the correlation between the two end-point genes is large enough and determined to be in one of the non-null components.
# @param x A vector of normalized correlation coefficients.
# @param max.it The maximum number of EM algorithm iterations (default=1000).
# @param tol The tolerance level to assess convergence of the EM algorithm (default=1e-12.)
# @return A list of the parameter estimates for the L2N model.
EM <- function(x, max.it=1000, tol=1e-12) {
  N   <- length(x)
  err <- 1
  # initialize the parameter values
  adjustMean <- mean(x) # centering the data around the mean
  x <- x - adjustMean
  # The parameters of the null ditribution, N(theta,tau) :
  theta <- mean(x)
  tau   <- 1
  # The location and scale parameters of the nonnull components:
  mu    <- abs(quantile(x,c(0.05,.95)))
  names(mu) <- c()
  sig   <- c(1, 1)
  # The initial probabilities of the three components:
  p0 <- 0.98
  p1 <- 0.01
  p2 <- 0.01
  # Set the initial component indicator variables:
  b1 <- rep(0,N)
  b2 <- rep(0,N)
  m1 <- 0
  m2 <- 0
  ct <- 0
  # Run the EM algorithm until the mixture fits the empirical
  # density well (total squared errors < tol)
  while (err > tol) {
    adjustMean <- adjustMean + theta
    x   <- x - theta # iteratively center the data, so that the mean of the
    # null component ends up being 0
    pos <- which(x > 0) # Fit the nonnull components according to the
    neg <- which(x < 0) #      sign of x

    d0 <- dnorm(x, theta, tau) # null component is normal
    d1 <- dlnorm(x, mu[1], sig[1])
    d2 <- dlnorm(-x, mu[2], sig[2])
    wtsm <- p0*d0 + p1*d1 + p2*d2    # The density of the mixture
    b1[-pos] <- 0
    b2[-neg] <- 0
    b1[pos] <- pmin(1,p1*d1[pos]/wtsm[pos])  # Posterior probabilities of the positive nonnull
    b2[neg] <- pmin(1,p2*d2[neg]/wtsm[neg])  # Posterior probabilities of the negative nonnull
    b0      <- 1 - (b1+b2)           # The posterior null probabilities
    # Update the component weights:
    p0 <- mean(b0)
    p1 <- mean(b1)
    p2 <- mean(b2)
    # Update the null component parameters:
    theta <- sum(b0*x)/sum(b0)
    tau   <- sqrt(sum(b0*(x-theta)^2)/sum(b0))
    d0    <- dnorm(x, theta, tau)
    # Update the nonnull (nonnull) components parameters:
    if (sum(b1[pos]) < 1e-2) {
      mu[1]  <- 0
      sig[1] <- 0
      d1     <- rep(0, N)
    } else {
      mu[1]  <- sum(b1[pos]*(log(x[pos])))/sum(b1[pos])
      sig[1] <- sqrt(sum(b1[pos]*(log(x[pos])-mu[1])^2)/sum(b1[pos]))
      d1     <- dlnorm(x, mu[1], sig[1])
    }

    if (sum(b2[neg]) < 1e-2) {
      mu[2]  <- 0
      sig[2] <- 0
      d2     <- rep(0, N)
    } else {
      mu[2]  <- sum(b2[neg]*(log(-x[neg])))/sum(b2[neg])
      sig[2] <- sqrt(sum(b2[neg]*(log(-x[neg])-mu[2])^2)/sum(b2[neg]))
      d2     <- dlnorm(-x, mu[2], sig[2])
    }

    # Check convergence
    err <- sum((p0*d0 + p1*d1 + p2*d2 - wtsm)^2)
    ct <- ct + 1
    if(ct > max.it)
      break
  }
  b1[-pos] <- 0
  b2[-neg] <- 0
  b1[pos] <- pmin(1,p1*d1[pos]/wtsm[pos])  # Posterior probabilities of the positive nonnull
  b2[neg] <- pmin(1,p2*d2[neg]/wtsm[neg])  # Posterior probabilities of the negative nonnull
  b0      <- 1 - (b1+b2)           # The posterior null probabilities
  pvals <- 2*(1-pnorm(abs(x), mean=0, sd=tau))
  bh <- p.adjust(pvals, method="BH")
  list(x=x, adjustMean=adjustMean,
       theta=theta, tau=tau,
       mu1=mu[1], s1=sig[1],
       mu2=mu[2], s2=sig[2],
       b0=b0, b1=b1, b2=b2,
       p.val=pvals, bh=bh,
       err=err, its=ct)
}


# Calculate the log-odds ratios to determine for each gene, in which
# of the three components in the L2N model, it belongs
logoddsValues <- function(y,theta,tau,mu1,s1,mu2,s2,p1,p2,vals=1:30) {
  ret <- matrix(0,nrow=length(vals),ncol=6)
  ret[,1] <- vals
  p0 <- 1-p1-p2
  xs <- seq(min(y),max(y),length=10000)
  pxs <- seq(1e-6,max(y),length=10000)
  nxs <- seq(min(y),-1e-6,length=10000)
  i <- 0
  for (val in vals) {
    i <- i + 1
    if (p1 < 1/length(y)) {
      ret[i,2] <- Inf
    } else {
      f <- function(x) { log((p1*dlnorm(x, mu1, s1))/
                               (p0*dnorm(x,  theta, tau)))-log(val) }
      rt <- try(uniroot(f,  lower = 1e-6, upper = max(y)), silent = T)
      if (class(rt) == "try-error")
        ret[i,2] <- Inf
      else
        ret[i,2] <- rt$root
    }
    if (p2 < 1/length(y)) {
      ret[i,3] <- -Inf
    } else {
      f <- function(x) { log((p2*dlnorm(-x, mu2, s2))/
                               (p0*dnorm(x,  theta, tau)))-log(val) }
      rt <- try(uniroot(f,  lower = min(y), upper = -1e-6), silent = T)
      if (class(rt) == "try-error")
        ret[i,3] <- -Inf
      else
        ret[i,3] <- rt$root
    }
    # type I:
    ret[i,4] <- pnorm(ret[i,3], theta, tau) +
      1 - pnorm(ret[i,2], theta, tau)
    # "Power":
    ret[i,5] <- (p1*(1-plnorm(ret[i,2], mu1, s1)) +
                   p2*(1-plnorm(-ret[i,3], mu2, s2)))/(p1+p2)
    # FDR:
    ret[i,6] <- p0*ret[i,4]/(p0*ret[i,4]+ret[i,5]*(p1+p2))
  }
  colnames(ret) <- c("ppr","Right","Left","TypeI","Power","FDR")
  ret
}


# calculate the posterior L2N mixture model density of x, given the parameter
# estimates
posteriorDensityL2N <- function(fit.em, x) {
  p0 <- mean(fit.em$b0)
  p1 <- mean(fit.em$b1)
  p2 <- mean(fit.em$b2)
  adjustMean <- fit.em$adjustMean + fit.em$theta
  x   <- x - fit.em$adjustMean
  pos <- which(x > 0) # Fit the nonnull components according to the
  neg <- which(x < 0) #      sign of x
  d0 <- dnorm(x, fit.em$theta, fit.em$tau) # null component is normal
  d1 <- dlnorm(x, fit.em$mu1, fit.em$s1)
  d2 <- dlnorm(-x, fit.em$mu2, fit.em$s2)
  wtsm <- p0*d0 + p1*d1 + p2*d2
  b1 <- rep(0, length(x))
  b2 <- rep(0, length(x))
  b1[pos] <- pmin(1,p1*d1[pos]/wtsm[pos])  # Posterior probabilities of the positive nonnull
  b2[neg] <- pmin(1,p2*d2[neg]/wtsm[neg])  # Posterior probabilities of the negative nonnull
  b0      <- 1 - (b1+b2)           # The posterior null probabilities
  list(b0=b0,b1=b1,b2=b2)
}


# Return the estimated density function of the mixture
mixtureDensityL2N <- function(fit.em, x) {
  mean(fit.em$b0)*dnorm(x, fit.em$theta, fit.em$tau) +
    mean(fit.em$b1)*dlnorm(x, fit.em$mu1, fit.em$s1) +
    mean(fit.em$b2)*dlnorm(-x, fit.em$mu2, fit.em$s2)
}


# Calculate the root mean squared error of the fitted mixture
GoodnessOfFit <- function(fit.em, mixturemodel="L2N") {
  x <- sort(fit.em$x)
  if(length(x) > 10000)
    x <- x[seq(1,length(x), length=10000)]
  diffs <- x[-1] - x[-length(x)]
  dnsfn <- approxfun(density(x,bw="SJ"))
  return(sqrt(sum((diffs* (dnsfn(x[-1])-mixtureDensityL2N(fit.em,x[-1])) ^2))))
}


#' Find clusters, and return node characteristics.
#'
#' Take an adjacency Matrix as input and find clusters. For each node, find the degree and clustering coefficient (CC). Then, calculate a centrality measure (type\*CC+1)\*deg. For type=0, it's just the degree. Note that setting type=1 we assign a higher value to nodes that not only have many neighbors, but the neighbors are highly interconnected. For example, suppose we have two components with k nodes, one has a star shape, and the other is a complete graph. With type=0 both graphs will get the same value, but with type=1 the complete graph will be picked by the algorithm first. Setting type to a negative value gives CC\*deg as the centrality measure.
#' @param A An adjacency Matrix(0/1).
#' @param minCtr The minimum centrality value to be considered for a cluster center (default=5).
#' @param type Determines how the centrality measure is computed.
#' @return A data frame with the following columns
#' \itemize{
#'  \item{labels} {Node label (e.g. gene names).}
#' \item{degree} {Node degree.}
#' \item{cc} {Node clustering coefficient.}
#' \item{ctr} {Node centrality measure: (type\*CC+1)\*deg, or CC\*deg if type is negative.}
#' \item{clustNo} {Cluster number.}
#' \item {iscenter} {1 for the node was chosen as the cluster's center, 0 otherwise.}
#' \item {intEdges} {Number of edges from the node to nodes in the same cluster.}
#' \item {extEdges} {Number of edges from the node to nodes NOT in the same cluster.}
#' \item {distCenter} {Standardized Manhattan distance to the central node.}
#' }
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    Sres <- edgefinder(SIM, ttl = "hub network")
#'    SimComp <- graphComponents(Sres$AdjMat)
#'    head(SimComp)
#' }
graphComponents <- function(A, minCtr=5, type=1) {
  stopifnot(grep("Matrix", class(A)) > 0)
  Vn <- ncol(A)
  ctrs <- rep(2*Vn, Vn)
  labels <- 1:Vn
  if(!is.null(rownames(A)))
    labels <- rownames(A)
  deg <- Matrix::rowSums(A)
  CC <- clusteringCoef(A)
  ctrs <- (type*CC+1)*deg
  if (type < 0)
    ctrs <- CC*deg
  clustersInfo <- data.frame(labels=labels, degree=deg, cc=CC, ctr=ctrs,
                             clustNo=rep(0,Vn), iscenter=rep(0,Vn),
                             intEdges=rep(0,Vn), extEdges=rep(0,Vn),
                             distCenter=rep(0,Vn))
  clustNo <- 1
  clustered <- which(deg < 1)
  while(length(clustered) < Vn) {
    notInCluster <- setdiff(1:Vn, clustered)
    if (max(ctrs[notInCluster]) < minCtr)
      return(clustersInfo)
    ctrnode <- notInCluster[which.max(ctrs[notInCluster])]
    # candidate cluster neighbors
    nbrs <- setdiff(sort(c(ctrnode, which(A[ctrnode,] != 0))), clustered)
    if(length(nbrs) > 1) {
      if (length(nbrs) > minCtr) {
        clustersInfo$iscenter[ctrnode] <- 1
        clustersInfo$clustNo[union(ctrnode,nbrs)] <- clustNo
        clustersInfo$intEdges[nbrs] <- Matrix::rowSums(A[nbrs,nbrs])
        if (length(nbrs) < ncol(A)) {
          clustersInfo$extEdges[nbrs] <- Matrix::rowSums(A[nbrs,-nbrs])
        } else {
          clustersInfo$extEdges[nbrs] <- 0
        }
        for (i in 1:length(nbrs)) {
          clustersInfo$distCenter[nbrs[i]] <- mean(xor(A[ctrnode,], A[nbrs[i],]))
        }
        clustNo <- clustNo + 1
      }  else {
        nbrs <- c()
      }
    } else {
      nbrs <- c()
    }
    clustered <- union(clustered, c(nbrs, ctrnode))
  }
  return(clustersInfo)
}



#' Show cluster characteristics.
#'
#' Takes an object obtained from graphComponents and prints and returns summary statistics.
#' @param clustersInfo Obtained from graphComponents.
#' @return A matrix with cluster number, number of nodes, and fivenum summaries for the degrees of nodes in the cluster, and the percentage of edges that are within the cluster.
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    WTComp <- graphComponents(WTres$AdjMat)
#'    (summtab <- summarizeClusters(WTComp))
#' }
summarizeClusters <- function(clustersInfo) {
  cat("Num of nodes:", nrow(clustersInfo),"\n")
  cat("Num of edges:", sum(clustersInfo$degree)/2,"\n")
  cat("Num of clusters:", max(clustersInfo$clustNo),"\n")
  cat("Num of unclustered nodes:", length(which(clustersInfo$clustNo == 0)),"\n")
  percentInCluster <- clustersInfo$intEdges/clustersInfo$degree
  percentInCluster[which(clustersInfo$degree == 0)] <- 0
  if (max(clustersInfo$clustNo) == 0)
    return(NULL)
  tab <- matrix(0,nrow=max(clustersInfo$clustNo),ncol=12)
  for (cnum in 1:max(clustersInfo$clustNo)) {
    tmpclusterInfo <- clustersInfo[which(clustersInfo$clustNo == cnum),]
    tab[cnum,] <- c(cnum,nrow(tmpclusterInfo), fivenum(tmpclusterInfo$degree),
                    fivenum(percentInCluster[which(clustersInfo$clustNo == cnum)]))
  }
  colnames(tab) <- c("Cluster","Nodes","degreeMin","degreeQ25","degreeMedian",
                     "degreeQ75","degreeMax","pctInClstMin","pctInClstQ25",
                     "pctInClstMedian", "pctInClstQ75","pctInClstMax")
  tab
}


#' Return an adjacency matrix after collapsing clusters into their central nodes.
#'
#' Takes an object obtained from graphComponents and prints summary statistics.
#' @param A An adjacency Matrix.
#' @param clustersInfo Obtained from graphComponents
#' @return A weighted adjacency matrix between clusters and unclustered nodes.
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    rownames(WTres$AdjMat) = rownames(WT)
#'    WTComp <- graphComponents(WTres$AdjMat)
#'    Adj1 <- collapsedGraph(WTres$AdjMat, WTComp) > 0
#'    plotBitmapCC(Adj1,showMinDegree = 2)
#' }
collapsedGraph <- function(A, clustersInfo) {
  collDim <- length(which(clustersInfo$clustNo == 0)) + max(clustersInfo$clustNo)
  collA <- Matrix::Matrix(0, ncol=collDim, nrow=collDim)
  inCluster <- which(clustersInfo$clustNo > 0)
  notInCluster <- which(clustersInfo$clustNo == 0)
  if (length(notInCluster) > 0) {
    collA[1:length(notInCluster), 1:length(notInCluster)] <- A[notInCluster, notInCluster]>0
  }
  if (length(rownames(A)) != nrow(A)) {
    rownames(A) <- 1:nrow(A)
  }
  rownames(collA) <- c(rownames(A)[notInCluster],
                       paste0("CLS",1:max(clustersInfo$clustNo)))
  for (i in 1:max(clustersInfo$clustNo)) {
    Ci <- which(clustersInfo$clustNo == i)
    if (length(notInCluster) > 0) {
      Atmp <- matrix(A[notInCluster,which(clustersInfo$clustNo==i)],
                     nrow=length(notInCluster), ncol=length(which(clustersInfo$clustNo==i)))
      collA[i+length(notInCluster),1:length(notInCluster)] <- Matrix::rowSums(Atmp)
    }
    if (i < max(clustersInfo$clustNo)) {
      for (j in (i+1):max(clustersInfo$clustNo)) {
        Cj <- which(clustersInfo$clustNo == j)
        collA[i+length(notInCluster),j+length(notInCluster)] <- sum(A[Ci,Cj])
      }
    }
  }
  collA + Matrix::t(collA)
}



#' Calculate the clustering coefficient of each node.
#'
#' @param A an adjacency Matrix (0/1).
#' @return A vector with the clustering coefficient of each node.
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    clusteringCoef(WTres$AdjMat)
#' }
#'
clusteringCoef <- function(A) {
  rsum <- Matrix::rowSums(A)
  cc <- rep(0,nrow(A))
  for (i in 1:nrow(A)) {
    if (rsum[i] <= 1)
      cc[i] <- 0
    else {
      nbrs <- which(A[i,] == 1)
      At <- A[nbrs, nbrs]
      cc[i] <- 0.5*sum(At)/choose(rsum[i],2)
    }
  }
  cc
}


#' Plot the histogram of the data and the fitted mixture distribution.
#'
#' The function is called by the edgefinder function.
#' @param fit.em The object (list) returned from the EM function with the parameter estimates for the L2N model.
#' @param gof The root mean-squared error of the fitted model (to appear in the title of the plot).
#' @param ttl The title of the plot (default="").
#' @param trim The proportion of extreme values on both sides of the distribution to eliminate from the plot (default=0.) This can be useful if a small number of values are so extreme, that the plot shows mostly the tails and a spike in the middle. Default=0.
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    plotMixture(WTres$fitted, WTres$rmse)
#' }
plotMixture <- function(fit.em, gof, ttl="", trim=0) {
  xlim <- quantile(fit.em$x, c(trim/2, 1-trim/2))
  brks <- min(80,floor(length(fit.em$x)/100))
  hist(fit.em$x, freq=FALSE, breaks=brks,
       main=sprintf("%s\nrMSE %2.2f",ttl, gof),
       xlim=xlim,xlab="x", border="white", col="wheat")
  xs <- seq(min(fit.em$x), max(fit.em$x), length=1000)
  p0 <- mean(fit.em$b0)
  p1 <- mean(fit.em$b1)
  p2 <- mean(fit.em$b2)
  lines(xs,  p0*dnorm(xs,  fit.em$theta, fit.em$tau), col=2, lwd=2)
  lines(xs,  p1*dlnorm(xs, fit.em$mu1,   fit.em$s1),  col=3, lwd=2)
  lines(-xs, p2*dlnorm(xs, fit.em$mu2,   fit.em$s2),  col=3, lwd=2)
  mxfit <- p0*dnorm(xs,fit.em$theta, fit.em$tau) +
    p1*dlnorm(xs, fit.em$mu1, fit.em$s1) +
    p2*dlnorm(-xs, fit.em$mu2, fit.em$s2)
  lines(xs, mxfit, lwd=3, col=4, lty=2)
}


#' Plot the degree of nodes versus the degree times the clustering coefficient.
#'
#' The x-axis represents the number of neighbors of each node, and the y-axis represents the proportion of neighbors which are connected to each other.
#' @param edgefinderobj The object (list) returned by edgefinder.
#' @param clusterInfo obtained from graphComponents. If not provided by the user, it will be computed on the fly.
#' @param highlightNodes A vector of node-numbers which will be shown in red. Default is NULL.
#' @export
#' @import stats graphics
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    WTComp <- graphComponents(WTres$AdjMat)
#'    plotDegCC(WTres,WTComp)
#' }
plotDegCC <- function(edgefinderobj, clusterInfo=NULL, highlightNodes=NULL) {
  if (is.null(clusterInfo))
    clusterInfo <-  graphComponents(edgefinderobj$AdjMat)
  cc0 <- clusterInfo$cc
  deg0 <- clusterInfo$degree
  plot(deg0, deg0*cc0,axes=F,xlim=c(0,max(deg0)),
       ylim=c(0,1.1*max(deg0*cc0)),main="",
       xlab=bquote("degree"),ylab=bquote("CC*degree"),
       col="thistle",pch=24,cex=0.5); axis(1); axis(2)
  grid(); abline(0,1,col="seagreen1", lwd=2)
  if (!is.null(highlightNodes))
    points(deg0[highlightNodes],(deg0*cc0)[highlightNodes],col=2,pch=24,cex=0.5)
}


#' Edge-indicator bitmap plot.
#'
#' Plot a bitmap in which a black dot corresponds to a pair of highly correlated genes (an edge in the graph).
#' The default is to show the nodes according to their order in the input.
#' By setting orderByDegree=T as below, it is possible to change the order and cluster them, and show them in increasing degree order (from left to right.)
#' @param AdjMat An adjacency Matrix (0/1).
#' @param clusterInfo obtained from graphComponents. If not provided by the user, it will be computed on the fly.
#' @param orderByCluster If false, show the bitmap is the original node order. If TRUE, show nodes by clusters, and sort by distance from the center of the cluster.
#' @param showMinDegree Non-negative integer indicating the minimum degree of nodes that should be displayed. Default=0 (all nodes).
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    WTComp <- graphComponents(WTres$AdjMat)
#'    plotBitmapCC(WTres$AdjMat)
#'    plotBitmapCC(WTres$AdjMat, WTComp, orderByCluster=TRUE)
#'    plotBitmapCC(WTres$AdjMat, WTComp, orderByCluster=TRUE, showMinDegree = 30)
#' }
plotBitmapCC <- function(AdjMat, clusterInfo=NULL, orderByCluster=FALSE, showMinDegree=0) {
  if(!is.null(clusterInfo))
    orderByCluster <- TRUE
  if (orderByCluster) {
    if (is.null(clusterInfo))
      clusterInfo <- graphComponents(AdjMat)
    nodeOrder <- order(clusterInfo$clustNo,clusterInfo$distCenter)
    AdjMat <- AdjMat[nodeOrder, nodeOrder]
  }
  showNodes <- which(Matrix::rowSums(AdjMat) >= showMinDegree)
  Matrix::image(AdjMat[showNodes, showNodes])
}


#' Plot cluster network
#'
#' Plot a cluster network with all the nodes and edges - the central node is marked by a black circle. The radius of each point corresponds to its degree. The opacity corresponds to the percentage of edges from the node that is in the cluster (the darker it is, the larger the percentage of edges is within the cluster.) The distance from the center corresponds to the relative dissimilarity with the central node. This is computed as the number of neighbors the node and the central node do not have in common.
#' @param AdjMat An adjacency Matrix (0/1).
#' @param clustNo The chosen cluster.
#' @param clusterInfo Obtained from graphComponents.
#' @param labels If set to TRUE, show node names (default=FALSE).
#' @param nodecol The color(s) of the nodes. Can be a single value or a vector of length equal to the number of rows in AdjMat
#' @param labelsize Text size of node labels.
#' @param figtitle The title of the plot (default=NULL).
#' @param negcor The pairs which are negatively correlated, to be drawn as red edges (default=NULL). If set to null, all edges will have the same color.
#' @param edgecols The colors to be used for edges. Default="grey88". If one value is given, all edges will be drawn using this color. If negcor is used to specify which edges are negatively correlated and edgecol contains two valid colors, the first is used for positive correlations, and the second for negative ones.
#' @export
#' @examples
#' \donttest{
#'    data(WT)
#'    WTres <- edgefinder(WT, ttl = "Wild Type")
#'    WTComp <- graphComponents(WTres$AdjMat)
#'    plotCluster(WTres$AdjMat, 5, WTComp)
#'    # plotCluster(WTres$AdjMat, 5, WTComp, negcor = WTres$lt, edgecols = c("grey88","orange"))
#' }
plotCluster <- function(AdjMat, clustNo, clusterInfo=NULL, labels=FALSE, nodecol="blue", labelsize=1, figtitle=NULL, negcor=NULL, edgecols="grey88") {
  if(is.null(clusterInfo))
    clusterInfo <- graphComponents(AdjMat)
  if(length(nodecol) < nrow(AdjMat))
    nodecol <- rep(nodecol[1],length=nrow(AdjMat))
  if (length(negcor) > 0) {
    tmpmat <- Matrix::Matrix(0,nrow(AdjMat), ncol(AdjMat))
    vec <- rep(0,choose(nrow(AdjMat), 2))
    vec[negcor] <- -2
    tmpmat[upper.tri(tmpmat)] <- vec
    AdjMat <- AdjMat + (tmpmat+Matrix::t(tmpmat))
  }

  ids <- which(clusterInfo$clustNo == clustNo)
  if (length(ids) > 0) {
    tmpA <- AdjMat[ids,ids]
    tmpclusterInfo <- clusterInfo[ids,]
    rads <- round(10*tmpclusterInfo$distCenter/max(tmpclusterInfo$distCenter))
    thetas <- rep(0,length(rads))
    intvls <- findInterval(rads,seq(1,10))
    for (intvl in unique(sort(intvls))) {
      pts <- which(intvls == intvl)
      thetas[pts] <- 3*intvl*pi/max(intvls)+seq(0,1.9*pi,length=length(pts))
    }
    sizes <- pmax(0.3,tmpclusterInfo$degree/max(tmpclusterInfo$degree))
    opacity <- 0.25+tmpclusterInfo$intEdges/tmpclusterInfo$degree
    opacity <- opacity/max(opacity)
    nodecol <- rgb(t(col2rgb(nodecol)/255),alpha=opacity)[ids]
    plot(rads*cos(thetas), rads*sin(thetas),cex=sizes*3, pch=19,axes=F,
         xlab="",ylab="",col=nodecol, main=figtitle)
    for (i in 1:ncol(tmpA)) {
      nbrs <- setdiff(which(abs(tmpA[i,]) == 1), 1:i)
      if(length(nbrs) > 0) {
        edgecol <- rep(edgecols[1], ncol(tmpA))
        edgecol[which(tmpA[i,nbrs] == -1)] <- "deepskyblue4"
        if (edgecols[2] %in% colours()) {
          edgecol[which(tmpA[i,nbrs] == -1)] <- edgecols[2]
        }
        for (j in nbrs) {
          lines(c(rads[i]*cos(thetas[i]), rads[j]*cos(thetas[j])),
                c(rads[i]*sin(thetas[i]), rads[j]*sin(thetas[j])),
                col=edgecol[j], lwd=0.5)
        }
      }
    }
    points(rads*cos(thetas), rads*sin(thetas),cex=sizes*3, pch=19, col=nodecol)
    if (labels)
      text(rads*cos(thetas), rads*sin(thetas), tmpclusterInfo$labels, pos=3, cex=labelsize)
    ctr <- which(tmpclusterInfo$iscenter==1)
    points(rads[ctr]*cos(thetas[ctr]), rads[ctr]*sin(thetas[ctr]),pch=21,
           cex=sizes[ctr]*3, col="black",lwd=2)
  } else {
    cat("Invalid cluster number\n")
  }
}


#' Return a Matrix with the shortest path distance between nodes (check up to numSteps.)
#'
#' return the adjacency matrix of expMat connecting neighbors up to numSteps away.
#' @param AdjMat An adjacency Matrix (0/1).
#' @param numSteps The maximum number of edges between pairs of nodes. If numSteps=0, returns the input matrix. numSteps=1 adds neighbors of direct neighbors, etc.
#' @return A Matrix containing the shortset paths between nodes i and j
#' @export
#' @examples
#' \donttest{
#'    data(SIM)
#'    Sres <- edgefinder(SIM, ttl = "hub network")
#'    AdjMat1 <- shortestPathDistance(Sres$AdjMat, numSteps=50)
#'    max(AdjMat1)
#'    Matrix::image(AdjMat1)
#' }
shortestPathDistance <- function(AdjMat, numSteps=0) {
  degs <- 1:ncol(AdjMat)
  if (numSteps == 0)
    return(AdjMat)
  An <- Ap <- minDist <- AdjMat
  for (i in 1:numSteps) {
    An <- Ap%*%AdjMat
    if (sum((An | Ap) - (An & Ap)) == 0)
      break
    minDist[(An > 0) & (Ap == 0) & (minDist == 0)] <- i
    Ap <- An
  }
  rownames(minDist) <- colnames(minDist) <- rownames(AdjMat)
  minDist
}


#' Gene Expression data for the WildType group
#'
#' WT is a matrix with normalized gene expression data containing 3454 differentially expressed genes (when compared with the duplication group) from 15 samples (columns) from the wild-type group.
#'
#' @docType data
#' @keywords datasets
#' @name WT
#' @usage data(WT)
#' @format A matrix with 3454 rows and 15 columns
#' @references \url{https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4430}
NULL


#' Gene Expression data for the Duplication group
#'
#' DUP is a matrix with normalized gene expression data containing 3454 differentially expressed genes (when compared with wild-type) from 12 samples (columns) from the duplication group.
#'
#' @docType data
#' @keywords datasets
#' @name DUP
#' @usage data(DUP)
#' @format A matrix with 3454 rows and 12 columns.
#' @references \url{https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4430}
NULL


#' Simulated gene Expression data using the huge package
#'
#' SIM is a a simulated dataset with a hub structure, consisting of 1000 nodes and 50 hubs
#'
#' @docType data
#' @keywords datasets
#' @name SIM
#' @usage data(SIM)
#' @format A 1000 by 200 matrix, representing 50 hubs
NULL


