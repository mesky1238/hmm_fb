forward <- function (init,A,B,ntimes=NULL,return.all=FALSE,
                     homogeneous=TRUE,useC=TRUE,na.allow=TRUE) {
  
  # Forward-Backward algorithm (used in Baum-Welch)
  # Returns alpha, beta, and full data likelihood
  
  # NOTE THE CHANGE IN FROM ROW TO COLUMN SUCH THAT TRANSPOSING A IS
  # NOT NECCESSARY ANYMORE IN COMPUTING ALPHA AND BETA BUT IS NOW
  # NECCESSARY IN COMPUTING XI A = T*K*K array with transition
  # probabilities, from row to column!!!!!!!  B = T*D*K matrix with
  # elements ab_{ij} = P(y_i|s_j) init = N*K vector with initial
  # probabilities
  
  # NOTE: to prevent underflow, alpha and beta are scaled, using sca
  
  # NOTE: xi[t,i,j] = P(S[t] = j & S[t+1] = i) !!!NOTE the order of i and j!!!
  
  # Note by simingz: this script is for log scale B!
  init <- log(init) # init to log scale
  A <- log(A)
  B <- cbind(B[,,1],B[,,2]) # simingz to save time: should be apply(B,c(1,3),sum)
  
  nt <- dim(B)[1]
  ns <- ncol(init)
  
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  
  alpha <- matrix(ncol=ns,nrow=nt)
  sca <- vector(length=nt)
  tA <- A[1,,]
  alpha <- t(alpha)
  for(case in 1:lt) {
    alpha[,bt[case]] <- init[case,] + B[bt[case],] # initialize
    sca[bt[case]] <- -log.sum.exp(alpha[,bt[case]])
    alpha[,bt[case]] <- alpha[,bt[case]] + sca[bt[case]]
    
    if(ntimes[case]>1) {
      ecase <- et[case]-1
      for(i in bt[case]:ecase) {
        alpha[,i+1] <- sweep_add_logexp(tA,alpha[,i]) + B[i+1,]
        sca[i+1] <- -log.sum.exp(alpha[,i+1])
        alpha[,i+1] <- sca[i+1] + alpha[,i+1]
      }
    }
  }  
  
  like <- -sum(sca)
  
  if(return.all) {
    res <- list(alpha=t(alpha),sca=sca,logLike=like)
  } else {
    res <- list(logLike=like)
  }
  res
}

forward_old <- function (init,A,B,ntimes=NULL,return.all=FALSE,
                         homogeneous=TRUE,useC=TRUE,na.allow=TRUE) {
  
  # A = T*K*K array with transition probabilities, from row to
  # column!!!!!!!  B = T*D*K matrix with elements ab_{ij} = P(y_i|s_j)
  # init = N*K vector with initial probabilities
  
  # T = total number of time points
  # K = number of states
  # D = dimension of observations (D>1 is multivariate)
  # N = number of participants
  
  # NOTE: xi[t,i,j] = P(S[t] = j & S[t+1] = i) !!!NOTE the order of i
  # and j!!!
  
  # note by simingz: this script is for log scale B!
  
  init <- log(init) # init to log scale
  A <- log(A)

  # Should be apply(B,c(1,3),sum), but for my specific problem, I know
  # D=1 and K=2, so no need to sum
  B <- cbind(B[,,1],B[,,2]) 
  nt <- dim(B)[1]
  ns <- ncol(init)
  
  lt <- length(ntimes)
  et <- cumsum(ntimes)
  bt <- c(1,et[-lt]+1)
  
  alpha <- matrix(ncol=ns,nrow=nt)
  beta <- matrix(ncol=ns,nrow=nt)
  sca <- vector(length=nt)
  xi <- array(dim=c(nt,ns,ns))
  
  for(case in 1:lt) {
    alpha[bt[case],] <- init[case,] + B[bt[case],] # initialize
    sca[bt[case]] <- -log.sum.exp(alpha[bt[case],])
    alpha[bt[case],] <- alpha[bt[case],] + sca[bt[case]]
    
    if(ntimes[case]>1) {
      for(i in bt[case]:(et[case]-1)) {
        alpha[i+1,] <- apply(sweep(A[1,,], 2, alpha[i,], "+"), 1, log.sum.exp) + B[i+1,]
        sca[i+1] <- -log.sum.exp(alpha[i+1,])
        alpha[i+1,] <- sca[i+1] + alpha[i+1,]
      }
    }
  }  
  
  like <- -sum(sca)
  
  if(return.all) {
    res <- list(alpha=alpha,sca=sca,logLike=like)
  } else {
    res <- list(logLike=like)
  }
  res
}

log.sum.exp <- function (x) {
    
  # Computes log(sum(exp(x))
  #
  # Uses offset trick to avoid numeric overflow:
  # http://jblevins.org/notes/log-sum-exp
  if (max(abs(x)) > max(x))
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}


