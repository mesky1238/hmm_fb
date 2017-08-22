load("inst/old_package/data/hmm.toy.Rd")
test_data <- list(init=init,A=A,B=B,ntimes=ntimes)
devtools::use_data(test_data)
library(microbenchmark)
first_bench <- microbenchmark(Optimized=forward(test_data$init, test_data$A, test_data$B, test_data$ntimes, return.all=FALSE, homogeneous=TRUE,useC=TRUE,na.allow=TRUE),
                      Naive=forward_old(test_data$init, test_data$A, test_data$B, test_data$ntimes, return.all=FALSE, homogeneous=TRUE,useC=TRUE,na.allow=TRUE),
                      times=10,unit="relative")
devtools::use_data(first_bench)

lA <- log(test_data$A[1,,])
nB <- cbind(test_data$B[,,1],test_data$B[,,2])

second_bench <- microbenchmark(
  Optimized=forward(test_data$init, test_data$A, test_data$B, test_data$ntimes,return.all=T),
  All_Rcpp= forward_rcpp(linit = log(test_data$init),lA = lA,B = nB,ntimes =  test_data$ntimes,return_all = T),
  Naive=forward_old(test_data$init, test_data$A, test_data$B, test_data$ntimes, return.all=T),times=10,unit="relative")
devtools::use_data(second_bench)

profvis_dat <- profvis(forward_old(init,A,B,ntimes=ntimes,return.all = T))
devtools::use_data(profvis_dat)
