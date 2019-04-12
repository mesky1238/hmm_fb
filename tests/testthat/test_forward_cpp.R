context("Forward")

test_that("sweep_add works as expected",{
  
  A <- matrix(1:4,2,2)
  b <- 1:2
  expect_equal(sweep_add(A,b),sweep(A, 2, b, "+"))
})

test_that("sweep_add_logexp works as expected",{
  
  A <- matrix(1:4,2,2)
  b <- 1:2
  expect_equal(sweep_add_logexp(A,b),apply(sweep(A, 2, b, "+"), 1, log.sum.exp))

})

test_that("I'm overall the same as the old oversion",{
  
  data("test_data")
  expect_equal(forward(test_data$init, test_data$A, test_data$B, test_data$ntimes, return.all=FALSE, homogeneous=TRUE,useC=TRUE,na.allow=TRUE),
               forward_old(test_data$init, test_data$A, test_data$B, test_data$ntimes, return.all=FALSE, homogeneous=TRUE,useC=TRUE,na.allow=TRUE))
})

test_that("I can generate the starting conditions in C++",{
  data("test_data")
  init <- test_data$init
  init <- log(init) # init to log scale
  A <- log(test_data$A)
  B <- cbind(test_data$B[,,1],test_data$B[,,2]) # simingz to save time: should be apply(B,c(1,3),sum)
  
  nt <- dim(B)[1]
  ns <- ncol(init)
  
  lt <- length(test_data$ntimes)
  et <- cumsum(test_data$ntimes)
  bt <- c(1,et[-lt]+1)
  cpp_bt <- gen_bt(linit = init,lA = A[1,,],B = B,ntimes = test_data$ntimes)
  expect_equal(bt,cpp_bt+1)
})

test_that("Entirely C++ gives identical results to optimized R+c++ code",{
  data("test_data")
  
  lA <- log(test_data$A[1,,])
  nB <- cbind(test_data$B[,,1],test_data$B[,,2]) # simingz to save time: should be apply(B,c(1,3),sum)
  of_res <- forward(test_data$init, test_data$A, test_data$B, test_data$ntimes, return.all=T, homogeneous=TRUE,useC=TRUE,na.allow=TRUE)
  rcpp_res <- forward_rcpp(linit = log(test_data$init),lA = lA,B = nB,ntimes =  test_data$ntimes,return_all = T)
  expect_equal(of_res,rcpp_res)
})
