# hmm_fb

To convert a R implementation of forward algorithm to C implementation.

The R function that we are planning to convert to C is called
`forward()`. The R script is [here](code/HMM-func-forward.R). To run
this function please see
[here](http://htmlpreview.github.com/?https://github.com/simingz/hmm_fb/blob/master/test/run_forward_probablity_algorithm.html).

This R function is adapted from the
[depmixS4 package](https://cran.r-project.org/web/packages/depmixS4/index.html). We
revised it because of the original one have scaling problems when the
emmission density is very small. We applied the log sum exponetial
trick to avoid this.

The original R implementation is [here](code/fb.R), the accompanied c
implementation is [here](code/fb.h) and [here](code/fb.c) .
