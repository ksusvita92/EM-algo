R scripts:
  - **sanity check.R**: there are 2 functions defined here, *rmixGamma* and *EMsanityCheck*. The former is used to generate random sample from a mixture Gamma model. The inputs are:
      - *n*: #of generated random samples
      - *alpha*: shape parameter
      - *beta*: scale parameter
      - *weight*: weight vector for each component.
     
     The later is used to perform EM-algorithm sanity test, giving linear regression test for true parameters vs estimates. The inputs are:
       - *nrep*: #of repetitions. Each repetition generates 100 data points.
       - *k*: #of components in the mixture model. Value ranges from 1 to 4.
       - *discretization*: default is F. 
       - *rounding*: default is F. If TRUE, the generated data points will be rounded to represent real SI.
       
  - bbb
