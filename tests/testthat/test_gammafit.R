

y = c(1090,2041,2725,1431,889,1199,1694,2310,1201,1011,995,601,1457,1942,1987)

set.seed(1) #the seed throwed an unexpected error in fitgammamodel2
th1 = fitgammamodel2(y,niter = 100)
th2 = c(752.4333,0.5435)
expect_equal(th1,th2,tol=1e-4)
