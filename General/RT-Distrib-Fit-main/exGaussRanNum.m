% Generate ex-Gaussian distributed random numbers
function r = exGaussRanNum(tau,mu,sigma,n,m)

r = normrnd(mu,sigma,n,m)+exprnd(tau,n,m);