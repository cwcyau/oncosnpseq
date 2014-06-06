function loglik = binopdfln(k, d, p)
%
% loglik = binopdfln(k, d, p)
%
loglik = gammaln(d+1) - gammaln(k+1) - gammaln(d-k+1) + k.*log(p) + (d-k).*log(1-p); 
