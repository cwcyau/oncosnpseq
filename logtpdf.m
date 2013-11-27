function [y, u] = logtpdf(x, m, S, v)
%
% one-dimensional log student t pdf
%
[p, n] = size(x);

z = bsxfun(@minus, x, m);

d = (1/S)*(z.*z);

y = -0.5*log(S) - 0.5*(v+p)*log(1 + d/v);
u = (v + p)./(v + d);
