function z = logsumexp(a, dim)

if nargin < 2
  	dim = 1;
end

ndim = ndims(a);

[amax, i] = max(a, [], dim);

dims = ones(1, ndim);
dims(dim) = size(a, dim);

a = bsxfun(@minus, a, amax);

z = amax + log(sum(exp(a), dim));

z(~isfinite(amax)) = amax(~isfinite(amax));

