function [vpath, phi_max] = viterbi(log_nu, log_obslik, log_T)

[S, T] = size(log_obslik);

phi = zeros(S, T);
delta = zeros(S, T, 'uint32');

for s = 1 : S
	phi(s, 1) = log_nu(s) + log_obslik(s, 1);
end

for t = 2 : T

	for s = 1 : S
		tmp = zeros(S, 1);
		for s_prev = 1 : S
			tmp(s_prev) = phi(s_prev, t-1) + log_obslik(s, t) + log_T(s_prev, s);
		end
		[phi(s, t), delta(s, t)] = max(tmp);
	end

end

vpath = zeros(1, T, 'uint32');
[phi_max, vpath(T)] = max(phi(:, T));
for t = T-1:-1:1
	vpath(t) = delta(vpath(t+1), t+1);
end

