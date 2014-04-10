function [log_nu, logT, log_nu_U, logT_U] = gentransmatS(phi, tumourState, params)

S = size(tumourState, 1);
U = params.U;

log_nu = zeros(S, 1);
logT = zeros(S, S);

for si = 1 : S

	for sj = 1 : S

		if si == sj
			logT(si, sj) = 0;
		end
		
		if si ~= sj
			logT(si, sj) = phi;
		end

	end

end


log_nu_U = zeros(U, 1);
logT_U = zeros(U, U);

for ui = 1 : U

	for uj = 1 : U

		if ui == uj
			logT_U(si, sj) = 0;
		end
		
		if ui ~= uj
			logT_U(si, sj) = phi;
		end

	end

end

