function [arrayind, log_nu, logT] = gentransmat(phi, delta, tumourState, uu)

S = size(tumourState, 1);
U = length(uu);

cn =  tumourState(:, 4);
loh = tumourState(:, 5);

log_nu = zeros(S*U, 1);
logT = zeros(S*U, S*U);

for si = 1 : S

	for ui = 1 : U
		
		ind = (si-1)*U + ui;
		
		arrayind(ind, 1) = si;
		arrayind(ind, 2) = ui;

		for sj = 1 : S

			for uj = 1 : U
			
				ind2 = (sj-1)*U + uj;
			
				if si == sj & ui == uj
					logT(ind, ind2) = 0;
				end
				
				if si ~= sj & ui == uj		
					logT(ind, ind2) = phi;
				end

				if si == sj & ui ~= uj 
					logT(ind, ind2) = delta;
				end

				if si ~= sj & ui ~= uj					
					logT(ind, ind2) = phi + delta;
				end				

			end

		end
		
	end

end

logT = logT - repmat(logsumexp(logT, 2), [1 S*U]);

