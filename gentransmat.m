function [arrayind, log_nu, logT] = gentransmat(phi, delta, tumourState, uu)

S = size(tumourState, 1);
U = length(uu);

cn =  tumourState(:, 4);
loh = tumourState(:, 5);

%arrayind = [];
%Z = 0;
%loh = tumourState(:, 5);
%for si = 1 : S
%	if loh(si) == 2 | loh(si) == -1
%		Z = Z + 1;
%		arrayind = [ arrayind; si 1 ];
%	end
%	if loh(si) == 0 | loh(si) == 1
%		Z = Z + U;
%		arrayind = [ arrayind; si*ones(U, 1) [1:U]' ];
%	end
%end

%log_nu = zeros(Z, 1);
%logT = zeros(Z, Z);

%for zi = 1 : Z

%	si = arrayind(zi, 1);
%	ui = arrayind(zi, 2);

%	for zj = 1 : Z

%		sj = arrayind(zj, 1);
%		uj = arrayind(zj, 2);

%		if si == sj & ui == uj
%			logT(zi, zj) = 0;
%		end
%		
%		if si ~= sj & ui == uj
%%			if cn(si) == cn(sj)
%%				if ( loh(si) == 0 & loh(sj) == 2 ) | ( loh(si) == 2 & loh(sj) == 0 )
%%					logT(zi, zj) = -1;	
%%				end
%%			else					
%				logT(zi, zj) = phi;
%%			end
%		end

%		if si == sj & ui ~= uj 
%			logT(zi, zj) = delta;
%		end

%		if si ~= sj & ui ~= uj
%%			if cn(si) == cn(sj)
%%				if ( loh(si) == 0 & loh(sj) == 2 ) | ( loh(si) == 2 & loh(sj) == 0 )
%%					logT(zi, zj) = -1 + delta;	
%%				end
%%			else					
%				logT(zi, zj) = phi + delta;
%%			end
%		end				

%	end

%end


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

