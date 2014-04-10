function seg = findsegments(chr, arm, pos, v, x, u, loglik, patients, options, params)

tumourState = options.tumourState;

seg = [];

N = length(x);
S = params.S;
U = params.U;

nseg = 0;

for chrNo = options.chrRange

	for armNo = 1 : 2

		chrloc = find( chr == chrNo & arm == armNo );
		n_chr = length(chrloc);

		pos_chr = pos(chrloc);
		v_chr = v(chrloc);
		x_chr = x(chrloc);
		u_chr = u(chrloc);
		loglik_chr = loglik(:, chrloc);
		
		startInd = 1;
		for i = 2 : n_chr

			cn_prev = tumourState(x_chr(i-1), 4);
			loh_prev = tumourState(x_chr(i-1), 5);
			majorcn_prev = tumourState(x_chr(i-1), 3);
			minorcn_prev = tumourState(x_chr(i-1), 2);			

			cn = tumourState(x_chr(i), 4);
			loh = tumourState(x_chr(i), 5);

			if ( x_chr(i) ~= x_chr(i-1) ) | i == n_chr

				nseg = nseg + 1;

				endInd = i-1;

				seg{nseg}.chromosome = chrNo;
				seg{nseg}.startInd = startInd;
				seg{nseg}.endInd = endInd;
				seg{nseg}.startPos = pos_chr(startInd);
				seg{nseg}.endPos = pos_chr(endInd);
				seg{nseg}.cn = cn_prev;
				seg{nseg}.loh = loh_prev;
				seg{nseg}.nAlt = 0;
				seg{nseg}.nprobes = endInd-startInd+1;
				seg{nseg}.ts = x_chr(i-1);
				seg{nseg}.u = u_chr(i-1);
				seg{nseg}.majorcn = majorcn_prev;
				seg{nseg}.minorcn = minorcn_prev;	
				seg{nseg}.loglik = zeros(1, S);
				seg{nseg}.u_alt = zeros(1, S); 
				
				range = startInd:endInd;
				
				% compute segmental likelihood for all tumour states
				for si = 1 : S
					ind = (si-1)*U + [1:U]; % select all rows relevant to this tumour state
					[ seg{nseg}.loglik(si), u_ind ] = max( sum(loglik_chr(ind, range), 2) ); % compute sum over this region and find max. 
					seg{nseg}.u_alt(si) = params.u0 + (1-params.u0)*params.u_range(u_ind); 
				end
				
				startInd = i;

			end

		end

	end

end



