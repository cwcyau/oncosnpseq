function seg = findsegments(chr, arm, pos, v, xd, x, xprev, u, loglik, patients, options, params)

tumourState = options.tumourState;

seg = [];

N = length(x);
S = params.S;
U = params.U;

nseg = 1;

for chrNo = options.chrRange

	for armNo = 1 : 2

		chrloc = find( chr == chrNo & arm == armNo );
		n_chr = length(chrloc);
		if n_chr == 0
			continue;
		end

		pos_chr = pos(chrloc);
		xd_chr = xd(chrloc);
		x_chr = x(chrloc);
		v_chr = v(chrloc);
		u_chr = u(chrloc);
		xprev_chr = xprev(chrloc);
		loglik_chr = loglik(:, chrloc);

		if xd_chr(1) == 1
			startInd = 1;
		else
			startInd = 0;
		end

		for i = 2 : n_chr

			cn = tumourState(x_chr(i-1), 4);
			loh = tumourState(x_chr(i-1), 5);
			majorcn = tumourState(x_chr(i-1), 3);
			minorcn = tumourState(x_chr(i-1), 2);			

			if xd_chr(i-1) == 0 & xd_chr(i) == 1
				startInd = i;
			end

			if ( xd_chr(i) == 0 & xd_chr(i-1) == 1 ) | ( ( startInd > 0 ) & ( i == n_chr ) )

				endInd = i-1;

				seg{nseg}.chromosome = chrNo;
				seg{nseg}.startInd = startInd;
				seg{nseg}.endInd = endInd;
				seg{nseg}.startPos = pos_chr(startInd);
				seg{nseg}.endPos = pos_chr(endInd);
				seg{nseg}.cn = cn;
				seg{nseg}.loh = loh;
				seg{nseg}.nAlt = 0;
				seg{nseg}.nprobes = endInd-startInd+1;
				seg{nseg}.ts = x_chr(i-1);
				seg{nseg}.u = u_chr(i-1);
				seg{nseg}.majorcn = majorcn;
				seg{nseg}.minorcn = minorcn;	
				seg{nseg}.loglik = zeros(1, S);

				range = startInd:endInd;
				for si = 1 : S
					ind = (si-1)*U + [1:U];
					[ seg{nseg}.loglik(si), u_ind ] = max( sum(loglik_chr(ind, range), 2) ); % compute sum over this region and find max. 
					seg{nseg}.u_alt(si) = params.u0 + (1-params.u0)*params.u_range(u_ind); 
				end

				startInd = 0;
				nseg = nseg + 1;

			end

		end

	end

end



