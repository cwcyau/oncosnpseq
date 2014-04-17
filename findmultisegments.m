function seg = findsegments(chr, arm, pos, v, vd, x, xprev, u, loglik, options, params)

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
		vd_chr = vd(chrloc);
		x_chr = x(chrloc);
		v_chr = v(chrloc);
		u_chr = u(chrloc);
		xprev_chr = xprev(chrloc);
		loglik_chr = loglik(:, chrloc);

		if vd_chr(1) == 1
			startInd = 1;
		else
			startInd = 0;
		end

		for i = 2 : n_chr

			cn = tumourState(x_chr(i-1), 4);
			loh = tumourState(x_chr(i-1), 5);
			majorcn = tumourState(x_chr(i-1), 3);
			minorcn = tumourState(x_chr(i-1), 2);			

			% if we are not in the middle of a segment
			if startInd == 0
				% start a new segment ...
				% if the current segmentation switches from being the same as the previous segmentation to different OR 
				% if the current segmentation differs from the previous and there is a change in state
				if ( vd_chr(i-1) == 0 & vd_chr(i) == 1 ) | ( vd_chr(i-1) == 1 & vd_chr(i) == 1 & v_chr(i-1) ~= v_chr(i) ) 
					startInd = i;
				end
			end

			% if we are in the middle of a segment
			if ( startInd > 0 )
				
				% end the segment if ...
				% the current segmentation switches from being different to the same as the previous segmentation OR
				% the current segmentation is different to the previous segmentation and there is a change in state OR
				% this is the end of the chromosome arm
				if ( vd_chr(i) == 0 & vd_chr(i-1) == 1 ) | ( vd_chr(i) == 1 & vd_chr(i-1) == 1 & v_chr(i) ~= v_chr(i-1) ) | ( i == n_chr )

					endInd = i-1;

					seg{nseg}.chromosome = chrNo;
					seg{nseg}.startInd = startInd;
					seg{nseg}.endInd = endInd;
					seg{nseg}.startPos = pos_chr(startInd);
					seg{nseg}.endPos = pos_chr(endInd);
					seg{nseg}.cn = cn;
					seg{nseg}.loh = loh;
					seg{nseg}.nAlt = 0;
					seg{nseg}.nprobes = endInd - startInd + 1;
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

					% if the segment was ended because of a change in state then restart the segment
					if ( vd_chr(i) == 1 & vd_chr(i-1) == 1 & v_chr(i) ~= v_chr(i-1) )
						startInd = i;
					else
						startInd = 0;
					end
					nseg = nseg + 1;

				end
			
			end

		end

	end

end



