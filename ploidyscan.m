function params = ploidyscan(chr, arm, k, d, dd, log_pr_gg, params, options)

%
% if fast mode
%
oldparams = params;
if options.fastmode
	params.u_range = 1e-3;
	params.U = 1;
	params.p_u = ones(params.S, 1);
end



N = length(k);
S = params.S;
U0 = params.U0;
u0_range = params.u0_range;
p_u0 = params.p_u0;

tumourState = options.tumourState;
read_depth_range = options.read_depth_range;
chrRange = options.chrRange;
lambda_1_range = options.lambda_1_range;

[ arrayind, log_prior_vec, logtransMat ] = gentransmat(params.phi, params.delta, tumourState, params.u_range);

n_ri = length(read_depth_range);
n_lev = length(lambda_1_range);


%
% grid search possible base read depth/normal contamination values
%
fprintf('Grid searching possible base read depth/normal contamination values: \n');
x_ri = zeros(1, N);
u_ri = zeros(1, N);
cn_ave = zeros(n_ri, U0);
u0Cost = -Inf + zeros(n_ri, U0);
x_cost = zeros(length(chrRange), 2);

% for each base read depth value
for ri = 1 : n_ri

	params.read_depth = read_depth_range(ri);
				
	fprintf('%d ->', params.read_depth);

	% for each normal contamination level
	for u0i = 1 : U0
	
		params.u0 = u0_range(u0i);

		% compute observation likelihood
		fprintf(' %1.1f', params.u0);
		log_pr_s = calclikelihoodLite(k, d, dd, log_pr_gg, params, options);

		% segment
		for chrNo = chrRange
			for armNo = 1 : 2
				chrloc = find( chr == chrNo & arm == armNo );
				n_chr = length(chrloc);
				if n_chr > 0
					[vpath, loglik_v] = viterbimex(log_prior_vec, log_pr_s(:, chrloc), logtransMat);
					x_ri(chrloc) = arrayind(vpath, 1);			
					x_cost(chrNo, armNo) = loglik_v; 
				end
			end
		end
		
		% compute average copy number and likelihood for this base read depth/normal contamination level
		cn_ave(ri, u0i) = mean(tumourState(x_ri, 4));
		u0Cost(ri, u0i) = log(p_u0(u0i)) + sum(x_cost(:));

	end

	fprintf('\n');
	
end
fprintf('\n');

%
% write scan results to file
%
disp(['Writing ploidy scan results to: ' options.outfile_scan]);
fid = fopen(options.outfile_scan, 'wt');
fprintf(fid, 'Haploid Read Depth\tNormal Fraction\tAverage Copy Number\tLog-likelihood\n');
for ri = 1 : n_ri
	for u0i = 1 : U0
		fprintf(fid, '%d\t%1.1f\t%1.1f\t%f\n', read_depth_range(ri), u0_range(u0i), cn_ave(ri, u0i), u0Cost(ri, u0i));
	end
end
fclose(fid);


% find peaks in log-likelihood surface
ploidyvals = [];
ploidycost = [];
ploidynormal = [];
ploidycn = [];

u0Cost

maxPts = peakfinder(u0Cost);
nPeaks = size(maxPts, 1);
for i = 1 : nPeaks
	% if there is only one peak, store details else only store if the peak satisfies the max./min. ploidy requirements
	if ( cn_ave(maxPts(i, 1), maxPts(i, 2)) >= options.minploidy & cn_ave(maxPts(i, 1), maxPts(i, 2)) <= options.maxploidy ) | ( nPeaks == 1 )
		ploidyvals = [ ploidyvals read_depth_range(maxPts(i, 1)) ];
		ploidycost = [ ploidycost u0Cost(maxPts(i, 1), maxPts(i, 2)) ];	
		ploidynormal = [ ploidynormal u0_range(maxPts(i, 2)) ];
		ploidycn = [ ploidycn cn_ave(maxPts(i, 1), maxPts(i, 2)) ];	
	end
end

% sort the peaks in order of likelihood
[ ploidycost, I ] = sort(ploidycost, 2, 'descend');
ploidyvals = ploidyvals(I);
ploidynormal = ploidynormal(I);
ploidycn = ploidycn(I);

disp([ 'Haploid read depths: ' num2str(ploidyvals, '%2.0f ') ]);
disp([ 'Copy number: ' num2str(ploidycn, '%2.1f ') ]);
disp([ 'Normal fraction: ' num2str(ploidynormal, '%1.1f ') ]);
disp([ 'Log-likelihood: ' num2str(ploidycost, '%f ') ]);

n_ploidy = length(ploidyvals);

params.ploidyvals = ploidyvals;
params.ploidycost = ploidycost;
params.ploidynormal = ploidynormal;
params.ploidycn = ploidycn;
params.n_ploidy = n_ploidy;


%
% plot heatmaps
%

% heatmap of log-likelihood
figure(101); clf;

subplot(2, 2, [1 2]);
hold on;
imagesc(u0_range, read_depth_range, u0Cost);
for i = 1 : n_ploidy 
	text( ploidynormal(i), ploidyvals(i), num2str(i) );	
end
axis( [ 0.9*min( u0_range ) 1.1*max( u0_range ) min(read_depth_range)-0.5 max( read_depth_range )+0.5  ] );
ylabel(gca, 'Haploid Coverage', 'Color', 'k');
xlabel(gca, 'Normal Contamination', 'Color', 'k');
title('(a) Log-likelihood');
hnd = colorbar;
ylabel(hnd, 'Log-Likelihood');

% heatmap of average copy number
subplot(2, 2, [3 4]);
hold on;
imagesc(u0_range, read_depth_range, cn_ave);
for i = 1 : n_ploidy 
	text( ploidynormal(i), ploidyvals(i), num2str(i) );	
end
axis( [ 0.9*min( u0_range ) 1.1*max( u0_range ) min(read_depth_range)-0.5 max( read_depth_range )+0.5  ] );
ylabel(gca, 'Haploid Coverage', 'Color', 'k');
xlabel(gca, 'Normal Contamination', 'Color', 'k');
title('(b) Ploidy');
hnd = colorbar;
ylabel(hnd, 'Average Copy Number');

print('-dpsc2', '-r300', options.outfile_cost);

save(options.matfile, 'u0_range', 'read_depth_range', 'u0Cost', 'cn_ave', 'options', 'params');

if options.fastmode
	params.U = oldparams.U;
	params.u_range = oldparams.u_range;
	params.p_u = oldparams.p_u;
end


