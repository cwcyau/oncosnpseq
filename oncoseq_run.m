function oncoseq_run(options)

% set up parameters
[params, options] = setup(options);

% load data
switch options.seqtype
	case 'illumina'
		disp('Using Illumina data mode.');
		[chr, arm, pos, k, d, dd, log_pr_gg] = loaddata(options, params);
	case 'cg'
		disp('Using Complete Genomics data mode.');
		[chr, arm, pos, k, d, dd, log_pr_gg] = loaddataCG(options, params);
	case 'lfr'
		disp('Using LFR data mode.');
		[chr, arm, pos, k, d, dd, log_pr_gg] = loaddataLFR(options, params);
end

% compile a random subset of training data
N_train = options.N_train;
tumourState = options.tumourState;

disp(['Number of tumour states: ' num2str(size(tumourState, 1))]);

CNmax = max(tumourState(:, 4));
d_s = nanmoving_average(d, 30);

params.lambda_s = mad( d - d_s, 1 );
params.lambda_s = std( d - d_s );

S = params.S;
U = params.U;
N = length(d);
loc = randperm(N);
N_train = min(N_train, N);
loc = loc(1:N_train);

loc = sort(loc);
chr_train = chr(loc);
arm_train = arm(loc);
pos_train = pos(loc);
d_train = d(loc);
dd_train = dd(loc);
k_train = k(loc);
log_pr_gg_train = log_pr_gg(:, loc);

if options.training == 1
	N = N_train;
	chr = chr_train;
	pos = pos_train;
	arm = arm_train;
	d = d_train;
	dd = dd_train;
	k = k_train;
	log_pr_gg = log_pr_gg_train;
end

% scan different ploidy/contamination values on training dataset
if options.fixed_range == 0
	options.read_depth_range = floor(median(dd)/options.maxploidy):1:ceil(median(dd)/options.minploidy);
end
params = ploidyscan(chr_train, arm_train, k_train, d_train, dd_train, log_pr_gg_train, params, options);

% output QC metrics
outputQC(params, options);

% do segmentation for each ploidy/contamination configuration
for i = 1 : params.n_ploidy

	% set parameters for this configuration
	for li = 1 : length(options.lambda_1_range)
		x{li} = zeros(1, N);
		u{li} = zeros(1, N);
		seg{li} = {};
	end
	segall = [];
	params.read_depth = params.ploidyvals(i);
	params.u0 = params.ploidynormal(i);

	options.outfile_plot = [ options.outdir '/' options.samplename '.' num2str(i) '.ps' ];

	% scan each chromosome and arm
	fprintf('Segmenting chromosome (Configuration %d): ', i);
	
	log_pr_su = zeros(S*U, N);
	log_pr_s = zeros(S, N);
	
	for chrNo = options.chrRange

		fprintf('%d ', chrNo);
	
		for armNo = 1 : 2

			chrloc = find( chr == chrNo & arm == armNo );
			n_chr = length(chrloc);
	
			if n_chr > 0
	
				% calculate observation likelihood
				[ log_pr_su(:, chrloc), log_pr_s(:, chrloc) ] = calclikelihoodLite(k(chrloc), d(chrloc), dd(chrloc), log_pr_gg(:, chrloc), params, options);	
				
				% do segmentation
				[ x_chr, seg_chr, u_chr, segall_chr ] = segment(chr(chrloc), arm(chrloc), pos(chrloc), k(chrloc), d(chrloc), dd(chrloc), log_pr_gg(:, chrloc), log_pr_su(:, chrloc), log_pr_s(:, chrloc), params, options);
				
				% compile segmentation results into array
				for li = 1 : length(x_chr)
					x{li}(chrloc) = x_chr{li};
					u{li}(chrloc) = u_chr{li};
					if ~isempty(seg_chr{li})
						seg{li} = { seg{li}{:} seg_chr{li}{:} };
					end
				end

				if ~isempty(segall_chr)
					if ~isempty(segall) 
						segall = { segall{:} segall_chr{:} };
					else
						segall = { segall_chr{:} };
					end
				end

			end
	
		end
	
	end
	fprintf('\n');

	% write segmentation results to files
	output2file(i, seg, segall, params, options);
	
	% plot output
	plotoutput(chr, arm, pos, k, d, dd, x, u, params, options);

	% dump diagnostic info
	if ~isdeployed
		if options.diagnostics
			if i == 1
				disp('Saving diagnostics information...');
				save(options.outfile_diagnostics, '-v7.3', 'chr', 'arm', 'pos', 'k', 'd', 'dd', 'x', 'u', 'seg', 'segall', 'log_pr_su', 'log_pr_s', 'params', 'options');
			end
		end
	end
	
end


