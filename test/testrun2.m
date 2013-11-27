addpath('../');

mex -outdir ../ ../lshmm_func.c 
mex -outdir ../ ../lshmm_multi_func.c 

%% default settings
options.maxCopy = 8; %% maximum copy number
options.read_depth_range = [10:40]; % range of haploid read coverages to scan
options.chrRange = [1:22]; % number of chromosomes to process
options.N_train = 15000; % number of training data points to use
options.maxploidy = 4.5; % maximum ploidy (max average copy number)
options.minploidy = 1.5; % minimum ploidy (min average copy number)
options.normalcontamination = 1;
options.tumourheterogeneity = 0;
options.u_levels = 6;
options.maxnormalcontamination = 0.5;
options.lambda_1_range = [ 200 100 50 30 10 ];
options.lambda_2 = 1;
options.training = 0;
options.read_error = 0.01;
options.seq_error = 0.01;
options.u_alpha = 2;
options.u0_alpha = 0.5;
options.u0_beta = 5;

options.lambda = linspace(50, 10, 5);
options.psi = 2;
options.eta = 3;
options.nu = linspace(50, 10, 5);
options.delta = 1;
options.karyolibfile = '/data/suzaku/yau/cancercopy/data/tcga_ludwig.dat';

options.tumourStateTable = '../config/tumourStatesOncoSNPSEQ.txt';
options.gcdir = '../gc/b37/';
options.mapdir = [];
options.hgtables = '../config/hgTables_b37.txt';

options.usekaryolib = 1;
options.outdir = '/data/suzaku/yau/cancercopy/output/withlib/';

for fi = 1 : 10

	for vi = 1 : 6

		options.samplename = [ 'testex_' num2str(fi) '_' num2str(vi) ];
		disp(options.samplename);

		% set up parameters
		[params, options] = setup(options);

		load(['/data/suzaku/yau/cancercopy/data/testex_' num2str(fi) '_' num2str(vi) '.mat' ]);

		options.outfile = [ options.outdir '/' options.samplename '.cnvs' ];
		options.outfile_all = [ options.outdir '/' options.samplename '.cnvall' ];
		options.outfile_cost = [ options.outdir '/' options.samplename '.cost.ps' ];	
		options.matfile = [ options.outdir '/' options.samplename '.mat' ];	
		options.outfile_qc = [ options.outdir '/' options.samplename '.qc' ];	

		chr2 = [];
		pos2 = [];
		k2 = [];
		d2 = [];
		for chrNo = options.chrRange
			loc = find(chr == chrNo);
			chr2 = [ chr2 chrNo*ones(1, length(loc)) ];
			pos2 = [ pos2 pos(loc) ];
			k2 = [ k2; k(loc) ];
			d2 = [ d2; d(loc) ];	
		end
		chr = chr2;
		pos = pos2;
		k = k2';
		d = d2';
		k = min(k, d);
		k = min(k, d-k);

		arm = ones(size(chr));
		dd = d;

		% read in database
		karyolib = readkaryolib(chr, pos, options);

		% compile a random subset of training data
		N_train = options.N_train;
		tumourState = options.tumourState;

		CNmax = max(tumourState(:, 4));
		d_s = nanmoving_average(d, 30);

		params.lambda_s = mad( d - d_s );

		S = params.S;
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

		if options.training == 1
			N = N_train;
			chr = chr_train;
			pos = pos_train;
			arm = arm_train;
			d = d_train;
			dd = dd_train;
			k = k_train;
		end


		% scan different ploidy/contamination values on training dataset
		params = ploidyscan(chr_train, arm_train, k_train, d_train, dd_train, params, options);

		% output QC metrics
		outputQC(params, options);


		%save testdat options tumourState chr_train pos_train arm_train k_train d_train params;

		% do segmentation for each ploidy/contamination configuration
		for i = 1 : params.n_ploidy

			% set parameters for this configuration
			for li = 1 : length(options.lambda_1_range)
				x{li} = zeros(1, N);
				u{li} = zeros(1, N);
				seg{li} = {};
			end
			segall = {};
			params.read_depth = params.ploidyvals(i);
			params.u0 = params.ploidynormal(i);

			options.outfile_plot = [ options.outdir '/' options.samplename '.' num2str(i) '.ps' ];
	
			if options.usekaryolib

				disp('Calculating observation likelihood ...');
				tic
				log_pr_s = calclikelihoodLite(k, d, dd, params, options);			
				toc

				disp('Segmenting ...');

				[ x, seg, u, segall, patients ] = segment_karyolib(chr, arm, pos, k, d, dd, log_pr_s, karyolib, params, options);

			else

				% scan each chromosome and arm
				fprintf('(%d): ', i);
				log_pr_s = zeros(S, N);
				for chrNo = options.chrRange

					fprintf('%d ', chrNo);
		
					for armNo = 1 : 2
	
						chrloc = find( chr == chrNo & arm == armNo );
						n_chr = length(chrloc);
		
						if n_chr > 0
		
							% calculate observation likelihood
							log_pr_s(:, chrloc)  = calclikelihoodLite(k(chrloc), d(chrloc), dd(chrloc), params, options);			
			
							% do segmentation
							[ x_chr, seg_chr, u_chr, segall_chr ] = segment(chr(chrloc), arm(chrloc), pos(chrloc), k(chrloc), d(chrloc), dd(chrloc), log_pr_s(:, chrloc), params, options);
			
							% compile segmentation results into array
							for li = 1 : length(x_chr)
								x{li}(chrloc) = x_chr{li};
								u{li}(chrloc) = u_chr{li};
								if ~isempty(seg_chr{li})
									seg{li} = { seg{li}{:} seg_chr{li}{:} };
								end
							end

							segall = { segall{:} segall_chr{:} };

						end
		
					end
		
				end
				fprintf('\n');

			end

			% write segmentation results to files
			output2file(i, seg, segall, params, options);
	
			% plot output
			plotoutput(chr, arm, pos, k, d, dd, x, u, params, options);
	
		end

	end

end

