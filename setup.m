function [ params, options ] = setup(options)
%
% set up parameter array and tumour states
%

if isempty(options.tumourStateTable)
	disp('Generating default tumour state table.');
	gentumourstatetable;
else
	disp([ 'Using user-specified tumour state table in: ' options.tumourStateTable ]);
	if exist(options.tumourStateTable, 'file')
		[ allele1, allele2, allele3, allele4, loh ] = textread(options.tumourStateTable, '%n %n %n %n %n', 'headerlines', 1);
		tumourState = zeros(length(loh), 5);
		tumourState(:, 1) = allele1;
		tumourState(:, 2) = allele2;
		tumourState(:, 3) = allele3;
		tumourState(:, 4) = allele4;
		tumourState(:, 5) = loh;

		% if using paired analysis, can remove germline loh states
		if options.paired
			loc = find( loh ~= 2 );
			tumourState = tumourState(loc, :);
		end
	else
		disp(['Error. The file: ' options.tumourStateTable ' does not exist.']);
		return;
	end
end
S = size(tumourState, 1);

for si = 1: S
	if tumourState(si, 5) == 0 % non-loh
		pg(si, :) = [ 1/3 1/6 1/6 1/3 ];
		normalState(si, :) = [ 0 1 1 2 ];
		genotypeState(si, :) = [ 1 2 3 4 ];
	end
	if tumourState(si, 5) == 1 % somatic loh
		normalState(si, :) = [ 0 1 1 2 ];
		genotypeState(si, :) = [ 1 2 3 4 ];
		pg(si, :) = [ 1/3 1/6 1/6 1/3 ];
	end
	if tumourState(si, 5) == 2 % germline loh
		normalState(si, :) = [ 0 0 2 2 ];
		genotypeState(si, :) = [ 1 1 4 4 ];
		pg(si, :) = [ 1/4 1/4 1/4 1/4 ];
	end	
end


K = 1;
G = 4;
U0 = 1;
U = 1;

if options.normalcontamination == 1
	u0_range = max(1e-3, options.u0_levels);
else
	u0_range = 1e-3;
end
if options.tumourheterogeneity == 1
	u_range = max(1e-3, options.u_levels);	
else
	u_range = 1e-3;
end
U0 = length(u0_range);
U = length(u_range);

% priors on tumour heterogeneity
p_u(1) = options.u_alpha;
p_u(2:U) = 1;
p_u = betapdf(u_range, options.u_alpha, options.u_beta);
p_u = p_u./sum(p_u);
p_u = repmat(p_u, [S 1]);

% priors on normal contamination
p_u0 = betapdf(u0_range, options.u0_alpha, options.u0_beta);
p_u0 = p_u0./sum(p_u0);

params.read_error = options.read_error;
params.seq_error = options.seq_error;
params.read_depth = 20;
params.u0 = 1e-3;
params.nu = 4;
params.u_range = u_range;
params.p_u = p_u;
params.S = S;
params.K = K;
params.G = G;
params.U = U;
params.U0 = U0;
params.u0_range = u0_range;
params.p_u0 = p_u0;
params.phi = -10;
params.delta = -20;

options.tumourState = tumourState;
options.normalState = normalState;
options.genotypeState = genotypeState;
options.pg = pg;

