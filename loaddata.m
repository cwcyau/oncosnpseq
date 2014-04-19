function [ chr, arm, pos, k, d, dd, log_pr_gg ] = loaddata(options, params)

disp(['Reading data file: ' options.infile]);
[pathstr, name, ext] = fileparts(options.infile);

tmpfile = options.infile;
if strfind(ext, 'gz')
	disp(['Found gzipped file: ' options.infile]);
	disp(['Gunzipping to: ' options.outdir]);
	gunzip(options.infile, options.outdir);		
	tmpfile = fullfile(options.outdir, name);
end
if strfind(ext, 'zip')
	disp(['Found zipped file: ' options.infile]);
	disp(['Unzipping to: ' options.outdir]);
	unzip(options.infile, options.outdir);		
	tmpfile = fullfile(options.outdir, name );
end
[chr, pos, ref, a, c, g, t] = textread(tmpfile, '%n %n %n %n %n %n %n %*[^\n]', 'delimiter', ',', 'headerlines', 1);

if strfind(ext, 'gz')
	disp(['Removing temporary file: ' tmpfile]);
	delete(tmpfile);
end
if strfind(ext, 'zip')
	disp(['Removing temporary file: ' tmpfile]);
	delete(tmpfile);
end

% if doing paired analysis, load normal data
if options.paired

	disp(['Reading data file: ' options.normalfile]);
	[pathstr, name, ext] = fileparts(options.normalfile);

	tmpfile = options.infile;
	if strfind(ext, 'gz')
		disp(['Found gzipped file: ' options.normalfile]);
		disp(['Gunzipping to: ' options.outdir]);
		gunzip(options.normalfile, options.outdir);		
		tmpfile = fullfile(options.outdir, name);
	end
	if strfind(ext, 'zip')
		disp(['Found zipped file: ' options.normalfile]);
		disp(['Unzipping to: ' options.outdir]);
		unzip(options.normalfile, options.outdir);		
		tmpfile = fullfile(options.outdir, name );
	end
	[chr_n, pos_n, ref_n, a_n, c_n, g_n, t_n] = textread(tmpfile, '%n %n %n %n %n %n %n %*[^\n]', 'delimiter', ',', 'headerlines', 1);

	if strfind(ext, 'gz')
		disp(['Removing temporary file: ' tmpfile]);
		delete(tmpfile);
	end
	if strfind(ext, 'zip')
		disp(['Removing temporary file: ' tmpfile]);
		delete(tmpfile);
	end

else

	chr_n = chr;
	pos_n = pos;
	ref_n = 0*ref;
	a_n = 0*a;
	c_n = 0*c;
	g_n = 0*g;
	t_n = 0*t;

end

[C, IA, IB] = intersect( [chr pos], [chr_n pos_n], 'rows' );

chr = chr(IA);
pos = pos(IA);
a = a(IA);
c = c(IA);
g = g(IA);
t = t(IA);
ref = ref(IA);

chr_n = chr_n(IB);
pos_n = pos_n(IB);
a_n = a_n(IB);
c_n = c_n(IB);
g_n = g_n(IB);
t_n = t_n(IB);
ref_n = ref_n(IB);

X = [ a c g t ];
k = max(X, [], 2);
d = ref + k;
dd = ref + k;
k = min(k, d-k);

Xn = [ a_n c_n g_n t_n ];
kn = max(Xn, [], 2);
dn = ref_n + kn;
kn = min(kn, dn-kn);

N = length(d);

loc = find( rand(1, N) < 0.5 );
k(loc) = d(loc)-k(loc);
kn(loc) = dn(loc)-kn(loc);


chrloc = [];
for chrNo = options.chrRange
	loc = find( chr == chrNo );
	chrloc = [ chrloc; loc ];
end

chr = chr(chrloc);
pos = pos(chrloc);
d = d(chrloc);
dd = dd(chrloc);
k = k(chrloc);
kn = kn(chrloc);
dn = dn(chrloc);
N = length(d);
arm = zeros(N, 1);


% assign arm numbers (1 -p, 2 - q)
disp(['Reading hg tables file: ' options.hgtables]);
[ chrArm, chrStart, chrEnd, armNo ] = textread(options.hgtables, '%n %n %n %n', 'headerlines', 1);
for ci = 1 : length(chrArm)

	loc = find( chr == chrArm(ci) & pos >= chrStart(ci) & pos < chrEnd(ci) );
	
	if ~isempty(loc)
		arm(loc) = armNo(ci);
	end

end

%
% GC and MAPPABILITY CORRECTION
%

if ~isempty(options.gcdir) | ~isempty(options.mapdir)

	fprintf('Doing Local GC content and Mappability Correction: ');
	for chrNo = options.chrRange

		fprintf('%g ', chrNo);

		chrloc = find( chr == chrNo );
		if isempty(chrloc)
			continue;
		end

		dd_chr = dd(chrloc);		
		if options.paired
			dn_chr = dn(chrloc) - mean(dn(chrloc));
		else
			dn_chr = ones(size(chrloc));
		end
		k_chr = k(chrloc);		
		pos_chr = pos(chrloc);

		% load gc content data
		if ~isempty(options.gcdir)

			if chrNo == 23
				gcfile = [ options.gcdir '/X_1k.txt' ];
			end
			if chrNo == 24
				gcfile = [ options.gcdir '/Y_1k.txt' ];
			end
			if chrNo < 23
				gcfile = [ options.gcdir '/' num2str(chrNo) '_1k.txt' ];
			end
			[gc_start, gc_end, gc_content] = textread(gcfile, '%n %n %n', 'headerlines', 0);
			gc_pos = 0.5*(gc_start + gc_end);

			% sort by genomic position
			[gc_pos, I] = sort(gc_pos);
			gc_content = gc_content(I);
			gc_content = (gc_content - nanmean(gc_content))./nanstd(gc_content);

			% interpolate	
			gc_chr = interp1(gc_pos, gc_content, pos_chr, 'nearest', 'extrap')';
			gc_chr = reshape(gc_chr, size(dd_chr));   

		else

			gc_chr = ones(size(d_chr));

		end 
	
		% load mappability data
		if ~isempty(options.mapdir)

			if chrNo == 23
				mapfile = [ options.mapdir '/X.txt' ];
			end
			if chrNo == 24
				mapfile = [ options.mapdir '/Y.txt' ];
			end
			if chrNo < 23
				mapfile = [ options.mapdir '/' num2str(chrNo) '.txt' ];
			end

			[map_chr, map_pos, map] = textread( mapfile,  '%n %n %n', 'headerlines', 1);

			% sort by genomic position
			[map_pos, I] = sort(map_pos);
			map_content = map_content(I);
			map_content = (map_content - nanmean(map_content))./nanstd(map_content);

			% interpolate		
			map_chr = interp1(map_pos, map_content, pos_chr, 'nearest', 'extrap')';
			map_chr = reshape(map_chr, size(d_chr));  

		else
		
			map_chr = ones(size(dd_chr));

		end

		nonzeroloc = find( dd_chr > 0 );

		if ~isempty(options.gcdir) & ~isempty(options.mapdir)
			betas = robustfit([ gc_chr map_chr dn_chr ], dd_chr);
			dd_chr = dd_chr - betas(2).*gc_chr - betas(3).*map_chr - betas(4).*dn_chr;
		end
		if ~isempty(options.gcdir) & isempty(options.mapdir)
			betas = robustfit([ gc_chr dn_chr ], dd_chr);
			dd_chr = dd_chr - betas(2).*gc_chr - betas(3).*dn_chr;
		end
		if isempty(options.gcdir) & ~isempty(options.mapdir)
			betas = robustfit([ map_chr dn_chr ], dd_chr);
			dd_chr = dd_chr - betas(2).*map_chr - betas(3).*dn_chr;
		end

		dd(chrloc) = dd_chr;

	end
	fprintf('\n');

end

k(k < 0) = 0;
d(d < 0) = 0;
dd(dd < 0) = 0;
dd = round(dd);
kn(kn < 0) = 0;
dn(dn < 0) = 0;

k = reshape(k, [1 N]);
d = reshape(d, [1 N]);
dd = reshape(dd, [1 N]);
kn = reshape(kn, [1 N]);
dn = reshape(dn, [1 N]);

% do germline genotyping
disp('Performing germline genotyping ...');
log_pr_gg = zeros(params.G, N);
if options.paired
	alpha = options.alpha;
	beta = options.beta;
	for gi = 1 : params.G
		log_pr_gg(gi, :) = gammaln(dn+1) - gammaln(kn+1) - gammaln(dn-kn+1) + gammaln(kn+alpha(gi)) + gammaln(dn-kn+beta(gi)) - gammaln(dn+alpha(gi)+beta(gi)) + gammaln(alpha(gi)+beta(gi)) - gammaln(alpha(gi))-gammaln(beta(gi));
	end
end

