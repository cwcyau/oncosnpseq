function [chr, arm, pos, k, d, dd, log_pr_gg] = loaddataCG(options)

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

[ chr, pos, var1, var2, Ta, Td, Na, Nd ] = textread(tmpfile, '%n %n %n %n %n %n %n %n', 'headerlines', 1);

if strfind(ext, 'gz')
	disp(['Removing temporary file: ' tmpfile]);
	delete(tmpfile);
end
if strfind(ext, 'zip')
	disp(['Removing temporary file: ' tmpfile]);
	delete(tmpfile);
end

disp('Removing low-quality loci ...');

goodloc = find(var1 > 0 & var2 > 0 & ( Td > Ta ) & ( Nd > Na ) );
chr = chr(goodloc);
pos = pos(goodloc);
var1 = var1(goodloc);
var2 = var2(goodloc);
Ta = Ta(goodloc);
Td = Td(goodloc);
Na = Na(goodloc);
Nd = Nd(goodloc);

dd = max(0, Td); 
d = max(0, Td);
k = max(0, Ta);
kn = max(0, Na);
dn = max(0, Nd);

N = length(d);
%loc = find( rand(1, N) < 0.5 );
%k(loc) = d(loc)-k(loc);
k = min(k, d-k);
%kn(loc) = dn(loc)-kn(loc);
kn = min(kn, dn-kn);

chrloc = [];
for chrNo = options.chrRange
	loc = find( chr == chrNo );
	chrloc = [ chrloc; loc ];
end

chr = chr(chrloc);
pos = pos(chrloc);
d = d(chrloc);
dd = dd(chrloc);
dn = dn(chrloc);
kn = kn(chrloc);
k = k(chrloc);
N = length(d);
arm = ones(N, 1);

% assign arm numbers (1 -p, 2 - q)
disp(['Reading hg tables file: ' options.hgtables]);
[ chrArm, chrStart, chrEnd, armNo ] = textread(options.hgtables, '%n %n %n %n', 'headerlines', 1);
for ci = 1 : length(chrArm)

	if ci == 1
		loc = find( chr == chrArm(ci) & pos <= chrEnd(ci) );
	else
		loc = find( chr == chrArm(ci) & pos >= chrStart(ci) );
	end

	if ~isempty(loc)
		arm(loc) = armNo(ci);
	end

end

if ~isempty(options.gcdir)

	fprintf('Doing Local GC content: ');
	for chrNo = options.chrRange

		fprintf('%g ', chrNo);

		chrloc = find( chr == chrNo );
		n_snps = length(chrloc);
	
		if n_snps > 0

			k_chr = k(chrloc);
			d_chr = d(chrloc);
			dn_chr = dn(chrloc); 
			dn_chr = dn_chr - mean(dn_chr);
			pos_chr = pos(chrloc);
	
			if chrNo < 23
				gcfile = [ options.gcdir '/' num2str(chrNo) '_1k.txt' ];
			end
			if chrNo == 23
				gcfile = [ options.gcdir '/X_1k.txt' ];
			end
			if chrNo == 24
				gcfile = [ options.gcdir '/Y_1k.txt' ];
			end
	
			[gcstart, gcend, gc] = textread(gcfile, '%n %n %n');
			gcpos = 0.5*(gcstart + gcend);
			gc = gc - mean(gc);
	
			gc_chr = interp1(gcpos, gc, pos_chr, 'linear', 'extrap');
		
			betas = robustfit([ gc_chr dn_chr], d_chr);

			k_chr = k_chr - (k_chr./d_chr).*( betas(2).*gc_chr + betas(3).*dn_chr );			
			k(chrloc) = k_chr;

			d_chr = d_chr - betas(2).*gc_chr - betas(3).*dn_chr;			
			d(chrloc) = d_chr;
	
		end
	
	end	

	fprintf('\n');

end

k = max(0, k);
d = max(0, d);
kn = max(0, kn);
dn = max(0, dn);


k = reshape(k, [1 N]);
d = reshape(d, [1 N]);
dd = d;
kn = reshape(kn, [1 N]);
dn = reshape(dn, [1 N]);

% do germline genotyping
disp('Performing germline genotyping ...');
alpha = options.alpha;
beta = options.beta;

log_pr_gg = zeros(2, N);
for gi = 1 : 2
	log_pr_gg(gi, :) = gammaln(dn+1) - gammaln(kn+1) - gammaln(dn-kn+1) + gammaln(kn+alpha(gi)) + gammaln(dn-kn+beta(gi)) - gammaln(dn+alpha(gi)+beta(gi)) + gammaln(alpha(gi)+beta(gi)) - gammaln(alpha(gi))-gammaln(beta(gi));
end

%b = kn./dn;

%b([1:15]+1000)
%log_pr_gg(:, [1:15]+1000)

%figure(100); clf;

%set(gcf, 'Renderer', 'Painters');

%col{1} = 'r';
%col{2} = 'g';
%col{3} = 'y';
%col{4} = 'b';

%n = length(b);
%ind = randperm(n);
%ind = ind(1:30000);
%ind = sort(ind);

%mysubplot(2, 1, 1);
%hold on;
%plot( kn(ind)./dn(ind), 'k.', 'markersize', 1);

%mysubplot(2, 1, 2);
%hold on;
%plot( b(ind), 'k.', 'markersize', 1);
%[ mxval, mxloc ] = max(log_pr_gg(:, ind), [], 1);
%for gi = 1 : 4
%	plot(find(mxloc==gi), b(ind(mxloc==gi)), '.', 'color', col{gi}, 'markersize', 2);
%end
%ylim([-0.05 1.05]);


%print -dpsc2 -r150 gg.ps;
%gzip gg.ps;
%delete gg.ps;

