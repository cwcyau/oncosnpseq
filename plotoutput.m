function plotoutput(chr, arm, pos, k, d, dd, x, u, params, options)

mxPlotPts = 30000;
fontSz = 6;
markSz = 0.5;


chrRange = options.chrRange;
outfile_plot = options.outfile_plot;
outfile_dat = options.outfile_dat;
n_lev = length(options.lambda_1_range);

tumourState = options.tumourState;
normalState = options.normalState;

read_depth = params.read_depth(1);

CNmax = max(tumourState(:, 4));
N = length(d);

winSz = N/1000;
d_s = nanmoving_average(dd, winSz);
b = k./d;

loc = find(rand(1, N) < 0.5);
b(loc) = 1 - b(loc);


XTickLab = [];
XTickLoc = [];
ind = 1;
for chrNo = chrRange
	for armNo = 1 : 2
		chrloc = find( chr == chrNo & arm == armNo );
		if ~isempty(chrloc)
			if armNo == 1
				XTickLoc(ind) = mean(chrloc);
				XTickLab{ind} = [ num2str(chrNo) 'p' ];
			else
				XTickLoc(ind) = mean(chrloc);
				XTickLab{ind} = [ num2str(chrNo) 'q' ];
			end
			ind = ind + 1;
		end
	end
end

% plot summary figure
disp(['Drawing summary figure: ' outfile_plot]);

figure(1); clf;

mplot = 2 + n_lev;
fontSz = 4;
set(gcf, 'Renderer', 'painters', 'PaperUnits', 'inches', 'PaperSize', [8.5 11.69 ]);
orient portrait;

pos_min = 0;
pos_max = N;

I = randperm(N);
I = I(1:min(N, mxPlotPts));
I = sort(I);

mysubplot(mplot, 1, 1);
hold on;
plot(I, dd(I), 'k.', 'markersize', 1);
plot(I, d_s(I), 'r.', 'markersize', 1);
for ci = 1 : CNmax
	line([0 N], ci*read_depth*[1 1], 'Color', 'g', 'LineWidth', 0.5);
end
set(gca, 'YTick', [0:25:200]);
set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTick', XTickLoc, 'XTickLabel', [], 'TickLength', [0.005 0.0125]);
ylabel('Read Depth', 'FontSize', fontSz);
ylim([-5 quantile(dd(I), 0.99)]);
xlim([pos_min pos_max]);
ax = axis;

for chrNo = chrRange
	for armNo = 1 : 2
		chrloc = find( chr == chrNo & arm == armNo );
		if ~isempty(chrloc)
			if armNo == 1
				col = [ 0 0 0 ];
			else
				col = 0.5*[ 1 1 1 ];
			end
			line( chrloc(1)*[1 1], [ax(3) ax(4)], 'color', col);
		end
	end
end

mysubplot(mplot, 1, 2);
hold on;
plot(I, b(I), 'k.', 'markersize', 1);
for chrNo = chrRange
	for armNo = 1 : 2
		chrloc = find( chr == chrNo & arm == armNo );
		if ~isempty(chrloc)
			if armNo == 1
				col = [ 0 0 0 ];
			else
				col = 0.5*[ 1 1 1 ];
			end		
			line( chrloc(1)*[1 1], [-CNmax CNmax], 'color', col);
		end
	end
end
ylim([-0.05 1.05]);
xlim([pos_min pos_max]);
ylabel('Lesser Allele Freq', 'FontSize', fontSz);
set(gca, 'FontSize', fontSz, 'Box', 'On', 'YTick', [0:0.2:1], 'XTick', XTickLoc, 'XTickLabel', [], 'TickLength', [0.005 0.0125]);


for lev = 1 : n_lev

	cn_vec = tumourState(x{lev}, 4);
	loh_vec = tumourState(x{lev}, 5);
	u_vec = u{lev};

	loh2_vec = loh_vec;

	loh_vec(loh_vec ~= 1) = NaN;
	loh_vec(loh_vec == 1) = -1;
	loh2_vec(loh2_vec ~= 2) = NaN;
	loh2_vec(loh2_vec == 2) = -1;

	mysubplot(mplot, 1, 2+lev);
	hold on;
	hnd = plot(I, cn_vec(I), 'ko', 'markersize', 1);
	plot(I, loh_vec(I), 'r.', 'markersize', 1);
	plot(I, loh2_vec(I), 'g.', 'markersize', 1);
	for chrNo = chrRange
		for armNo = 1 : 2
			chrloc = find( chr == chrNo & arm == armNo );
			if ~isempty(chrloc)
				if armNo == 1
					col = [ 0 0 0 ];
				else
					col = 0.5*[ 1 1 1 ];
				end			
				line( chrloc(1)*[1 1], [-CNmax CNmax], 'color', col);
			end
		end
	end
	plot(I, CNmax*u_vec(I), 'b.', 'markersize', 1);
	xlim([pos_min pos_max]);
	ylim([-1.75 CNmax+0.5]);
	ylabel('C', 'FontSize', fontSz);
	set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTick', XTickLoc, 'YTick', [-1 0:1:CNmax], 'TickLength', [0.005 0.0125]);
	if lev < n_lev
		set(gca, 'XTickLabel', []);
	end
	if lev == n_lev
		set(gca, 'XTick', XTickLoc, 'XTickLabel', XTickLab);		
		xlabel('Chromosome', 'FontSize', fontSz);
	end

end

hnd = suptitle('Overview');
set(hnd, 'FontWeight', 'bold', 'FontSize', 12);

fillPage(gcf);

print('-dpsc2', '-r300', outfile_plot);

fprintf('Drawing chromosome figure: ');
for chrNo = options.chrRange

	% plot per-chromosome figure
	fprintf('%g ', chrNo);

	figure(1); clf;

	mplot = 2 + n_lev;
	fontSz = 4;
	set(gcf, 'Renderer', 'painters', 'PaperUnits', 'inches', 'PaperSize', [8.5 11.69 ]);
	orient portrait;

	mxPlotPts = 30000;

	chrloc = find( chr == chrNo );	
	n_chr = length(chrloc);
	if n_chr == 0
		continue;
	end

	d_chr_s = d_s(chrloc);
	dd_chr = dd(chrloc);
	b_chr = b(chrloc);
	pos_chr = pos(chrloc)/1e6;

	pos_min = min(pos_chr);
	pos_max = max(pos_chr);

	I = randperm(n_chr);
	I = I(1:min(n_chr, mxPlotPts));
	I = sort(I);

	mysubplot(mplot, 1, 1);
	hold on;
	plot(pos_chr(I), dd_chr(I), 'k.', 'markersize', 1);
	plot(pos_chr(I), d_chr_s(I), 'r.', 'markersize', 1);
	for ci = 1 : CNmax
		line([pos_min pos_max], ci*read_depth*[1 1], 'Color', 'g', 'LineWidth', 0.5);
	end
	set(gca, 'YTick', [0:25:200]);
	set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTickLabel', [], 'TickLength', [0.005 0.0125]);
	ylabel('Read Depth', 'FontSize', fontSz);
	ylim([-5 quantile(dd, 0.99)]);
	xlim([pos_min pos_max]);
	ax = axis;


	mysubplot(mplot, 1, 2);
	hold on;
	plot(pos_chr(I), b_chr(I), 'k.', 'markersize', 1);
	ylim([-0.05 1.05]);
	xlim([pos_min pos_max]);
	ylabel('Lesser Allele Freq', 'FontSize', fontSz);
	set(gca, 'FontSize', fontSz, 'Box', 'On', 'YTick', [0:0.2:1], 'XTickLabel', [], 'TickLength', [0.005 0.0125]);


	for lev = 1 : n_lev

		cn_chr_vec = tumourState(x{lev}(chrloc), 4);
		loh_chr_vec = tumourState(x{lev}(chrloc), 5);
		loh2_chr_vec = loh_chr_vec;

		loh_chr_vec(loh_chr_vec ~= 1) = NaN;
		loh_chr_vec(loh_chr_vec == 1) = -1;
		loh2_chr_vec(loh2_chr_vec ~= 2) = NaN;
		loh2_chr_vec(loh2_chr_vec == 2) = -1;

		u_chr_vec = u{lev}(chrloc);

		mysubplot(mplot, 1, 2+lev);
		hold on;
		hnd = plot(pos_chr(I), cn_chr_vec(I), 'ko', 'markersize', 1);
		plot(pos_chr(I), loh_chr_vec(I), 'r.', 'markersize', 1);
		plot(pos_chr(I), loh2_chr_vec(I), 'g.', 'markersize', 1);
		plot(pos_chr(I), CNmax*u_chr_vec(I), 'b.', 'markersize', 1);
		xlim([pos_min pos_max]);
		ylim([-1.5 CNmax+0.5]);
		
		ylabel('C', 'FontSize', fontSz);
		set(gca, 'FontSize', fontSz, 'Box', 'On', 'YTick', [-1 0:1:CNmax], 'TickLength', [0.005 0.0125]);
		if lev < n_lev
			set(gca, 'XTickLabel', []);
		end
		if lev == n_lev	
			xlabel('Position / Mb', 'FontSize', fontSz);
		end

	end

	hnd = suptitle(['Chromosome: ' num2str(chrNo) ]);
	set(hnd, 'FontWeight', 'bold', 'FontSize', 12);

	fillPage(gcf);

	print('-dpsc2', '-r300', '-append', outfile_plot);

end
fprintf('\n');

%
% plot residuals
%
figure(1); clf;

mplot = 2*n_lev;
fontSz = 4;
set(gcf, 'Renderer', 'painters', 'PaperUnits', 'inches', 'PaperSize', [8.5 11.69 ]);
orient portrait;

pos_min = 0;
pos_max = N;

mxPlotPts = 30000;

I = randperm(N);
I = I(1:min(N, mxPlotPts));
I = sort(I);

seqdata = []; % data object to store residuals for later plotting
seqdata.chr = chr;
seqdata.arm = arm;
seqdata.pos = pos;
seqdata.coverage = dd;
seqdata.lesserAlleleFrac = b;

for lev = 1 : n_lev

	cn_vec = tumourState(x{lev}(:), 4);
	loh_vec = tumourState(x{lev}(:), 5);
	u_vec = u{lev}(:);

	d_s = u_vec.*normalState(x{lev}(:), 4)*read_depth + (1-u_vec).*tumourState(x{lev}(:), 4)*read_depth;
	g_s = zeros(N, 2);
	g_s(:, 1) = ( u_vec.*normalState(x{lev}(:), 1) + (1-u_vec).*tumourState(x{lev}(:), 1) )./( u_vec.*normalState(x{lev}(:), 4) + (1-u_vec).*tumourState(x{lev}(:), 4) );
	g_s(:, 2) = ( u_vec.*normalState(x{lev}(:), 2) + (1-u_vec).*tumourState(x{lev}(:), 2) )./( u_vec.*normalState(x{lev}(:), 4) + (1-u_vec).*tumourState(x{lev}(:), 4) );
	
	mysubplot(mplot, 1, (lev-1)*2+1);
	hold on;
	plot(I, dd(I), 'k.', 'markersize', 1);
	plot(I, d_s(I), 'ro', 'markersize', 1);
	for ci = 1 : CNmax
		line([0 N], ci*read_depth*[1 1], 'Color', 'g', 'LineWidth', 0.5);
	end
	set(gca, 'YTick', [0:25:200]);
	set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTick', XTickLoc, 'XTickLabel', [], 'TickLength', [0.005 0.0125]);
	ylabel('Read Depth', 'FontSize', fontSz);
	ylim([-5 quantile(dd(I), 0.99)]);
	xlim([pos_min pos_max]);
	ax = axis;
	for chrNo = chrRange
		for armNo = 1 : 2
			chrloc = find( chr == chrNo & arm == armNo );
			if ~isempty(chrloc)
				if armNo == 1
					col = [ 0 0 1 ];
				else
					col = 0.5*[ 1 1 1 ];				
				end
				line( chrloc(1)*[1 1], [ax(3) ax(4)], 'color', col);
			end
		end
	end

	mysubplot(mplot, 1, (lev-1)*2+2);
	hold on;
	plot(I, b(I), 'k.', 'markersize', 1);
	plot(I, g_s(I, 1), 'ro', 'markersize', 1);
	plot(I, 1-g_s(I, 1), 'ro', 'markersize', 1);
	plot(I, g_s(I, 2), 'ro', 'markersize', 1);	
	plot(I, 1-g_s(I, 2), 'ro', 'markersize', 1);
	for chrNo = chrRange
		for armNo = 1 : 2
			chrloc = find( chr == chrNo & arm == armNo );
			if ~isempty(chrloc)
				if armNo == 1
					col = [ 0 0 0 ];
				else
					col = 0.5*[ 1 1 1 ];				
				end			
				line( chrloc(1)*[1 1], [-10 10], 'color', col);
			end
		end
	end
	ylim([-0.05 1.05]);
	xlim([pos_min pos_max]);
	ylabel('Lesser Allele Freq', 'FontSize', fontSz);
	set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTick', XTickLoc, 'YTick', [0:0.2:1], 'XTickLabel', [], 'TickLength', [0.005 0.0125]);
	if lev == n_lev
		set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTick', XTickLoc, 'YTick', [0:0.2:1], 'XTickLabel', XTickLab, 'TickLength', [0.005 0.0125]);
	end

	%
	% store residuals for plotting
	%
	seqdata.coverage_residuals{lev} = d_s;
	seqdata.lesserAlleleFrac_residuals{lev} = g_s;
	seqdata.cn_vec{lev} = cn_vec;
	seqdata.loh_vec{lev} = loh_vec;
	seqdata.u_vec{lev} = u_vec;

end

hnd = suptitle('Model fit and residuals');
set(hnd, 'FontWeight', 'bold', 'FontSize', 12);

fillPage(gcf);

print('-dpsc2', '-r300', '-append', outfile_plot);

gzip(outfile_plot);
delete(outfile_plot);

save(outfile_dat, 'seqdata', '-v7.3');

