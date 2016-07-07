clear all;
close all;
clc;

addpath('../');
addpath('../external/');

rand('state', 1);
randn('state', 1);

options.hgtables = '../config/hgTables_b37.txt'; % human genome annotation table
options.gcdir = '/data/cyau/software/gc/b37/'; % directory of local GC content files
options.outdir = '/data/cyau/ahmed/output/'; % output directory
options.tumourStateTable = '../config/tumourStatesSimple.txt';

options.seqtype = 'lfr'; % sequencing type 'cg' or 'illumina'

seqerror = 0.05;
readerror = 0.05;

options.samplename = 'primary-all'; % sample name
options.infile = '/data/cyau/ahmed/summary/primary-all.consensus.oncoseq'; % input data file

%options.samplename = 'OP1019-consensus'; % sample name
%options.infile = '/data/cyau/ahmed/OP1019/OP1019-consensus.oncoseq'; % input data file


%options.samplename = 'OP1019-preB'; % sample name
%options.infile = '/data/cyau/ahmed/OP1019/OP1019-preB.oncoseq'; % input data file

%options.samplename = 'OP1019-postB'; % sample name
%options.infile = '/data/cyau/ahmed/OP1019/OP1019-postB.oncoseq'; % input data file


tic
oncoseq( 	'--read_depth_range', '[5:15]', ...
			'--chr_range', '[1:23]', ...
			'--n_train', '60000', ...
			'--maxploidy', '3.5', ...
			'--minploidy', '1.5', ...
			'--normalcontamination', ...
			'--maxnormalcontamination', '0.5', ...
			'--lambda1', '[ 5000 1000 500 100 50 ]', ...
			'--seqerror', num2str(seqerror), ...
			'--readerror', num2str(readerror), ...
			'--tumourstatestable', options.tumourStateTable, ...
			'--hgtable', options.hgtables, ...
			'--gcdir', options.gcdir, ...
			'--seqtype', options.seqtype, ...
			'--samplename', options.samplename, ...
			'--infile', options.infile, ...
			'--outdir', options.outdir ...
		);
toc




