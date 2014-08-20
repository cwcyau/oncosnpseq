clear all;
close all;
clc;

addpath('../');
addpath('../external/');

rand('state', 1);
randn('state', 1);

options.hgtables = '../config/hgTables_b37.txt'; % human genome annotation table
options.gcdir = '/data/cyau/software/gc/b37/'; % directory of local GC content files
options.outdir = '/data/cyau/temp/'; % output directory
options.tumourStateTable = '../config/tumourStates.txt';

options.samplename = 'maria';
options.infile = '/home/cyau/Dropbox/maria/input.output.forCYau/snp_SS6003119'; % input data file
options.normalfile = '/home/cyau/Dropbox/maria/input.output.forCYau/snp_SS6003118'; % input data file

[chr, pos, ref, a, c, g, t] = textread(options.infile, '%n %n %n %n %n %n %n %*[^\n]', 'delimiter', ',', 'headerlines', 1);
[chr_n, pos_n, ref_n, a_n, c_n, g_n, t_n] = textread(options.normalfile, '%n %n %n %n %n %n %n %*[^\n]', 'delimiter', ',', 'headerlines', 1);


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


dd = ref + max([a c g t], [], 2);
dn = ref_n + max([a_n c_n g_n t_n], [], 2);

figure(1); clf;
set(gcf, 'Renderer', 'Painters');
plot(dd, dn, 'k.');

print -dpsc2 -r300 temp2.ps
gzip temp2.ps
delete temp2.ps


