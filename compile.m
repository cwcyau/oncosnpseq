function compile()


rootdir = pwd;

destdir = 'executables/';

mkdir(destdir);

mcc -v -R -nodisplay -m oncoseq.m -d executables/ -o oncoseq
