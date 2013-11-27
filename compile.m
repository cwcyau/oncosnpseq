function compile()

rootdir = pwd;
archType = getenv('ARCH');

destdir = 'executables/';

mkdir(destdir);

file{1} = 'oncoseq.m'
file{2} = 'oncoseq_run.m'

count = 3;

extdir = [];
extdir{1} = 'external/';

indir{1} = '.';
ndir = length(indir);


for di = 1 : ndir
	dirContents = dir(indir{di});
	nFiles = length(dirContents);
	for fi = 3 : nFiles
		if ( dirContents(fi).isdir == 0 )
			filename = dirContents(fi).name;
			if ~isempty(strfind(filename(end-1:end), '.m')) || ~isempty(strfind(filename(end-5:end), '.mexa64'))
 				file{count} = [ rootdir '/' indir{di} '/' filename ];
 				count = count + 1;
			end
		end
	end
end

nextdir = length(extdir);
for di = 1 : nextdir
	dirContents = dir(extdir{di});
	nFiles = length(dirContents);
	for fi = 3 : nFiles
		if ( dirContents(fi).isdir == 0 )
			filename = dirContents(fi).name;
			if ~isempty(strfind(filename(end-1:end), '.m')) || ~isempty(strfind(filename(end-5:end), '.mexa64'))
 				file{count} = [ rootdir '/' extdir{di} '/' filename ];
 				count = count + 1;
			end
		end
	end
end


str = [ 'mcc -v -R -nodisplay -m' ];
for i = 1 : count-1
     str = [ str ' ' file{i} ];
end
str = [ str ' -d ' destdir ' -o oncoseq ' ];

str
eval(str);



