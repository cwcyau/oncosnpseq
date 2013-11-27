% generate tumour table
tumourState(1, :) = [ 0 0 0 0 1 ];
tumourState(2, :) = [ 0 0 1 1 1 ];
tumourState(3, :) = [ 0 1 1 2 0 ];
tumourState(4, :) = [ 0 0 2 2 1 ];
for ci = 3 : options.maxCopy

	for bi = 0 : ci
		if bi <= ci-bi
			if bi == 0 
				loh = 1;
			else
				loh = 0;
			end
			tumourState = [ tumourState; [ 0 bi ci-bi ci loh ] ];
		end
	end	
	
end
