function maxPts = peakfinder(X)

[m, n] = size(X);

maxPts = [];
neigb = [ -1 1; 0 1; 1 1; -1 0; 1 0; -1 -1; 0 -1; 1 -1 ];

for i = 1 : m

	for j = 1 : n
	
		% Calculate the neighbor coordinate
		ii = i + neigb(:, 2); % row 
		jj = j + neigb(:, 1); % column

		% Check if neighbor is inside image
		ins = (ii>=1) & (jj>=1) & (jj<=n) & (ii<=m); 	

		loc = find( ins == 1 );
		ind = sub2ind([m n], ii(loc), jj(loc));

		if X(i, j) >= max(X(ind))
			maxPts = [ maxPts; [i j] ];
		end   

	end

end

if isempty(maxPts)
	[ ii, jj ] = find( X == max(X(:)) );
	maxPts(1, 1) = ii(1);
	maxPts(1, 2) = jj(1);
end

