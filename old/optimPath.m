function [bestPath, totalCost, phi, delta] = optimPath(obslike, cost)

if nargin < 2
	cost = 10;
end

[nStates, nSnps] = size(obslike);

if isscalar(cost) == 1
	costMatrix = cost*ones(nStates, nStates);
	for i = 1 : nStates
		costMatrix(i, i) = 0;
	end
else
	costMatrix = cost;
end

[bestPath, phi, delta] = optimPathC(obslike, nStates, nSnps, costMatrix);
bestPath = bestPath + 1;
totalCost = max(phi(:, end));
