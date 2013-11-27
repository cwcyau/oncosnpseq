function [bestPath, totalCost, phi, delta] = optimMultiPath(obslike, cost, x, cost2)

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

[bestPath, phi, delta] = optimMultiPathC(obslike, nStates, nSnps, costMatrix, x, cost2);
bestPath = bestPath + 1;
totalCost = max(phi(:, end));
