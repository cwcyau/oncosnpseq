function [ vpath, loglik ] = viterbimex(log_nu, log_obslike, log_transMat)
%
% vpath = viterbimex(log_nu, log_obslike, log_transMat)
%
[nStates, nSnps] = size(log_obslike);

assert(length(log_nu) == nStates)
assert(size(log_transMat, 1) == nStates)
assert(size(log_transMat, 2) == nStates)

[vpath, loglik] = viterbi_liteC(log_nu, log_transMat, log_obslike, nStates, nSnps);

vpath = vpath + 1;

