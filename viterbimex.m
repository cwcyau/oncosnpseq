function [ vpath, loglik ] = viterbimex(log_nu, log_obslike, log_transMat)
%
% vpath = viterbimex(log_nu, log_obslike, log_transMat)
%
[nStates, nSnps] = size(log_obslike);

[vpath, loglik] = viterbi_liteC(log_nu, log_transMat, log_obslike, nStates, nSnps);

vpath = vpath + 1;

