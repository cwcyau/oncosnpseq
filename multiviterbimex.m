function [ vpath, loglik ] = multiviterbimex(log_nu, log_obslike, log_transMat, prev_path, penalty)
%
% vpath = viterbimex(log_nu, log_obslike, log_transMat)
%
[nStates, nSnps] = size(log_obslike);

[vpath, loglik] = multiviterbi_liteC(log_nu, log_transMat, log_obslike, prev_path, penalty, nStates, nSnps);

vpath = vpath + 1;

