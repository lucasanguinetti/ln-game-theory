% efficiencyFunction  Compute the throughput given SINR and number of
%                     information bits
%
% inputs:        in (SINR)
%                L (number of information bits per packet)
%
% outputs:       out (throughput)
%
function out=efficiencyFunction(in, L)
out=(1-exp(-in)).^L;
