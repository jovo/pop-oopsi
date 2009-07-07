function [patt,freq]=NETFIT_compress_cols(H)
%Function compresses firing patterns into pattern x frequencies pairs.
% [patterns,frequencies]=NETFIT_compress_cols(H)
% H is a NxT matrix of firing patterns recorded on N neurons in T time-bins;
% pattern is a NxK matrix of unique firing patterns and frequencies is a
% 1xK vector of frequencies of these patterns in H.

B    = sortrows(H');
ind  = [find(max(diff(B),[],2));size(B,1)];
patt = B(ind,:)';
freq = diff([0;ind(:)])';