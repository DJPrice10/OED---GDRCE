function [C,number, idx1] = howmanyunique(A)

[B, idx1] = sortrows(A);

[C, ib, ~] = unique(B, 'rows');

number = diff([ib; size(B,1)+1]);