function [Zc, mu] = centerRows(Z)

% Compute data mean
mu = mean(Z,2);

% Subtract mean
Zc = Z - mu;
