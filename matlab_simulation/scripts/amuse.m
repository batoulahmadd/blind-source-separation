function [y] = amuse(X, tau)

    if nargin < 2
        tau = 1;  % Default time delay
    end

    % Get the size of the input data
    [m, N] = size(X);  % m = number of observed signals, N = number of samples

    % Step 1: Centering - Make data zero mean
    [Xc, mu] = centerRows(X);
    X=Xc;
    
    % Step 2: Pre-whitening using Singular Value Decomposition (SVD)
    [X_white, T] = whitenRows(X);
    
    % Step 3: Time-delayed covariance matrix
    X_tau = [zeros(m, tau), X_white(:, 1:N-tau)];  % Time-delayed version of X_white
    R_tau = (X_white * X_tau') / (N-1);            % Time-delayed covariance matrix

    % Step 4: Eigenvalue Decomposition (EVD) of R_tau
    [B, ~] = eig(R_tau);  % B contains the separating matrix

    % Step 5: Recover the source signals
    y = B' * X_white;  % Separated source signals
  
end
