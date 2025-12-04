function [Zica, W, T, mu] = fastICA(Z, r, type, flag)

% Constants
TOL = 1e-9;          % Convergence tolerance
MAX_ITERS = 500;     % Maximum number of iterations

% Parse inputs
if nargin < 4
    flag = 1;
end
if nargin < 3
    type = 'kurtosis';
end
n = size(Z, 2);

% Set algorithm type
if strcmpi(type, 'kurtosis')
    nonlinearity = @(x) x.^3;
    nonlinearity_deriv = @(x) 3 * x.^2;
elseif strcmpi(type, 'negentropy')
    nonlinearity = @(x) x .* exp(-0.5 * x.^2);
    nonlinearity_deriv = @(x) (1 - x.^2) .* exp(-0.5 * x.^2);
elseif strcmpi(type, 'tanh')
    nonlinearity = @(x) tanh(x);
    nonlinearity_deriv = @(x) 1 - tanh(x).^2;
else
    error('Unsupported type ''%s''', type);
end

% Center and whiten data
[Zc, mu] = centerRows(Z);
[Zcw, T] = whitenRows(Zc);

% Normalize rows to unit norm
normRows = @(X) bsxfun(@rdivide, X, sqrt(sum(X.^2, 2)));

% Initialize weights
W = normRows(randn(r, size(Z, 1)));  % Use normal distribution for better initialization
k = 0;
delta = inf;

% Prepare status updates
if flag
    fmt = sprintf('%%0%dd', ceil(log10(MAX_ITERS + 1)));
    str = sprintf('Iter %s: max(1 - |<w%s, w%s>|) = %% .4g\\n', fmt, fmt, fmt);
    fprintf('***** Fast ICA (%s) *****\n', type);
end

% Perform Fast ICA
while delta > TOL && k < MAX_ITERS
    k = k + 1;
    
    % Save last weights
    Wlast = W;
    
    % Compute activations
    Sk = W * Zcw;
    
    % Compute non-Gaussianity functions and their derivatives
    G = nonlinearity(Sk);
    Gp = nonlinearity_deriv(Sk);
    
    % Update weights
    W = (G * Zcw') / n - bsxfun(@times, mean(Gp, 2), W); % multiply
    W = normRows(W);
    
    % Decorrelate weights using Gram-Schmidt (faster alternative to SVD)
    for i = 1:r
        for j = 1:i-1
            W(i,:) = W(i,:) - (W(i,:) * W(j,:)') * W(j,:);
        end
        W(i,:) = W(i,:) / norm(W(i,:));
    end
    
    % Update convergence criteria
    delta = max(1 - abs(dot(W, Wlast, 2)));
    
    % Print status update
    if flag
        fprintf(str, k, k, k - 1, delta);
    end
end
if flag
    fprintf('\n');
end

% Independent components
Zica = W * Zcw;

end

