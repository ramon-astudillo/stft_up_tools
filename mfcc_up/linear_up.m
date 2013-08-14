% function [mu_y,Sigma_y] = linear_up(mu_x,Sigma_x,W,diagcovflag)
%
% Computes the mean and variance of the random variable transformed through
% the linear transformation W.
%
% diagcovflag = 1 ignores the covariances
%
% Ramon F. Astudillo

function [mu_y,Sigma_y] = linear_up(mu_x,Sigma_x,W,diagcovflag)

% ARGUMENT HANDLING
if nargin < 3
    %defaults
    diagcovflag = 1;  % Only diagonal covariances considered
end

% Get sizes
L     = size(mu_x,2);
[J,K] = size(W);
% Use sparse matrix to avoid NaNs in calculations
W = sparse(W);

% MEANS
mu_y = W*mu_x;

% Quick exit for deterministic case
if all(Sigma_x(:)==0)
    Sigma_y = zeros(size(mu_y));
    return
end

% VARS
if diagcovflag
    Sigma_y = (W.^2)*Sigma_x;
else
    % TODO: UNROLL THIS LOOP
    Sigma_y = zeros(J,J,L);
    
    % Initial Covariance Diagonal
    if length(size(Sigma_x)) == 2
        % For each frame
        for l = 1:L
            Sigma_y(:,:,l) = W * diag(Sigma_x(:,l)) * W';
        end
    else
        % For each frame
        for l = 1:L
            Sigma_y(:,:,l) = W * Sigma_x(:,:,l) * W';
        end
    end
end