% function [mu_y,Sigma_y] = log_up_logN(mu_x,Sigma_x,diagcov_flag)
%
% Computes the propagation of a log-normal random variable through the 
% logarithm
% 
% Input: mu_x, Sigma_x  Mean and covariance
%        diagcov-flag   Forces ignoring covariances
%
% Please refer to Chapter 6 of 
%
%   [Astudillo2010] R. F. Astudillo, "Integration of short-time fourier domain speech enhancement and observation uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische Universitaet Berlin, 2010.
% 
% For details. This is available online under: 
%
%   http://opus.kobv.de/tuberlin/volltexte/2010/2676/pdf/astudillo_ramon.pdf 
%
% Ramon F. Astudillo

function [mu_y,Sigma_y] = log_up_logN(mu_x,Sigma_x,diagcov_flag)

if nargin < 3
    diagcov_flag = 1;
end

% Quick exit for deterministic case
if all(Sigma_x==0)
    mu_y    = log(mu_x);
    Sigma_y = zeros(size(mu_y));
    return
end

% Get size
[J,L] = size(mu_x);

% CHECK FOR DIAGONAL MATRIX
if diagcov_flag

    % Log-Normal approximation    
    Sigma_y  = log(Sigma_x./(mu_x.^2) + 1);
    mu_y     = log(mu_x) - 0.5*Sigma_y;
    
else
    
    % Initialization
    Sigma_y  = zeros(J,J,L);
    mu_y     = zeros(J,L);
    % For each frame
    for l=1:L
        Sigma_y(:,:,l) = log(Sigma_x(:,:,l)./(mu_x(:,l)*mu_x(:,l)') + 1);
        mu_y (:,l)     = log(mu_x(:,l)) - 0.5*diag(Sigma_y(:,:,l));
    end

end
