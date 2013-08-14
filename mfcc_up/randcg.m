% Sample from a complex Gaussian distribution of mean mu_X and variance Sigma_X
%
% Ramon F. Astudillo

function sample_X = randcg(mu_X,Sigma_X,n_samples)

% Get sizes
[K,L] = size(mu_X);

% COMPLEX GAUSSIAN 
sample_X = (real(mu_X) + sqrt(Sigma_X/2) .* randn(K,L*n_samples))...
         + (imag(mu_X) + sqrt(Sigma_X/2) .* randn(K,L*n_samples))*sqrt(-1);
