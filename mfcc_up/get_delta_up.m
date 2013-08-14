% function [mu_y,Sigma_y] = get_delta_up(mu_x,Sigma_x,window,simplediffs)
%
% Compute delta coefficients for a Gaussian uncorrelated uncertain features
%
% Input:  mu_x,Sigma_x  Matrix of means and variances of the Gaussian 
%                       variable
%
% Output: mu_y,Sigma_y
%         window        Half the size of the accelerations window
%         simplediffs   Scale the deltas and accelerations
%
% Ramon F. Astudillo

function [mu_y,Sigma_y] = get_delta_up(mu_x,Sigma_x,window,simplediffs)

[I,L]   = size(mu_x);
denom   = 2*sum([1:window].^2);

mu_y    = zeros(I,L);
Sigma_y = zeros(I,L);

for l=1:L
    num    = 0;
    varacc = 0;
    for theta=1:window
        % Check first and last frame bounds
        ind1 = min(l+theta,L);
        ind2 = max(l-theta,1);
        
        % Compute simple or standard diffs
        if simplediffs=='T'
            num    = num    + mu_x(:,ind1)-mu_x(:,ind2);
            varacc = varacc + Sigma_x(:,ind1)+Sigma_x(:,ind2);

        else
            num    = num    + theta      * (mu_x(:,ind1)    - mu_x(:,ind2));
            varacc = varacc + (theta.^2) * (Sigma_x(:,ind1) + Sigma_x(:,ind2));
        end
    end
    mu_y(:,l)    = num/denom;
    Sigma_y(:,l) = varacc/(denom.^2);
end