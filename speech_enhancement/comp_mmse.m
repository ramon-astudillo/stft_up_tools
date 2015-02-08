function [hat_X, Lambda] = comp_mmse(hat_x, Lambda_x, alpha, method)

% Propagate the posterior to attain non-linear estimators
switch method
    
    case 'Wiener'
        
        % MMSE estimate is directly in complex domain
        hat_X  = hat_x;
        Lambda = Lambda_x;

    case 'MFCC'
        
        % MMSE estimate is directly in complex domain
        hat_X  = hat_x;
        Lambda = Lambda_x;
        
    case 'LSA'
        
        % Retrieve amplitude estimate, use mixture phase
        hat_X  = exp(hat_x).*alpha;
        Lambda = zeros(size(hat_X));
        
    case 'STSA'
        
        % Retrieve amplitude estimate, use mixture phase
        hat_X  = hat_x.*alpha;
        Lambda = zeros(size(hat_X));
        
    case 'STSP'
        
        % Retrieve amplitude estimate, use mixture phase
        hat_X  = sqrt(hat_x).*alpha;
        Lambda = zeros(size(hat_X));
        
end
