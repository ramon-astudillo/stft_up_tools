% function[mu_y,Sigma_y] = append_deltas_up(mu_x,Sigma_x,targetkind,deltawindow,accwindow,simplediffs)
%
% Computes and appends delta and acceleration coefficients follwoing HTKs 
% configuration  
%
% Input: mu_x,Sigma_x  Matrix of means and variances of the Gaussian 
%                      variable
%              
%        HTK parameters 
%
%        targetkind    Feature extraction name  
%        deltawindow   Half the size of the delta window
%        accwindow     Half the size of the accelerations window
%        simplediffs   Scale the deltas and accelerations
%
% Output: mu_y,Sigma_y Matrix of means and variances.
%
% Ramon F. Astudillo

function[mu_y,Sigma_y] = append_deltas_up(mu_x,Sigma_x,targetkind,deltawindow,accwindow,simplediffs)

% BY DEFAULT, EQUAL TO INPUT
mu_y    = mu_x;
Sigma_y = Sigma_x;

% DELTA
if ~isempty(strfind(targetkind,'_D')) || ~isempty(strfind(targetkind,'_A'))
    
    [mu_delta,Sigma_delta] = get_delta_up(mu_x,Sigma_x,deltawindow,simplediffs);
    
    % If deltas and accelerations solicited
    if ~isempty(strfind(targetkind,'_D'))
        mu_y                   = [ mu_x  ; mu_delta ];
        Sigma_y                = [Sigma_x; Sigma_delta ];
    else
    % If only accelerations solicited    
        mu_y                   = mu_x;
        Sigma_y                = Sigma_x;
    end
        
end

% ACCELERATION
if ~isempty(strfind(targetkind,'_A'))
    [mu_acc,Sigma_acc] = get_delta_up(mu_delta,Sigma_delta,accwindow,simplediffs);
    mu_y               = [ mu_y   ; mu_acc ];
    Sigma_y            = [ Sigma_y; Sigma_acc ];
end