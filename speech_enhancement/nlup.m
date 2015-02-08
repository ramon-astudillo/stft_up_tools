% function [hat_x,Lambda_x] = nlup(mu_XcY,Sigma_XcY,method)        
%
% Obtain non-linear estimators from a Wiener posterior by propagation
%
% Input: mu_XcY      [K, L] Mean STFT of the Wiener filter of K bins and L frames
%
% Input: Sigma_XcY   [K, L] Variance of the Wiener filter (minimum MSE)
%
% Input: method      ('Wiener', 'LSA', 'STSA', 'PSD') Type of estimator in the
%                    respective domain
%
% Output: hat_x      MMSE estimate in the respective domain. Note: This is not
%                    the spectral estimate!. You need to invert the 
%                    non-linearity to account for this.   
%
% Output: Lambda_x   Variance of the estimate, only computable for Wiener, STSA
%                    and PSD
%
% Ramón F. Astudillo May 2014         

function [hat_x,Lambda_x] = nlup(mu_XcY,Sigma_XcY,method)        

% Propagate the posterior to attain non linear estimators
switch method
    case 'Wiener' 
        %
        hat_x    = mu_XcY;
        Lambda_x = Sigma_XcY;
            
    case 'MFCC' 
        %
        hat_x    = mu_XcY;
        Lambda_x = Sigma_XcY;
 
    case 'LSA'
        %
        nu       = (abs(mu_XcY).^2)./Sigma_XcY;        % SNR of Wiener posterior (for the noise Fourier coeff.)
        hat_x    = log(abs(mu_XcY)) + .5*expint(nu);   % MMSE-LSA estimate
        Lambda_x = zeros(size(hat_x));
        
    case 'STSA'
        %
        nu        = (abs(mu_XcY).^2)./Sigma_XcY;                      % SNR of Wiener posterior (for the noise Fourier coeff.)
        hat_x     = gamma(1.5)*sqrt(Sigma_XcY).*exp(-nu/2).*( (1+nu).*besseli(0,nu/2) + nu.*besseli(1,nu/2)); 
        Lambda_x  = Sigma_XcY;
        % APPROXIMATION FOR HIGH NU VALUES: BESSELI has sometimes a
        % little overflow problem due to NU being to high. In this
        % cases we use the approximation  GH1 = mu_XcY./Ymat
        % [Ephraim1985].
        hat_x(nu>1300) = abs(mu_XcY(nu>1300));  
        
     case 'PSD'    
        % Estimate noise amplitude with MMSE-PSD
        hat_x     = abs(mu_XcY).^2 + Sigma_XcY;  
        Lambda_x  = 2*Sigma_XcY.*abs(mu_XcY).^2 + Sigma_XcY.^2;

end
       
