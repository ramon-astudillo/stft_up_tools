% function[mu_y,Sigma_y] = abs_up_Xi(mu_X,Sigma_X,usepower,min_var,Chik)
%
% This function performs uncertainty propagation of an uncertain STFT 
% through the magnitude or squared magnitude non-linearities of an 
% uncertain STFT.
%
% Input: mu_X    is a matrix of complex valued random variables
%
%        Sigma_X are the correspoding variances
%
%        usepower = 1 forces computing PSD, STSA otherwise
%
%        min_var  is the variance floor
%
% Two modes are possible
%
%        Chik = 2 STFT is assumed to be a matrix of uncorrelated circularly
%                 symmetric complex Gaussian variables. amplitude is 
%                 therefore RICE and PSD scaled Chi-square with two df
%  
%        Chik = 1 Amplitude of the STFT is assumed to be Chi with one 
%                 degree of freedom. PSD with two.
%
% Output: mu_y    Mean fo the Rice or Chi distribution
%         Sigma_y Variance of the Rice or Chi distribution
%
% Please refer to Chapter 5 of [Astudillo 2010] for details 
%
%   [Astudillo 2010] R. F. Astudillo, "Integration of short-time fourier domain speech enhancement and observation uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische Universitaet Berlin, 2010.
% 
% For details. This is available online under: 
%
%   http://opus.kobv.de/tuberlin/volltexte/2010/2676/pdf/astudillo_ramon.pdf 
%    
% Ramon F. Astudilo

function [mu_y,Sigma_y] = abs_up_Xi(mu_X,Sigma_X,usepower,min_var,Chik)

% ARGUMENT HANDLING
if nargin < 5
    Chik = 2;          % Uncertainty model used (default RICE)
end

if nargin < 4
    min_var = 1e-6;    % Treshold to consider a random variable deterministic
end

% INITIALIZE
k3_y   = zeros(size(mu_X));
k4_y   = zeros(size(mu_X));
skew_y = zeros(size(mu_X));
kurt_y = zeros(size(mu_X));
            
% IF SPECTRUM IS UNCERTAIN 
% non-null Variances present (In fact at least one element over treshold)
if any(Sigma_X(:) > min_var)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Non-central Chi square one degree of freedom 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This is the quivalent closed form solution to the approximation used 
    % in 
    %
    % D. Kolossa, A. Klimas, R. Orglmeister "Separation and Robust Recognition of Noisy, Convolutive Speech Mixtures using Time-Frequency Masking and Missing Data Techniques", in Proceedings of the WASPAA 2005, pp. 82-85, New Paltz, NY, USA, October 16-19, 2005.
    %
    % See [Astudillo 2010, Ch. 5]
    
    if Chik == 1
        % ABS FROM RANDOM VARIABLE
        if usepower=='T'
            % PSD => NON-CENTRAL CHI SQUARE 1 DOF
            % Chi^2 parameters
            %
            % k=1;                                 % Degrees of freedom
            % l=(abs(mu).^2)./(lambda/2);          % The other parameter
            %
            % Momments

            mu_y     = Sigma_X/2 + abs(mu_X).^2;
            Sigma_y  = Sigma_X .* (Sigma_X/2 + 2*abs(mu_X).^2);
            
            % Resolve deterministic variables
            mu_y(Sigma_X<min_var)    = abs(mu_X(Sigma_X<min_var)).^2;
            Sigma_y(Sigma_X<min_var) = 0;
            k3_y(Sigma_X<min_var)    = 0;
            k4_y(Sigma_X<min_var)    = 0;
            skew_y(Sigma_X<min_var)  = 0;
            kurt_y(Sigma_X<min_var)  = 0;
            
        else
            % ABS => NON-CENTRAL CHI 1 DOF (Folded Normal unit variance)
            % Chi parameters
            %
            % k=1;                                 % Degrees of freedom
            % l=sqrt((abs(mu).^2)./(lambda/2));    % The other parameter
            %
            % Moments (Overflow robust reformulation)

            warning off
            mu_y    = sqrt(Sigma_X/2) .* sqrt(2/pi) .* exp(-(abs(mu_X).^2)./Sigma_X) + abs(mu_X).*(1-2.*normcdf(-abs(mu_X)./sqrt(Sigma_X/2)));
            Sigma_y = (abs(mu_X).^2 + Sigma_X/2) - mu_y.^2;
            warning on
            
            % Resolve deterministic variables
            mu_y(Sigma_X<min_var)    = abs(mu_X(Sigma_X<min_var));
            Sigma_y(Sigma_X<min_var) = 0;
            k3_y(Sigma_X<min_var)    = 0;
            k4_y(Sigma_X<min_var)    = 0;
            skew_y(Sigma_X<min_var)  = 0;
            kurt_y(Sigma_X<min_var)  = 0;
            
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Non-central Chi square one degree of freedom (Rice Uncertainty Model)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif  Chik==2
        % ABS FROM RANDOM VARIABLE
        if usepower=='T'
            % PSD => NON-CENTRAL CHI SQUARE 2 DOF
            % Chi^2 parameters
            %
            % k=2;                                 % Degrees of freedom
            % l=(abs(mu).^2)./(lambda/2);          % The other parameter
            %
            % Momments
            
            warning off
            mu_y                     = Sigma_X + abs(mu_X).^2;
            Sigma_y                  = Sigma_X .* (Sigma_X + 2*abs(mu_X).^2);
            k3_y                     = 2*Sigma_X.^2 .* (Sigma_X + 3*abs(mu_X).^2);
            k4_y                     = 6* Sigma_X.^3 .* (Sigma_X + 4*abs(mu_X).^2); 
            skew_y                   = (2*(Sigma_X.^(1/2)).*(Sigma_X + 3*abs(mu_X).^2))./((Sigma_X + 2*abs(mu_X).^2).^(3/2));
            kurt_y                   = (6*Sigma_X.*(Sigma_X + 4*abs(mu_X).^2))./((Sigma_X + 2*abs(mu_X).^2).^2)+3;
            warning on
            
            % Resolve deterministic variables
            mu_y(Sigma_X<min_var)    = abs(mu_X(Sigma_X<min_var)).^2;
            Sigma_y(Sigma_X<min_var) = 0;
            k3_y(Sigma_X<min_var)    = 0;
            k4_y(Sigma_X<min_var)    = 0;
            skew_y(Sigma_X<min_var)  = 0;
            kurt_y(Sigma_X<min_var)  = 0;
            
        else
            % ABS => NON-CENTRAL CHI 2 DOF   (Rice unit variance)
            % Chi parameters
            %
            % k=2;                                 % Degrees of freedom
            % l=sqrt((abs(mu).^2)./(lambda/2));    % The other parameter
            %
            % Momments (Overflow robust reformulation)
 
            % Generalized Laguerre polinomial, that evil monster
            warning off
            x                        = -(abs(mu_X).^2)./Sigma_X; 
            % Laguerre 12
            L12                      = exp(x/2).*((1-x).*besseli(0,-x/2) - x.*besseli(1,-x/2));  
            % Correction for large rice SNR ratios that cause binaries 
            % malfunction. It uses approximation of Bessel function for 
            % high values 
            L12(Sigma_X<min_var)     = 0;                                                       
            % Laguerre 32
            L32                      = exp(x/2).*( (3 - 6*x + 2*x.^2).*besseli(0,-x/2) + (2*x.^2 - 4*x).*besseli(1,-x/2) );
            % Correction for large evidence to uncertainty ratios that 
            % cause binaries malfunction. It uses approximation of Bessel 
            % function for high values 
            % Correction for too low values
            L32(Sigma_X<min_var)     = 0;                                                      
            warning on 
            mu_y                     = sqrt(Sigma_X/2).*sqrt(pi/2).*L12;
            Sigma_y                  = Sigma_X + abs(mu_X).^2 - mu_y.^2;
            
            % Rice -> Gaussian approximation for NaNs
            idx_NaN                  = ~isfinite(L12); 
            mu_y(idx_NaN)            = sqrt((abs(mu_X(idx_NaN).^2) + Sigma_X(idx_NaN)/2));
            Sigma_y(idx_NaN)         = Sigma_X(idx_NaN)/2;
            
            % Higher order momments
            E3                       = (Sigma_X/2).^1.5 .* sqrt(pi/2).*L32;
            k3_y                     = E3 - 3*Sigma_y.*mu_y - mu_y.^3;
            k4_y                     = Sigma_X .* (Sigma_X + 2*abs(mu_X).^2);

            % Resolve deterministic variables
            mu_y(Sigma_X<min_var)    = abs(mu_X(Sigma_X<min_var));
            Sigma_y(Sigma_X<min_var) = 0;
            k3_y(Sigma_X<min_var)    = 0;
            k4_y(Sigma_X<min_var)    = 0;
            skew_y(Sigma_X<min_var)  = 0;
            kurt_y(Sigma_X<min_var)  = 0;
            
        end
    end
else
    % ABS FROM DETERMINISTIC SIGNAL
    if usepower=='T'
        mu_y = abs(mu_X).^2;
    else
        mu_y = abs(mu_X);
    end
    Sigma_y= zeros(size(mu_X));
end
