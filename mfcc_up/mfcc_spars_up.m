% function [mu_x,Sigma_x] = mfcc_up(X,p,cf) 
%
% This function perfoms uncertainty propagation for MFCC based feature 
% extractions. It transforms an uncertain STFT into MFCC domain. The input
% STFT is assumed to be a matrix of uncorrelated scaled Bernoulli random 
% variables. The output is assumed to be a matrix of Gaussian variables.
%
% Input:  X        Observed STFT (one or two chanels)
%         p        Source presence probability for each time frequency bin
%
% Output: mu_x     matrix of real values means of the Gaussian variable
%         Sigma_x  matrix of variances
%
% Please refer to  
%
%        F. Nesta and M. Matassoni and R. F. Astudillo, "A flexible spatial blind source extraction framework for robust speech recognition in noisy environments", 2nd International Workshop on Machine Listening in Multisource Environments CHiME, pp. 33-38, June 2013.
%   
% For details. This is available online under: 
%
%   spandh.dcs.shef.ac.uk/chime_workshop/papers/pP3_nesta.pdf
%
% Ramon F. Astudillo

function [mu_x,Sigma_x] = mfcc_spars_up(X,p,cf) 

% CONFIGURATION INITIALIZATION
if ~isfield(cf,'W') || ~isfield(cf,'T')
	[cf.W,cf.T] = init_mfcc_HTK(cf);
end

% SHORT TIME SPECTRAL AMPLITUDE OR PERIODOGRAM
% Single uncertain channel
if size(X,3) == 1

    if cf.usepower=='T'
        % The expected value and variance of the source are those of a bernoulli event
        mu_x    = p.*abs(X).^2;
        Sigma_x = p.*(1-p).*abs(X).^4;
    else
        mu_x    = p.*abs(X);	
        Sigma_x = p.*(1-p).*abs(X).^2;
    end

% Two uncertain channels, delay and sum pointing in the frontal-direction
elseif size(X,3) == 2
     
    % Input are two bernulli spectrums 
    if cf.usepower=='T'
        % The expected value and variance of the source are those of a bernoulli event
        mu_x    = p(:,:,1).*abs(X(:,:,1)).^2 + p(:,:,2).*abs(X(:,:,2)).^2;	
        Sigma_x = p(:,:,1).*(1-p(:,:,1)).*abs(X(:,:,1)).^4 + p(:,:,2).*(1-p(:,:,2)).*abs(X(:,:,2)).^4;
    else
        mu_x    = p(:,:,1).*abs(X(:,:,1)) + p(:,:,2).*abs(X(:,:,2));	
        Sigma_x = p(:,:,1).*(1-p(:,:,1)).*abs(X(:,:,1)).^2 + p(:,:,2).*(1-p(:,:,2)).*abs(X(:,:,2)).^2;
    end   
    
else
       error('Not implemented yet for more than two channels')
end

% If this is the targetkind solicited exit
if strfind(cf.targetkind,'SPECTROGRAM'); return; end

% MEL FILTERBANK
[mu_x,Sigma_x] = linear_up(mu_x,Sigma_x,cf.W,cf.diagcov_flag);

% MEL-FLOOR
% This is, in fact, ignored in the propagation, it better not be very high
mu_x(mu_x<cf.melfloor) = cf.melfloor;

% If this is the targetkind solicited exit
if strfind(cf.targetkind,'MELSPEC'); return; end

% LOGARITHM
if strcmp(cf.log_prop,'LOGUT')    
    % Unscented Transform
    [mu_x,Sigma_x] = log_up_ut(mu_x,Sigma_x,cf.diagcov_flag);
elseif strcmp(cf.log_prop,'LOGN') 
    % Log-Normal assumption
    [mu_x,Sigma_x] = log_up_logN(mu_x,Sigma_x,cf.diagcov_flag);
end

% If this is the targetkind solicited exit
if strfind(cf.targetkind,'FBANK'); return; end

% DISCRETE COSINE TRANSFORM
[mu_x,Sigma_x] = linear_up(mu_x,Sigma_x,cf.T,cf.diagcov_flag);

% APPEND ENERGY AS LAST TERM  (ONLY IF _0 used?)
mu_x    = [ mu_x(2:end,:)    ; mu_x(1,:) ];
% Fix to get only the diagonal from full covariances
if ~cf.diagcov_flag   
    L = size(Sigma_x,3);
    Sigma_tmp = zeros(size(Sigma_x,2),L);
    for l = 1:L
        Sigma_tmp(:,l) = diag(Sigma_x(:,:,l));
    end
    Sigma_x = Sigma_tmp;
end
Sigma_x = [ Sigma_x(2:end,:) ; Sigma_x(1,:) ];
