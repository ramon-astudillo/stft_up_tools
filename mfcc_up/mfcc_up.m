% function [mu_x,Sigma_x] = mfcc_up(mu_X,Sigma_X,cf) 
%
% This function perfoms uncertainty propagation for MFCC based feature 
% extractions. It transforms an uncertain STFT into MFCC domain. The input
% STFT is assumed to be a matrix of uncorrelated circularly symmetric 
% complex Gaussian random variable. The output is assumed to be a matrix
% of Gaussian variables.
%
% Input:  mu_X     matrix of complex means of the complex Gaussian
%         Sigma_X  matrix of variances
%
% Output: mu_x     matrix of real values means of the Gaussian variable
%         Sigma_x  matrix of variances
%
% Please refer to Chapter 6 of 
%
%   [Astudillo 2010] R. F. Astudillo, "Integration of short-time fourier domain speech enhancement and observation uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische Universitaet Berlin, 2010.
% 
% For details. This is available online under: 
%
%   http://opus.kobv.de/tuberlin/volltexte/2010/2676/pdf/astudillo_ramon.pdf 
%
% Ramon F. Astudillo

function [mu_x,Sigma_x] = mfcc_up(mu_X,Sigma_X,cf) 

% CONFIGURATION INITIALIZATION
if ~isfield(cf,'W') || ~isfield(cf,'T')
	[cf.W,cf.T] = init_mfcc_HTK(cf);
end

% SHORT TIME SPECTRAL AMPLITUDE OR PERIODOGRAM
[mu_x,Sigma_x] = abs_up_Xi(mu_X,Sigma_X,cf.usepower,cf.min_var,cf.Chik);

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
    [mu_x,Sigma_x] = log_up_ut(mu_x,Sigma_x,cf.melfloor,cf.diagcov_flag);
elseif strcmp(cf.log_prop,'LOGN') 
    % Log-Normal assumption
    [mu_x,Sigma_x] = log_up_logN(mu_x,Sigma_x,cf.diagcov_flag);
end

% If this is the targetkind solicited exit
if strfind(cf.targetkind,'FBANK'); 
    if ~cf.diagcov_flag & any(Sigma_x(:))
        Sigma_x = get_diagonal(Sigma_x, cf);
    end
    return; 
end

% DISCRETE COSINE TRANSFORM
[mu_x,Sigma_x] = linear_up(mu_x,Sigma_x,cf.T,cf.diagcov_flag);

% APPEND ENERGY AS LAST TERM  (ONLY IF _0 used?)
mu_x    = [ mu_x(2:end,:)    ; mu_x(1,:) ];
% Fix to get only the diagonal from full covariances
if ~cf.diagcov_flag & any(Sigma_x(:))  
    Sigma_x = get_diagonal(Sigma_x, cf);
end
Sigma_x = [ Sigma_x(2:end,:) ; Sigma_x(1,:) ];




%
% subfunctions
%

function Sigma_x = get_diagonal(Sigma_x, cf)

% Fix to get only the diagonal from full covariances
L         = size(Sigma_x,3);
Sigma_tmp = zeros(size(Sigma_x,2),L);
for l = 1:L
    Sigma_tmp(:,l) = diag(Sigma_x(:,:,l));
end
Sigma_x = Sigma_tmp;
