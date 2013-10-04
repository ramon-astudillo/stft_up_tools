% A small example of how to use STFT-UP to propagate a scaled Bernoulli  
% posterior in STFT domain into MFCC domain. The speech probabilities 
% are here generated with a kind-of-Wiener filter (see below)
%
% For the propagation formulas please refer to 
%
% 	[Nesta 2013] F. Nesta, M. Matassoni, R. F. Astudillo, "A FLEXIBLE SPATIAL BLIND SOURCE EXTRACTION FRAMEWORK FOR ROBUST SPEECH RECOGNITION IN NOISY ENVIRONMENTS", In 2nd International Workshop on Machine Listening in Multisource Environments (CHiME), pages 33-38, June 2013
% 
% available online under: 
%
%	spandh.dcs.shef.ac.uk/chime_workshop/papers/pP3_nesta.pdf
%
% Note that a kind-of-Wiener filter is just a particular way of 
% generating speech probabilities. 
%
% Ramon F. Astudillo, last revision Aug 2013

% Clean the house
clear,clc

%%%%%%%%%%%%%%%%%%%%
%  INITIALIZATION
%%%%%%%%%%%%%%%%%%%%

% Add needed tools
addpath('stft')
addpath('mfcc_up')  

% HTK CONFIGURATION MFCC_0_D_A_Z
%config.sourceformat   = 'NOHEAD';       % Does not apply
%config.sourcekind     = 'WAVEFORM';     % Does not apply
config.sourcerate     = 625;
config.targetkind     = 'MFCC_0_D_A_Z'; 
config.targetrate     = 100000;
%config.savecompressed = 'T';            % Does not apply
%config.savewithcrc    = 'T';            % Does not apply
config.windowsize     = 250000;
config.usehamming     = 'T';
config.preemcoef      = 0.9700;
config.numchans       = 26;
config.ceplifter      = 22;
config.numceps        = 12;
config.enormalise     = 'T';
config.zmeansource    = 'T';
%config.byteorder      = 'VAX';          % Does not apply    
%config.hifreq         = -1;             % Not supported 
%config.lofreq         = -1;             % Not supported          
config.usepower       = 'T';                                  
config.deltawindow    = 2;              
config.accwindow      = 2;              
config.simplediffs    = 'F';  
config.melfloor       = exp(-10);
% Translates HTK STFT parameters into MATLAB ones
config                = init_stft_HTK(config);

% Initialize MFCCs (compute Mel-fiterbank and DCT matrices)
[config.W,config.T]   = init_mfcc_HTK(config);

% UNCERTAINTY PROPAGATION CONFIGURATION
% Note that 'LOGUT' approx. can not be used with config.usepower = 'T'; and 
% might break due to non positive definite matirces when 
config.log_prop           = 'LOGN';   % 'LOGN': Log-normal/CGF approximation, 
                                      % 'LOGUT': Dies with spasity
                                      % uncertainty (DONT USE IT!)
config.diagcov_flag       = 0;        % 0 = Full covariance after Mel-filterbank considered
config.min_var            = 1e-6;     % Floor for the uncertainty
config.Chik               = 2;        % Actually not needed here, kept for compatibility


%%%%%%%%%%%%%%%%%%%%
% SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%

% SIMULATED NOISY SIGNAL
% We will use an example from the PASCAL-CHiME 2011. You can download this from
% http://spandh.dcs.shef.ac.uk/projects/chime/PCC/datasets.html
[y_t,fs] = wavread('DATA/s29_pbiz6p.wav');
[d_t,fs] = wavread('DATA/s29_pbiz6p_enhanced_target.wav');
[x_t,fs] = wavread('DATA/s29_pbiz6p_enhanced_noise.wav');
% STFT
Y        = stft_HTK(y_t,config);
D        = stft_HTK(d_t,config);
X        = stft_HTK(x_t,config);

% We use the ratio of amplitudes to attain an estimate of speech activity. This is 
% not the only way it could be done.
p   = abs(X)./(abs(X) + abs(D));

% Get sizes
[K,L,n]  = size(Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROPAGATION THROUGH FEATURE EXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ANALYTIC SOLUTION USING STFT-UP
% STFT-UP solution for MFCCs                                  
[mu_x,Sigma_x] = mfcc_spars_up(Y,p,config);
% Deltas and Accelerations
[mu_x,Sigma_x] = append_deltas_up(mu_x,Sigma_x,config.targetkind,config.deltawindow,config.accwindow,config.simplediffs);
% Cepstral Mean Subtraction
[mu_x,Sigma_x] = cms_up(mu_x,Sigma_x,config.targetkind);

% MONTE CARLO SOLUTION
% Not the fastest implementation, but spares some code. 1e3-1e4 samples 
% should give reasonable execution time
max_samples = 1e3;   
% This will store first and second order statistics
mu_x_MC     = zeros(size(mu_x));
mu_x2_MC    = zeros(size(mu_x));
% Draw samples and accumulate statistics
for n = 1:max_samples
    % Draw from scaled Bernoulli STFT
    sample_X_1 = Y(:,:,1).*(rand(K,L) < p(:,:,1));
    sample_X_2 = Y(:,:,2).*(rand(K,L) < p(:,:,2));
    % Beamformer
    sample_X   = sample_X_1 + sample_X_2;
    % Transform through conventional mfccs (zero variance)                             
    sample_x = mfcc_up(sample_X,zeros(size(sample_X)),config);
    % Deltas and Accelerations
    sample_x = append_deltas_up(sample_x,zeros(size(sample_x)),config.targetkind,config.deltawindow,config.accwindow,config.simplediffs);
    % Cepstral Mean Subtraction
    sample_x = cms_up(sample_x,zeros(size(sample_x)),config.targetkind);
    % Accumulate statistics
    mu_x_MC  = ((n-1)*mu_x_MC  + sample_x)/n;
    mu_x2_MC = ((n-1)*mu_x2_MC + sample_x.^2)/n;
    % To see if there is time for a coffe, inform each time 1% done 
    if ~mod(n,ceil(0.01*max_samples))
        fprintf('\rProcessed %3d%%',ceil(100*n/max_samples))
    end
end 
fprintf('\n')
% Variances
Sigma_x_MC = mu_x2_MC - mu_x_MC.^2;

% PLOT
figure
subplot(2,2,1)
imagesc(mu_x),colorbar
xlabel('frames')
ylabel('Mean of the MFCCs')
title('Analytic solution (STFT-UP)')
subplot(2,2,2)
imagesc(Sigma_x),colorbar
title('Variance of the MFCCs')
xlabel('frames')
ylabel('Variance of the MFCCs')
title('Analytic solution (STFT-UP)')
subplot(2,2,3)
imagesc(mu_x_MC),colorbar
xlabel('frames')
ylabel('Mean of the MFCCs')
title('Monte Carlo solution')
subplot(2,2,4)
imagesc(Sigma_x_MC),colorbar
xlabel('frames')
ylabel('Variance of the MFCCs')
title('Monte Carlo solution')
