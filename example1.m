% A small example of how to use STFT-UP to propagate a complex Gaussian 
% posterior in STFT domain into MFCC domain. The STFT posterior is here 
% generated with a Wiener filter
%
% For the propagation formulas please refer to Chapters 5 and 6 of 
%
%   [Astudillo2010] R. F. Astudillo, "Integration of short-time fourier domain speech enhancement and observation uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische Universitaet Berlin, 2010.
% 
% available online under: 
%
%   http://opus.kobv.de/tuberlin/volltexte/2010/2676/pdf/astudillo_ramon.pdf 
%
% Note that a Wiener filter is just a particular way of generating an 
% uncertain description of the STFT of the clean signal. Refer to
%
%   [Astudillo2013] R. F. Astudillo, R. Orglmeister, "Computing MMSE Estimates and Residual Uncertainty directly in the Feature Domain of ASR using STFT Domain Speech Distortion Models", IEEE Transactions on Audio, Speech and Language Processing, Vol. 21 (5), pp 1023-1034, 2013
% 
% for details on Wiener posterior as a measure of uncertaity also 
%
%   http://www.astudillo.com/ramon/research/
%
% for other forms of computing uncertainty
%
% Note that this example is thought for the tutorial of Interspeech 2012. 
% In coherence with the slides, noise is denoted with N rather than the 
% usual D.
%
% Ramon F. Astudillo, last revision Jun 2013

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
% config.diagcov_flag = 0;
config.log_prop           = 'LOGN';   % 'LOGN': Log-normal/CGF approximation, 
                                      % 'LOGUT': Unscented transform for the logarithm propagation. 
config.diagcov_flag       = 0;        % 0 = Full covariance after Mel-filterbank considered
config.min_var            = 1e-6;     % Floor for the uncertainty
config.Chik               = 2;        % Use Chi with one or two degrees of freedom, see [Astudillo2010, Ch. 5]

%%%%%%%%%%%%%%%%%%%%
% SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%

% SIMULATED NOISY SIGNAL
% We will use an example from the GRID audiovisual corpus. You can download this from
% http://spandh.dcs.shef.ac.uk/gridcorpus/examples/s2_swwp2s.wav 
[x,fs] = wavread('s2_swwp2s.wav');
% Downsample to 16KHz
x = resample(x,16,25);
% We will artificially add white noise as an example
n      = 0.25*mean(abs(x))*randn(size(x));
y      = x + n;
% Compute STFT
Y      = stft_HTK(y,config);
% Compute STFT of noise
N      = stft_HTK(n,config);

% Get sizes
[K,L]  = size(Y);

% WIENER FILTER WITH SIMULATED A PRIORI KNOWLEGDE
% Simulate ideal noise variance knowlegde
Lambda_N       = repmat(mean(abs(N).^2,2),1,L);
% Get speech variance by power subtraction (in reality decision directed 
% approaches are better)
% Stablish a power floor at 1e-6
Lambda_X       = max(abs(Y).^2 - Lambda_N,1e-6);
% Wiener Estimate
hat_X          = Lambda_X./(Lambda_X+Lambda_N).*Y;
% Residual MSE (Wiener posterior variance)
Lambda         = Lambda_N.*Lambda_X./(Lambda_X+Lambda_N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROPAGATION THROUGH FEATURE EXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ANALYTIC SOLUTION USING STFT-UP
% STFT-UP solution for MFCCs                                  
[mu_x,Sigma_x] = mfcc_up(hat_X,Lambda,config);
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
    % Draw from the complex Gaussian Wiener posterior
    sample_X = randcg(hat_X,Lambda,1);
    % Transform through conventional mfccs (zero variance)                             
    sample_x = mfcc_up(sample_X,zeros(size(hat_X)),config);
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