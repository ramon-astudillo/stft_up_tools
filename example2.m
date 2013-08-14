% A small example of a realistic use of STFT-UP to produce MMSE-MFCC 
% estimates as explained in
%
%   [1] R. F. Astudillo, R. Orglmeister, "Computing MMSE Estimates and Residual Uncertainty directly in the Feature Domain of ASR using STFT Domain Speech Distortion Models", IEEE Transactions on Audio, Speech and Language Processing, Vol. 21 (5), pp 1023-1034, 2013
% 
% To reproduce the results you will need the AURORA4 corpus, which you can 
% purchase from LDC. Training followed Vertanne's HTK training recipe for 
% WSJ
%
%   http://www.keithv.com/software/htk/
%
% The HTK modified binaries to perform uncertainty decoding and modified 
% imputation can be compiled after aplying the patches found on 
%
%   http://www.astudillo.com/ramon/research/stft-up/
%
% In order to estimate noise variance we used the default IMCRA setting 
% 
%   [2] Israel Cohen, Noise Spectrum estimation in Adverse Environments:Improved Minima Controlled Recursive Averaging. IEEE. Trans. Acoust. Speech Signal Process. VOL. 11, NO. 5, Sep 2003.
%
% The code is also provided (see this script for a running example)
%
% Ramon F. Astudillo, last revision Jun 2013

% Clean the house
clear,clc

% Add needed tools
addpath('stft')
addpath('mfcc_up')
addpath('speech_enhancement')
% addpath('voicebox')   % Add this to write files in HTK format (see end of this file)

%%%%%%%%%%%%%%%%%%%%
%  INITIALIZATION
%%%%%%%%%%%%%%%%%%%%

% HTK CONFIGURATION MFCC_0_D_A
%config.sourceformat   = 'NOHEAD';       % Does not apply
%config.sourcekind     = 'WAVEFORM';     % Does not apply
config.sourcerate     = 625;
config.targetkind     = 'MFCC_0_D_A';  
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
config              = init_stft_HTK(config);

% Initialize MFCCs (compute Mel-fiterbank and DCT matrices)
[config.W,config.T] = init_mfcc_HTK(config);

% SPEECH ENHANCEMENT CONFIGURATION
alpha     = 0.92;                     % Decision directed a priori SNR estimation parameter. 
dB_xi_min = -25;                      % Limits minimal a priori SNR 
imcra     = init_IMCRA(config.nfft/2+1);

% UNCERTAINTY PROPAGATION CONFIGURATION
% Note that 'LOGUT' approx. can not be used with config.usepower = 'T'; and 
% might break due to non positive definiteness when config.diagcov_flag = 0
config.log_prop           = 'LOGN';   % 'LOGN': Log-normal/CGF approximation, 
                                      % 'LOGUT': Unscented transform for the logarithm propagation. 
config.diagcov_flag       = 0;        % 0 = Full covariance after Mel-filterbank considered
config.min_var            = 1e-6;     % Floor for the uncertainty
config.Chik               = 2;        % Use Chi with one or two degrees of freedom, see [Astudillo2010, Ch. 5]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL PROCESSING + FEATURE EXTRACTION FOR EACH FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of files we are going to use (read from an scp for example)
file_list = {'s2_swwp2s.wav'};

for i=1:length(file_list)
	
	%
	% SIMULATED NOISY SIGNAL
	%
	
	% We will use an example from the GRID audiovisual corpus. You can download this from
	% http://spandh.dcs.shef.ac.uk/gridcorpus/examples/s2_swwp2s.wav 
	[x,fs] = wavread(file_list{i});
    % Downsample to 16KHz
    x = resample(x,16,25);
	% We will artificially add white noise as an example
	d      = 0.25*mean(abs(x))*randn(size(x));
	y      = x + d;
	% Compute STFT
	Y      = stft_HTK(y,config);
	% Compute STFT of noise
	D      = stft_HTK(d,config);
	% Get sizes
	[K,L]  = size(Y);

    %
    % SPEECH ENHANCEMENT + FEATURE EXTRACTION
    %
    
    % As show in [1], any non-linear estimator can be computed in two
    % stages. 
    %
    %    1) We derive the posterior associated to the Wiener filter.
    % 
    %    2) We transform this posterior through the corresponding non
    %       linearity. 
    %
    % If this non-linearity is a feature extraction, we have have a MMSE 
    % estimator in that feature domain. STFT-UP can be the used for that 
    % purpose.
    
    % 1) Compute Wiener posterior
    % This will hold the Wiener estimated clean speech
    hat_X_W        = zeros(size(Y));
    % This will hold the residual estimation uncertainty, in other words
    % the variance of the Wiener posterior
    Lambda         = zeros(size(Y));
    % Initialize noise power
    imcra.Lambda_D = abs(Y(:,1)).^2;
    % Initialize Gain and a posteriori SNR
    GH1            = ones(size(Y(:,1)));
    Gamma          = GH1;
    % Loop over frames
    for l=1:L 
        
        % SNR ESTIMATION (II)
        % A posteriori SNR
        new_Gamma = (abs(Y(:,l)).^2)./imcra.Lambda_D;                      % [2, eq.3]
        % Decision directed a priori SNR estimation, with lower bound
        xi        = alpha*(GH1.^2).*Gamma + (1-alpha)*max(new_Gamma-1,0);  % [2, eq.32]
        xi        = max(xi,10^(dB_xi_min/20));                             
        % Update Gamma
        Gamma     = new_Gamma;

        % WIENER Posterior
        % Mean (Wiener filter)
        hat_X_W(:,l) = xi./(1+xi).*Y(:,l);
        % Variance (residual MSE)
        Lambda(:,l)  = xi./(1+xi).*imcra.Lambda_D;
        % Get the gain as well
        GH1          = xi./(1+xi);
        
        % SNR ESTIMATION (I), yes it is done in this order
        % IMCRA estimation of noise variance
        imcra = IMCRA(imcra,Y(:,l),Gamma,xi); 
    end
    
    % 2) Propagate that posterior through the non-linearities to attain the
    % non-linear estimators
    
    % MMSE-LSA + MFCC
    % We propagate into log(|X|) domain and invert the point estimate. We
    % also use the noisy phase, yielding
    hat_X_LSA       = hat_X_W.*exp(.5*expint((abs(hat_X_W).^2)./Lambda));   
    % Compute MFCCs
    hat_x_LSA       = mfcc_up(hat_X_LSA,zeros(size(hat_X_W)),config);
    % Deltas and Accelerations
    hat_x_LSA       = append_deltas_up(hat_x_LSA,zeros(size(hat_x_LSA)),config.targetkind,config.deltawindow,config.accwindow,config.simplediffs);

    % MMSE-MFCC
    % The point estimate is produced directly in MFCC domain
    [mu_x,Sigma_x] = mfcc_up(hat_X_W,Lambda,config);
    % Deltas and Accelerations
    [hat_x,Sigma_x] = append_deltas_up(mu_x,Sigma_x,config.targetkind,config.deltawindow,config.accwindow,config.simplediffs);

    % 
    % WRITING INTO A HTK-FORMAT FILE
    %
    
    % In order for this setup to be used with e.g. the AURORA or CHiME
    % datasest we need to write the features in HTK format. This can be
    % attained with the writehtk.m function of the voicebox toolbox
    % the targetkind used is USER (9), see help writehtk for details.
    % Note the transposed feature matrix
    
    % writehtk('s2_swwp2s.mfcc',hat_x',(config.windowsize-config.overlap)/config.fs,9)
    
end
