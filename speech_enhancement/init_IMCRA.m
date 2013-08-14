% cf = init_IMCRA(K,varargin)
%
% Initialization function of Improved Minima Controlled Recursive Averaging
% (IMCRA) after
%
%  [3] Israel Cohen, Noise Spectrum estimation in Adverse Environments:
%  Improved Minima Controlled Recursive Averaging. IEEE. Trans. Acoust.
%  Speech Signal Process. VOL. 11, NO. 5, Sep 2003.
%
% This is to be used with the function IMCRA.m, see this function for more
% info.
%
%  Input: K is the number of frequency bins below the niquist frequency
%
%  Any parameter of the default config can be changed through the varargin
%  e.g. 
%
%  imcf = init_IMCRA(256,'alpha_s',0.85,'IS',20)
%
% Ramón F. Astudillo

function cf = init_IMCRA(K,varargin)
 
% IMCRA DEFAULT CONFIGURATION
% Consult [3] on how to set the paremeters
cf.IS                 = 10;               % Number of frames used to initialize noise
cf.w                  = 1;                % How many adjacent bins (before and after bin k) are used to smooth spectrogram in frecuency (S_f(k,l))
cf.alpha_s            = 0.9;              % Adaptation rate for the spectrogram smoothing in time
cf.U                  = 8;                % U spectrogram values are stored to perform minimum tracking
cf.V                  = 15;               % Each V frames minimum tracking will take place
%D                    = 120;              % Number of frames before the actual frame considered in minimum traking. (Irrelevant its included in Bmin)
cf.Bmin               = 1.66;             % Bias of minimum noise estimate
cf.Gamma0             = 4.6;              % A priori signal to noise ratio treshold for rough speech absence I (other conditions must also be met)
cf.Gamma1             = 3;                % A priori signal to noise ratio treshold for speech absence q computation (other conditions must also be met)
cf.zeta0              = 1.67;             % Threshold for "something similar to minimun estimated a priori signal to noise ratio but with spectrogram instead of signal [eq 18]"
cf.alpha_d            = 0.85;             % Adaptation rate for speech probability dependent time smoothing paramater
cf.beta               = 1.47;             % Bias compensation when speech is absent for noise variance estimation
cf.normalized_hanning = 1;                % Allows hanning window to be normalized

% Overload defaults specified through varargin
% Sanity Check: Even number of arguments
if mod(length(varargin),2)
    error(sprintf('The number of arguments in VARARGIN must be even (pairs of variable_name variable_value).'))
end
for i=1:2:length(varargin)
    % Sanity Check: Argument not being a string or known config parameter
    if ~ischar(varargin{i}) || ~isfield(cf,varargin{i})
        error(sprintf('%dth argument in VARARGIN is invalid',i))
    end
    % Overload parameter
    cf.(varargin{i}) = varargin{i+1};
end  

% Get number of frequency bins
cf.K = K;
% Set counter for number of frames processed
cf.l = 0;
% Set counter for the buffer update
cf.j = 0; 

% To do the smoothing in frequency we use a matrix of indices on each frame
% and a corresponding matrix of row-stacked hanning windows
if cf.normalized_hanning
    cf.smth_win = repmat((hanning(2*cf.w+1)/sum(hanning(2*cf.w+1)))',cf.K,1);
else
    cf.smth_win = repmat(hanning(2*cf.w+1)',K,1);
end
% Create coresponding matrix of indices
cf.smth_mat = (repmat(1:cf.K,2*cf.w+1,1) + repmat(-cf.w:cf.w,cf.K,1)')';
% Ignore indices out of bounds
idx              = find(cf.smth_mat<=0 | cf.smth_mat > cf.K);
cf.smth_mat(idx) = 1;
cf.smth_win(idx) = 0;

% Smoothed spectrograms
cf.S             = [];    % Smoothed Spectrogram first iteration
cf.tilde_S       = cf.S;  % Smoothed Spectrogram second iteration
cf.Smin          = cf.S;  % Smoothed Spectrogram minimum first iteration
cf.tilde_Smin    = cf.S;  % Smoothed Spectrogram minimum first second iteration
cf.Smin_sw       = cf.S;  % Smoothed Spectrogram minimum running minimum
cf.tilde_Smin_sw = cf.S;  % Second smoothed Spectrogram minimum running minimum
% Store buffers
cf.Storing       = [];    % Smoothed Spectrogram minimum first iteration store buffer
cf.tilde_Storing = [];    % Smoothed Spectrogram minimum second iteration store buffer
% Other parameters
cf.ov_Lambda_D   = [];    % Biased noise variance estimate
cf.Lambda_D      = [];    % Unbiased noise variance estimate 
cf.p             = 1;     % A posteriori speech presence probability