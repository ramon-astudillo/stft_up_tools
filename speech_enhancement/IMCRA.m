% function cf = IMCRA(cf,Y_l,Gamma,xi)
%
% Implementation of Improved Minima Controlled Recursive Averaging(IMCRA) 
% after
%
%  [1] Israel Cohen, Noise Spectrum estimation in Adverse Environments:
%  Improved Minima Controlled Recursive Averaging. IEEE. Trans. Acoust.
%  Speech Signal Process. VOL. 11, NO. 5, Sep 2003.
%
% Input: cf     IMCRA configuration from previous iteration or obtained from
%               init_IMCRA.m
%
%        Y_l    Current frame of the STFT of size [K, 1]
%
%        Gamma  A posteriori signal to noise ratio    
%
%           xi  A priori signal to noise ratio 
%
% Output:   
%
%         cf.Lambda_D is the noise variance estimate 
%         cf.p        is the a posteriori speech probability
%
%
% Example of typical use for noisy STFT matrix Y
%    
%   % This can be done once for all files and all frames of each file
%   icf      = init_IMCRA(size(Y,1)); 
%   % File loop
%   for f=1:n_files
%       % Initialize noise variance for this file
%       icf.Lambda_D  = abs(Y(:,1)).^2;
%       % Frame loop
%       for l=1:n_frames       
%           % Here speech enhancement code e.g. Wiener + Decision directed
%           % You may use icf.Lambda_D or icf.p (initialized to 1)
%           ... 
%           % Update IMCRA parameters
%           icf      = IMCRA(icf,Y(:,l),Gamma,xi);       
%
% Ramón F. Astudillo

function cf = IMCRA(cf,Y_l,Gamma,xi)

% Increase frame counter
cf.l = cf.l + 1;

% If in first frame, initialize buffers
if cf.l == 1    
    % Smoothed spectrograms
    cf.S             = sum(cf.smth_win.*abs(Y_l(cf.smth_mat)).^2,2); % Smoothed Spectrogram first iteration
    cf.tilde_S       = cf.S;                                         % Smoothed Spectrogram second iteration
    cf.Smin          = cf.S;                                         % Smoothed Spectrogram minimum first iteration
    cf.tilde_Smin    = cf.S;                                         % Smoothed Spectrogram minimum first second iteration
    cf.Smin_sw       = cf.S;                                         % Smoothed Spectrogram minimum running minimum
    cf.tilde_Smin_sw = cf.S;                                         % Second smoothed Spectrogram minimum running minimum
    % Store buffers
    cf.Storing       = [];                                           % Smoothed Spectrogram minimum first iteration store buffer
    cf.tilde_Storing = [];                                           % Smoothed Spectrogram minimum second iteration store buffer
    % Other parameters
    cf.ov_Lambda_D   = abs(Y_l).^2;                                  % Biased noise variance estimate
    cf.Lambda_D      = cf.ov_Lambda_D;                               % Unbiased noise variance estimate
    cf.p             = ones(size(Y_l));                              % A posteriori speech presence probability
end
    
% If in initialization segment, update noise stats only
% Note: Keep in mind that IS might be zero
if cf.l <= cf.IS
    % Compute minimum statistics smoothed spectrograms
    Sf         = sum(cf.smth_win.*abs(Y_l(cf.smth_mat)).^2,2);             % [3,eq.14]
    % Time smoothing
    cf.S       = cf.alpha_s*cf.S + (1-cf.alpha_s)*Sf;                      % [3,eq.15]
    % Update running minimum
    cf.Smin    = min(cf.Smin,cf.S);
    cf.Smin_sw = min(cf.Smin_sw,cf.S);
    % Compute smoothed spectrogram
    cf.Lambda_D = cf.alpha_d*cf.Lambda_D + (1-cf.alpha_d)*abs(Y_l).^2;
    % Set a posteriori speech probability to zero
    cf.p           = zeros(size(Y_l));
    
else
    
    % FIRST MINIMA CONTROLLED VAD
    % This provides a rough VAD to eliminate relatively strong speech
    % components towards the second power spectrum estimation
    Sf         = sum(cf.smth_win.*(abs(Y_l(cf.smth_mat)).^2),2);           % [3,eq.14]
    % Time smoothing
    cf.S       = cf.alpha_s*cf.S+(1-cf.alpha_s)*Sf;                        % [3,eq.15]
    % update running minimum
    cf.Smin    = min(cf.Smin,cf.S);
    cf.Smin_sw = min(cf.Smin_sw,cf.S);
    % Indicator function for VAD
    Gamma_min     = (abs(Y_l).^2)./(cf.Bmin*cf.Smin);                      % [3,eq.18]
    zeta          = cf.S./(cf.Bmin*cf.Smin);                               % [3,eq.21]
    I             = zeros(size(Y_l));
    I((Gamma_min < cf.Gamma0 ) & (zeta < cf.zeta0)) = 1;                   % [3,eq.21]
    
    % SECOND MINIMA CONTROLLED VAD
    % This provides the speech probability needed to compute the final
    % noise estimation. The hard VAD index I, computed in the previous
    % estimation, is here used to exclude strong speech components.
    idx              = find(sum(I(cf.smth_mat),2) == 0);                   % [3,eq.26]
    warning off
    cf.tilde_Sf      = sum(cf.smth_win.*I(cf.smth_mat).*(abs(Y_l(cf.smth_mat)).^2),2)./sum(cf.smth_win.*I(cf.smth_mat),2);
    warning on
    cf.tilde_Sf(idx) = cf.tilde_S(idx);
    % Time smoothing
    cf.tilde_S       = cf.alpha_s*cf.tilde_S+(1-cf.alpha_s)*cf.tilde_Sf;   % [3,eq.27]
    % Update running minimum
    cf.tilde_Smin    = min(cf.tilde_Smin,cf.tilde_S);                      % [3,eq.26]
    cf.tilde_Smin_sw = min(cf.tilde_Smin_sw,cf.tilde_S);                   % [3,eq.27]
    
    % A PRIORI SPEECH ABSENCE
    tilde_Gamma_min  = (abs(Y_l).^2)./(cf.Bmin*cf.tilde_Smin);
    tilde_zeta       = cf.S./(cf.Bmin*cf.tilde_Smin);                      % [3,eq.28]
    % Speech absence
    q                = zeros(size(Y_l));
    idx_1            = find((tilde_Gamma_min <= 1) & (tilde_zeta < cf.zeta0));          % [3,eq.29]
    q(idx_1)         = ones(size(q(idx_1)));
    idx_2            = find((1 < tilde_Gamma_min) & (tilde_Gamma_min < cf.Gamma1) & (tilde_zeta < cf.zeta0));
    q(idx_2)         = (repmat(cf.Gamma1,size(q(idx_2)))-tilde_Gamma_min(idx_2))./(repmat(cf.Gamma1,size(q(idx_2)))-ones(size(q(idx_2))));
    
    % A POSTERIORI SPEECH PROBABILITY
    nu               = Gamma.*xi./(1+xi);
    cf.p             = zeros(size(Y_l));
    cf.p(q < 1)      = (1+(q(q < 1)./(1-q(q < 1))).*(1+xi(q < 1)).*exp(-nu(q < 1))).^(-1); % [3,eq.7]
    
    % PROBABILITY DRIVEN RECURSIVE SMOOTHING
    % Smoothing parameter
    tilde_alpha_d    = cf.alpha_d+(1-cf.alpha_d)*cf.p;                                     % [3,eq.11]
    % UPDATE NOISE SPECTRUM ESTIMATE
    cf.ov_Lambda_D   = tilde_alpha_d.*cf.ov_Lambda_D + (1-tilde_alpha_d).*abs(Y_l).^2;     % [3,eq.10]
    % Bias correction
    cf.Lambda_D      = cf.beta*cf.ov_Lambda_D;                                             % [3,eq.12]
end

% UPDATE MINIMUM TRACKING
cf.j = cf.j+1;
if cf.j == cf.V
    % Minimum traking for the first estimation
    % store Smin_sw
    if size(cf.Storing,2) <= cf.U
        cf.Storing = [cf.Storing cf.Smin_sw];
    else
        cf.Storing = [cf.Storing(:,2:end) cf.Smin_sw];
    end
    % Set Smin to minimum
    cf.Smin = min(cf.Storing,[],2);
    % Let Smin_sw = S
    cf.Smin_sw = cf.S;
    % Minimum traking for the second estimation
    % store Smin_sw
    if size(cf.tilde_Storing,2) <= cf.U
        cf.tilde_Storing = [cf.tilde_Storing cf.tilde_Smin_sw];
    else
        cf.tilde_Storing = [cf.tilde_Storing(:,2:end) cf.tilde_Smin_sw];
    end
    % Set Smin to minimum
    cf.tilde_Smin    = min(cf.tilde_Storing,[],2);
    % Let Smin_sw = tilde_S
    cf.tilde_Smin_sw = cf.tilde_S;
    % reset counter
    cf.j=0;
end