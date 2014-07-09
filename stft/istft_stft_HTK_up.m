% function MSE2 = istft_stft_HTK_up(MSE,config1,config)
%
% Propagation of the variances associated to a matrix of non-identically
% distributed circularly symmetric complex Gaussian random variables through
% inverse STFT and STFT transformations with different configurations. The
% approximation ignores the correlations induced by overlaping frames on the
% acros time and frequency on the final MSE, see
%
% [1] R. F. Astudillo and S. Braun and E. A. P. Habets "A Multichannel Feature 
% Compensation Approach for Robust ASR in Noisy and Reverberant Environments",
%  REVERB workshop 2014. 
%
% Input: MSE        [K, L] Variances of the STFT estimate of K freq bins and 
%                   L frames           
% 
% Input: config1    structure containing the STFT configuration of the inverse
%                   STFT, see stft_HTK..
%
% Input: config     structure containing the STFT configuration of the direct 
%                   STFT
%
% Output:MSE        
%
% Ramón F. Astudillo


function MSE2 = istft_stft_HTK_up(MSE,config1,config)

% Get parameters of STFT into a suitable format
ftsize = config1.nfft;
h      = config1.windowsize - config1.overlap; 

% Get sizes
[K1,L1] = size(MSE); 

% Get parameters of 2nd STFT
xlen   = ftsize + L1 * h;
shift2 = config.windowsize - config.overlap;
L2     = fix((xlen-config.windowsize)/shift2) + 1;
K2     = config.nfft/2+1;

% SLOW NON-BUILT IN iFFT
N      = ftsize;
% This is a bit convoluted ...
iW     = ifft(eye(N));
% Zero-padding
iW     = iW(1:config1.windowsize,:);
win    = get_window(config1);
iW     = diag(win)*iW;
% Slow non built-in DFT
W      = fft(eye(config.nfft));
% Zero padding
W      = W(:,1:config.windowsize)*diag(get_window(config1));

% Extend to symmetric spectrum
ft     = [MSE; conj(MSE((ftsize/2):-1:2,:))];
% This will store the 2nd STFT
MSE2 = zeros(2*K2-2,L2);

% Scan all possible positions of 1st STFT frames overlapping with 2nd STFT
% frames
t1 = ((1:L1)-1)*h;
t2 = ((1:L2)-1)*shift2;
dist = (repmat(t1',1,L2) + config1.windowsize) - repmat(t2,L1,1);
% Add each overlapping component for each frame shift delta
for delta = 1:config.windowsize + config1.windowsize

    % For each 2nd STFT frame start, find 1st STFT frame end at that delta
    [idx_l1,idx_l2] = find(dist == delta);
    
    % If at least one 1st STFT frame found at that delta, add it in the
    % corresponding 2nd STFT frame
    if any(idx_l1)
        % Depending on the overlap, the matrix is different.
        if delta < config.windowsize    
            %              delta
            % t1          |------>
            %  |-----------------|
            %             |---------|
            %            t2
            R = W(:,1:delta)*iW(end-delta+1:end,:);
        elseif (config.windowsize <= delta) && (delta <= config1.windowsize)
            %             delta   
            %  t1    |------------->
            %    |-----------------|
            %        |---------|
            %       t2
            R = W*iW(end-delta+1:end-delta+config.windowsize,:);
        else
            %           delta
            % t1|--------------------->
            %       |-----------------|
            %   |---------|
            %  t2
            overlap = config1.windowsize-(delta-config.windowsize);
            R       = W(:,end-overlap+1:end)*iW(1:overlap,:);
        end
        % Add those frames to total
        MSE2(:,idx_l2) = MSE2(:,idx_l2) + (abs(R).^2)*ft(:,idx_l1);
    end
    
%    if ~mod(delta,50)
%        fprintf('\rISTFT-UP %2.0f%%',100*delta/(config.windowsize + config1.windowsize-1))
%    end
end

% Get only lower part
MSE2 = MSE2(1:end/2+1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function window = get_window(cf)

% SELECT ADECUATE WINDOW
if ~isfield(cf,'windowtype')
    cf.windowtype = 'hamming';
end

switch cf.windowtype
    case 'hanning'
        window = hanning(cf.windowsize);
    case 'hamming'
        window = hamming(cf.windowsize);
    case 'kaiser'
        window = kaiser(cf.windowsize);
    case 'rectangular'
        window = rectwin(cf.windowsize);
    case 'sinus'
        a=0.54; b=-0.46; Phi=pi/cf.windowsize;
        win_rec(1:cf.windowsize)=sqrt(cf.overlap)/sqrt(cf.windowsize);
        win_sin=2*win_rec./(sqrt(4*a^2+2*b^2)).*(a+b*cos(2*pi*(1:cf.windowsize)./cf.windowsize + Phi));
        window=win_sin.';
    otherwise
        error('Unknown windowtype %s',windowtype)
end
