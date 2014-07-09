% function x = istft(X,windowsize,shift,nfft)
%
% Computes the inverse of the short time Fourier domain transform of an 
% incoming signal
%
% Input: X  Complex values STFT of K freq bins and L frames 
%            [K, L, nMics, nSignals]
%
% Output:  x  Real valued vectors of T samples [T, nSignals, nMics]
%            e.g. nMics = 2 for stereo wavs.  
%
%         cf Structure containing configuration, see init_stft_HTK.m
%
% Ramon F. Astudillo


function x = istft_HTK(X,cf)

% Get windowsize, number of frames, number of microphones and signals 
[K L nMics nSignals] = size(X);

% Initialize windowed time domain signal
x_framed             = zeros(cf.windowsize,L,nMics,nSignals);

% For each signal and microphone
for l = 1:nSignals
    for m = 1:nMics
              
        % Complete complex conjugate part
        X_tmp             = [ X(:,:,m,l) ; conj(flipud(X(2:K-1,:,m,l)))];
        
        % IDFT
        tmp_framed        = real(ifft(X_tmp,cf.nfft));
        
        % zero padding
        x_framed(:,:,m,l) = tmp_framed(1:cf.windowsize,:);
         
    end
end

% FRAMED TIME DOMAIN SIGNAL PROCESSING
window_func =  get_window(cf);
% Apply window function
x_framed = x_framed .* repmat(window_func, [1, size(x_framed,2) , nMics , nSignals]);

% Compute inverse framing
x        = iframing(x_framed,cf.windowsize,cf.windowsize - cf.noverlap);

%PRE-EMPHASIS (optional)
if isfield(cf,'preemcoef')
    x=filter(1,[1 -cf.preemcoef],x);
end


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
