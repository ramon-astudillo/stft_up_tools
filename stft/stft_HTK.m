% function X = stft_HTK(x,cf)
%
% Computes the short time Fourier domain transform (STFT) for a time domain
% signal  x in similar form to HTK.
%
% Input:  x  Real valued vectors of T samples [T, nSignals, nMics]
%            e.g. nMics = 2 for stereo wavs.  
%
%         cf Structure containing configuration, see init_stft_HTK.m
%
% Output: X  Complex values STFT of K freq bins and L frames 
%            [K, L, nMics, nSignals]
%
% Ramon F. Astudillo


function X = stft_HTK(x,cf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME DOMAIN PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OFFSET COMPENSATION FOR ADVANCED FRONT-END
if isfield(cf,'offset_comp') & cf.offset_comp 
    % OPTIONAL OFFSET COMPENSATION
    x = OffComp(x);
end
    
% PRE-EMPHASIS 
if isfield(cf,'preemcoef') & cf.preemcoef 
    x = filter([1 -cf.preemcoef],1,x);
end

% Compute framing
x_framed = framing(x,cf.windowsize,cf.windowsize - cf.noverlap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRAMED TIME DOMAIN PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get windowsize, number of frames, number of microphones and signals 
[N, L, nMics, nSignals] = size(x_framed);

% OPTIONAL REMOVING OF MEAN PER WINDOW
if isfield(cf,'zmeansource') & strcmp(cf.zmeansource,'T')
    x_framed = zmeansubs(x_framed);
end

% COMPUTE LOG ENERGY FOR ADVANCED FRONT END [ES 202 050 p.14]
if isfield(cf,'htkenergy') & ~cf.htkenergy
    
    cf.rawte = LogE(x_framed);
      
elseif isfield(cf,'htkenergy') & cf.htkenergy
    
    cf.rawte = LogE_HTK(x_framed);
    
end

% FRAMED TIME DOMAIN SIGNAL PROCESSING
window_func =  get_window(cf);
% Apply window function
x_framed = x_framed .* repmat(window_func,[1 size(x_framed,2) nMics nSignals]);

% Initialize STFT
X     = zeros(cf.nfft/2+1,L,nMics,nSignals);

% For each signal and microphone
for l = 1:nSignals    
    for m = 1:nMics
        
        % Compute DFT
        X_tmp      = fft(x_framed(:,:,m,l),cf.nfft);
        
        % Get inferior part of the spectrum
        X(:,:,m,l) = X_tmp(1:cf.nfft/2+1,:);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function[s_of]=OffComp(s_in)
%
% Performs Offset Compensation as described in the ETSI standard

function s_of = OffComp(s_in)

% Get number of samples, number of microphones and signals 
[T, nMics, nSignals] = size(s_in);

% For all microphone signals compute Offseted Version
s_of=zeros(size(s_in));
for l = 1:nSignals 
    for m = 1:nMics
        s_of(1,m,l)=s_in(1,m,l);
        for i=2:length(s_in(:,m,l))
            s_of(i,m,l)=s_in(i,m,l)-s_in(i-1,m,l)+0.999*s_of(i-1,m,l);
        end
    end
end

% function  x_framed = zmeansubs(x_framed)
%
% Per window removing of tools

function  x_framed = zmeansubs(x_framed)

% Get number of samples, number of microphones and signals 
[N, L ,nMics, nSignals] = size(x_framed);
 
% For all microphone signals compute Offseted Version
for l = 1:nSignals 
    for m = 1:nMics
        x_framed(:,:,m,l) = x_framed(:,:,m,l) - repmat(mean(x_framed(:,:,m,l),1),size(x_framed,1),1);
    end
end


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

% rawte = LogE(x_framed)
%
% Log-Energy

function rawte = LogE(x_framed)

% Get number of samples, number of microphones and signals 
[N, L ,nMics, nSignals] = size(x_framed);
 
% Initialize
rawte = zeros(L,nMics,nSignals);

for l = 1:nSignals 
    for m = 1:nMics
        
        % LOG ENERGY OF FRAME
        rawte(:,m,l) = 0.5+(16/log(2))*log((64+sum(x_framed(:,:,m,l).^2))/64);
        
        % Floor energy
        idx = find(cf.rawte<-50);
        if ~isempty(idx);  
            rawte(idx,m,l) = -50; 
        end
    end
end

% rawte = LogE_HTK(x_framed)
%
% Log-Energy according to HTK

function rawte = LogE_HTK(x_framed)

% Get number of samples, number of microphones and signals 
[N, L ,nMics, nSignals] = size(x_framed);

% Initialize
rawte = zeros(L,nMics,nSignals);

for l = 1:nSignals 
    for m = 1:nMics

        % ELSE HTKs LOG ENERGY  [HTKBook321 p.65 / ES 201 108 p.10]
        % LOG ENERGY OF FRAME
        rawte(:,m,l)=log(sum(tmp_framed.^2));
        
        % Floor energy
        if ~isfield(cf,'enormalise') | strcmp(cf.enormalise,'T')
            rawte = cf.rawte - repmat(max(rawte)+1,1,length(rawte));
        end

    end
end