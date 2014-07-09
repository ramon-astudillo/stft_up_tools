% function x = iframing(x_framed,windowsize,shift)
% 
% Recovers time domain signal x from framed signal x_framed by overlap and add
% with windows of length WINDOWSIZE and shifting SHIFT samples each time.
%
% Ramon F. Astudillo 

function x = iframing(x_framed,windowsize,shift)

% Get length of file, number of microphones and signals 
[N, L, nMics, nSignals] = size(x_framed);

% Initialize time domain signal
x = zeros(L*shift + windowsize,nMics,nSignals,1);

% COMPUTE FRAMES
% For all microphone signals 
for s = 1:nSignals
    for m = 1:nMics 
        % For each frame
        for l=1:L
             
           % Overlap and add 
           x((l-1)*shift + 1:(l-1)*shift + 1 + windowsize - 1,m,s) = x((l-1)*shift + 1:(l-1)*shift + 1 + windowsize-1,m,s) + x_framed(:,l,m,s);
        
        end
    end
end

% Squeeze empty dimensions
x = squeeze(x);
