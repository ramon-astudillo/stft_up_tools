% function x_framed = framing(x,windowsize,shift)
% 
% Extracts overlapping windows from a time domain signal x with windows of
% length WINDOWSIZE and shifting SHIFT samples each time.
%
% Ramon F. Astudillo 

function x_framed = framing(x,windowsize,shift)

% Get length of file, number of microphones and signals 
[T, nMics, nSignals] = size(x);

% Compute corresponding number of frames  
L = fix((T-windowsize)/shift) + 1;           % Number of frames. Must satisfy (L-1)*(shift) + windowsize < T
    
% Now we create index for each window by cloning the index for the first     
% window N times and adding the offset for each window.    
Index    = (repmat(1:windowsize,L,1) + repmat((0:(L-1))'*(shift),1,windowsize))';

% Initialize framed time domain signal
x_framed = zeros(windowsize,L,nMics,nSignals);

% COMPUTE FRAMES
% For all microphone signals 
for l = 1:nSignals    
    for m = 1:nMics    
        
        % DIVIDE SIGNAL INTO A MATRIX OF OVERLAPPING WINDOWS 
        % We will need this dummy var
        dummy_var  = x(:,m,l);
        
        % FRAMING MATRIX
        tmp_framed = dummy_var(Index);
              
        % WINDOWING 
        x_framed(:,:,m,l) = tmp_framed;        
    end
end
% Eliminate empty dimensions   
x_framed = squeeze(x_framed);  
