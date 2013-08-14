% function cf = init_stft_HTK(cf)
%
% Convert HTK's STFT configuration parameters for HCopy into MATLAB STFT 
% parameters to be used with mfcc_up. Admits also stft-up-tools-compatible  
% parameters as input. If there is a conflict it issues an error
%
% Input:  Structure containing HTK cofiguration file parameters as fields
%
% Output: Structure with translated HTK fields, suitable for stft-up-tools
%
% Ramon F. Astudillo

function cf = init_stft_HTK(cf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK FOR SAMPLING FREQUENCEY OR SOURCERATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF SAMPLING FREQUENCY SPECIFIED
if  isfield(cf,'fs') 
    % IF ALSO SOURCERATE SPECIFIED
    if isfield(cf,'sourcerate')
        % If different
        if cf.fs ~= 1e7/cf.sourcerate         
            % Issue error 
            error('Conflicting configuration fields, FS is %d Hz while SOURCERATE is %d (%d Hz)',cf.fs,cf.sourcerate,1e7/cf.sourcerate)
        % If same
        else
            % Issue warning
            warning('Redundant configuration fields, FS is the same as 1e7/SOURCERATE')
        end
    end
    
% If SOURCERATE provided instead   
elseif isfield(cf,'sourcerate')
    % Convert it
    cf.fs = 1e7/cf.sourcerate;
    
% ELSE ERROR    
else
    error('At least FS or SOURCERATE fields must be provided')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK FOR WINDOWSIZE AND TARGETRATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF WINDOWSIZE AND TARGETRATE SPECIFIED
if isfield(cf,'windowsize') && isfield(cf,'targetrate')
    
    % IF ALSO NFFT SPECIFIED
    if isfield(cf,'nfft'); 
        % If different
        if cf.nfft ~= fix(cf.fs*cf.windowsize/1e7)
            % Issue error
            error('Conflicting configuration fields, NFFT is %d samples while WINDOWSIZE is %d (%d samples)',cf.nfft,cf.windowsize,fix(cf.fs*cf.windowsize/1e7))
            
        % If the same    
        else
            % Issue warning
            warning('Redundant configuration fields, NFFT is the same as fix(cf.fs*cf.windowsize/1e7)')
        end
    end
    
    % IF ALSO NOVERLAP SPECIFIED
    if isfield(cf,'noverlap'); 
        % If different
        if cf.nfft ~= fix(cf.fs*(cf.windowsize-cf.targetrate)/1e7)
            % Issue error
            error('Conflicting configuration fields, NOVERLAP is not coherent with WINDOWSIZE and TARGETRATE')
            
        % If the same    
        else
            % Issue warning
            warning('Redundant configuration fields, NOVERLAP or WINDOWSIZE and TARGETRATE')
        end
    end
    
    % IF HAMMING WINDOW SPECIFIED
    if strcmp(cf.usehamming,'T')
        cf.windowtype = 'hamming';
    end
    
    % Compute the parameters anyway
    cf.nfft       = fix(cf.fs*cf.windowsize/1e7);
    cf.noverlap   = fix(cf.fs*(cf.windowsize-cf.targetrate)/1e7);
    
    % CHAPUZA for compatibility with EMSP's code
    cf.overlap    = cf.noverlap;
    cf.windowsize = cf.nfft;     % Note that this is in samples not like HTKs!!  
   
    % REMOVE HTK FIELDS THAT ARE NOT NEEDED
    if isfield(cf,'targetrate')
        cf=rmfield(cf,'targetrate'); 
    end
    if isfield(cf,'usehamming')
        cf=rmfield(cf,'usehamming');
    end

% IF NO NFTT OR OVERLAP SPECIFIED
elseif ~isfield(cf,'nfft') || (~isfield(cf,'overlap') && ~isfield(cf,'noverlap'))  
    
    % Issue error
    error('At least NFFT, OVERLAP (or NOVERLAP) fields [MATLAB] OR WINDOWSIZE, TARGETRATE fileds [HTK] must be provided')
    
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY OTHER HTK OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If specified in the HTK configuration use twice the required FFT size 
if isfield(cf,'doublefft') && cf.doublefft == 'T'				
	cf.nfft=cf.nfft*2;
end

% If NFFT was not even, make it even          
cf.nfft = 2^nextpow2(cf.nfft); 