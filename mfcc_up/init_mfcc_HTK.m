% function [W,T] =init_mfcc_HTK(cf)
%
% Calculates the matrices for the Mel-filterbank and DFT transformations 
% used on the Mel-frequency Cepstral coefficients (MFCCs)
%
% Input:  Structure with HTK configuration translated to stft-up-tools 
%         format, see init_stft_HTK.m
%             
% Output: W Mel-filterbank matrix, T DCT matrix
%
% This version is directly based on the HTKBooks explanation and slightly
% differ from the one derived from the MASVs original code. Except for the
% fact that HTK book uses Melf = 2595*log10(1+f/700).
%
% Ramon F. Astudillo

function [W,T] = init_mfcc_HTK(cf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEL-FILTERBANK WEIGTHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EQUALY SPACED BINS IN MEL-DOMAIN (freq bins under fs/2) 
cent_mf = linspace(0,1127*log(1 + 0.5*cf.fs/700),cf.numchans+2);
% Transform them back to Fourier domain
cent_f  = (exp(cent_mf/1127)-1)*700;

% CREATE WEIGHT MATRIX OF FILTERBANK
W = zeros(cf.numchans,cf.nfft/2+1);
% For each filter center 
for m = 2:length(cent_f)-1
    % For each frequency bin
    for k = 1:(cf.nfft/2+1)
        % FIND WEIGHT         
        % Frecuency corresponding to bin k
        frec = (k-1)*cf.fs/cf.nfft;
        % If it belongs to ascending ramp of the filter
        if (frec <= cent_f(m)) && (frec >= cent_f(m-1))
            W(m-1,k) = (frec-cent_f(m-1))/(cent_f(m)-cent_f(m-1));
        % if it belongs to descending ramp of the filter
        elseif (frec <= cent_f(m+1)) && (frec >= cent_f(m))
            W(m-1,k) = 1-((frec-cent_f(m))/(cent_f(m+1)-cent_f(m)));
        else
        % Does not belong to filter    
            W(m-1,k) = 0;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCT WEIGTHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NORMAL DCT FROM SIGNAL
% NOTE: This is not exactly the same using Matlabs DCT()
T = sqrt(2.0/cf.numchans) * cos( pi/cf.numchans * (0:cf.numceps)' * ([1:cf.numchans]-0.5));

% CEPLIFTER
if cf.ceplifter ~= 0
    lift_matr = 1+cf.ceplifter/2 * sin([0:cf.numceps]'*pi/cf.ceplifter);
    % Multiply by ceplifter matrix
    T =  diag(lift_matr)*T;
end
