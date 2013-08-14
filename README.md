stft_up_tools
=============

This is part of the code I used for speech enhancement and Short-Time Fourier Transform Uncertainty Propagation (STFT-UP) in my thesis and related papers. I developed the first versions while doing my Ph.D. at the Technische Universtitaet Berlin. The code therefore benefited from the previous works of Prof. Dorothea Kolossa and her students in particular Diep Huyn and Alexander Klimas. Latest version or a link to an online repository should be found here

        http://www.astudillo.com/ramon/research/stft-up/ 

Feel free to use it entirely or partially for your own code. Just reference the corresponding papers (and the URL above for e.g. the IMCRA implementation). For STFT-UP in general you can cite 

        [1] R. F. Astudillo, "Integration of short-time Fourier domain speech enhancement and observation
        uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische 
        Universitaet Berlin, 2010.

for the MMSE-MFCC and other MMSE estimators attained by propagating the Wiener posterior using STFT-UP you can cite

        [2] R. F. Astudillo, R. Orglmeister, "Computing MMSE Estimates and Residual Uncertainty directly in
        the Feature Domain of ASR using STFT Domain Speech Distortion Models", IEEE Transactions on Audio, 
        Speech and Language Processing, Vol. 21 (5), pp 1023-1034, 2013

To reproduce experiments using HTK, you will also need the voicebox toolbox   

        http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html 

see also

        http://www.astudillo.com/ramon/research/

for addittional references


Current Version contains code for:

**Single Channel Speech Enhancement** 

This includes the IMCRA noise variance estimator 

        [3] I. Cohen, "Noise Spectrum Estimation in Adverse Environments: Improved Minima Controlled Recursive
        Averaging", in IEEE Trans. on Speech and Audio Processing, Vol 11 (5), pp 1063-6676, 2003  

along with MMSE-LSA and MMSE-MFCC (attained through STFT-UP) estimators.

**MFCC feature extraction with uncertainty propagation**

This includes the STFT-UP implementation for a MFCC feature extraction with delta, acceleration features and cepstral mean subtraction (see example1.m). With aid of the voicebox toolbox the code also allows to replace HCopy with a custom matlab function (see example2.m). This can be used to reproduce the experiments in [1] 

Ramón F. Astudillo, last revision Aug 2013
