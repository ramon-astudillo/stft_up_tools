stft_up_tools
=============

This is part of the code I used for speech enhancement and Short-Time Fourier Transform Uncertainty Propagation (STFT-UP) in my thesis and related papers. I developed the first versions while doing my Ph.D. at the Technische Universtitaet Berlin. The code therefore benefited from the previous works of Prof. Dorothea Kolossa and her students in particular Diep Huyn and Alexander Klimas. Latest version or a link to an online repository should be found here

        http://www.astudillo.com/ramon/research/stft-up/ 

Feel free to use it entirely or partially for your own code. Just reference the corresponding papers (and the repo URL for e.g. the IMCRA implementation). 

Current Version contains code for:

**Improved Minima Controlled Recursive Averaging noise variance estimator (IMCRA)** 

        example2.m

This is a noise variance estimator based on minimum statistics and though for non-stationary additive noises. It also provides a probabilistic voice activity detection, see 

        [1] I. Cohen, "Noise Spectrum Estimation in Adverse Environments: Improved Minima Controlled 
        Recursive Averaging", in IEEE Trans. on Speech and Audio Processing, Vol 11 (5), pp 1063-6676,
        2003  

**Uncertainty Propagation for the MFCC Features**

        example1.m

This propagates an uncertain description of the STFT of a signal into MFCC domain. The model used for STFT uncertainties is a complex Gaussian distribution. The example compares the analytic solutions for STFT-UP using the Log-normal approximation or the Uncented Transform with the Monte Carlo solution. To cite Uncertainty Propagation in general please use

        [2] R. F. Astudillo, "Integration of short-time Fourier domain speech enhancement and observation
        uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische 
        Universitaet Berlin, 2010.

**A MMSE-MFCC Estimator derived with STFT-UP**

        example2.m

Here the previous STFT-UP apprximation is used to propagate the residual mean square error (MSE) of the Wiener filter, thus attaining a minimum MSE (MMSE) estimator in MFCC domain (MMSE-MFCC), see

        [3] R. F. Astudillo, R. Orglmeister, "Computing MMSE Estimates and Residual Uncertainty directly
        in the Feature Domain of ASR using STFT Domain Speech Distortion Models", IEEE Transactions on
        Audio, Speech and Language Processing, Vol. 21 (5), pp 1023-1034, 2013

for details. This matlab script is also though to process a batch of files and thus allow to reproduce the front-end used in this paper. To reproduce experiments using HTK, you will also need the voicebox toolbox   

        http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html 

and the patches to modify HTK to perform Uncertainty Decoding or Modified Imputation, see

        http://www.astudillo.com/ramon/research/stft-up

**Sparsity Based Uncertainty Model and Propagation**

        example1b.m

The approach in [2] models the uncertainty over the value of each Fourier coefficient under and additivity assumption. Here we model uncertainty over either source or noise being active under an sparsity assumption. This leads to a scaled Bernoulli uncertainty model than can also be (approximately) propagated. This is described in

        [4] Francesco Nesta, Marco Matassoni, Ramon Fernandez Astudillo, A FLEXIBLE SPATIAL BLIND SOURCE 
        EXTRACTION FRAMEWORK FOR ROBUST SPEECH RECOGNITION IN NOISY ENVIRONMENTS, In 2nd International 
        Workshop on Machine Listening in Multisource Environments (CHiME), pages 33-38, June 2013

Ramón F. Astudillo, last revision Oct 2013
