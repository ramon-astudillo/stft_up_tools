% function [mu_y,Sigma_y] = log_up_ut(mu_x,Sigma_x,melfloor,diagcov_flag,kappa)
%
% Computes the propagation through the logarithm using the unscented
% transform
% 
% Input: mu_x, Sigma_x  Mean and covariance
%        melfloor       Lower bound for the Mel-filterbank 
%        diagcov-flag   Forces ignoring covariances
%        kappa          Scaling factor for the Uncented Transform 
%
% Ouput: mu_y, Sigma_y  Mean and variance of the logarithm of the input
%                       variable
%        
% Please refer to Chapter 6 of 
%
%   [Astudillo2010] R. F. Astudillo, "Integration of short-time fourier domain speech enhancement and observation uncertainty techniques for robust automatic speech recognition", Ph.D. dissertation, Technische Universitaet Berlin, 2010.
% 
% For details. This is available online under: 
%
%   http://opus.kobv.de/tuberlin/volltexte/2010/2676/pdf/astudillo_ramon.pdf 
%
% Ramon F. Astudillo

function [mu_y,Sigma_y] = log_up_ut(mu_x,Sigma_x,melfloor,diagcov_flag,kappa)

% Get size
[J,L] = size(mu_x);

% ARGUMENT HANDLING
if nargin < 5
	kappa        = 3-J;
end

if nargin < 4
	kappa        = 3-J;
	diagcov_flag = 1;
end

% CHECK FOR DIAGONAL MATRIX
if diagcov_flag
    
    % INITIALIZE
    mu_y    = zeros(size(mu_x));
    Sigma_y = zeros(size(mu_x));
    
    for l=1:L
        % Get Sigma points
        [S_UT,W_UT]  = get_UT_pdf(mu_x(:,l),Sigma_x(:,l),kappa);
        % replicate
        W_UT         = repmat(W_UT,size(mu_x,1),1);
        % Transform them
        LS_UT        = log(max(S_UT,melfloor));
        % Estimate Mean
        mu_y(:,l)    = sum(W_UT.*LS_UT,2);
        % Estimate Var
        Sigma_y(:,l) = sum(W_UT.*(LS_UT -  repmat(mu_y(:,l),1,size(LS_UT,2))).^2,2);
    end    
    
else

    % INITIALIZE
    mu_y    = zeros(size(mu_x));
    Sigma_y = zeros(size(mu_x,1),size(mu_x,1),size(mu_x,2));
    
    for l=1:L
        % Get Sigma points
        % Note: It might be the case that the covariance matrix is not
        % positive definite. This can be patched but is here left untouched
        % as it is kind of a CHAPUZA.
        [S_UT,W_UT] = get_UT_pdf(mu_x(:,l),Sigma_x(:,:,l),kappa);
        % replicate
        W_UT        = repmat(W_UT,size(mu_x,1),1);
        % Transform them through the non-linearity
        LS_UT       = log(max(S_UT,melfloor));
        % Estimate Mean
        mu_y(:,l)   = sum(W_UT.*LS_UT,2);
        % Estimate Cov
        delta       =  LS_UT -  repmat(mu_y(:,l),1,size(LS_UT,2));
        % For each sigma point
        for i=1:size(LS_UT,2)
            Sigma_y(:,:,l) =  Sigma_y(:,:,l) + W_UT(1,i)*(delta(:,i) * delta(:,i)');
        end
        
    end  

end


function [S_UT,W_UT] = get_UT_pdf(mu,Sigma,kappa,model_independently)

% GET SIZE OF INPUT
N = size(mu,1);

% ARGUMENT HANDLING
if nargin < 4
    % Default is to model a [N 1] mu vector as a 2*N+1-points distribution
    % even if the covariance is diagonal. Now it is also possible to solve
    % each expectation independently* in this case having N 3-points
    % distributions. Set this flag to 1 for this.
    %
    % (*) But only for point-wise non-linearities!
    model_independently = 0;
end

if nargin < 3
    model_independently = 0;
    % Default is Gaussian
    kappa = 3-N ;
end

% ARGUMENT HANDLING
if isempty(kappa)
    % Default is Gaussian
    kappa = 3-N;
end

% SECURITY CHECK (FOR SQUARE MATRICES)
if size(Sigma,1) ==  size(Sigma,2)
    [dum,not_positive_definite] = chol(Sigma);
    if not_positive_definite
        error('The Covariance matrix is not positive definite')
    end
end

% IF THE INPUT DISTRIBUTION IS ONE-DIMENSIONAL, IT JUST NEEDS 3 SIGMA
% POINTS
if N == 1 
    % DETERMINE WEIGTHS 
    W_UT = [kappa/(N + kappa) repmat(0.5/(N + kappa),1,2*N)];
    % DETERMINE SIGMA POINTS
    S_UT = [mu mu + sqrt((N + kappa)*diag(Sigma)) mu - sqrt((N + kappa)*diag(Sigma))]; 
    
% ELSE IF THE INPUT DISTRIBUTION IS N-DIMENSIONAL BUT STATISTICALLY
% INDEPENDENT AND WE WISH TO MODEL EACH OF THE N VARIABLES WITH ITS OWN SET
% OF 3 SIGMA POINTS
elseif model_independently
    % FORCE N TO BE 1
    N = 1; 

    % DETERMINE WEIGTHS
    W_UT = repmat([kappa/(N + kappa) repmat(0.5/(N + kappa),1,2*N)],size(mu,1),1);
    % DETERMINE SIGMA POINTS
    S_UT = [mu mu + sqrt((N + kappa)*diag(Sigma)) mu - sqrt((N + kappa)*diag(Sigma))];

% ELSE IF THE INPUT DISTRIBUTION IS N-DIMENSIONAL, WE WISH TO MODEL IT AS
% A SINGLE DISTRIBUTION OF 2*N+1 SIGMA POINTS AND IT IS A FULL COVARIANCE!!
elseif size(Sigma,1) ==  size(Sigma,2) 
    % DETERMINE WEIGTHS
    W_UT = [kappa/(N + kappa) repmat(0.5/(N + kappa),1,2*N)];
    % DETERMINE SIGMA POINTS
    % We need to compute the matrix square root 
    MSR  = sqrtm((N+kappa)*Sigma);
    S_UT = [      mu        repmat(mu,1,N) + MSR'   repmat(mu,1,N) - MSR'];  
else
    % DETERMINE WEIGTHS
    W_UT = [kappa/(N + kappa) repmat(0.5/(N + kappa),1,2*N)];
    % DETERMINE SIGMA POINTS
    % We need to compute the matrix square root 
    MSR  = diag(sqrt((N+kappa)*Sigma));
    S_UT = [      mu        repmat(mu,1,N) + MSR'   repmat(mu,1,N) - MSR'];  
end