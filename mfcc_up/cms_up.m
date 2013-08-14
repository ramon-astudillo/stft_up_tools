% function [mu_y,Sigma_y] = cms_up(mu_x,Sigma_x,targetkind)
%
% Input:  mu_x,Sigma_x  Matrix of means and variances of the Gaussian 
%                       variable
%
% Output: mu_y,Sigma_y  Matrix of means and variances of the output 
%                       Gaussian variable
%
% Ramon F. Astudillo

function [mu_y,Sigma_y] = cms_up(mu_x,Sigma_x,targetkind)

%
% The extra targetkind parameters is for compatibility with HTK, DO NOT use
% it if you are not having HTK compatible TARGETKIND definitions
%

% Default is we perform CMS 
if nargin < 3
    targetkind = '_Z';
end

% By default, no action
mu_y    = mu_x;
Sigma_y = Sigma_x;

% APPEND CEPSTRAL MEAN SUBSTRACTION
% Calculated after the rest
% Its local to this file not loaded from any file
if ~isempty(strfind(targetkind,'_Z'))
    % Appliying that if
    % x -> N(mu1,lambda1)
    % y -> N(mu2,lambda2)
    % Then
    % a*x -> N(a*mu1,a^2*lambda1)
    % x+y -> N(mu1+mu2,lambda1+lambda2)
    % x-y -> N(mu1-mu2,lambda1+lambda2)
    %
    % BUT! If x and y correlated then
    % Var(x-y)= Varx(x) - Var(y) - 2 Cov(xy)
    % In this case the only correlation is with itself
    % Cov(x (z + ... + x + ...+ k)/N)= Var(x)/N
    %
    % number of frames
    L = size(mu_x,2);
    % MEANS
    mu_CM    = mean(mu_x,2);
    mu_y     = mu_x - repmat(mu_CM,1,L);
    % VARS
    if any(Sigma_x)
        Sigma_CM = mean(Sigma_x,2)/L;
        Sigma_y  = Sigma_x + repmat(Sigma_CM,1,L) - 2*Sigma_x/L;
    else
        Sigma_y = zeros(size(mu_y));
    end
end