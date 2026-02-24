function S = gh_inv_map(L_p,k0)
%% bijective map from log-Cholesky factor to SPSD with k0 rank
%  written by Yuanyao Tan
%  updated on 2024/03/26
%
% Input:
%   L_p: log-Cholesky factor of S
%   k0: rank of this matrix
% Output:
%   S: m x m symmetric semi-definite matrix with fixed rank k0
% 
%%

l = length(L_p);
m = (2 * l / k0 + k0 - 1) / 2;

% Map from L' to S
l0 = sum(sum(tril(true(m))));
vech_L1 = zeros(l0,1);
vech_L1(1:l) = L_p;
L1 = zeros(m,m);
L1(tril(true(m))) = vech_L1;
L = L1;
diag_L = diag(L);
L(1:(m+1):end) = exp(diag_L);
L = L(:,1:k0);
S = L * L';


end