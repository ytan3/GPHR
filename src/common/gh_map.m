function L_p = gh_map(S, k0)
%% bijective map from SPSD with k0 rank to its log-Cholesky factor
%  written by Yuanyao Tan
%  updated on 2024/03/27
% 
% Input:
%   S: m x m symmetric semi-definite matrix with fixed rank k0
%   k0: rank of this matrix
% Output:
%   L_p: log-Cholesky factor of S
%%

m = size(S,1);
l = k0*(m+m-k0+1)/2;

[R0, ~] = cholcov(S); % R0: k0 x m, S = R0 * R0^T
[Q, R] = qr(R0);
L_m = R(1:k0,:)';
% S_m = L_m*L_m';
neg = -(diag(L_m) < 0);
neg(neg==0) = 1;
L_m = tril(L_m);
L1_m = L_m * diag(neg); % correspondg to L
L1_m(1:(m+1):end) = log(diag(L1_m));
vech_L1 = L1_m(tril(true(size(L1_m)))); % elements from and below the main diagonal
L_p = vech_L1;


% % Map S to L'
% [U, d] = eig(S, 'vector'); % A == U*diag(d)*U', return eigenvalues in a column vector 
% d(d < 0) = 0;  % Set negative eigenvalues of A to zero
% 
% % Q*R == sqrt(D)*U', so A == U*diag(d)*U'== R'*Q'*Q*R == R'*R == L*L'
% [~, R] = qr(diag(sqrt(d))*U');
% L = R';
% 
% neg = -(diag(L) < 0);
% neg(neg==0) = 1;
% L1 = L * diag(neg);
% L1(1:(m+1):end) = log(diag(L1));
% vech_L1 = vech(L1); % elements from and below the main diagonal
% L_p = vech_L1(1:l);

end