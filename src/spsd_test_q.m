function [res1, res2, CIs, B, q, BIC_vals] = spsd_test_q(X,cmat,ord,q_set,k0,S,alpha)
% two step testing + boostrap multiple testing
% BIC model selection to choose number of factor q
% Input:
    % X : covariates
    % cmat : connectivity matrices
    % ord : index for reordering
    % q_set: candidate set of number of factors
    % k0 : rank of the connectivity matrices
    % S: boostrap times
    % alpha: significance level
% Output:
    % res1: significant covariate with connectivity in step 1
         % each row is [covariate index, ROI-1 index (1), p value]
    % res2: significant covariate correlated to the connectivity in step 2
         % each row is [covariate index, ROI-k index, lower CI, upper CI]
    % CIs: (p x m) matirx of pvalues
    % B: estimated coeffecient matrix
    % q: selected number of factors
    % BIC_vals: BIC val for q_set

%%
[m,~,n] = size(cmat);
p = size(X,2);

%% Reorder
disp('Reorder')
l = k0*(m+m-k0+1)/2;
Y = zeros(n,l);

% mapping
for i = 1:n
    A = cmat(:,:,i);
    
    %reordering
    Pe = eye(m);
    Pe = Pe(ord,:);  % shuffle the rows of P according to the permutation Perm
    
    A = Pe*A*Pe';
    Y(i,:) = gh_map(A,k0);
end

%% Estimation
disp('Estimation')
%Centered X and Y
mu_X = mean(X);
mu_Y = mean(Y);
X_n = X - mu_X;
Y_n = Y - mu_Y;
XX_inv = inv(X_n'*X_n);
B = XX_inv*(X_n'*Y_n); %p x l
% B = (X_n'*X_n)\(X_n'*Y_n); %p x l
% B0 = mu_Y - mu_X*B;
Res = Y_n - X_n * B;
Sigma = (Y_n-X_n*B)'*(Y_n-X_n*B)/n;

[U, d] = eig(Sigma, 'vector');
% resort d with descending
[d,ind] = sort(d,'descend');
U = U(:,ind);

% Candidate values for q
% q_set = 5:5:50; % set maximum number of factors
BIC_vals = zeros(length(q_set), 1);
logL_vals = zeros(length(q_set), 1);

for r = 1:length(q_set)
    q = q_set(r); % number of factors

    % Estimate sigma and construct V_q
    sigma_q = mean(d(q+1:end));
    Lq = diag(d(1:q));
    Uq = U(:,1:q);
    V_q = Uq * sqrtm(Lq - sigma_q * eye(q));
    Sigma_q = V_q * V_q' + sigma_q * eye(l);

    % Compute log-likelihood
    logdet_Sigma = sum(log(eig(Sigma_q) + 1e-10)); % stability
    Mahalanobis = trace((Res / Sigma_q) * Res'); 
    logL = -n * l / 2 * log(2*pi) - n/2 * logdet_Sigma - 0.5 * Mahalanobis;

    % Degrees of freedom = q*l (loadings) + 1 (sigma)
    df = p*l + q * l - q*(q-1) + 1;
    BIC = -2 * logL + log(n) * df;

    logL_vals(r) = logL;
    BIC_vals(r) = BIC;
end

% Select optimal q
[~, idx_opt] = min(BIC_vals);
q = q_set(idx_opt);
fprintf('Selected number of latent factors: q = %d\n', q);


sigma = sum(d((q+1):l))/(l-q);
Ur = U(:,1:q); %top q eigenvectors
Lr = diag(d(1:q)); %largest eigen values
V = Ur*sqrtm(Lr-sigma*diag(ones(q,1)));
C = V*V' + sigma * diag(ones(l,1));

U = Y_n - X_n*B; % residuals
R = U'*U/(n-p);
b0 = 0;

%% 1 Test first beta beta_{j1}
disp('Step 1: T-test for beta_{j1}')
T1 = zeros(p,1);
p1 = zeros(p,1);
for j = 1:p
    b = B(j,1);
    se1 = sqrt(R(1,1)*XX_inv(j,j));
    T1(j) = (b - b0) / se1;
    % Compute the two-tailed p-value for the T-test
    p1(j) = 2 * (1 - tcdf(abs(T1(j)), n - p));
end
ind1 = find(p1 < alpha);
res1 = [ind1, ones(length(ind1),1), p1(ind1)];

%% 2 Bootstrap Testing MaxT
disp('Step 2: Boostrap-based Multiple Testing')
rng("default")

% set the number of bootstrap samples
% S = 500;
e_indices = randi(n,[n, S]);

% Initialize array to store max values
M_stars = zeros(S, 1);
% alpha = 0.1;

CIs = zeros(m,3,p);

for j = 1:p
    % j = 2
    se = sqrt(diag(R(2:m,2:m))*XX_inv(j,j));
    Ts = zeros(m,S);
    % j

    % Bootstrap iterations
    for s = 1:S
        % s
        % Step 1: Sample e_i* with replacement from residuals
        idx = e_indices(:,s);
        e_star = U(idx,:);
        X_n_s = X_n(idx,:);

        XX_inv_s = inv(X_n_s'*X_n_s);
        B_e = XX_inv_s*(X_n_s'*e_star);
        R_e = (e_star-X_n_s*B_e)'*(e_star-X_n_s*B_e)/(n-p);

        % Step 2: Compute max |T(j,k)(e_i*)| for each ROI-k
        for k = 2:m
            b = B_e(j,k);
            se_e = sqrt(R_e(k,k)*XX_inv_s(j,j));
            Ts(k,s) = abs(b - b0)/se_e;
        end

        max_T = max(Ts(:,s)); 

        % Store max value
        M_stars(s) = max_T;
    end
    
    % Step 4: Compute empirical 1-alpha quantile and construct simultaneous CI
    c_alpha_hat = quantile(M_stars, 1 - alpha);
    CIs(2:m,:,j) = [B(j,2:m)', B(j,2:m)' - c_alpha_hat*se, B(j,2:m)' + c_alpha_hat*se];
end


covers_zero = squeeze(CIs(2:m, 2, : ) <= 0 & CIs(2:m, 3, :) >= 0);
sum(~covers_zero)

[ind_ROI,ind_c] = find(~covers_zero); % ind_ROI is in ord = [hub_index1, others]'
CI_re = reshape(permute(CIs(2:m,:,:), [1,3,2]),[(size(CIs,1)-1)*size(CIs,3),size(CIs,2)]);
ci_sig = CI_re(sub2ind(size(covers_zero), ind_ROI, ind_c),:);

if isempty(ind_c)
    res2 = [];
else
    res2 = [ind_c, ord(ind_ROI+1), ci_sig]; 
    % ind_ROI didn't count the first hub index
end

% res = [res1;res2];

end