function [] = RDA_ABIDE505_rcc_q(alpha)

% Apply Restricted Connectivity-on-scalar Regression and Detect Hub-centered Edges in ABIDE505 dataset, with rank selection for the latent factor model
% Input
%   alpha: significance level of the test

%% load data
home_dir = '/gpfs/home/yt20ba/git/Res-CosReg';
addpath(genpath(home_dir))

disp("load data...")
dat = load('Cov_Mat_505_aal_harm_gh_rank67.mat');

% load the connectivity matrix
info0 = dat.info;
cmat_aal = dat.Cov_Mat_h;
% Mat_ts = dat.Mat_ts;

m = size(cmat_aal,1);
n0 = size(cmat_aal,3);
TOL = 1e-6;

% QR for rank
rs0 = zeros(n0,1);
rrs0 = zeros(n0,1);
k0 = 67;

for i = 1:n0
%     i
    A = cmat_aal(:,:,i);
    rs0(i) = rank(A);
    rrs0(i) = rank(A(:,1:k0)); 
end
ind1 = (rs0==k0)&(rrs0==k0);
sprintf([num2str(sum(ind1)),' connectivity matrices with rank = ', ...
                    num2str(k0), ' and the first ', num2str(k0),...
                 ' columns are linear independent'])


cmat = cmat_aal(:,:,ind1);
info = info0(ind1,:);

% load covariates
dx = info(:,2) - 1;
age0 = info(:, 3);
age = (age0 - mean(age0))/std(age0); % standardized age
sex = info(:,4) - 1; % 0 = Male; 1 = Female
X = [dx, age, sex];
p = size(X,2);
[m,~,n] = size(cmat);

%% Two-step testing
% consider hubs
% hub_index = [19 20 3 4 31 34 43 44 29 30 65 51 52 32 33 67 68];
hub_index = [68];
all_roi = 1:m;
others = setdiff(all_roi, hub_index);

% set up parameters
S = 1000;

for j = 1:length(hub_index)
    j
    h_roi = hub_index(j);
    hub_index1 = [h_roi,setdiff(hub_index,h_roi)];
    ord = [hub_index1, others]';
    q_set = 1:12;
    [res1, res2, CIs, B, q, BIC_vals] = spsd_test_q(X,cmat,ord,q_set,k0,S,alpha);
    
    Res{j} = res2;
    CI{j} = CIs;
end

pack_res = [];
for j = 1:length(hub_index)
    res_j = [repmat(hub_index(j),size(Res{j},1),1), Res{j}];
    pack_res = [pack_res;res_j];
end

csv_filename = sprintf("%s/results/RDA/pack_res505_%d_%.2f_q.csv", home_dir, length(hub_index),alpha);
mat_filename = sprintf("%s/results/RDA/Res505_%d_%.2f_q.mat", home_dir, length(hub_index),alpha);

% csvwrite(csv_filename, pack_res);
% save(mat_filename,'Res','CI','hub_index','q_set','BIC_vals')


end

