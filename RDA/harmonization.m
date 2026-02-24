
%% Loading
home_dir = '/Users/yuanyaotan/Library/CloudStorage/OneDrive-FloridaStateUniversity/Git/Res-CosReg';
data_dir = '/Users/yuanyaotan/Library/CloudStorage/OneDrive-FloridaStateUniversity/Resource/Data/ABIDE';
scrip_dir = '/Users/yuanyaotan/Library/CloudStorage/OneDrive-FloridaStateUniversity/Research/Spsd-Reg/Code/ComBatHarmonization/Matlab/scripts';
addpath(genpath(home_dir))
addpath(scrip_dir)

%% Harmonization
% The fourth argument, 1, indicates parametric adjustements; 
% 0, non-parametric adjustments
% data_harmonized = combat(dat, batch, mod, 0);

load([data_dir,'/connec_mat/Cov_Mat_aal_505.mat']) % 505 subjects in 8 sites

Cov_Mat = Cov_Mat_approx;
atlas = "aal";
n = size(Cov_Mat,3);
m = size(Cov_Mat,1);
k0 = 67;

l = k0 * (m+m-k0+1)/2;
Y = zeros(l, n);

for i = 1:n
    c_mat = Cov_Mat(:,:,i);
    Y(:,i) = gh_map(c_mat, k0);
end

% only age, sex
mod = [info(:,3), info(:,4)-1];
% 1 for parametric way, 0 nonparametric way
[data_harmonized, std_harmonized] = combat(Y, center_ids, mod, 1);
% data_harmonized = combat(Y, center_ids, mod, 0);


% Initialize Cov_Mat
Cov_Mat_h = zeros(m, m, n);
Cov_Mat_h_std = zeros(m, m, n);
ranks = zeros(1, n);

% Loop over each column of Y
for i = 1:n
    y = data_harmonized(:, i);
    y_std = std_harmonized(:, i);
    Cov_Mat_h(:, :, i) = gh_inv_map(y,k0);    
    Cov_Mat_h_std(:, :, i) = gh_inv_map(y_std,k0); 
    ranks(i) = rank(Cov_Mat_h(:, :, i));  % Compute rank of the i-th matrix
end

% output_path = sprintf('%s/connec_mat/Cov_Mat_505_%s_harm_gh_rank67.mat', data_dir, atlas);
% save(output_path,'Cov_Mat_h','centers','center_ids','nameMap','sub_ids','info','mod')


%% Check the harmonization effect (KW test)
data_noharm = load([data_dir,'/connec_mat/Cov_Mat_aal_505.mat'], 'Cov_Mat_approx','center_ids','centers','info'); % 505 subjects in 8 sites
data_harm = load([data_dir,'/connec_mat/Cov_Mat_505_aal_harm_gh_rank67.mat'], 'Cov_Mat_h', 'Cov_Mat_h_std', 'center_ids'); % 505 subjects in 8 sites

Cov_Mat_approx = data_noharm.Cov_Mat_approx;
Cov_Mat_h = data_harm.Cov_Mat_h;
Cov_Mat_h_std = data_harm.Cov_Mat_h_std;
center_ids = data_noharm.center_ids;
centers = data_noharm.centers;

[m, ~, n] = size(Cov_Mat_approx);
l0 = m * (m-1) / 2;

% vech FC
Y = zeros(l0, n);
Y_h = zeros(l0, n);
Y_h_std = zeros(l0, n);
for i = 1:n
    c_mat = Cov_Mat_approx(:, :, i);
    tril_c = tril(c_mat, -1); 
    Y(:,i) = tril_c(tril_c~=0);

    c_mat_h = Cov_Mat_h(:, :, i);
    tril_c = tril(c_mat_h, -1); 
    Y_h(:,i) = tril_c(tril_c~=0);

    c_mat_std = Cov_Mat_h_std(:, :, i);
    tril_c = tril(c_mat_std, -1); 
    Y_h_std(:,i) = tril_c(tril_c~=0);
end


% kruskalwallis before harmonization
pvals = zeros(l0,1);
for j = 1:l0
    dat_j = Y(j,:); % j-th element over n subjects
    [p, tbl, stats] = kruskalwallis(dat_j, centers, "off");
    % [p, tbl, stats] = kruskalwallis(dat_j, centers, "on");

    pvals(j) = p;
end
% Apply FDR correction
pvals_fdr = mafdr(pvals, 'BHFDR', true);
sig_noharm = sum(pvals_fdr < 0.05) / l0

% kruskalwallis after harmonization
pvals_h = zeros(l0,1);
for j = 1:l0
    dat_j = Y_h(j,:); % j-th element over n subjects
    [p, ~, ~] = kruskalwallis(dat_j, center_ids, "off");
    pvals_h(j) = p;
end
pvals_h_fdr = mafdr(pvals_h, 'BHFDR', true);
sig_harm = sum(pvals_h_fdr < 0.05) / l0


%% Density Plot for each site
v = 433;
dat = Y(v, :);
dat_h = Y_h(v, :);
dat_h_std = Y_h_std(v, :);
[p, tbl, stats] = kruskalwallis(dat, centers, "on");
[p_h, ~, ~] = kruskalwallis(dat_h, centers, "on");


% Site names
sites = unique(centers);
% get dataset for each site
dat_site = cell(length(sites),1);
dat_site_h = cell(length(sites),1);
dat_site_h_std = cell(length(sites),1);

for t = 1:length(sites)

    dat_t = Y(v, strcmp(centers, sites{t}));
    dat_site{t} = dat_t;

    dat_t_h = Y_h(v, strcmp(centers, sites{t}));
    dat_site_h{t} = dat_t_h;

    dat_t_h_std = Y_h_std(v, strcmp(centers, sites{t}));
    dat_site_h_std{t} = dat_t_h_std;
end


colors = [0 0 1;        % Leuven - Blue
          0 1 1;        % MaxMun - Cyan
          0 1 0;        % NYU - Green
          0.5 1 0;      % SDSU - Light Green
          1 0.6 0;      % Stanford - Orange
          1 0 0;        % Yale - Red
          0.5 0.6 0;    % UCLA - Olive Green
          0.9 0 0.9];   % KKI - Magenta


% Initialize arrays to store maximum density values for each site
num_sites = length(sites);
max_density = zeros(num_sites, 1);
min_density = zeros(num_sites, 1);
max_x = -inf;
min_x = inf;

% Calculate density ranges and x-limits
for t = 1:num_sites
    % Compute density for original data
    [f_before, x_before] = ksdensity(dat_site{t});
    [f_after, x_after] = ksdensity(dat_site_h{t});
    [f_after_std, x_after_std] = ksdensity(dat_site_h_std{t});

    % Store overall density range
    max_density(t) = max(max(f_before), max(f_after));
    min_density(t) = min(min(f_before), min(f_after));
    
    % Update global x-limits
    min_x = min([min_x, min(x_before), min(x_after)]);
    max_x = max([max_x, max(x_before), max(x_after)]);
end



figure('Position', [0 100 600 1000]);
for t = 1:length(sites)
    subplot(8,1,t)
    dat = dat_site{t};
    [f,xi] = ksdensity(dat);
    
    % Plot the density as a filled curve
    fill([xi fliplr(xi)], [f zeros(size(f))], ...
        colors(t, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    xlim([min_x, max_x])  % Apply consistent x-limits

    % Remove y-axis tick marks
    set(gca, 'YTick', []);

    hold on
    set(gca, 'YTick', [])    
    set(gca, 'Color', 'none')  % Transparent background
    ylim([0, max_density(t)]); % Use site-specific ylim
    box off                 % Turn off the box to eliminate unnecessary borders
    ylabel(sites{t}, 'FontSize', 14)
end
xlabel('Before Harmonization', 'FontSize', 14)
filename = sprintf('%s/results/RDA/real_505/density_v%s_beforeHarm.png', home_dir, int2str(v)); % Replace variable_name with your actual variable
% exportgraphics(gcf, filename, 'Resolution', 300);


figure('Position', [600 100 600 1000]);
for t = 1:length(sites)
    subplot(8,1,t)
    dat = dat_site_h{t};
    [f,xi] = ksdensity(dat);
    % Plot the density as a filled curve
    fill([xi fliplr(xi)], [f zeros(size(f))], ...
        colors(t, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    xlim([min_x, max_x])  % Apply consistent x-limits

    set(gca, 'YTick', [])
    set(gca, 'Color', 'none')  % Transparent background
    ylim([0, max_density(t)]); % Use site-specific ylim
    box off                 % Turn off the box to eliminate unnecessary borders
    ylabel(sites{t}, 'FontSize', 14)
    hold on;
end

xlabel('After Harmonization', 'FontSize', 14)
filename = sprintf('%s/results/RDA/real_505/density_v%s_afterHarm.png', home_dir, int2str(v)); % Replace variable_name with your actual variable
% exportgraphics(gcf, filename, 'Resolution', 300);


%% Density Plot for each site overlap

figure('Position', [0 100 600 400]);
legend_handles = gobjects(length(sites), 1); % Preallocate legend handles

for t = 1:length(sites)
    % subplot(8,1,t)
    dat = dat_site{t};
    [f,xi] = ksdensity(dat);

    % Plot the density as a filled curve
    fill([xi fliplr(xi)], [f zeros(size(f))], ...
        colors(t, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on

    % % Plot the density as a filled curve
    % fill([xi fliplr(xi)], [f zeros(size(f))], ...
    %     colors(t, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    legend_handles(t) = plot(nan, nan, 'Color', colors(t, :), 'LineWidth', 2);
end

% Formatting
set(gca, 'Color', 'none'); % make background transparent
box off % Remove border
xlim([min_x, max_x])  % Apply consistent x-limits
ylim([min_density(t), max_density(t)]); % Use site-specific ylim
xlabel('Before Harmonization', 'FontSize', 14)
ylabel('Density', 'FontSize', 14);
legend(legend_handles, sites, 'Location', 'northeast', 'FontSize', 10);
filename = sprintf('%s/results/RDA/real_505/density_overlap_v%s_beforeHarm.png', home_dir, int2str(v)); % Replace variable_name with your actual variable
% exportgraphics(gcf, filename, 'Resolution', 300);


figure('Position', [600 100 600 400]);
for t = 1:length(sites)
    % subplot(8,1,t)
    dat = dat_site_h{t};
    [f,xi] = ksdensity(dat);

    % Plot the density as a filled curve
    fill([xi fliplr(xi)], [f zeros(size(f))], ...
        colors(t, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on

   legend_handles(t) = plot(nan, nan, 'Color', colors(t, :), 'LineWidth', 2);
end

set(gca, 'Color', 'none');
box off % Remove border
xlim([min_x, max_x])  % Apply consistent x-limits
ylim([min_density(t), max_density(t)]); % Use site-specific ylim
xlabel('After Harmonization', 'FontSize', 14)
ylabel('Density', 'FontSize', 14);
legend(legend_handles, sites, 'Location', 'northeast', 'FontSize', 10);
filename = sprintf('%s/results/RDA/real_505/density_overlap_v%s_afterHarm.png', home_dir, int2str(v)); % Replace variable_name with your actual variable
exportgraphics(gcf, filename, 'Resolution', 300);

