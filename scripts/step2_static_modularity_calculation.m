%% Static global modularity calculation

%Side calculation of static whole brain modularity on a given set of regions.
%Useful measure of whole brain flex for correlations.

%Reproducing the following since code not available.
%Daws, R.E., Timmermann, C., Giribaldi, B. et al. Increased global integration in the brain after psilocybin therapy for depression. Nat Med 28, 844–851 (2022).
%https://doi.org/10.1038/s41591-022-01744-z


%Uses BCT and genlouvain toolbox

%% Init
close all
clear all

%% Path
%main repo
paths.home = (cd);
addpath(genpath(paths.home))
%adds bct and genlouvain tools


paths.out = fullfile(paths.home,'results','dynamic_flex');
if ~exist(paths.out)
    mkdir(paths.out);
end


%% Load Data
load(fullfile(paths.home,'data','ROI_timeseries_PSC_noGSR_scrub.mat'));
%same as step 1 .mat
time_series_data = roi_timeseries.schaefer232;

n_sub = size(time_series_data, 1);
n_ses = size(time_series_data, 2);
n_workers = 10;

%Turn timeseries structure into transformed roi x roi structure
M = {n_sub,n_ses};
for sub = 1:n_sub
    for ses = 1:n_ses
        M{sub,ses} = corr(time_series_data{sub,ses});
    end
end

%% Params
n_rep = 100;
n_null = 100;
n_rew = 1;
q = nan(n_sub, n_ses);
q_null_mean = nan(n_sub, n_ses);

%% Static modularity calculation
for sub = 1:n_sub
    sub
    for ses = 1:n_ses
            Qb = 0;
            for rep = 1:n_rep
                A = M{sub, ses}(:, :);
                A = A .* (A > 0); %set all negative fc to  0
                [~, Qt] = community_louvain(A, 1, []);
                if Qt > Qb
                    Qb = Qt;
                end
                q(sub, ses) = Qb;
            end
     end
end
save(fullfile(paths.out,'mean_static_modularity.mat'), 'q');

%% Randomized modularity calculation
workerpool = parpool(n_workers); %LONG
parfor sub = 1:n_sub
    q_null = zeros(1, n_null);
    sub
    for ses = 1:n_ses
            A = M{sub, ses}(:, :);
            A = A .* (A > 0);
            if sum(isnan(A(:))) > 0
                q_null_mean{sub, ses} = 0;
            else
                q_null = zeros(1, n_null);
                for null = 1:n_null
                    B = randmio_und(A, n_rew);
                    Qb = 0;
                    for rep = 1:n_rep
                        [~, Qt] = community_louvain(B, 1, []);
                        if Qt > Qb
                            Qb = Qt;
                        end
                    end
                    q_null(null) = Qb;
                end
                q_null_mean(sub, ses) = mean(q_null);
            end
     end
end
delete(workerpool)
save(fullfile(paths.out,'mean_null_modularity.mat'), 'q_null_mean');

%% Normalized modularity calculation
q_norm = nan(n_sub, n_ses);
for sub = 1:n_sub
    for ses = 1:n_ses
            q_norm(sub, ses) = q(sub, ses) / q_null_mean(sub, ses);
    end
end
save(fullfile(paths.out,'mean_normalized_modularity.mat'), 'q_norm');
