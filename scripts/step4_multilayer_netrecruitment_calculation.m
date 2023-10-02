%% Network integration/segregation

%Extracts recruitment and integration values using a
%permutation approach per network.

%requires the helper functions
%normalize_networks_mean

%Reproducing the following since code not available.
%Daws, R.E., Timmermann, C., Giribaldi, B. et al. Increased global integration in the brain after psilocybin therapy for depression. Nat Med 28, 844–851 (2022).
%https://doi.org/10.1038/s41591-022-01744-z


%% Init
close all
clear all

%% Path
paths.home = (cd);
addpath(genpath(paths.home))

paths.out = fullfile(paths.home,'results','dynamic_ica_flex');
if ~exist(paths.out)
    mkdir(paths.out);
end

%% Load data and labels
load(fullfile(paths.out,'dynamic_multilayer_allegiance.mat'));

%No variable network roi labels per scan (fixed here since we use single parc.
netfile = xlsread(fullfile(paths.home,'mni_atlas','parc_networks','schaefer18networks_idx.xlsx'))
yeoROI = netfile(:,1);
yeoID = netfile(:,2);
roi_idx = yeoID;
%but this can be fixed just ask

%% Params
n_sub = size(opt_mean_allegiance_mat, 1);
n_ses = size(opt_mean_allegiance_mat, 2);
n_nets = unique(roi_idx);
n_perms = 1000;
n_workers = 8;

%% Begin extraction
for sub = 1:n_sub
    for ses = 1:n_ses
        %get allegiance matrices averages across window x optimisations
        in_dat = opt_mean_allegiance_mat{sub,ses};
        %get roi labels per subject (variable since ica)
        in_labels = roi_idx;
        %extract mean across perms and windows
        norm_mean_allegiance{sub,ses} =  normalize_networks_mean(in_dat,in_labels,n_perms); %across all windows

    end
end
%% Now to check stats
save(fullfile(paths.out,"dynamic_netrecruitment.mat"),"norm_mean_allegiance");
