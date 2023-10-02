%% Allegiance matrix generation from multilayer (windowed) modularity

%Calculates changes in network allegiance across windows and optimisaton
%runs

%from this network integration can be devised.

%requires the helper functions
% allegiance_matrix
% allegiance_matrix_opti
% all_window_allegiance_mean
% single_window_allegiance
% single_window_allegiance_mean


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

%% Load Louvain data 
load(fullfile(paths.out,'dynamic_multilayer_modularity.mat'));

%% Params
n_sub = size(modules, 1);
n_ses = size(modules, 2);
n_win = size(modules{1}, 2);
n_opt = size(modules{1}, 3);
n_workers = 10;


%% Main extraction
workerpool = parpool(n_workers); %LONG
parfor sub = 1:n_sub
    for ses = 1:n_ses
        %restructure to opt x roi x window
        in_dat = permute(modules{sub,ses},[3, 1, 2]);
        %extract mean across perms and windows
        opt_mean_allegiance_mat{sub,ses} = allegiance_matrix_opti(in_dat); %across all windows

        %extract mean across perms per window
        opt_meanwin_allegiance_mat{sub,ses} = all_window_allegiance_mean(in_dat);
    end
end
delete(workerpool)

%% Save
save(fullfile(paths.out,'dynamic_multilayer_allegiance.mat'),'opt_mean_allegiance_mat','opt_meanwin_allegiance_mat');
