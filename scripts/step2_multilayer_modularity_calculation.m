%% Multilayer community detection algoritm

%Again, extends to variable roi x rois if need be.
%Reproducing since code not available.
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

%% Data
load(fullfile(paths.home,'results','dynamic_ica_flex','sliding_window_fc.mat'));
M = dynamic_correlation_matrices;
n_sub = size(M, 1);
n_ses = size(M, 2);
n_win = size(M{1}, 3);

%% Define operation parameters
gamma = 1;
omega = 1;
n_rep = 100;
n_workers = 10; %for parfor

%% Initialize cells to store data
modularity_mean = cell(n_sub, n_ses);
modules = cell(n_sub, n_ses);

workerpool = parpool(n_workers); %VERY SLOW
parfor sub = 1:n_sub
 %for sub = 1:n_sub
    for ses = 1:n_ses
        %Prepare objects
        n_roi = size(M{sub, ses, 1}, 1); % Determine the number of ROIs for this scan
        A = cell(1, n_win);
        B = spalloc(n_roi * n_win, n_roi * n_win, (n_roi + n_win) * n_roi * n_win);
        twomu = 0;

        %Generate null models
        for win = 1:n_win
            %Copy network with positive weights thresholding
            M_win= M{sub,ses};
            A{win} =  M_win(:,:,win) .* (M_win(:,:,win) > 0);
            k = sum(A{win}); % node degree
            twom = sum(k); % mean network degree
            twomu = twomu + twom; % increment
            indx = (1:n_roi) + (win - 1) * n_roi; % find indices
            B(indx, indx) = A{1,win} - gamma * (k' * k) / twom; % fill B matrix
        end
        twomu = twomu + 2 * omega * n_roi * (n_win - 1);

        B = B + omega / 2 * spdiags(ones(n_roi * n_win, 2), [-n_roi, n_roi], n_roi * n_win, n_roi * n_win);
        B = B + omega * spdiags(ones(n_roi * n_win, 2), [-2 * n_roi, 2 * n_roi], n_roi * n_win, n_roi * n_win);

        %Calculate multilayer modules
        for rep = 1:n_rep
            clc;
            fprintf('Subject = %i, Session = %i N_roi = %i ',  sub, ses, n_roi );
            [S, Q] = genlouvain(B);
            Q = Q / twomu;

            S = reshape(S, n_roi, n_win);

            modularity_mean{sub, ses}(:, rep) = Q;
            modules{sub, ses}(:,:,rep) = S;
        end
    end
end
delete(workerpool)

%% Save
save(fullfile(paths.out,'dynamic_multilayer_modularity.mat'), 'modularity_mean','modules');
