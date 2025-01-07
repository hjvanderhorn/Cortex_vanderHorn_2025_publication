%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PROJECT LAPTOP1 v3 and LAPTOP 2020 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script uses centroids from previous analyses and runs LEiDA on new
% data (using AAL1 parcellation) and then computes state measures. It does
% so by computing the cityblock distance between time points (first
% eigenvector per time point) and the previous centroids and labeling them
% accordingly based on smallest distance. In other word, each time point is
% assigned to a state (i.e., a previous centroid) based on smallest distance.
%
% Static FC is also computed between AAL regions (i.e., Fisher z transformed corr).
%
% A list of full paths for AAL csv files is needed. This can be made using
% the CreateList scripts.
% 
% We run separately for LAPTOP1 visit 3 and for LAPTOP2020 and then
% merge the data later (also with data from previous analyses of visit 1
% and 2 of LAPTOP1).
%

clc
clear

%% Variables to load and make
% so you run this for LAPTOP1 v3 and LAPTOP2020 separately

% Adjust the following by commenting/uncommenting:
%load('ListRESTfilesLAPTOP2020', 'List'); 
%OutputFileName = 'Data_LAPTOP2020.mat' ; 
load('ListRESTfilesLAPTOP1_v3', 'List'); 
OutputFileName = 'Data_LAPTOP1_v3.mat'; 

%% Add LEIDA, gift and spm to path
addpath('T:/research/analysis/human/amayer/shared/MAYER_ALL/andy/Hans/LEIDA/LEiDA-master/LEiDA-master/');
addpath(genpath('T:/research/analysis/human/amayer/shared/MAYER_ALL/andy/Hans/gift/gift-master/GroupICATv4.0c/'));
addpath('T:/research/analysis/human/amayer/shared/apps/matlab/spm12_Hans');
addpath('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\scripts');

%% Extract the centroids for k is 10 from previous analyses to be used as starting point
cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');
load('LEIDA_060523_old_laptop_data.mat', 'CTot'); % centroids from previous analyses
PreviousCentroids = CTot{10}; % we need the ones from k is 10.

%% Load limits meanFD for exclusion
% from previous analyses of LAPTOP1 v1 and v2
load('LEIDA_060523_old_laptop_data.mat', 'Limit_meanFD_HC', 'Limit_meanFD_mTBI');

%% Extract dynamic and static FC data

% preallocate for some speed:
TotalMemberships = zeros(647, size(List,1)); % for F, MDT and NT; 647 (Tmax) is number of time points after kicking out first and last time point post-Hilbert
TotalStaticFC = zeros(size(List,1), (116*115/2)); % for static fc between 116 AAL regions; number of connections is 116*115/2

% loop to extract membership indices (i.e., to what state observations/rows
% belong) for each subject in List (a character array with paths to csv aal
% tc files), and to compute all state measures we need
for i = 1:size(List,1)
    disp(['Working on file ', num2str(i), ' out of ', num2str(size(List,1))])
    
    % this function computes everything we need:
    [SubjectMembershipIndex, V1Data, F, MDT, NT, EasyForStatsProbVector, StaticFC, ...
    Mean_iFC_state9_vector] = ExtractDataLEIDAandStatic(List(i,:), PreviousCentroids); 
    
    % save k means membership index for subject in total Tmax * Subj matrix
    TotalMemberships(:,i) = SubjectMembershipIndex;
    
    % save V1 data for subject in total Subj*Tmax*AAL matrix
    TotV1Data(i,:,:) = V1Data;
    
    % save mean phase coherences for state 9 in total matrix, in case we
    % want to look at connectivity (i.e., phase coherence) between AAL
    % regions within this state specifically
    TotMean_iFC_state9(i,:) = Mean_iFC_state9_vector;
    
    % save state measures in Subj*Measures matrices
    fracttime(i,:) = F + 1;
    dwelltime(i,:) = MDT + 1;
    NumTransitions(i,:) = NT + 1;
    Probabilities(i,:) = EasyForStatsProbVector + 1;
    % NB: we add 1 to make it compatible with gamma. Make sure you subtract
    % this 1 when you are plotting for example boxplots.
    
    TotalStaticFC(i,:) = StaticFC; % add static FC for subject to total matrix
    
    clear SubjectMembershipIndex V1Data F MDT NT EasyForStatsProbVector StaticFC Mean_iFC_state9_vector
    
end

% combine state measures data into one big Subject*Measures matrix, so we
% can convert this into a table
TotalStateMeasures = [fracttime dwelltime NumTransitions Probabilities];

%% Create a state measures table as well as a static FC table

% This requires getting indices for columns first (10 refers to k = 10)
IndicesWorkDataF = 1:10;
IndicesWorkDataMDT = 10+1:10*2;
IndicesWorkDataNT = 10*2+1;
IndicesWorkDataProbs = 10*2+2:10*2+1+10*10;

% And now we can make the variable names
for tt = 1:size(IndicesWorkDataF,2)
    StateMeasuresVariableNames{IndicesWorkDataF(tt)} = ['F', num2str(tt)];
end

for tt = 1:size(IndicesWorkDataMDT,2)
    StateMeasuresVariableNames{IndicesWorkDataMDT(tt)} = ['MDT', num2str(tt)];
end

for tt = 1:size(IndicesWorkDataNT,2)
    StateMeasuresVariableNames{IndicesWorkDataNT(tt)} = 'NT';
end

for tt = 1:size(IndicesWorkDataProbs,2)
    StateMeasuresVariableNames{IndicesWorkDataProbs(tt)} = ['P', num2str(tt)];
end

% Create variable names for static FC connections
AALlabels = readtable('T:/research/analysis/human/amayer/shared/MAYER_ALL/andy/Hans/LEIDA/AAL/aal.nii.txt');
ROI_Names = AALlabels(:,2) ; % ROI names
ROI_Names = table2cell(ROI_Names); % put it in cell array, easier to use later.

ComponentPairs = nchoosek([1:116],2); % get the possible combinations for the connections (these thus match one triangular part of the conn matrix)

% now loop so that you construct the connection names
for i=1:size(ComponentPairs,1) % put it all together in strings
    ComponentPairsNamed{:,i} = [ROI_Names{ComponentPairs(i,1)}, 'with', ...
        ROI_Names{ComponentPairs(i,2)}];
end

% make table state measures
TableTotalStateMeasures = array2table(TotalStateMeasures, 'VariableNames',...
    StateMeasuresVariableNames);

% make table static FC
TableTotalStaticFC = array2table(TotalStaticFC, 'VariableNames', ... 
    ComponentPairsNamed);
 
save(OutputFileName);  % save it all.

% - for LAPTOP1 visit 3 now run: CoupleLAPTOP1withClinical.m
% - for LAPTOP 2020 now run CoupleLAPTOP2020withClinical.m 
% - afer you ran both, run StatsTotal.m to merge LAPTOP1 and 20202 and perform stats



