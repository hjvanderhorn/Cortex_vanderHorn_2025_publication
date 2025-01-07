
% Script to perform stats on connectivity measures and to create plots

clc
clear

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');

%% load the data you need
load('DataLAPTOP1_v3_coupled.mat', 'MergedStateMeasuresLAPTOP1')
load('DataLAPTOP1_v3_coupled.mat', 'MergedStaticFCMeasuresLAPTOP1')
load('DataLAPTOP2020_coupled.mat', 'LAPTOP2020_TableTotalStateMeasuresNew_after_motion_excl',...
    'LAPTOP2020_TableTotalStaticFC_New_after_motion_excl')

%for plotting errors in shades:
addpath(genpath('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\scripts\raacampbell-shadedErrorBar-19cf3fe'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Analyses dynamic FC (F9) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine LAPTOP1 (v1/2 merged with v3) with LAPTOP2020

TotalDataStateMeasures = [MergedStateMeasuresLAPTOP1 ; ... 
    LAPTOP2020_TableTotalStateMeasuresNew_after_motion_excl];

% since we may need other covariates, we want to put clin data in same
% order (now for the entire dataset). We need to do this again since we are
% now merging with data from previous analyses:
URSIandVisitREST = table2array(TotalDataStateMeasures(:,1:2));
ClinicalData = readtable('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\scripts\LAPTOP1_2020_combined_clean_concise_age_fine_long.xlsx');
URSIandVisitClinicalData = table2array(ClinicalData(:,[2,3]));

% NB: I now reversed the order and match Clinical Data to REST data (does
% not really matter which order, but it feels better to keep all merged
% REST data as is)
[LocInRESTData, LocInClindata] = ismember(URSIandVisitREST, URSIandVisitClinicalData, 'rows'); 
ClinicalDataNew = ClinicalData(LocInClindata,:);

%% sanity checks:
isequal(sum(LocInRESTData), size(URSIandVisitREST,1)) % this should be equal to size first dimension of data
isequal(table2array(ClinicalDataNew(:,2:3)), URSIandVisitREST)

% little extra sanity check, unique URSI with visit?
UniqueRowsURSIVisitStateMeasures = unique(URSIandVisitREST,'rows');
isequal(size(UniqueRowsURSIVisitStateMeasures,1), size(URSIandVisitREST,1))


%% adjustments to table to make it ready for glme

% add scanner (i.e., TRIOTIM vs PRISMA) as variable
TotalDataStateMeasures.Scanner = [zeros(size(MergedStateMeasuresLAPTOP1,1),1); ...
    ones(size(LAPTOP2020_TableTotalStateMeasuresNew_after_motion_excl,1),1)];

% % write for merging purposes in SPSS
% T = TotalDataStateMeasures(:,[1:4, 13, 126]); % URSI, visit, DX, meanFD,
% F9, Scanner.
% writetable(T,'F9.xlsx');
% clear T;

% specify variables as categorical
TotalDataStateMeasures.Scanner = categorical(TotalDataStateMeasures.Scanner);
TotalDataStateMeasures.DX = categorical(TotalDataStateMeasures.DX);
TotalDataStateMeasures.Visit = categorical(TotalDataStateMeasures.Visit);

% add age and sex:
TotalDataStateMeasures = [TotalDataStateMeasures  ClinicalDataNew(:,[6,7,8])]; % AgeEnrolled, Sex, Age_fine
TotalDataStateMeasures.Sex= categorical(TotalDataStateMeasures.Sex); % define sex as categorical

% add a variable for if we want to compare previous paper data (laptop1 visit 1/2) with added
% data (i.e., visit 3 laptop1 and laptop2020 all visits)
TotalDataStateMeasures.OldvsNewData = [ones(221,1);zeros(681,1); ...
    ones(size(LAPTOP2020_TableTotalStateMeasuresNew_after_motion_excl,1),1)];
TotalDataStateMeasures.OldvsNewData = categorical(TotalDataStateMeasures.OldvsNewData);


% Now let's sort the table based on URSI and visit
TotalDataStateMeasures = sortrows(TotalDataStateMeasures,[1,2]);

% Now exclude a few subjects that were included, but had braces at V1, and
% some subjects that had mental problems or new mTBI ocurring during study (after v1):
% (17385,91936 braces at V1; 54210, 55185, 57605: new psych/neuro between
% v1 and v2). For 30252 V2 had to be removed since it was accidentily
% processed using V1 data (so V1 can stay); we remove it at this stage
% since it was already included in the analyses for the previous paper (it
% is laptop1).

URSI_visit_tmp = [table2array(TotalDataStateMeasures(:,1)) ...
    double(table2array(TotalDataStateMeasures(:,2)))];
AdditionalExclIndex = [find(URSI_visit_tmp(:,1)==17385); ...
    find(URSI_visit_tmp (:,1)==54210 & URSI_visit_tmp (:,2)==3);
    find(URSI_visit_tmp (:,1)==55185 & URSI_visit_tmp (:,2)==3);
    find(URSI_visit_tmp (:,1)==57605 & URSI_visit_tmp (:,2)==2);
    find(URSI_visit_tmp(:,1)==91936);
    find(URSI_visit_tmp (:,1)==30252 & URSI_visit_tmp (:,2)==2)];

TotalDataStateMeasures(AdditionalExclIndex,:) = [];

%% run glme

model_F9_age_sex = fitglme(TotalDataStateMeasures, ...
    'F9 ~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

% to get the type III effects for the entire term (F-test)
AnovaModelF9 = anova(model_F9_age_sex)

% % now with age fine
% model_F9_agefine_sex = fitglme(TotalDataStateMeasures, ...
%     'F9 ~ 1 + DX * Visit + Scanner + Age_fine + Sex + meanFD + (1|URSI)', 'Distribution', 'gamma', ...
%     'Link', 'log', 'DummyVarCoding', 'effects', 'CovariancePattern', 'CompSymm')
% 
% % to get the type III effects for the entire term (F-test)
% anova(model_F9_agefine_sex)

% some varaibles that are needed for below
visit = double(table2array(TotalDataStateMeasures(:,2)));
group = grp2idx(table2array(TotalDataStateMeasures(:,3)));
age = double(table2array(TotalDataStateMeasures(:,127)));


% do some post-hoc tests for Group x Visit effect
DataDynV1 = TotalDataStateMeasures(visit==1,:);
DataDynV2 = TotalDataStateMeasures(visit==2,:);
DataDynV3 = TotalDataStateMeasures(visit==3,:);

model_dyn_v1 = fitglm(DataDynV1, ...
    'F9 ~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects')

model_dyn_v2 = fitglm(DataDynV2, ...
    'F9 ~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects')

model_dyn_v3 = fitglm(DataDynV3, ...
    'F9 ~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects')

%% Plots

F9 = table2array(TotalDataStateMeasures(:,13));
F9_orig = F9-1; % since I added 1 during computation to make it compatible with gamma
TotalDataStateMeasures.F9_orig = F9_orig;

make_figure_FC(F9_orig, group, visit, 'Fractional Occupancy State 9', ...
    'Fractional Occupancy', 'northeast', 'none');

% to make figure narrower I choose 'none' in above and save manually after
% resize
cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures');
print('Fractional Occupancy State 9', '-dtiff', '-r300');
cd ..

% 130th column is now F9 orig
make_figures_age(TotalDataStateMeasures(:,[1,2,131,127]),'F9_orig',...
    'Fractional Occupancy State 9','Visit', 'Fractional Occupancy',  '-dtiff');

cd ..

% make_figures_age_fine(TotalDataStateMeasures(:,[1,2,130,129]),'F9_orig',...
% 'F9','Visit', 'F9');

% if you want plot residuals with respect to all factors except age and
% then plot relationship with age:
model_F9_sex = fitglme(TotalDataStateMeasures, ...
    'F9 ~ 1 + DX * Visit + Scanner + Sex + meanFD + (1|URSI)', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky');
Residuals = residuals(model_F9_sex);
TotalDataStateMeasures.Residuals = Residuals;

make_figures_age(TotalDataStateMeasures(:,[1,2,132,127]),'Residuals',...
    'Fractional Occupancy State 9','Visit', 'Residuals Fractional Occupancy', '-dtiff');
cd ..

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Analyses static FC %%%%% %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%% Combine data
TotalDataStatic = [MergedStaticFCMeasuresLAPTOP1 ; ...
    LAPTOP2020_TableTotalStaticFC_New_after_motion_excl];

% sanity checks
RowsURSIVisitStaticFC = table2array(TotalDataStatic(:,1:2));
UniqueRowsURSIVisitStaticFC = unique(RowsURSIVisitStaticFC,'rows');
isequal(size(UniqueRowsURSIVisitStaticFC,1), size(RowsURSIVisitStaticFC,1))

%% Prepare table for stats
% add scanner variable
TotalDataStatic.Scanner = [zeros(size(MergedStaticFCMeasuresLAPTOP1,1),1); ...
    ones(size(LAPTOP2020_TableTotalStaticFC_New_after_motion_excl,1),1)];

% % write for merging purposes in spss
% T = TotalDataStatic(:,[1:4, 2105, 5450, 5895, 6675, 6676,6677 ]);
% writetable(T,'Static.xlsx');

% define variables as categorical
TotalDataStatic.DX = categorical(TotalDataStatic.DX);
TotalDataStatic.Visit = categorical(TotalDataStatic.Visit);
TotalDataStatic.Scanner = categorical(TotalDataStatic.Scanner);

% add age and sex:
TotalDataStatic = [TotalDataStatic ClinicalDataNew(:,[6,7,8])]; % order is same as for state measures data, so concat is allowed
TotalDataStatic.Sex= categorical(TotalDataStatic.Sex); % define sex as categorical

% add a variable for if we want to compare previous paper data (laptop1 visit 1/2) with added
% data (i.e., visit 3 laptop1 and laptop2020 all visits)
TotalDataStatic.OldvsNewData = [ones(221,1);zeros(681,1); ...
    ones(size(LAPTOP2020_TableTotalStateMeasuresNew_after_motion_excl,1),1)];
TotalDataStatic.OldvsNewData = categorical(TotalDataStatic.OldvsNewData);

% now sort based on URSI and visit
TotalDataStatic = sortrows(TotalDataStatic,[1,2]);

% Now exclude a few subjects that were included, but had braces at V1, and
% new mental problems of mTBI occurring during study (after v1):
TotalDataStatic(AdditionalExclIndex,:) = [];

% Little sanity check to see if tables for dFC and sFC are same:
isequal(TotalDataStateMeasures(:,[1:2]), TotalDataStatic(:,[1:2])) % ursi-visit part is equal

% extract names connections used for making formula lme; indexes (2101 etc)
% are obtained from previous data, and are the significant connections. We
% add 4 because connections start at column 5 in the table.
rSMA_CBL_10L_name = char(TotalDataStatic.Properties.VariableNames(2101+4))
L_R_Precuneus_name = char(TotalDataStatic.Properties.VariableNames(5446+4))
L_R_Thalamus_name = char(TotalDataStatic.Properties.VariableNames(5891+4))

%% run models
model_rSMA_CBL_10L = fitlme(TotalDataStatic, ...
    [rSMA_CBL_10L_name, '~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

model_L_R_Precuneus = fitlme(TotalDataStatic,...
    [L_R_Precuneus_name, '~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

model_L_R_Thalamus = fitlme(TotalDataStatic,...
    [L_R_Thalamus_name, '~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

% to get the type III effects for the entire term
AnovaModel_rSMA_CBL_10L = anova(model_rSMA_CBL_10L)
AnovaModel_L_R_Precuneus = anova(model_L_R_Precuneus)
AnovaModel_L_R_Thalamus = anova(model_L_R_Thalamus)

% post hoc for group x Visit for precuneus
% do some post-hoc tests for Group x Visit effect
DataStaticV1 = TotalDataStatic(visit==1,:);
DataStaticV2 = TotalDataStatic(visit==2,:);
DataStaticV3 = TotalDataStatic(visit==3,:);

model_precuneus_v1 = fitglm(DataStaticV1, ...
    [L_R_Precuneus_name, '~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD'], 'DummyVarCoding', 'effects')

model_precuneus_v2 = fitglm(DataStaticV2, ...
    [L_R_Precuneus_name, '~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD'], 'DummyVarCoding', 'effects')

model_precuneus_v3 = fitglm(DataStaticV3 , ...
    [L_R_Precuneus_name, '~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD'], 'DummyVarCoding', 'effects')

%% now plot

% NB: we have set outputformat to 'none', so we can manipulate image before
% printing (size, legends etc)

rSMA_CBL_10L = table2array(TotalDataStatic(:,2101+4));

subplot(3,1,1); 
make_figure_FC(rSMA_CBL_10L, group, visit,'Right SMA with Left Cerebellum 10',...
    'Fisher z score', 'northwest', 'none');

L_R_Precuneus = table2array(TotalDataStatic(:,5446+4));

subplot(3,1,2); 
make_figure_FC(L_R_Precuneus, group, visit,'Left with Right Precuneus',...
    'Fisher z score', 'northwest', 'none');

L_R_Thalamus = table2array(TotalDataStatic(:,5891+4));

subplot(3,1,3); 
make_figure_FC(L_R_Thalamus, group, visit,'Left with Right Thalamus',...
    'Fisher z score', 'northwest', 'none');

% now you can adjust size, remove second and third legend, make first
% legend horizontal and put above plots
% and then save using below.
cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures');
print('staticFC_narrow', '-dtiff', '-r300');
cd ..

% age effects

make_figures_age(TotalDataStatic(:,[1,2,2101+4,6676]),rSMA_CBL_10L_name,...
    'Right SMA with left cerebellum','Visit', 'Fisher z score', '-dtiff');
cd ..

make_figures_age(TotalDataStatic(:,[1,2,5891+4,6676]),L_R_Thalamus_name,...
    'Left with right thalamus','Visit', 'Fisher z score', '-dtiff');
cd ..

% if you want plot residuals with respect to all factors except age and
% then plot relationship with age:
model_rSMA_lCBL_sex = fitlme(TotalDataStatic, ...
    [rSMA_CBL_10L_name, '~ 1 + DX * Visit + Scanner + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky');

ResidualsRightSMAleftCBL = residuals(model_rSMA_lCBL_sex);
TotalDataStatic.ResidualsRightSMAleftCBL = ResidualsRightSMAleftCBL;

make_figures_age(TotalDataStatic(:,[1,2,6680,6676]),'ResidualsRightSMAleftCBL',...
    'Right SMA with Left Cerebellum 10','Visit', ...
    'Residuals Fisher z score', 'none'); % none so I can adjust the size before print

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures');
print('Age_Residuals_rSMA_lCBL', '-dtiff', '-r300');
cd ..

model_L_R_Thalamus_sex = fitlme(TotalDataStatic, ...
    [L_R_Thalamus_name, '~ 1 + DX * Visit + Scanner + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky');

ResidualsL_R_Thalamus = residuals(model_L_R_Thalamus_sex);
TotalDataStatic.ResidualsL_R_Thalamus = ResidualsL_R_Thalamus;

make_figures_age(TotalDataStatic(:,[1,2,6681,6676]),'ResidualsL_R_Thalamus',...
    'Left with Right Thalamus','Visit', ...
    'Residuals Fisher z score', 'none'); % none so I can adjust the size before print

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures');
print('Age_Residuals_l_r_thalamus', '-dtiff', '-r300');
cd ..


% %% Age fine with static 
% 
% % run models
% model_rSMA_CBL_10L_agefine = fitlme(TotalDataStatic, ...
%     [rSMA_CBL_10L_name, '~ 1 + DX * Visit + Scanner + Age_fine + Sex + meanFD + (1|URSI)'],...
%     'DummyVarCoding', 'effects', 'CovariancePattern', 'CompSymm')
% 
% model_L_R_Precuneus_agefine = fitlme(TotalDataStatic,...
%     [L_R_Precuneus_name, '~ 1 + DX * Visit + Scanner + Age_fine + Sex + meanFD + (1|URSI)'],...
%     'DummyVarCoding', 'effects', 'CovariancePattern', 'CompSymm')
% 
% model_L_R_Thalamus_agefine = fitlme(TotalDataStatic,...
%     [L_R_Thalamus_name, '~ 1 + DX * Visit + Scanner + Age_fine + Sex + meanFD + (1|URSI)'],...
%     'DummyVarCoding', 'effects', 'CovariancePattern', 'CompSymm')
% 
% % to get the type III effects for the entire term
% anova(model_rSMA_CBL_10L_agefine)
% anova(model_L_R_Precuneus_agefine)
% anova(model_L_R_Thalamus_agefine)


%% show sample sizes for groups and visits
% to compare with sample sizes for making flowchart in spss

% N HC
disp('sample sizes for v1, v2 and v3 for HC');
N_HC_v1 = sum(group==1 & visit ==1)
N_HC_v2 = sum(group==1 & visit ==2)
N_HC_v3 = sum(group==1 & visit ==3)

%NTBI
disp('sample sizes for v1, v2 and v3 for TBI');
N_mTBI_v1 = sum(group==2 & visit ==1)
N_mTBI_v2 = sum(group==2 & visit ==2)
N_mTBI_v3 = sum(group==2 & visit ==3)

%% save it all
save('StatsResults_01222023');


% Also make and write table for merging with spss clinical measures.
% So you can import this xlsx file into spss, save as spss file and use
% this for merger.
TableForMergerSPSS = [TotalDataStateMeasures(:,1:4) TotalDataStateMeasures(:,13) TotalDataStateMeasures(:,[126,129:130]) ...
    TotalDataStatic(:,2101+4) TotalDataStatic(:,5446+4) TotalDataStatic(:,5891+4) ];
 writetable(TableForMergerSPSS , 'TableForMergerSPSS.xlsx');

 %% To conduct analyses in new data only (so without the laptop1 visit 1/2 that we already published)
 
NewStateData = TotalDataStateMeasures;
NewStateData = NewStateData(double(Test.OldvsNewData) == 2,:);
model_F9_age_sex_new_data = fitglme(NewStateData, ...
'F9 ~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)', 'Distribution', 'gamma', ...
'Link', 'log', 'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')
anova(model_F9_age_sex_new_data)

NewStaticData = TotalDataStatic;
NewStaticData  = NewStaticData (double(Test.OldvsNewData) == 2,:);

model_rSMA_CBL_10L_new_data = fitlme(NewStaticData, ...
    [rSMA_CBL_10L_name, '~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

model_L_R_Precuneus_new_data = fitlme(NewStaticData,...
    [L_R_Precuneus_name, '~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

model_L_R_Thalamus_new_data = fitlme(NewStaticData,...
    [L_R_Thalamus_name, '~ 1 + DX * Visit + Scanner + AgeEnrolled + Sex + meanFD + (1|URSI)'],...
    'DummyVarCoding', 'effects', 'CovariancePattern', 'FullCholesky')

% to get the type III effects for the entire term
anova(model_rSMA_CBL_10L_new_data)
anova(model_L_R_Precuneus_new_data)
anova(model_L_R_Thalamus_new_data)