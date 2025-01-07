
% script to couple rest visit 3 data for LAPTOP1 with clinical data, and
% merge this with previous LAPTOP1 data (visit 1 and 2).

%% load data
clc 
clear

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');
load('Data_LAPTOP1_v3.mat', 'List')
load('Data_LAPTOP1_v3.mat', 'TableTotalStateMeasures')
load('Data_LAPTOP1_v3.mat', 'TableTotalStaticFC')

load('LEIDA_060523_old_laptop_data.mat', 'Limit_meanFD_HC')
load('LEIDA_060523_old_laptop_data.mat', 'Limit_meanFD_mTBI')

%% Extract URSI and visit from the filelist of visit 3 data
% we need this info to put clinical and REST data is same order

disp('Extract URSI and visit from filelist...')

URSIandVisitREST = zeros(size(List,1),2);

% let's fill this matrix 
for i=1:size(List,1)
    disp(['Working on file', num2str(i), ' of ', num2str(size(List,1))]);
    
    URSIandVisitREST(i,1)=str2num(char(extractBetween(List(i,:), 'M871', [filesep 'visit']))); %NB: on Linux the slash needs to be forward
    URSIandVisitREST(i,2)=str2num(char(extractBetween(List(i,:), 'visit', [filesep 'REST']))); %NB: on Linux the slash needs to be forward
    
end

%% Load clinical data  

ClinicalData = readtable('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\scripts\LAPTOP1_2020_combined_clean_concise_age_fine_long.xlsx');
URSIandVisitClinicalData = table2array(ClinicalData(:,[2,3])); 

%% put clinical and rest data in same order 
[LocInClinData, LocInRESTdata] = ismember(URSIandVisitClinicalData,URSIandVisitREST, 'rows'); % where are subjects from SPSS in total filelist rest
% LocInClinData is logical indicating if clinical data is in REST data
% LocInRESTdata is index for where clinical data is in REST data (with
% zeros indicating that it is not in REST data at all)

LocInRESTdata(LocInRESTdata==0)=[]; % remove subjects (zeros) in vector that are apparently not in REST data
% this results in a vector with indices in REST data that correspond with the
% logical ones in LocInClinData (alternative with same solution:
% LocInRESTdata = LocInRESTdata(LocInClinData) );

% make selection of clinical data that are also in REST data
ClinicalDataNew = ClinicalData(LocInClinData,:); % also make new list with variables in clin data that are also in REST
URSIandVisitClinicalDataNew  = URSIandVisitClinicalData(LocInClinData,:); % and new URSI variable

TableTotalStateMeasuresNew = TableTotalStateMeasures(LocInRESTdata,:); % make new table with subjects for which clinical data is available, and in that order
URSIandVisitREST_New = URSIandVisitREST(LocInRESTdata,:); % new URSI and visit variable as well
TableTotalStaticFC_New = TableTotalStaticFC(LocInRESTdata,:); % and make new static fc table as well

% is data really in same order?
CheckOrder1 = isequal(URSIandVisitREST_New, URSIandVisitClinicalDataNew)
CheckOrder2 = isequal(URSIandVisitREST_New, table2array(ClinicalDataNew(:,[2,3])))

%% Make table of clinical data

TableClinicalMeasures = array2table(table2array(ClinicalDataNew(:,[2:5])), 'VariableNames',...
    {'URSI', 'Visit', 'DX', 'meanFD'});
    
%% Identify subjects to exclude based on mean FD 
% limits are from previous LAPTOP1 visit 1/2 study

group = table2array(TableClinicalMeasures(:,3)); % group variable
meanFD = table2array(TableClinicalMeasures(:,4)); % meanFD
ExcludeMotionHC = find(meanFD>Limit_meanFD_HC & group == 0); % find subjects that are above limit for HC
ExcludeMotionTBI = find(meanFD>Limit_meanFD_mTBI & group == 1); % and for mTBI
ExcludeTotal = [ExcludeMotionHC; ExcludeMotionTBI]; % combine these subjects into Exclude variable

%% Combine FC with clinical data and exclude motion outliers

% combine and prename final state measures table
TableTotalStateMeasuresNew_after_motion_excl = [TableClinicalMeasures...
    TableTotalStateMeasuresNew];

TableTotalStateMeasuresNew_after_motion_excl(ExcludeTotal,:) = []; % exclude motion outliers

% combine and prename final static fc table
TableTotalStaticFC_New_after_motion_excl = [TableClinicalMeasures...
    TableTotalStaticFC_New];

TableTotalStaticFC_New_after_motion_excl(ExcludeTotal,:) = []; % exclude motion outliers
    

%% Now merge visit 3 data with visit 1 and 2 data

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\AnalysesResults_LAPTOP_PreprocessedOldWay\doublecheck_110823');

% merge dFC data
Visit_1_2_data_dFC = readtable('Table_for_k_is_010.xlsx');
MergedStateMeasuresLAPTOP1 = [TableTotalStateMeasuresNew_after_motion_excl; Visit_1_2_data_dFC];

% merge static FC data
Visit_1_2_data_staticFC = readtable('TableStaticConnAAL_old_laptop_data_doublecheck_081123.csv');
MergedStaticFCMeasuresLAPTOP1 = [TableTotalStaticFC_New_after_motion_excl; Visit_1_2_data_staticFC];

%% save it all.
cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');
save('DataLAPTOP1_v3_coupled'); 

% now you can merge with LAPTOP2020 data (which happens in StatsTotal.m).
