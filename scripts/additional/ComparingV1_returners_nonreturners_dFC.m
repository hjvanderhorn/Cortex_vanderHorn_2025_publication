clc
clear

load('StatsResults_01222023', 'TotalDataStateMeasures');

% variables you need to divide/extract data etc
URSI_base = TotalDataStateMeasures.URSI;
Visit = double(TotalDataStateMeasures.Visit); % need to get back to double since you defined as cat

% URSI's per visit
URSI_visit1 = URSI_base(Visit == 1);
URSI_visit2 = URSI_base(Visit == 2);
URSI_visit3 = URSI_base(Visit == 3);

% find URSI for subjects who have V1 and V2
Present_Visit_1_and_2 = ismember(URSI_visit1, URSI_visit2);
URSI_present_Visit_1_and_2 = URSI_visit1(Present_Visit_1_and_2);

% find URSI for subjects who have V1 and V2 who also have V3
Present_Visit_1_2_and_3 = ismember(URSI_present_Visit_1_and_2, URSI_visit3);
URSI_present_visit_1_2_and_3 = URSI_present_Visit_1_and_2(Present_Visit_1_2_and_3);

% now build an index variable for subjects that have all visits
Index_all_visits = [];
for i = 1:size(URSI_present_visit_1_2_and_3,1)
Index_all_visits = [Index_all_visits; find(URSI_base == URSI_present_visit_1_2_and_3(i))];
end

% get the data for subject who have all visits
TotalDataStateMeasures_all_visits = TotalDataStateMeasures(Index_all_visits,:);
IndexV1_for_all_returners = double(TotalDataStateMeasures_all_visits.Visit) == 1;
TotalDataStateMeasures_all_visits_V1 = TotalDataStateMeasures_all_visits(IndexV1_for_all_returners,:);

% get data for subjects that do not have all visits (ergo, subjects with
% either only V1, or V1 and V2, or V1 and V3, or V2 and V3, or only V2 or
% only V3:
TotalDataStateMeasures_not_all_visits = TotalDataStateMeasures;
TotalDataStateMeasures_not_all_visits(Index_all_visits,:)=[];

% now obtain the duplicates here, which means you identify subjects who
% have V1 and V2, or V1 and V3, or V2 and V3.
URSI_not_all_visits = TotalDataStateMeasures_not_all_visits.URSI;
[~, UniqueIndex] = unique(URSI_not_all_visits, 'stable');
Repeat_index = setdiff(1:numel(URSI_not_all_visits), UniqueIndex); % Note that we know that subjects have at maximum 2 visits, otherwise they would have all visits

% extract data for subjects that only have V1
URSI_only_one_visit = URSI_not_all_visits;
URSI_only_one_visit([Repeat_index Repeat_index-1],:)=[]; % remove the ones with two visits (we do minus 1 to also remove the one that the duplicate belongs too, the first)
[~ , PresentOnlyOneVisit] = ismember(URSI_only_one_visit, URSI_base); % find only the ones that have one visit (this is mostly V1, but can also be a V2 or V3)
TotalDataStateMeasures_onlyV1 = TotalDataStateMeasures(PresentOnlyOneVisit,:); % get this data, and call it V1, but below we will remove the only V2s and only V3s
Only_Visit2 = double(TotalDataStateMeasures_onlyV1.Visit) == 2; % now identify the ones that have V2 as single visit
TotalDataStateMeasures_onlyV1(Only_Visit2,:) = []; % and remove
Only_Visit3 = double(TotalDataStateMeasures_onlyV1.Visit) == 3; % do same for single V3
TotalDataStateMeasures_onlyV1(Only_Visit3,:) = [];

% all right, let's do some stats. You now only need glm since you got rid
% of visit term.

model_v1_non_returners = fitglm(TotalDataStateMeasures_onlyV1, ...
    'F9 ~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects')

model_v1_returners = fitglm(TotalDataStateMeasures_all_visits_V1, ...
    'F9 ~ 1 + DX + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects')


% alternatively, after discussion with Andy, we simply compare R vs NR,
% irrespective of DX

TotalDataStateMeasures_all_visits_V1.RvsNR = zeros(size(TotalDataStateMeasures_all_visits_V1,1),1);
TotalDataStateMeasures_onlyV1.RvsNR = ones(size(TotalDataStateMeasures_onlyV1,1),1);

TotalTable = [TotalDataStateMeasures_all_visits_V1; TotalDataStateMeasures_onlyV1];
TotalTable.RvsNR = categorical(TotalTable.RvsNR);

model_returners_vs_nonreturners = fitglm(TotalTable, ...
    'F9 ~ 1 + RvsNR + DX + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects'); % Note that I added DX too, since for NR the N for HC is way less

% to perform these analyses per group separately
group_return_vs_nonreturn = double(TotalTable.DX);
TotalTableHC = TotalTable(group_return_vs_nonreturn==1,:);
TotalTableTBI = TotalTable(group_return_vs_nonreturn==2,:);

model_returners_vs_nonreturners_HC = fitglm(TotalTableHC, ...
    'F9 ~ 1 + RvsNR + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects') 

model_returners_vs_nonreturners_TBI = fitglm(TotalTableTBI, ...
    'F9 ~ 1 + RvsNR + Scanner + AgeEnrolled + Sex + meanFD', 'Distribution', 'gamma', ...
    'Link', 'log', 'DummyVarCoding', 'effects') 