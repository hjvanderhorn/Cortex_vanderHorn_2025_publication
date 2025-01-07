
% script to make a list with paths to aal tc files for visit 3 laptop1
% this list can be used to perform leida using ProjectVisit3.m

clear
clc

addpath('T:/research/analysis/human/amayer/shared/apps/matlab/spm12_Hans');

cd('T:/research/analysis/human/amayer/laptop_20019/data/derivatives/'); % rootfolder

% this gives same list:
folderList = dir(fullfile(pwd, 'sub-*'));

% remove everything that is not a folder
for i=1:size(folderList,1)
IsFolder(i,1) = isfolder(folderList(i).name);
end

folderList = folderList(IsFolder);
% now remove the last entry, since this is not a subject:
folderList(size(folderList,1)) = [];


disp('Creating list with tc files...');
List = []; % make empty List to which you add all paths to csv files

% let's make the list
for i = 1:size(folderList,1)
    disp(['Working on folder', num2str(i), ' of ', num2str(size(folderList,1))]);
    cd(folderList(i).name);
    
    % search for aal tc file within the REST folder
    tmp = spm_select('FPListRec', pwd, '^.*.aal_tc.csv');
    
    if ~isempty(tmp) % there could be no rest at all
        
    List = [List; tmp]; % add filepath to list
    
    end
    
    cd .. % go back to rootfolder data
    clear tmp;
    
end

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');
save('ListRESTfilesLAPTOP2020', 'List');