
% script to make a list with paths to aal tc files for visit 3 laptop1
% this list can be used to perform leida using ProjectVisit3.m

clear
clc

addpath('T:/research/analysis/human/amayer/shared/apps/matlab/spm12_Hans');

cd('T:/research/analysis/human/amayer/laptop_20753/data'); % rootfolder

% this gives same list:
folderList = dir(fullfile(pwd, 'M871*'));
folderList(440)=[]; % remove the one that is not a subjectfolder we want to use (with ~Template in name)

disp('Creating list with tc files...');
List = []; % make empty List to which you add all paths to csv files

% let's make the list
for i = 1:size(folderList,1)
    disp(['Working on folder', num2str(i), ' of ', num2str(size(folderList,1))]);
    cd(folderList(i).name);
    
    % check if there is a visit 3 folder
    if isfolder('visit3')
        cd('visit3');
    else 
        cd .. % go back up one step to get back in the rootfolder
        continue % if not present, continue with next iteration
    end
    
    % check if there is a REST folder
    if isfolder('REST1')
        cd('REST1');
    else % if not present continue with next iteration
        cd ../.. % go back up two steps to get back in the rootfolder
        continue
    end
    
    % search for aal tc file within the REST folder
    tmp = spm_select('FPListRec', pwd, '^.*.aal_tc.csv');
    
    if ~isempty(tmp) % extra safety step, there might be no csv file in the folder for whatever reason
        
    List = [List; tmp]; % add filepath to list
    
    end
    
    cd ../../.. % go back to rootfolder data
    clear tmp;
    
end

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');
save('ListRESTfilesLAPTOP1_v3', 'List');

