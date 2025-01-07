%%% code to make inputs for BrainNet Viewer (Xia et al, PlosOne, 2013).
% this script is an adaptation from the original MakeInputBrainNetViewer.m
% script, and now we only make a positive community for state 9

clc
clear

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\results');
load('LEIDA_060523_old_laptop_data.mat', 'CTot')
States = CTot{10};


test = States(9,:);

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\scripts\Analyses_LAPTOP_PreprocessedOldWay\BrainNetViewer_20191031\Data\ExampleFiles\AAL90\test_states')
T = readtable('Node_AAL116.txt');


index_pos = find(test>0) ;% positive element indices in centroid; 
T2 = T;
T2{index_pos,4} = 5; %5 will become red in jet in BrainNet if you set threshold in BrainNet low for blue and high for red 
index_neg = 1:116;
index_neg(index_pos)=[]; % now you have a vector with the negative element indices

% give nodes size based on absolute value in Vc
for i=1:116
    T2{i,5} = abs(test(i));
end

T2(index_neg,:)= []; % remove negative community

cd('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\figures\BrainNetFigures');
writetable(T2, 'NodesState9_only_positive.txt', 'Delimiter', '\t', 'WriteVariableNames', 0); 

% edges; simply binarize all edges between pos elements
M = ones(1,size(T2,1));
M2 = M.*M'; % outer product now results in edges of 1 between pos values, and 0 between the rest
M2(1:size(M2,1)+1:end)=0; % !!! added this after making state 2 fig, don't think it matters, but maybe test in the end
writematrix(M2, 'EdgeBinaryState9_only_positive.txt', 'Delimiter','tab'); % adjust number

clear T2 M M2 test index_pos index_neg


% now use these in BrainNet viewer