    
function [SubjectMembershipIndex, V1Data, F, MDT, NT, EasyForStatsProbVector, StaticFC, ...
    Mean_iFC_state9_vector] = ExtractDataLEIDAandStatic(InputData, PreviousCentroids)

N = 116; % number of AAL regions
Tmax = 647; % number of time points after cutting off first and last
Phases=zeros(N,(Tmax+2)); % I put +2 here, because first it needs all time points (649), later after Hilbert I remove first and last (Tmax = 647)

%% load tc data
FileTCs = readtable(InputData); % read data as table
signaldata = table2array(FileTCs(:,3:end)); % extract time courses and put in array

%% Get the Phase of BOLD data
for seed=1:N
    timeseriedata=demean(detrend(signaldata(:,seed)));
    Phases(seed,:) = angle(hilbert(timeseriedata));
end

Phases = Phases(:,2:648); % kick out first and last time point to account for edge effects (and that makes it 647)

%% Get the iFC leading eigenvector at each time point

V1Data = zeros(Tmax, N); % empty matrix to save all V1's
Tot_iFC = zeros(Tmax, N, N); % empty matrix to save all iFCs

for t=1:Tmax % NB: Solonso-Martinez just selects everything except for first and last here, multiple ways to Rome
    
    iFC=zeros(N);
    
    for n=1:N
        for p=1:N
            iFC(n,p)=cos(adif(Phases(n,t),Phases(p,t))); % adif not really needed, similar results as below
            %iFC(n,p)=cos(Phases(n,t)-Phases(p,t)); % but then we have
            %to do minus instead of , (because that's an argument
            %within adif)
        end
    end
    
    % make total iFC matrix
    Tot_iFC(t,:,:) = iFC;
    
    % old code Cabral, which I used so far
    [eVec,eigVal]=eig(iFC);
    eVal=diag(eigVal);
    [val1, i_vec_1] = max(eVal);
    V1 = eVec(:,i_vec_1);
    
    %         % which can be replaced by, let's see if it makes a difference
    %         [V1,~] = eigs(iFC,1);
    %         % I did read some posts that eigs is tricky, cause it uses random
    %         % starting points, so it may produce different results from run to
    %         % run. So maybe best to stick to eig. I tested eig vs eigs, no
    %         % differences with respect to LEiDA.
    
    % below is needed to ensure most elements have neg values (by
    % convention):
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end
    % I believe just if sum(V1)>0 then V1=-V1 would also be sufficient,
    % as in paper S. Alonso-Martinez
    
    V1Data(t,:)=V1;
    %         % Calculate the variance explained by the leading eigenvector
    %         Var_Eig(s,t)=val1/sum(eVal);
    clear V1 iFC 
end

%% run kmeans clustering using previous centroids, that's why you can leave
% k []; after discussion with Matlab community, better to use pdist2

% compute cityblock distance between observations for subj and the previous
% centroids, and get the corresponding labels (membership indexes)
[~, SubjectMembershipIndex] = pdist2(PreviousCentroids,V1Data,'cityblock','Smallest',1);

SubjectMembershipIndex = SubjectMembershipIndex'; % transpose since it's a row vector

%SubjectMembershipIndex = kmeans(V1Data, [], 'Distance','cityblock',...
%    'MaxIter', 5000, 'Start', PreviousCentroids);

%% extract mean iFC matrix for state 9
% we could decide to look at static phase coherence connections for this
% state as well

tmpIndex9 = sum(SubjectMembershipIndex==9);

if tmpIndex9 > 1
    Mean_iFC_state9 = squeeze(mean(Tot_iFC(SubjectMembershipIndex==9,:,:)));
elseif tmpIndex9 == 1
    Mean_iFC_state9 = Tot_iFC(SubjectMembershipIndex==9,:,:);
else
    Mean_iFC_state9 = nan(N,N);
end

% extract lower triangular portion of matrix and vectorize (connections
% correspond with nchoosek(1:116,2).
Mean_iFC_state9_vector = icatb_mat2vec(Mean_iFC_state9);

clear Tot_iFC Mean_iFC_state9

%% compute state summary measures
[F, ~, MDT, NT] = icatb_dfnc_statevector_stats(SubjectMembershipIndex, 10); % k=10 states

%% compute switch probabilities
indSwitch=[1:10];
SwitchMatrix=zeros(1,10,10);
Ctime=zeros(1,Tmax);

T=1:Tmax;

%check in which state the person starts
Ctime(1)=find(indSwitch==SubjectMembershipIndex(T(1)));

%for all other timepoints
for t=2:length(T)
    Ctime(t)=find(indSwitch==SubjectMembershipIndex(T(t)));
    % create switchmatrix (containing absolute numbers of switches
    % between states)
    SwitchMatrix(1,Ctime(t-1),Ctime(t))=SwitchMatrix(1,Ctime(t-1),Ctime(t))+1;
end

% compute probabilities of swichting from one state to another,
% this is done relative to all switching possibilities per state,
% so per column
for c=1:10
    SwitchProb(1,:,c)=SwitchMatrix(1,c,:)/sum(SwitchMatrix(1,c,:));
end

tmp = SwitchProb(1,:,:);
EasyForStatsProbVector = tmp(:);
% this vector contains all probabilities columnwise extracted from the
% SwitchProb matrix (so is case of 10 states, first value is prob stat 1 to
% state 1, second is state 1 to state 2 ...... columnwise to state 10 to
% state 10). For the end results I disegard the self transitions, but I do
% output these to the table.

%% compute static FC
c = icatb_corr(signaldata); % compute correlations between time courses
c(1:size(c, 1) + 1:end) = 0; % put 0 on diagonal (not really needed, since we extract off diag)
c = icatb_mat2vec(c); % extract one triangular part and put it in vector.
StaticFC = icatb_r_to_z(c); % Fisher z transform correlations.


end
