
% script to compute effect sizes for connectivity findings

% first we do it for the Group effects, by using the wide format file (that
% I had made already for some testing after meeting with Erik), and
% computing the average across visits (marginal mean).
T = readtable('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\PipelineClinicalAssociations\tests_after_stats_meeting\LAPTOP1+2_complete_base_01222024_long_format_merged_back_to_wide.xlsx');
Group = T.DX;

Thal_average = T.ThalamicFC_average;
d_thal_average = computeCohen_d(Thal_average(Group==0), Thal_average(Group==1))

SMACBL_average = T.rSMA_leftCBL_average;
d_SMACBL_average = computeCohen_d(SMACBL_average(Group==0), SMACBL_average(Group==1))

F9_average = T.F9_average;
d_F9_average = computeCohen_d(F9_average(Group==0), F9_average(Group==1))

Precun_average = T.Precun_average;
d_precun_average = computeCohen_d(Precun_average(Group==0), Precun_average(Group==1))

% for Group x Visit, we compute d per visit, so we need long formatted file
T2 = readtable('T:\research\analysis\human\amayer\shared\MAYER_ALL\andy\Hans\LEIDA\Analyses_restFMRI_LAPTOP_total\PipelineClinicalAssociations\TableForMergerSPSS.xlsx');

Visit = str2num(char(T2.Visit));
Group = str2num(char(T2.DX));

F9 = T2.F9

d_v1_F9 = computeCohen_d(F9(Group==0 & Visit ==1), F9(Group==1 & Visit ==1))
d_v2_F9 = computeCohen_d(F9(Group==0 & Visit ==2), F9(Group==1 & Visit ==2))
d_v3_F9 = computeCohen_d(F9(Group==0 & Visit ==3), F9(Group==1 & Visit ==3))

Precun = T2.Precuneus_LwithPrecuneus_R;

d_v1_Prec = computeCohen_d(Precun(Group==0 & Visit ==1), Precun(Group==1 & Visit ==1))
d_v2_Prec = computeCohen_d(Precun(Group==0 & Visit ==2), Precun(Group==1 & Visit ==2))
d_v3_Prec = computeCohen_d(Precun(Group==0 & Visit ==3), Precun(Group==1 & Visit ==3))