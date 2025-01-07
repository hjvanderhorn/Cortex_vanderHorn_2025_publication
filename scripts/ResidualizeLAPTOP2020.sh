#!/bin/bash

##################################################################################################
################### Script to residualize LAPTOP2020 data ########################################
##################################################################################################


# NB: filesize becomes a lot bigger after running this, while I don't adjust the output grid at all, not quite sure why this is.


cd /export/research/analysis/human/amayer/laptop_20019/data/derivatives # the rootfolder for data

file=list_sorted.txt # file with paths to func folders

for i in `cat $file` # loop over paths in List

do

	cd ${i}

	echo ${i}

	rest_input_file=$(find ~+ -type f -name "*_task-rest_preproc-dc.nii.gz")  # distortion corrected input

	rest_tx_epi2T1_file=$(find ~+ -type f -name "*_task-rest-to-T1.aff12.1D")  # warp rest (epi) to T1

	BaseName=$(echo ${rest_tx_epi2T1_file} | cut -d'_' -f1-3) # get basename for output

	cat_matvec ${rest_tx_epi2T1_file} -I > ${BaseName}_T1_2_epi.1D # inverse warp epi to T1 gives us T1 to epi

	# movement parameter files
	rest_movement=$(find ~+ -type f -name "*_task-rest_movement_regressor-mc.1D") 
	rest_movement_deriv=$(find ~+ -type f -name "*_task-rest_movement_regressor-mc-deriv.1D") 


	cd ../anat

	# tissue files
	spm_wm_file=$(find ~+ -type f -name "*_spm-wm_mask.nii.gz") 
	spm_csf_file=$(find ~+ -type f -name "*_spm-csf_mask.nii.gz") 

	# register wm and csf to epi
	3dAllineate -overwrite -master ${rest_input_file} -1Dmatrix_apply ${BaseName}_T1_2_epi.1D -input ${spm_wm_file} -final NN -prefix ${BaseName}_spm-wm_mask_Reg.nii.gz
	3dAllineate -overwrite -master ${rest_input_file} -1Dmatrix_apply ${BaseName}_T1_2_epi.1D -input ${spm_csf_file} -final NN -prefix ${BaseName}_spm-csf_mask_Reg.nii.gz

	# extract tissue signals from epi
	3dROIstats -overwrite -quiet -mask ${BaseName}_spm-wm_mask_Reg.nii.gz ${rest_input_file} > ${BaseName}_Tissue.Regressor.wm.1D
	3dROIstats -overwrite -quiet -mask ${BaseName}_spm-csf_mask_Reg.nii.gz ${rest_input_file} > ${BaseName}_Tissue.Regressor.csf.1D
	

	# regress out tissue signals and motion 
	3dDeconvolve \
		-input ${rest_input_file} \
		-polort 2 \
		-allzero_OK \
		-GOFORIT 10 \
		-num_stimts 14 \
		-stim_file 1 ${BaseName}_Tissue.Regressor.csf.1D\
		-stim_label 1 'CSF'\
		-stim_minlag 1 0\
		-stim_maxlag 1 0\
		-stim_file 2 ${BaseName}_Tissue.Regressor.wm.1D\
		-stim_label 2 'WM'\
		-stim_minlag 2 0\
		-stim_maxlag 2 0\
		-stim_file 3 ${rest_movement}'[1]'\
		-stim_label 3 'ROLL'\
		-stim_minlag 3 0\
		-stim_maxlag 3 0\
		-stim_file 4 ${rest_movement}'[2]'\
		-stim_label 4 'PITCH'\
		-stim_minlag 4 0\
		-stim_maxlag 4 0\
		-stim_file 5 ${rest_movement}'[3]'\
		-stim_label 5 'YAW'\
		-stim_minlag 5 0\
		-stim_maxlag 5 0\
		-stim_file 6 ${rest_movement}'[4]'\
		-stim_label 6 'dS'\
		-stim_minlag 6 0\
		-stim_maxlag 6 0\
		-stim_file 7 ${rest_movement}'[5]'\
		-stim_label 7 'dL'\
		-stim_minlag 7 0\
		-stim_maxlag 7 0\
		-stim_file 8 ${rest_movement}'[6]'\
		-stim_label 8 'dP'\
		-stim_minlag 8 0\
		-stim_maxlag 8 0\
		-stim_file 9 ${rest_movement_deriv}'[1]'\
		-stim_label 9 'ROLL_deriv'\
		-stim_minlag 9 0\
		-stim_maxlag 9 0\
		-stim_file 10 ${rest_movement_deriv}'[2]'\
		-stim_label 10 'PITCH_deriv'\
		-stim_minlag 10 0\
		-stim_maxlag 10 0\
		-stim_file 11 ${rest_movement_deriv}'[3]'\
		-stim_label 11 'YAW_deriv'\
		-stim_minlag 11 0\
		-stim_maxlag 11 0\
		-stim_file 12 ${rest_movement_deriv}'[4]'\
		-stim_label 12 'dS_deriv'\
		-stim_minlag 12 0\
		-stim_maxlag 12 0\
		-stim_file 13 ${rest_movement_deriv}'[5]'\
		-stim_label 13 'dL_deriv'\
		-stim_minlag 13 0\
		-stim_maxlag 13 0\
		-stim_file 14 ${rest_movement_deriv}'[6]'\
		-stim_label 14 'dP_deriv'\
		-stim_minlag 14 0\
		-stim_maxlag 14 0\
		-errts ${BaseName}_task-rest_preproc-dc_resid.nii.gz\
		-jobs 10 \
		-xsave\
		-float\
		-nobucket\
		-overwrite

done
