#!/bin/bash

#####################################################################################################################################################################
######## this script runs bandpass filters residualized REST data, and then warps it into tlrc space using pre-existent warps, and extract AAL time courses. ########
#####################################################################################################################################################################

roi=/export/research/analysis/human/amayer/shared/MAYER_ALL/andy/Hans/LEIDA/AAL/AAL_tta_resampled_3x3x3+tlrc.HEAD # load roi (which is AAL1 atlas converted to tlrc space and resliced to data dimensions)

rootfolder=/export/research/analysis/human/amayer/laptop_20019/data/derivatives # the rootfolder for data

cd ${rootfolder}

file=list_sorted.txt # file with paths to func folders


########### BP filter ############

for i in `cat $file` # loop over paths in List

do
	echo "$i"
	
	cd ${i} # go into the subject/visit folder
	
	FilteredFile=$(find ~+ -type f -name "*_task-rest_preproc-dc_resid.nii.gz") # find the residualized data  (with respect to WM, CSF and motion) ~+ (is called tilde expansion) and gives you the pwd, so this will be pasted before filename to give you the full path
	
	OutputFileName=$(echo ${FilteredFile} | cut -d'.' -f1) # make base output file name 
	
	#3dBandpass -prefix ${OutputFileName}_bp_0.01_0.1 0.01 0.1 ${FilteredFile} # old version, they (the royal AFNI they) recommend 3dTproject now.
	3dTproject -prefix ${OutputFileName}_bp_0.01_0.1.nii.gz -passband 0.01 0.1 -input ${FilteredFile} -overwrite
	
done


########## warping to tlrc and extracting aal tc ############

rootfolder=/export/research/analysis/human/amayer/laptop_20019/data/derivatives # the rootfolder for data

cd ${rootfolder}

file=list_sorted.txt


for i in `cat $file` # loop over paths in List

do
	echo "$i"
	
	cd ${i} # go into the subject/visit folder
	
	# warp into tlrc space
	
	RestREG=$(find ~+ -type f -name "*_task-rest-to-TLRC.NWARP.nii.gz") # find concatenated matrix that contains warp from REST to tlrc
	
	SourceData=$(find ~+ -type f -name "*_bp_0.01_0.1.nii.gz") # find residualized and BP filtered rest data
	
	cd ../anat
	
	Anat_qwarp=$(find ~+ -type f -name "*_T1w_SKSP.TLRC.nii.gz") # non-linear warp (of T1) to TT_N27 (so T1 warped into tlrc space)
	
	OutputFileName2=$(echo ${SourceData} | cut -d'.' -f1-3) # base output file name
	
	3dNwarpApply -overwrite -source ${SourceData} -nwarp ${RestREG} -master ${Anat_qwarp} -newgrid 3 -interp wsinc5 -prefix ${OutputFileName2}_warped_3x3x3.nii.gz
	
	# extract aal tc
	
	cd ${i} # go back into func folder
	
	WarpedDataPresent=$(find ~+ -type f -name "*_warped_3x3x3.nii.gz")
	
	OutputFileName3=$(echo ${WarpedDataPresent} | cut -d'_' -f1-3) # base output file name 
	
	if [ -z "$WarpedDataPresent" ]; then
		echo "No warped data present" # because if there is no warped data we don't want any output (empty aal_tc csv), which will happen if we don't use this if statement
		
		else
		3droistats -overwrite -mask ${roi} ${WarpedDataPresent} > ${OutputFileName3}_aal_tc.csv # extract time courses for rois
	fi
	
done
