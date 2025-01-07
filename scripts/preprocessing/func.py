#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@package func
Preprocessing module
Author: jling@mrn.org

VZ added physio correction

Created on Mon Mar 29 16:28:37 2021

@author: jling
"""


from mayerlab.preprocessing import utils
import os




class pipelines:
    
    def standard_physio(subj):
        
        """
        Parameters
        ----------
        subj : TYPE
            fully process each task with physio correction
            if physio recordings files are available
            if not, uses the standard processing;
            names of physio corrected datasets include pc-
            
            despike(subj,task,run)
            physio_prep(subj,task,run)
            physio_corr(subj,task,run,icard,iresp)
            physio_plot(subj,task,run)
            
            preprocess both uncorrected and corrected data for each run
            iproc=0: uncorrected, iproc=1: physio-corrected
               timeshift(subj,task,run,iproc)
               imReg_2d(subj,task,run,iproc)
               imReg_3d(subj,task,run,iproc)
            estimate_run_motion(subj,task,run)
            
            estimate_field_distortion(subj,task)
            apply_distortion_correction_sbref(subj,task)
            generate_task_mask(subj,task)
            
            process both uncorrected and corrected data for each run
               apply_distortion_correction_run(subj,task,run,iproc)
               create_motion_regressors_run(subj,task,run,iproc)
               
            concatenate_dc_runs(subj,task,run_list,iproc)
            concatenate_motion_regressors(subj,task,run_list,iproc)
            
            spatial_normalization(subj,task)
            
        Returns
        -------
        None.
    
        """
            
        try:
            
            # all get tasks and runs from protocol
            proc_map = subj.get_subj_task_map()
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing" )
            
            for task, run_number in proc_map.items():

                if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
    
                    # pull sbref for this task
                    pull_sbref(subj, task)
                    
                    run_list = range(1, run_number + 1)
                    pc_list = list()
                
                    for run in run_list:
 
                        despike(subj, task, run)
                        physio_prep(subj, task, run)
                        ret = physio_corr(subj, task, run, icard=1, iresp=1)
                        if ret > 0:
                            proc_list = [0,1]
                        else:
                            proc_list = [0]
                        pc_list.append(ret)  
                        physio_plot(subj, task, run)
                                
                        for iproc in proc_list:
                            timeshift(subj, task, run, iproc)
                            imReg_2d(subj, task, run, iproc)
                            imReg_3d(subj, task, run, iproc)
                                
                        estimate_run_motion(subj, task, run)
                    
                    # prepare data for distortion correction 
                    estimate_field_distortion(subj, task)
                    apply_distortion_correction_sbref(subj, task)
                    generate_task_mask(subj, task)
                    
                    # perform distortion correction separately for each fMRI run
                    # because not all runs may have had successful physio correction
                    # also, different CVR runs may have had different CO2 levels
                    
                    for run in run_list:
                        apply_distortion_correction_run(subj, task, run, iproc=0)
                        create_motion_regressors_run(subj, task, run, iproc=0)
                        if pc_list[run-1] == 1:
                            apply_distortion_correction_run(subj, task, run, iproc=1)
                            create_motion_regressors_run(subj, task, run, iproc=1)
                            
                    # concatenate physio-uncorrected datasets
                    concatenate_dc_runs(subj, task, run_list, iproc=0)
                    concatenate_motion_regressors(subj, task, run_list, iproc=0)
                    
                    # concatenate physio-corrected datasets, only if all exist
                    if sum(pc_list) == run_number:
                        concatenate_dc_runs(subj, task, run_list, iproc=1)
                        concatenate_motion_regressors(subj, task, run_list, iproc=1)
                        
                    qc.snapshot_mask(subj, task)       
                    proc_cleanup(subj, task, 0)
                    proc_cleanup(subj, task, 1)
                   
                # for each task regardless of 'preproc-dc.nii.gz' status
                spatial_normalization(subj, task)
                qc.snapshot_spatial_normalization(subj, task)
                    
            # after all tasks complete
            qc.aggregate_qc_pdf(subj)
            qc.run_mriqc(subj)
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing Complete")
            
        except:
            raise
 

 
            
    def standard(subj):
        
        """
        
        Parameters
        ----------
        subj : TYPE
            fully process each task.
            
            despike(subj,study,run)
            timeshift(subj,study,run)
            imReg_2d(subj,study,run)
            imReg_3d(subj,study,run)
            estimate_run_motion(subj,study,run)
                
            #    finalize each run type
            #    condense each run into a final preproc file
            #    ready for level1
            #
            #    combine runs into single volume
            #    apply field distortion correction to unwarp data
            #    mask for run_type (from SBRef)
            #    affine alignment. Keep separate so that user can redo just this step
            #    full concatenated alignment path through affine and non-linear
            
            
        Returns
        -------
        None.
    
        """
            
        try:
            
            # all get tasks and runs from protocol
            proc_map = subj.get_subj_task_map()
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing" )
            
            for task, run_number in proc_map.items():

                if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
                    
                    if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc.nii.gz'))):
    
                        # do this before each task is run
                        # pull sbref for this task
                        pull_sbref(subj, task)
                
                        for run in range(1, run_number + 1):
                            # do this for each run            
                            despike(subj, task, run)
                            timeshift(subj, task, run)
                            imReg_2d(subj, task, run)
                            imReg_3d(subj, task, run)
                            estimate_run_motion(subj, task, run)
                   
                    
                    # do this for the task after each run processed
                    concate_runs(subj, task, run_list=range(1, run_number + 1))
                    create_motion_regressors(subj, task, run_list=range(1, run_number + 1))
                    create_motion_deriv_file(subj, task, run_list=range(1, run_number + 1))
                    estimate_field_distortion(subj, task)
                    apply_distortion_correction_sbref(subj, task)
                    generate_task_mask(subj, task)
                    apply_distortion_correction_concat(subj, task)
                    qc.snapshot_mask(subj, task)
                    clean_proc(subj, task)
                    

                # for each task regardless of 'preproc-dc.nii.gz' status
                spatial_normalization(subj, task)
                qc.snapshot_spatial_normalization(subj, task)
                    
                
                
            # after all tasks complete
            qc.aggregate_qc_pdf(subj)
            qc.run_mriqc(subj)
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing Complete")
            
        except:
            raise            


            
        
    def laptop(subj):
        
        try:
            
            # for LAPTOP you need the whole process
            # get tasks and runs from protocol
            proc_map = subj.get_subj_task_map()
    
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing")
            
            
            for task, run_number in proc_map.items():
                
                if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
                    
                    if not os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc.nii.gz'))):
                        
                        # do this before each task is run
                        # pull sbref for this task
                        pull_sbref(subj, task)
                
                        for run in range(1, run_number + 1):
                            # do this for each run            
                            despike(subj, task, run)
                            timeshift(subj, task, run)
                            imReg_2d(subj, task, run)
                            imReg_3d(subj, task, run)
                            estimate_run_motion(subj, task, run)
                   
                    
                    # do this for the task after each run processed
                    concate_runs(subj, task, run_list=range(1, run_number + 1))
                    create_motion_regressors(subj, task, run_list=range(1, run_number + 1))
                    create_motion_deriv_file(subj, task, run_list=range(1, run_number + 1))
                    estimate_field_distortion(subj, task)
                    apply_distortion_correction_sbref(subj, task)
                    generate_task_mask(subj, task)
                    apply_distortion_correction_concat(subj, task)
                    spatial_normalization(subj, task, space='TLRC')
                    qc.snapshot_mask(subj, task)
                    qc.snapshot_spatial_normalization(subj, task, space='TLRC')

                    clean_proc(subj, task)

            qc.aggregate_qc_pdf(subj)
            qc.run_mriqc(subj)

            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : Preprocessing Complete" )
        
        except:
            raise        




    
def physio_prep(subj, task, run):
    
    """
    VZ prepare physio data for fMRI correction

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    Returns
    -------
    None.
    
    """
    
    from pathlib import Path
    import numpy as np
    import json
    
    # Siemens clock interval (2.5 ms)
    TIC = 2.5e-3
    
    # puls sampling interval: 2 * TIC
    # resp sampling interval: 8 * TIC
    
    func_path = subj.derivatives_path('func')  
    
    infolog = subj.bids_physio_path(bidsType='func', bidsLabel='physio', task=task, run=run, extension='_info.log')
    pulslog = subj.bids_physio_path(bidsType='func', bidsLabel='physio', task=task, run=run, extension='_puls.log')
    resplog = subj.bids_physio_path(bidsType='func', bidsLabel='physio', task=task, run=run, extension='_resp.log')
    
    if (infolog and os.path.isfile(infolog)) and (pulslog and os.path.isfile(pulslog)) and (resplog and os.path.isfile(resplog)):
        
        out_info = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_info.1D'))
        out_puls = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_puls.1D'))
        out_resp = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_resp.1D'))
        out_para = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_params.1D'))
        out_slof = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'slice_offsets.1D'))
        
        # additional files for plotting
        out_puls2 = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_puls2.1D'))
        out_resp2 = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_resp2.1D'))

        # extract sequence parameters from json file
        bold_json = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='json')
        with open(bold_json, 'r') as json_file:
            json_data = json.load(json_file)
        TR = json_data['RepetitionTime']
        ST = np.array(json_data['SliceTiming'], dtype=np.float32)
        NS = len(ST)
 
        # parsing info log file
        infolog = Path(infolog).resolve()
        lines = infolog.read_text().splitlines()
    
        for line in [line for line in lines if line]:

            # strip any leading and trailing whitespace and comments
            line = line.split('#')[0].strip()

            if '=' in line:
                varname, value = [item.strip() for item in line.split('=')]
                if varname == 'NumSlices':
                    nrslices = int(value)
                if varname == 'NumVolumes':
                    nrvolumes = int(value)
                if varname == 'FirstTime':
                    firsttime = int(value)
                if varname == 'LastTime':
                    lasttime = int(value)

        # default puls array size for 10 min run
        P = np.zeros((120000,2), dtype=np.int32)
        ind = int(0)
    
        # parsing puls log file
        pulslog = Path(pulslog).resolve()
        lines = pulslog.read_text().splitlines()
    
        for line in [line for line in lines if line]:

            # strip any leading and trailing whitespace and comments
            line = line.split('#')[0].strip()

            if '=' in line:
                varname, value = [item.strip() for item in line.split('=')]
                if varname == 'SampleTime':
                    puls_samp = int(value)
            else:
                # puls data
                dataitems = line.split()
                dataitems = [dataitems[n] if n < len(dataitems) else '0' for n in range(5)]

                # if the first column isn't numeric, it is probably the header
                if not dataitems[0].isdigit():
                    continue
            
                curtic = int(dataitems[0])
                curval = int(dataitems[2])
                
                if ind > 0:
                    # compare with previous point
                    tdiff = curtic - pretic
                    vdiff = curval - preval
                    steps = int(np.rint(tdiff / puls_samp))
                    if steps == 0:
                        steps = 1
                    add = int(np.rint(vdiff / steps))
                    # fill in missing points if any
                    #if (steps > 1):
                        #print("%5d %1d %8d %8d %4d %4d" % (ind, steps, pretic, curtic, preval, curval))
                    for j in range(1,steps):
                        P[ind,0] = pretic + j * puls_samp
                        P[ind,1] = preval + j * add
                        #print("%5d %8d %4d" % (ind, P[ind,0], P[ind,1]))
                        ind = ind + 1
                    P[ind,0] = curtic
                    P[ind,1] = curval
                    pretic = curtic
                    preval = curval
                    ind = ind + 1

                else:
                    P[ind,0] = curtic
                    P[ind,1] = curval
                    pretic = curtic
                    preval = curval
                    ind = ind + 1
            
        P = P[0:ind,:]
        Np = len(P)
        #print(Np)
              
        # default resp array size for 10 min run
        R = np.zeros((30000,2), dtype=np.int32)
        ind = int(0)
    
        # parsing resp log file
        resplog = Path(resplog).resolve()
        lines = resplog.read_text().splitlines()
    
        for line in [line for line in lines if line]:

            # strip any leading and trailing whitespace and comments
            line = line.split('#')[0].strip()

            if '=' in line:
                varname, value = [item.strip() for item in line.split('=')]
                if varname == 'SampleTime':
                    resp_samp = int(value)
            else:
                # resp data
                dataitems = line.split()
                dataitems = [dataitems[n] if n < len(dataitems) else '0' for n in range(5)]

                # if the first column isn't numeric, it is probably the header
                if not dataitems[0].isdigit():
                    continue
            
                curtic = int(dataitems[0])
                curval = int(dataitems[2])
                
                if ind > 0:
                    # compare with previous point
                    tdiff = curtic - pretic
                    vdiff = curval - preval
                    steps = int(np.rint(tdiff / resp_samp))
                    if steps == 0:
                        steps = 1
                    add = int(np.rint(vdiff / steps))
                    # fill in missing points if any
                    #if (steps > 1):
                        #print("%5d %1d %8d %8d %4d %4d" % (ind, steps, pretic, curtic, preval, curval))
                    for j in range(1,steps):
                        R[ind,0] = pretic + j * resp_samp
                        R[ind,1] = preval + j * add
                        #print("%5d %8d %4d" % (ind, R[ind,0], R[ind,1]))
                        ind = ind + 1
                    R[ind,0] = curtic
                    R[ind,1] = curval
                    pretic = curtic
                    preval = curval
                    ind = ind + 1

                else:
                    R[ind,0] = curtic
                    R[ind,1] = curval
                    pretic = curtic
                    preval = curval
                    ind = ind + 1
            
        Nr = ind
        
        # resp acquisition is shorter than puls acquisition
        # check if can add an extra point to make them closer
        
        lastpulstic = P[Np-1,0]
        lastresptic = R[Nr-1,0]
        
        if lastpulstic - lastresptic >= resp_samp:
            R[Nr,0] = lastresptic + resp_samp
            R[Nr,1] = R[Nr-1,1]
            R = R[0:Nr+1,:]
        else:
            R = R[0:Nr,:]
         
        Nr = len(R)
        #print(Nr)

        # number of resp data points for the sequence duration
        Npts = int(np.rint((TR * nrvolumes) / (TIC * resp_samp)))
        
        Pout = np.zeros((Npts,2), dtype=np.int32)
        Rout = np.zeros((Npts,2), dtype=np.int32)
        Rout = R[Nr-Npts:Nr,:]

        # find puls point closest to the last resp point
        lastresptic = Rout[Npts-1,0]
        exitflag = 0
        j = 1
        while exitflag == 0:
            pulsind = Np - j
            pulstic = P[pulsind,0]
            #print("%5d %8d" % (pulsind, pulstic))
            if abs(pulstic - lastresptic) < puls_samp:
                exitflag = 1
            else:
                j = j + 1
            
        # sub-sample puls data to resp sampling
        step = int(np.rint(resp_samp / puls_samp))
        N1 = pulsind - (Npts - 1) * step
        N2 = pulsind + step
        #print("%5d %8d %5d %5d" % (pulsind, pulstic, N1, N2))
        Pout = P[N1:N2:step,:]
        
        # save puls and resp values
        np.savetxt(out_puls, Pout[:,1], fmt='%4d', delimiter=' ')
        np.savetxt(out_resp, Rout[:,1], fmt='%4d', delimiter=' ')
        
        # save both tics and values for puls and resp
        finfo = open(out_info, "w")
        for i in range(Npts):
            tv = TIC * resp_samp * i
            pt = Pout[i,0]; pv = Pout[i,1]
            rt = Rout[i,0]; rv = Rout[i,1]
            print("%5d  %7.3f  %8d  %4d  %8d  %4d" % (i, tv, pt, pv, rt, rv), file=finfo)
        finfo.close()
        
        # save slice offsets
        fslof = open(out_slof, "w")
        for i in range(NS):
            print("%8.6f" % ST[i], file=fslof)
        fslof.close()    
        
        # save parameters needed for RetroTS
        fpara = open(out_para, "w")
        freq = int(np.rint(1.0e+0 / ( TIC * resp_samp )))
        print("%5.3f %d %d %d" % (TR, NS, nrvolumes, freq), file=fpara)
        fpara.close()   
        
        # prepare arrays for plotting physio traces in several subplots
        # sublot length in seconds
        Lsec = 60
        Lpts = Lsec * freq
        # number of subplots
        Nsub = int(np.ceil(Npts / Lpts))
        Pcol = np.zeros((Lpts*Nsub,1), dtype=np.int32)
        Rcol = np.zeros((Lpts*Nsub,1), dtype=np.int32)
        for i in range(Npts):
            Pcol[i,0] = Pout[i,1]
            Rcol[i,0] = Rout[i,1]
        Pmat = np.reshape(Pcol, (Lpts,Nsub), order='F')
        Rmat = np.reshape(Rcol, (Lpts,Nsub), order='F')
        
        # save reshaped puls and resp arrays
        np.savetxt(out_puls2, Pmat, fmt='%4d', delimiter=' ')
        np.savetxt(out_resp2, Rmat, fmt='%4d', delimiter=' ')
        
        print("Successfully prepared data for physio correction")
        
    else:
        print("One or more of required physio log files is not found")

        
        
        
def physio_corr(subj, task, run, icard=1, iresp=1):
    
    """
    VZ correct fMRI data using physio waveforms

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    icard : include correction of cardiac artifacts
            set to 0 to exclude from correction, e.g. if PULS recording is unreliable
    iresp : include correction of respiratiory artifacts
            set to 0 to exclude from correction, e.g. if RESP recording is unreliable
    Returns
    -------
    1 - if correction is successful or the corrected dataset exists
    0 - if correction is unsuccessful and no corrected dataset written
    
    """
    
    import math
    
    func_path = subj.derivatives_path('func')  
    
    puls1D = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_puls.1D'))
    resp1D = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_resp.1D'))
    para1D = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_params.1D'))
    slof1D = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'slice_offsets.1D'))
    
    if os.path.isfile(puls1D) and os.path.isfile(resp1D) and os.path.isfile(para1D) and os.path.isfile(slof1D):
        
        i_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'dspk.nii.gz'))
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-dspk.nii.gz'))
        
        if not os.path.isfile(o_file):
            utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : Physio Correction")
            
            # extract sequence parameters from para1D file
            for line in open(para1D):
                fields = line.split(" ")
                TR = float(fields[0])
                NS = int(fields[1])
                NV = int(fields[2])
                FR = int(fields[3])
                
            prefT = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio'))
            fileT = prefT + ".slibase.1D"
            fileD = prefT + "_sli_det.1D"
            cens  = os.path.join(func_path, subj.bids_file_join('cens.1D'))
            temp1 = os.path.join(func_path, subj.bids_file_join('temp1.1D'))
            temp2 = os.path.join(func_path, subj.bids_file_join('temp2.1D'))
            temp3 = os.path.join(func_path, subj.bids_file_join('temp3.1D'))
            temp4 = os.path.join(func_path, subj.bids_file_join('temp4.1D'))
            
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('temp*')))
            utils.purge_paths(fileT, fileD, cens)
            
            # define flags for customized output
            out_flags = ""
            if icard == 0:
                out_flags = out_flags + " -cardiac_out 0 "
            if iresp == 0:
                out_flags = out_flags + " -rvt_out 0 -respiration_out 0 "
            
            try:
                # run RetroTS.py to determine physio regressors for each slice
                utils.Cmd("RetroTS.py" + out_flags + " -c " + puls1D + " -r " + resp1D + " -p " + str(FR) + 
                      " -n " + str(NS) + " -v " + str(TR) + " -slice_order " + slof1D + " -prefix " + prefT).run()
            except:
                print("RetroTS.py program failed to complete")
                return 0
            
            # analyze RetroTS.py output
            if os.path.isfile(fileT):
                try:
                    fT = open(fileT, 'r')
                except:
                    print("Could not open the RetroTS.py output file:")
                    print(fileT)
                    return 0
                    
                dataT = fT.read()
                all_lines = dataT.splitlines()
                for lind in range(len(all_lines)):
                    full_line = all_lines[lind]
                    cind = full_line.find('ColumnLabels')
                    if cind >= 0:
                        n0 = full_line.count('s0')
                        n1 = full_line.count('RVT')
                        n2 = full_line.count('Resp')
                        n3 = full_line.count('Card')
                        break
                # number of regressor time series generated for each slice
                nr = math.floor((n1 + n2 + n3) / NS)
                fT.close()
                
                if nr == n0:
                    regnum = nr
                    #print(regnum)
                else:
                    print("The RetroTS.py output file has incorrect number of columns:")
                    print(fileT)
                    return 0                
            else:
                print("The RetroTS.py output file is not found:")
                print(fileT)
                return 0
            
            # define polort automatically as in 3dDeconvolve
            polort = 1 + math.floor(TR * NV / 150)
            #print(polort)
            
            # define censor for the first time point - the 1st volume will be excluded anyway
            utils.Cmd("1deval -overwrite -start 0 -num " + str(NV) + " -expr 'step(t)' > " + cens).run()
            
            # perform polynomial detrending 
            utils.Cmd("1dtranspose " + fileT + " " + temp1).run()
            utils.Cmd("3dTproject -cenmode NTRP -censor " + cens + " -polort " + str(polort) + " -input " + temp1 + " -prefix " + temp2).run()
            utils.Cmd("1dtranspose " + temp2 + " " + temp3).run()
            utils.Cmd("1dnorm -normx " + temp3 + " " + fileD).run()
                
            # prepare individual slices using one dspk volume
            for sl in range(NS):
                fcut_sl = "zcut%02d" % sl
                fcut_un = "zcut%02d_un" % sl
                prefS = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fcut_sl))
                prefU = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fcut_un))
                utils.Cmd("3dZcutup -overwrite -keep " + str(sl) + " " + str(sl) + " -prefix " + prefS + " " + i_file + "[1]").run()
                utils.Cmd("3dcalc -overwrite -expr 'notzero(a)+iszero(a)' -a " + prefS + "+orig -prefix " + prefU).run()

            # remove all existing regressors, in case different regressors are selected via out_flags
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_reg*')))
            
            # prepare 3D+time regressor for each ireg
            for ireg in range(regnum):
                
                freg = "physio_reg%02d.nii.gz" % ireg
                fileR = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), freg))
                
                for sl in range(NS):
                    offset =  regnum * sl
                    col = offset + ireg
                    fcut_un = "zcut%02d_un" % sl
                    fcut_ts = "zcut%02d_un_ts" % sl
                    prefU = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fcut_un))
                    prefM = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fcut_ts))
                    utils.Cmd("1deval -overwrite -expr 'a' -a " + fileD + "'['" + str(col) + "']' > " + temp4).run()
                    utils.Cmd("3dcalc -overwrite -expr 'a*b' -TR " + str(TR) + " -a " + prefU + "+orig -b " + temp4 + " -prefix " + prefM).run()

                fcat = "zcut??_un_ts+orig.HEAD"
                exprC = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fcat))
                utils.Cmd("3dZcat -overwrite -datum float -prefix " + fileR + " " + exprC).run()
                # prepare a string of physio regressor file names
                if ireg == 0:
                    reg_str = fileR
                else:
                    reg_str = reg_str + " " + fileR
                
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'zcut*')))
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('temp*')))
            utils.purge_paths(cens)
            
            # prepare a uniform volume using one dspk volume
            prefV = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'vol1'))
            utils.Cmd("3dcalc -overwrite -expr 'notzero(a)+iszero(a)' -a " + i_file + "[1] -prefix " + prefV).run()
            
            # prepare 3D+time regressors using Legendre polynomials up to polort
            for ipol in range(polort+1):
                
                fpol = "physio_pol%1d.nii.gz" % ipol
                fileP = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fpol))
                
                fp = "poly%1d.1D" % ipol
                pol = os.path.join(func_path, subj.bids_file_join(fp))
                exp = "'Pleg(" + str(ipol) + ", 2*t/" + str(NV-1) + "-1)'"
                utils.Cmd("1deval -overwrite -start 0 -num " + str(NV) + " -expr " + exp + " > " + pol).run()

                utils.Cmd("3dcalc -overwrite -expr 'a*b' -TR " + str(TR) + " -a " + prefV + "+orig -b " + pol + " -prefix " + fileP).run()
                # prepare a string of polynomial regressor file names
                if ipol == 0:
                    pol_str = fileP
                else:
                    pol_str = pol_str + " " + fileP
              
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'vol1*')))
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('poly*')))
            
            #print(reg_str)
            #print(pol_str)
            
            # run 3dTfitter program with physio regressors and polynomial terms
            fileB = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_bet.nii.gz'))
            fileF = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_fit.nii.gz'))   
            
            utils.Cmd("3dTfitter -overwrite -l2fit -RHS " + i_file + " -LHS " + reg_str + " " + pol_str + " -prefix " + fileB + " -fitts " + fileF).run()
            
            # compute polynomial baseline using beta coefficients
            fileS = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_sum.nii.gz'))
            utils.Cmd("3dcalc -overwrite -expr 'a*0' -a " + fileF + " -prefix " + fileS).run()
            for ipol in range(polort+1):
                ibet = regnum + ipol
                fpol = "physio_pol%1d.nii.gz" % ipol
                fileP = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), fpol))
                utils.Cmd("3dcalc -overwrite -expr 'a*b+c' -a " + fileP + " -b " + fileB + "[" + str(ibet) + "] -c " + fileS + " -prefix " + fileS).run()
                    
            # perform correction: subtract the fit and restore the polynomial baseline
            utils.Cmd("3dcalc -overwrite -expr 'posval(a-b+c)' -a " + i_file + " -b " + fileF + " -c " + fileS + " -prefix " + o_file).run()
            
            # compute temporal SNR before and after
            tsnr1 = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'tsnr_dspk.nii.gz'))
            tsnr2 = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'tsnr_pc-dspk.nii.gz'))
            tsnrC = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'tsnr_pc-dspk_pch.nii.gz'))
            stat1 = os.path.join(func_path, subj.bids_file_join('stat1.nii.gz'))
            stat2 = os.path.join(func_path, subj.bids_file_join('stat2.nii.gz'))
            
            utils.Cmd("3dTstat -overwrite -mean -stdev -prefix " + stat1 + " " + i_file).run()
            utils.Cmd("3dcalc -overwrite -expr 'step(a-100)*a/b' -a " + stat1 + "[0] -b " + stat1 + "[1] -prefix " + tsnr1).run()
            utils.Cmd("3dTstat -overwrite -mean -stdev -prefix " + stat2 + " " + o_file).run()
            utils.Cmd("3dcalc -overwrite -expr 'step(a-100)*a/b' -a " + stat2 + "[0] -b " + stat2 + "[1] -prefix " + tsnr2).run()
            utils.Cmd("3dcalc -overwrite -expr '100*step(b-1)*(a-b)/b' -a " + tsnr2 + " -b " + tsnr1 + " -prefix " + tsnrC).run()
            
            # clean up
            utils.purge_paths(fileF, fileS)
            utils.purge_paths(os.path.join(func_path, subj.bids_file_join('stat*')))
            
            print("The physio correction procedure is completed successfully")
            return 1
        
        else:
            print("The output physio corrected dataset already exists:")
            print(o_file)
            return 1
        
    else:
        print("One or more of required input 1D files is not found:")
        print(puls1D + "\n" + resp1D + "\n" + para1D + "\n" + slof1D)
        return 0




def physio_plot(subj, task, run):
    
    """
    VZ plot data related to physio correction

    """
    
    func_path = subj.derivatives_path('func')  
    
    puls = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_puls2.1D'))
    resp = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_resp2.1D'))
    regs = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_sli_det.1D'))
        
    puls_p = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_puls.png'))
    resp_p = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_resp.png'))
    regs_p = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'physio_regs_s0.png'))
    
    puls_t = str(subj.ursi.full) + "_v" + str(subj.visit) + "_" + task + "_run" + str(run) + "_cardiac"
    resp_t = str(subj.ursi.full) + "_v" + str(subj.visit) + "_" + task + "_run" + str(run) + "_respiratory"
    regs_t = str(subj.ursi.full) + "_v" + str(subj.visit) + "_" + task + "_run" + str(run) + "_regressors_for_slice0"
    
    phys_str = "1dplot.py -dpi 300 -xlabel 'each subplot is 60 sec' "
    regs_str = "1dplot.py -dpi 300 -xlabel 'volumes' "
    
    if os.path.isfile(puls):
        utils.Cmd(phys_str + "-colors red -title " + puls_t + " -prefix " + puls_p + " -infiles " + puls).run()
        
    if os.path.isfile(resp):
        utils.Cmd(phys_str + "-colors blue -title " + resp_t + " -prefix " + resp_p + " -infiles " + resp).run()
        
    if os.path.isfile(regs):
        # plot the first 13 time series
        utils.Cmd(regs_str + "-title " + regs_t + " -prefix " + regs_p + " -infiles " + regs + "'[0..12]'").run()




def pull_sbref(subj, task):
    """
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    import shutil
    func_path = subj.derivatives_path('func')
    utils.make_path(func_path)
    
    i_file = subj.bids_file_path(bidsType='func', bidsLabel='sbref', task=task, extension='nii.gz')
    if not os.path.isfile(i_file):
        raise ValueError("Cannot pull" + i_file)
        
    o_file = os.path.join(func_path, os.path.basename(i_file))
    if not os.path.isfile(o_file):
        shutil.copyfile(i_file, o_file)

 
 
    
def slice_timing_to_file(subj, i_file_json, o_file_txt):
    """
    Extract slice timeing from bids sidecar into text file

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    import json
        
    try:
        with open(i_file_json, 'r') as json_file:
            side_car_data = json.load(json_file)
        
        with open(o_file_txt, 'w') as timing_file:
            timing_file.write(" ".join(map(str, side_car_data['SliceTiming'])))
    
        return 
    except:
        raise
     



def despike(subj, task, run, ignore=1):
    """
    despike raw fMRI data

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    ignore : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """

    func_path = subj.derivatives_path('func')
    utils.make_path(func_path)
    i_file = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='nii.gz')
    o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'dspk.nii.gz'))
    
    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : Despike")

        # run despike
        utils.Cmd("3dDespike -overwrite -NEW -ignore " + str(ignore) + " -nomask -cut 4 4.5 -prefix "+o_file+" " + i_file).run()




def timeshift(subj, task, run, iproc=0, ignore=1):
    """
    timeshift despike data

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    iproc: TYPE, optional
        DESCRIPTION: 0 - use dspk dataset
                     1 - use pc-dspk dataset
    ignore : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    
    func_path = subj.derivatives_path('func')  
    
    if iproc > 0:
        i_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-dspk.nii.gz'))
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-tshft.nii.gz'))
    else:
        i_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'dspk.nii.gz'))
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'tshft.nii.gz'))        

    timing_json = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='json')
    timing_txt = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'slice_timing.1D'))
    slice_timing_to_file(subj, timing_json, timing_txt)

    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : Slice Timing Correction")
        utils.Cmd("3dTshift -overwrite -tpattern @"+timing_txt+" -ignore "+ str(ignore) + " -tzero 0 -prefix "+ o_file +" " + i_file).run()




def imReg_2d(subj, task, run, iproc=0):
    """
    2d registration of despike->tshift output to SBREF base file

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    iproc: TYPE, optional
        DESCRIPTION: 0 - use tshft dataset
                     1 - use pc-tshft dataset
    Returns
    -------
    None.

    """

    # registration target is SBREF
    base_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    # input is despike->tshift output, save 2d reg parameters
    
    if iproc > 0:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-tshft.nii.gz'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-2dreg.nii.gz'))
        reg_params_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-2dreg'))
    else:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'tshft.nii.gz'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '2dreg.nii.gz'))
        reg_params_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '2dreg'))        

    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : 2D Registration")
        utils.Cmd(f"2dImReg -overwrite -input {i_file} -basefile {base_file} -base 0 -prefix {o_file} -dprefix {reg_params_file}").run()

        
        

def imReg_3d(subj, task, run, iproc=0):
    """
    3d registration of despike->tshift->2dreg data to SBREF base file to correct motion
    This outputs parameters used as motion regressors. See create_run_motion_regressors()
    ".Movement.Regressor.mc.1D"):

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    iproc: TYPE, optional
        DESCRIPTION: 0 - use 2dreg dataset
                     1 - use pc-2dreg dataset
    Returns
    -------
    None.

    """

    # registration target is SBREF
    base_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    # despike->tshift->2dreg input
    
    if iproc > 0:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-2dreg.nii.gz'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-3dreg.nii.gz'))
        mat_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-mat.mc'))
        mot_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-log.mc'))
    else:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '2dreg.nii.gz'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '3dreg.nii.gz'))
        mat_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'mat.mc'))
        mot_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'log.mc'))        
    
    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : 3D Registration")
        utils.Cmd("3dvolreg -overwrite -dfile " + mot_file + " -1Dmatrix_save " + mat_file + " -prefix " + o_file + " -base " + base_file + " " + i_file ).run()




def estimate_run_motion(subj, task, run):
    """
    Estimate movement parameters
    Register the raw, unprocessed run data to SBRef to estimate motion
    This is used for motion quantification

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run : TYPE
        DESCRIPTION.
    Returns
    -------
    None.

    """
    
    i_file = subj.bids_file_path(bidsType='func', bidsLabel='bold', task=task, run=run, extension='nii.gz')
    base_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '3dreglog.raw'))
    md_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'maxdisp.raw'))
    
    tmp = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('tmp.1D'))
    tmpd = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('tmp_delt.1D'))
        
    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : Estimate raw motion parameters")
        utils.Cmd("3dvolreg -overwrite -dfile " + o_file + " -prefix NULL -base 0 -base " + base_file + " -maxdisp1D " + tmp + " " + i_file).run()
        utils.Cmd("1dcat -overwrite -form '%6.3f' " + tmp + " " + tmpd + " > " + md_file).run()
        utils.purge_paths(tmp, tmpd)




def concate_runs(subj, task, run_list, ignore=1):
    """
    

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    func_path = subj.derivatives_path('func')
    run_input_paths = []
    for run_number in run_list:
        run_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task, 'run-'+str(run_number), '3dreg.nii.gz')))

    o_file = os.path.join(func_path, subj.bids_file_join('task-' + task, 'preproc-uncorrected.nii.gz'))

    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : Concatenate Runs")
        utils.Cmd("3dTcat -prefix " + o_file + " " + " ".join([file + "'["+str(ignore)+"-$]'" for file in run_input_paths])).run()

        
        

def concatenate_dc_runs(subj, task, run_list, iproc=0):
    
    """
    VZ - concatenate distortion corrected fMRI runs
    initial volumes (specified by ignore) are already removed
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run_list : TYPE
        DESCRIPTION.
    iproc: TYPE, optional
        DESCRIPTION: 0 - use preproc-dc datasets
                     1 - use pc-preproc-dc datasets
    Returns
    -------
    None.

    """
    func_path = subj.derivatives_path('func')
    run_input_paths = []
    
    if iproc > 0:
        for run in run_list:
            run_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+str(run), 'pc-preproc-dc.nii.gz')))
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task, 'pc-preproc-dc.nii.gz'))
    else:
        for run in run_list:
            run_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+str(run), 'preproc-dc.nii.gz')))
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))

    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + " : Concatenate distortion corrected runs")
        utils.Cmd("3dTcat -overwrite -prefix " + o_file + " " + " ".join([file for file in run_input_paths])).run()
        


        
def concatenate_motion_regressors(subj, task, run_list, iproc=0):
    
    """
    VZ concatenate motion regressors and derivatives across runs
    also concatenate raw motion regressors and max displacements
    initial points (specified by ignore) are already removed
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.
    run_list : TYPE
        DESCRIPTION.
    iproc: TYPE, optional
        DESCRIPTION: 0 - use mreg-mc.1D and mreg-mc-deriv.1D
                     1 - use pc-mreg-mc.1D and pc-mreg-mc-deriv.1D
    Returns
    -------
    None.

    """
    func_path = subj.derivatives_path('func')
    mot_input_paths = []
    der_input_paths = []
    raw_input_paths = []
    max_input_paths = []
    
    if iproc > 0:
        for run in run_list:
            mot_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+ str(run), 'pc-mreg-mc.1D')))
            der_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+ str(run), 'pc-mreg-mc-deriv.1D')))
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'pc-movement_regressor-mc.1D'))
        d_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'pc-movement_regressor-mc-deriv.1D'))
    else:
        for run in run_list:
            mot_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+ str(run), 'mreg-mc.1D')))
            der_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+ str(run), 'mreg-mc-deriv.1D')))
            raw_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+ str(run), 'mreg-raw.1D'))) 
            max_input_paths.append(os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'run-'+ str(run), 'maxd-raw.1D'))) 
        o_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc.1D'))
        d_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc-deriv.1D'))
        r_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'movement_regressor.1D'))
        x_file = os.path.join(func_path, subj.bids_file_join('task-' + task.lower(), 'max_displacement.1D'))

    if not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + " : Concatenate movement regressors and derivatives")
        utils.Cmd("cat " + ' '.join(mot_input_paths) + " > " + o_file).run()
        utils.Cmd("cat " + ' '.join(der_input_paths) + " > " + d_file).run()
        
    if iproc <= 0:
        if not os.path.isfile(r_file):
            utils.Cmd("cat " + ' '.join(raw_input_paths) + " > " + r_file).run()
        if not os.path.isfile(x_file):
            utils.Cmd("cat " + ' '.join(max_input_paths) + " > " + x_file).run()
        


        
def create_motion_regressors_run(subj, task, run, iproc=0, ignore=1):
    
    """
    VZ movement regressors suitable for a single run level1 analysis
    also removes initial points specified by ignore
    also computes derivatives of movement regressors
    also prepares movement regressors for raw data
    
    iproc: TYPE, optional
        DESCRIPTION: 0 - use log.mc files
                     1 - use pc-log.mc files

    """
    
    if iproc > 0:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-log.mc'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-mreg-mc.1D'))
        d_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-mreg-mc-deriv.1D'))
    else:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'log.mc'))   
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'mreg-mc.1D')) 
        d_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'mreg-mc-deriv.1D'))
        
    r_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '3dreglog.raw'))
    x_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'maxdisp.raw'))
    mr_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'mreg-raw.1D'))
    mx_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'maxd-raw.1D'))
    
    if os.path.isfile(i_file) and not os.path.isfile(o_file):
        utils.prettyOut('task-' + task.lower() + ' run-' + str(run) + " : movement regressors and derivatives")
        utils.Cmd("tail -q -n +" + str(ignore+1) + " " + i_file + " > " + o_file).run()
        utils.Cmd("1d_tool.py -overwrite -infile " + o_file + " -derivative -write " + d_file).run()
        
    if iproc <= 0: 
        if not os.path.isfile(mr_file):
            utils.Cmd("tail -q -n +" + str(ignore+1) + " " + r_file + " > " + mr_file).run()
        if not os.path.isfile(mx_file):
            utils.Cmd("tail -q -n +" + str(ignore+1) + " " + x_file + " > " + mx_file).run()
        

    
        
        
def create_motion_regressors(subj, task, run_list, ignore=1):
    """
    creates file appropriate for level1

    """

    mot_quant_input_paths = []
    mot_reg_input_paths = []
    for run_number in run_list:
        mot_quant_input_paths.append(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run_number), '3dreglog.raw')))
        mot_reg_input_paths.append(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run_number), 'log.mc')))
        
    # original movement regressor created with raw epi input. This is used for motion quantification
    # create from available run data
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor.1D'))
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : movement_regressor.1D")
        utils.Cmd("tail -q -n +" + str(ignore+1) + " " + ' '.join(mot_quant_input_paths) +" > " + o_file).run()


    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc.1D'))
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : movement_regressor-mc.1D")
        utils.Cmd("tail -q -n +" + str(ignore+1) + " " + ' '.join(mot_reg_input_paths) +" >> " + o_file).run()


    
    
def create_motion_deriv_file(subj, task, run_list):

    """
    Create derivative of motion regressor file

    """

    i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc.1D'))    
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'movement_regressor-mc-deriv.1D'))

    run_length = []
    for run_number in run_list:
        with open(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run_number), '3dreglog.raw'))) as f:
            run_length.append(str(sum(1 for _ in f) - 1))
        
    #    generate motion derivative regessors if they don't exist
    if not os.path.isfile(o_file):
        utils.prettyOut(task.lower() + " : movement_regressor-mc-deriv.1D")
        utils.Cmd("1d_tool.py"
                " -infile " + i_file +
                " -set_run_lengths " +  ' '.join(run_length) +
                " -derivative"
                " -write " + o_file).run()
                

        
        
def estimate_field_distortion(subj, task):
    """
    Estimate distortion from fmaps

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    task : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # locate fmaps for task. if there is no fmap set specifically for task, use 'fmri'
    # approach is to give each task its own map set even if it is just duplication.
    # perhaps a better approach would be to map these in the visit settings or the BIDS sidecar attributed intended for such
    #
    
    # if there is neither a complete set of task specific or fmri general fmaps, issue a warning then skip distoration estimate
    ap_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-'+task+'_dir-ap_epi.nii.gz') )
    if not os.path.isfile(ap_file):
        ap_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-fmri_dir-ap_epi.nii.gz') )
        
    pa_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-'+task+'_dir-pa_epi.nii.gz') )
    if not os.path.isfile(pa_file):
        pa_file = os.path.join(subj.bids_path, 'fmap', subj.bids_file_join('acq-fmri_dir-pa_epi.nii.gz') )
                
    appa_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'fmap_appa.nii.gz'))
    topup_prefix = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'topup_results_appa'))
    topup_field = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'topup_field_appa'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'fmap_appa-dc.nii.gz'))
    
    # run topup
    if not os.path.isfile(topup_field + '.nii.gz'):
        
        utils.prettyOut("estimating field distortions for "+ task)
        
        # concatenate fmaps
        if not os.path.isfile(appa_file):
            utils.Cmd("fslmerge -t {} {} {}".format(appa_file, ap_file, pa_file)).run()

        # estimate distortion
        utils.Cmd("topup --verbose --imain=" + appa_file + " --datain=" + subj.study.paramsPath + "/epi_params.txt --config=" + subj.study.paramsPath+ "/b02b0.cnf --out=" + topup_prefix + " --fout=" + topup_field + " --iout=" + o_file).run()

    # clean up
    if os.path.isfile(o_file):
        utils.purge_paths(appa_file, o_file)





        
def apply_distortion_correction(subj, task, i_file, o_file, param_indx=1):
    """
    apply distortion correction to correct distortions for epi task    
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    i_file : TYPE
        DESCRIPTION.
    o_file : TYPE
        DESCRIPTION.
    param_indx : TYPE, optional
        DESCRIPTION. The default is 1.
        1,2 = AP; 3,4 = PA (in current params file)
        the index here is used to specify the acquisition encoding of the input and match it to the epi_params file.
        For our current setup, it will be '1' to specify that the entire thing was acquired A -> P
        
    Returns
    -------
    None.

    """
    utils.prettyOut("Applying distortion correction")
    if not os.path.isfile(o_file):
        
        topup_field = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'topup_field_appa'))
        
        if not os.path.isfile(o_file):
            utils.Cmd("applytopup " + 
                      " --verbose " +
                      " --imain=" + i_file +
                      " --inindex=" + str(param_indx) + 
                      " --topup=" + topup_field +
                      " --datain=" + subj.study.paramsPath  + "/epi_params.txt"
                      " --method=jac " +
                      " --interp=spline" + 
                      " --datatype=short" + 
                      " --out=" + o_file).run()
    
            # remove negative values around the brain as a result of jac interpolation from dewarped output
            # make sure datatype INT
            utils.Cmd("fslmaths {0} -abs {0} -odt short".format(o_file)).run()        



    
    
def apply_distortion_correction_sbref(subj, task):
    
    sbref_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc.nii.gz'))
    apply_distortion_correction(subj, task, i_file=sbref_file, o_file=o_file)




def apply_distortion_correction_concat(subj, task):
    
    i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-uncorrected.nii.gz'))
    o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))
    apply_distortion_correction(subj, task, i_file=i_file, o_file=o_file)



    
def apply_distortion_correction_run(subj, task, run, iproc=0, ignore=1):
    
    """
    VZ distortion correction for a single run
    also removes initial volumes specified by ignore
    iproc: TYPE, optional
        DESCRIPTION: 0 - use 3dreg for a given run
                     1 - use pc-3dreg for a given run
    """
    
    if iproc > 0:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-3dreg.nii.gz'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'pc-preproc-dc.nii.gz'))
    else:
        i_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), '3dreg.nii.gz'))
        o_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-' + str(run), 'preproc-dc.nii.gz'))        
      
    temp = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('temp.nii.gz'))
    if os.path.isfile(temp):
        utils.purge_paths(temp)
    
    if not os.path.isfile(o_file):
        apply_distortion_correction(subj, task, i_file=i_file, o_file=temp)
        utils.Cmd("3dTcat -overwrite -prefix " + o_file + " " + temp + "'[" + str(ignore) + "-$]'").run()
        utils.purge_paths(temp)



    
def generate_task_mask(subj, task):
    
    """
    create mask on base(SBREF) image using dist corrected file

    """

    try:
        sbref_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc.nii.gz'))
        mask_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz'))
        
        if not os.path.isfile(mask_file):
            utils.prettyOut(task + " : Make mask")
            utils.Cmd("3dAutomask -dilate 2 -erode 1 -prefix " + mask_file + " " + sbref_file).run()
        
    except:
        raise
        
    # apply this mask to sbref
    apply_task_mask(subj, task, sbref_file, os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz')))
    
    
    
    
    
def apply_task_mask(subj, task, i_file, o_file):
    """ 
    Apply task mask to i_file naming it o_file
    """
    if not os.path.isfile(i_file):
        raise ValueError("Input file not found")
        
    if not os.path.isfile(o_file):
        mask_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz'))
        utils.Cmd( "fslmaths {} -mas {} {}".format(i_file, mask_file, o_file)).run()
            
            
            
            
    
def align_to_T1_BBR(subj, task, template_key='mni_icbm152_nlin_asym_2009c'):
    """
    Boundary-Based Registration works better than ANTs to MNI template
    
    """
    
    anat_path = subj.derivatives_path(bidsType='anat')
    sbref_masked_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
    t1_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
    t1_brain_path = os.path.join(anat_path, subj.bids_file_join('T1w_brain.nii.gz'))
    t1_wmseg_path = os.path.join(anat_path, subj.bids_file_join('spm-wm_prob.nii.gz'))
    tx_for_pref = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'to-T1.aff'))
    tx_inv_pref = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'from-T1.aff'))
    
    try:
        
        if not os.path.isfile(tx_for_pref + '.mat'):   
            
            try:
                # do a standard flirt pre-alignment
                utils.Cmd("flirt -ref {} -in {} -dof 12 -omat {}".format(t1_brain_path, sbref_masked_file, os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'init.mat')))).run()
        
                # now run the bbr
                utils.Cmd("flirt -ref {} -in {} -dof 12 -cost bbr -wmseg {} -init {} -omat {} -schedule {}/etc/flirtsch/bbr.sch".format(
                    t1_brain_path, sbref_masked_file, t1_wmseg_path,
                    os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'init.mat')),
                    tx_for_pref + '.mat',
                    utils.Env.fsl_path
                    )).run()
            
            except:
                raise
                
            finally:
                # remove init mat
                utils.purge_paths(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'init.mat')))
        
    
    
        # also save out inverse transform
        if not os.path.isfile(tx_inv_pref + '.mat'):
            utils.Cmd("convert_xfm -omat {} -inverse {}".format(tx_inv_pref + '.mat', tx_for_pref + '.mat')).run()
            
        
        # convert FSL transform to ITK format for ANTs
        if not os.path.isfile(tx_for_pref + '.mat.itk.txt'):
            utils.Cmd(utils.Env.c3d_path + "/c3d_affine_tool -ref {} -src {} {} -fsl2ras -oitk {}".format(t1_path, sbref_masked_file, tx_for_pref + '.mat', tx_for_pref + '.mat.itk.txt' )).run()
          
        
        # collapse the transformations to a final displacement field
        composite_warp_path = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'to-'+template_key+'_Composite.warp.nii.gz'))
        composite_in_path = os.path.join(anat_path, subj.bids_file_join('T1w-to-'+template_key+'_Composite.h5'))
        if not os.path.isfile(composite_warp_path):
            utils.Cmd("antsApplyTransforms -d 3 -o '[{},1]' -t {} -t {} -r {}".format(composite_warp_path, composite_in_path, tx_for_pref + '.mat.itk.txt', utils.Env.template_lib[template_key]['brain'])).run()

    
    except:
        raise
        
    finally:
        utils.purge_paths(tx_for_pref + '_fast_*', tx_for_pref + '_init*',tx_for_pref + '.nii.gz')
    
    
    



def align_to_T1_AFNI(subj, task, cost_func="lpc", move='giant', cmass="-cmass cmass", opts="", purge=False):

    """
    Performs affine registration and outputs transform dwi2tlrc.aff12.1D

    Expects masked b0 image will already be available.

    :param subj: utils.Subj
    :param study: utils.Study
    :param cost_func: valid cost functions for 3dAllineate
    :param move:
    :param cmass:
    :param opts: Any other valid 3dAllineate flags to pass along
    :param purge: By default this function will not overwrite transforms as those may have been checked. Use this if you don't want to manually delete the previous transform.

    """
    anat_path = subj.derivatives_path(bidsType='anat')
    func_path = subj.derivatives_path(bidsType='func')
    tx_o_path = os.path.join(func_path, subj.bids_file_join('task-' + task.lower() + '-to-T1.aff12.1D'))
    sbref_masked_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
    
    
    # if alignment matrix exists and not told to blow it away
    if purge or not os.path.isfile(tx_o_path):
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : "+task.lower()+" : Spatial Normalization" )

        # just so user understands
        if purge:
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " :Purging previous transform matrix")
            utils.purge_paths(tx_o_path)    


        #    create alignment workspace
        workspace_path = os.path.join(func_path, "regTemp")
        utils.make_path(workspace_path)

        # check for inputs
        if not os.path.isfile(sbref_masked_file):
            utils.mask_data(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc.nii.gz')),
                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz')),
                            sbref_masked_file)


        #    go in there cuz afni script is too dumb to use full paths correctly
        current_path = os.getcwd()
        os.chdir(workspace_path)

        #    do alignment, sadly afni script not made callable from python
        if (move) and (move in ['big', 'giant', 'ginormous']):
             move_str = "-"+move+"_move"
        else:
             move_str = ""

        utils.Cmd("align_epi_anat.py \
                -epi2anat \
                -anat " + os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.nii.gz')) + " \
                -epi " + sbref_masked_file + " \
                -epi_base 0 \
                -epi_strip None \
                -volreg off \
                -deoblique on \
                -anat_has_skull no \
                -tshift off \
                -align_centers yes \
                "+ move_str +" \
                -cost "+ cost_func +" \
                -Allineate_opts '-source_automask+4'").run()

        
        import shutil        
        shutil.copy( subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked_al_mat.aff12.1D'), tx_o_path)

        # clean up
        os.chdir(current_path)
        shutil.rmtree(workspace_path, ignore_errors=True)


    return tx_o_path






def spatial_normalization(subj, task, space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
    """
    
    Creates a final Nonlinear tramsform from functional to template.
    
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    try:
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : "+task.lower()+" : Spatial Normalization" )
        
        
        if space == 'TLRC':
            
            # make sure there is an affine transform
            # if target space is TLRC, use AFNI method for alignment to T1
            # if MNI, use BBR method            
            affine_tx = align_to_T1_AFNI(subj, task)
            
            anat_path = subj.derivatives_path(bidsType='anat')        
            o_path = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz'))
    
            if not os.path.isfile(o_path):
                utils.Cmd("3dNwarpCat -warp1 {} -warp2 {} -warp3 {} -prefix {} ".format(os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.WARP.nii.gz')),
                                                                                        os.path.join(anat_path, subj.bids_file_join('T1w-to-'+space+'.aff12.1D')),
                                                                                        affine_tx,
                                                                                        o_path)).run()
                
        else:
            affine_tx = align_to_T1_BBR(subj, task, template_key=template_key)
            

        
    except:
        raise








def apply_spatial_normalization(
        subj, task, i_file, o_path=None, o_file=None, master_path=None, interp='Linear', regrid=None, fwhm=None,
        space='MNI', template_key='mni_icbm152_nlin_asym_2009c'):
    """
    Spatially normalize i_file into o_path with option flags
    o_file will be i_file with appended flags options
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    :param template_key:
    :param space:
    :param fwhm:
    :param regrid:
    :param interp:
    :param master_path:
    :param subj:
    :param o_file:
    :param o_path:
    :param i_file:
    :param task:

    """
    try:

        ####################################
        # MNI SPACE
        ####################################
        if space == 'MNI':  
            
            if not os.path.isfile(i_file):
                raise ValueError("Input file does not exist")
                
            # split input into path and file
            (i_path, i_file_name) = os.path.split(i_file)
    
            # reduce input file to root
            i_file_root = i_file_name.replace('.nii', '').replace('.gz', '')
            
            # use i_path as o_path if user did not specify
            if o_path is None:
                if o_file is None:
                    o_path = i_path
                    utils.make_path(o_path)
                else:
                    (o_path, o_file_name) = os.path.split(o_file)
                    utils.make_path(o_path)
                    
                    
            # build output filename
            if o_file is None:
                o_file_name = i_file_root
                if regrid is not None:
                    o_file_name = o_file_name + '_{0}x{0}x{0}'.format(regrid)
        
                if fwhm is not None:
                    o_file_name = o_file_name + '_GB{0}'.format(fwhm)
                    
                o_file_name = o_file_name + '_'+space+'.nii.gz'
            
            else:
                (o_path, o_file_name) = os.path.split(o_file)
                
            
            if master_path is None:
                master_path = utils.Env.template_lib[template_key]['brain']                 

            #
            utils.Cmd(f"antsApplyTransforms --verbose --float -d 3 --input-image-type 3 "
                      f"-o {os.path.join(o_path, o_file_name)} "
                      f"-t {os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'to-'+template_key+'_Composite.warp.nii.gz'))} "
                      f"-r {master_path} -i {os.path.join(i_path, i_file_name)} "
                      f"--output-data-type short --interpolation {interp}").run()
            
            # update NIFTI header for template space correction for afni viewer
            from mayerlab.preprocessing.anat import afni
            afni.set_header_space(space='MNI', f_path=os.path.join(o_path, o_file_name))

            if regrid is not None:
                utils.Cmd(f"3dresample -overwrite -dxyz {str(regrid)} {str(regrid)} {str(regrid)} -input {os.path.join(o_path, o_file_name)} -prefix {os.path.join(o_path, o_file_name)}").run()

            if fwhm is not None:
                utils.Cmd(f"3dmerge -overwrite -1blur_fwhm {str(fwhm)} -doall -prefix {os.path.join(o_path, o_file_name)} {os.path.join(o_path, o_file_name)}").run()




        ####################################
        # TLRC SPACE
        ####################################
        if space == 'TLRC':
            
            if not os.path.isfile(i_file):
                raise ValueError("Input file does not exist")
                
            # split input into path and file
            (i_path, i_file_name) = os.path.split(i_file)
    
            # reduce input file to root
            i_file_root = i_file_name.replace('.nii', '').replace('.gz','')        
            
            # use i_path as o_path if user did not specify
            if o_path == None:
                if o_file == None:
                    o_path = i_path
                    utils.make_path(o_path)
                else:
                    (o_path, o_file_name) = os.path.split(o_file)
                    utils.make_path(o_path)
                
                
            # build output filename
            if o_file == None:
                o_file_name = i_file_root
                if regrid != None:
                    o_file_name = o_file_name + '_{0}x{0}x{0}'.format(regrid)
        
                if fwhm != None:
                    o_file_name = o_file_name + '_GB{0}'.format(fwhm)
                    
                o_file_name = o_file_name + '_'+space+'.nii.gz'
            
            else:
                (o_path, o_file_name) = os.path.split(o_file)
                
            
            
            
            anat_path = subj.derivatives_path(bidsType='anat')
            
            if master_path == None:
                if space == 'MNI':
                    master_path = os.path.join(anat_path, subj.bids_file_join('T1w_MNI.nii.gz'))
                else:
                    master_path = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.TLRC.nii.gz'))
            
            
            
            if fwhm == None:
                if regrid == None:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            interp,
                                                                                                            os.path.join(o_path, o_file_name)
                                                                                                            )).run()
                    
                else:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -newgrid {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            str(regrid),
                                                                                                            interp,
                                                                                                            os.path.join(o_path, o_file_name)
                                                                                                            )).run()
                    
                
            else:
                
                if regrid == None:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            interp,
                                                                                                            os.path.join(o_path, 'tmp-'+o_file_name)
                                                                                                            )).run()                
                else:
                    utils.Cmd("3dNwarpApply -source {} -nwarp {} -master {} -newgrid {} -interp {} -prefix {} ".format(os.path.join(i_path, i_file_name),
                                                                                                            os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '-to-'+space+'.NWARP.nii.gz')),
                                                                                                            master_path,
                                                                                                            str(regrid),
                                                                                                            interp,
                                                                                                            os.path.join(o_path, 'tmp-'+o_file_name)
                                                                                                            )).run()
    
                utils.Cmd("3dmerge -1blur_fwhm {} -doall -prefix {} {}".format(fwhm, os.path.join(o_path, o_file), os.path.join(o_path, 'tmp-'+o_file_name) )).run()
                utils.purge_paths(os.path.join(o_path, 'tmp-'+o_file_name))
                
            print("\nInput: " + os.path.join(i_path, i_file_name))    
            print("Output: " + os.path.join(o_path, o_file_name))
                
    except:
        raise








def apply_blur(subj, i_file, fwhm=None, o_path=None):
    """
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    try:

        
        if not os.path.isfile(i_file):
            raise ValueError("Input file does not exist")
            
        # split input into path and file
        (i_path, i_file) = os.path.split(i_file)

        # reduce input file to root
        i_file_root = i_file.replace('.nii', '').replace('.gz','')        
        
        # use i_path as o_path if user did not specify
        if o_path == None:
            o_path = i_path
            utils.make_path(o_path)
            
            
        # build output filename
        o_file = i_file_root
        if fwhm != None:
            o_file = o_file + '_GB{0}'.format(fwhm)
            
        o_file = o_file + '.nii.gz'
        
    
        utils.Cmd("3dmerge -1blur_fwhm {} -doall -prefix {} {}".format(fwhm, os.path.join(o_path, o_file), os.path.join(i_path, i_file) )).echo()
            
        print("\nOutput:" + os.path.join(o_path, o_file))
            
    except:
        raise



def proc_cleanup(subj, task, iproc=0):
    """
    VZ still need to decide which results to keep

    """
    
    if os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
        
        if iproc <= 0:
            # keep dspk and tshft for now
            try:
                utils.purge_paths( os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_2dreg*')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_3dreg.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_3dreglog.raw')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_maxdisp.raw')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_log.mc'))
                              )
            except:
                raise
    
        if iproc > 0:
            # keep pc-dspk and pc-tshft for now
            try:            
                utils.purge_paths( os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_pc-2dreg*')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_pc-3dreg.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_pc-log.mc'))
                              )                              
            except:
                raise




def clean_proc(subj, task):
    """
    
    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_sbref-dc.nii.gz')),
    
    try:
        
        if os.path.isfile(os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, 'preproc-dc.nii.gz'))):
            
            utils.purge_paths( os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_dspk.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_tshft.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_2dreg*')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_3dreg.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_3dreglog.raw')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_log.mc')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_sbref.nii.gz')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_run-?_slice_timing.1D')),
                               os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower() + '_preproc-uncorrected.nii.gz'))
                              )
        
    except:
        raise



def generateResiduals(subj, study, fwhm, run_label, base_suffix="dist_corr"):
    """
    generate residuals for resting state analyzes
    Deprecate?
    rsfc_create_residuals() can generate residuals with or without blurring prior
    The only difference in this function then is that bandpass filters and spatially normalizes the residuals

    """
    #    build subject_path
    subject_path = study.path + "/" + subj.ursi.full

    #    workaround fact that regressor files contain visit
    #    and any other differences between studies that have or have not a visit....
    # subject_visit_path = ""
    if subj.visit:
        # subject_visit_path = "." + subj.visit
        subject_path = subject_path + "/" + subj.visit

    #    run_label indicates where processing data is
    #    run_type indicates final processing outpute
    run_type = run_label.rstrip("0123456789")
    final_blur_fwhm = str(fwhm)

    #    check for final output only before running
    if not os.path.isfile(
            subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3+tlrc.HEAD"):

        #    deposit junk in /tmp
        os.chdir("/tmp")

        #    define inputs to check

        # allow null base suffix
        if len(base_suffix):
            base_suffix = "." + base_suffix
        rest_input_file = subject_path + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc" + base_suffix + "+orig"
        rest_tx_epi2T1_file = subj.get_transform_path(
            study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"
        rest_tx_T12epi_file = subj.get_transform_path(
            study) + study.label + "." + subj.ursi.short + "." + run_type + "_T12epi.aff12.1D"
        rest_movement_file = study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D"
        spm_wm_file = subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm.mask.60+orig"
        spm_csf_file = subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf.mask.60+orig"

        #    check inputs
        if not os.path.isfile(rest_input_file + ".HEAD"):
            raise utils.CmdProcError("\nMissing Input " + rest_input_file + "\n")

        if not os.path.isfile(rest_tx_T12epi_file):

            if not os.path.isfile(rest_tx_epi2T1_file):

                raise utils.CmdProcError("\nMissing Input " + rest_tx_epi2T1_file + "\n")

            else:

                cmd = utils.Cmd("cat_matvec -ONELINE " + rest_tx_epi2T1_file + " -I > " + rest_tx_T12epi_file)
                cmd.run()

        if not os.path.isfile(rest_movement_file):
            raise utils.CmdProcError("\nMissing Input " + rest_movement_file + "\n")

        if not os.path.isfile(spm_wm_file + ".HEAD"):
            raise utils.CmdProcError("\nMissing Input " + spm_wm_file + "\n")

        if not os.path.isfile(spm_csf_file + ".HEAD"):
            raise utils.CmdProcError("\nMissing Input " + spm_csf_file + "\n")

        #    now start creating the tissue regressor files
        utils.prettyOut("Resting State Residuals : Creating the tissue regressor files")
        if not os.path.isfile(
                subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D"):

            if not os.path.isfile(
                    subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg+orig.HEAD"):
                cmd = utils.Cmd("3dAllineate "
                                "-master " + rest_input_file + " "
                                                               "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
                                                                                                          "-input " + spm_wm_file + " "
                                                                                                                                    "-final NN "
                                                                                                                                    "-prefix " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg")

                cmd.run()

            cmd = utils.Cmd("3dROIstats -quiet "
                            "-mask " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg+orig "
                            + rest_input_file + " > " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D")
            cmd.run()

        if not os.path.isfile(
                subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D"):

            if not os.path.isfile(
                    subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg+orig.HEAD"):
                cmd = utils.Cmd("3dAllineate "
                                "-master " + rest_input_file + " "
                                                               "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
                                                                                                          "-input " + spm_csf_file + " "
                                                                                                                                     "-final NN "
                                                                                                                                     "-prefix " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg")
                cmd.run()

            cmd = utils.Cmd("3dROIstats -quiet "
                            "-mask " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg+orig "
                            + rest_input_file + " > " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D")
            cmd.run()

            # cmd = utils.Cmd("rm -f " +subject_path + "/anat/" + study.label + "." + subj.ursi.short+ ".anat.spm.csf." + run_type + "_Reg+orig* " )
            # return_code = subprocess.call(cmd, shell=True)
            # if return_code:
            # raise utils.CmdProcError(return_code)

        #    generate motion derivative regessors if they don't exist
        if not os.path.isfile(
                study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D"):
            cmd = utils.Cmd("1d_tool.py"
                            " -infile " + rest_movement_file +
                            " -set_nruns 1"
                            " -derivative"
                            " -write " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D")

            cmd.run()

        #    blur input to glm
        if not os.path.isfile(
                subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig.HEAD"):
            cmd = utils.Cmd("3dmerge"
                            " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm +
                            " -1blur_fwhm " + final_blur_fwhm +
                            " -doall "
                            + rest_input_file)

            cmd.run()

        #    run glm to generate residuals
        utils.prettyOut("Resting State Residuals : Run glm to generate residuals")
        if not os.path.isfile(
                subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + "+orig.HEAD"):
            cmd = utils.Cmd("3dDeconvolve"
                            " -overwrite"
                            " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig"
                                                                                                                                                                           " -polort 2"
                                                                                                                                                                           " -allzero_OK"
                                                                                                                                                                           " -GOFORIT 10"
                                                                                                                                                                           " -num_stimts 14"
                                                                                                                                                                           " -stim_file 1 " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D"
                                                                                                                                                                                                                                                                           " -stim_label 1 'CSF'"
                                                                                                                                                                                                                                                                           " -stim_minlag 1 0"
                                                                                                                                                                                                                                                                           " -stim_maxlag 1 0"
                                                                                                                                                                                                                                                                           " -stim_file 2 " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D"
                                                                                                                                                                                                                                                                                                                                                                           " -stim_label 2 'WM'"
                                                                                                                                                                                                                                                                                                                                                                           " -stim_minlag 2 0"
                                                                                                                                                                                                                                                                                                                                                                           " -stim_maxlag 2 0"
                                                                                                                                                                                                                                                                                                                                                                           " -stim_file 3 " + rest_movement_file + "'[1]'"
                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_label 3 'ROLL'"
                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_minlag 3 0"
                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_maxlag 3 0"
                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_file 4 " + rest_movement_file + "'[2]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_label 4 'PITCH'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_minlag 4 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_maxlag 4 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_file 5 " + rest_movement_file + "'[3]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_label 5 'YAW'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_minlag 5 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_maxlag 5 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_file 6 " + rest_movement_file + "'[4]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_label 6 'dS'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_minlag 6 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_maxlag 6 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_file 7 " + rest_movement_file + "'[5]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_label 7 'dL'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_minlag 7 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_maxlag 7 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_file 8 " + rest_movement_file + "'[6]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_label 8 'dP'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_minlag 8 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_maxlag 8 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           " -stim_file 9 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[1]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                " -stim_label 9 'ROLL_deriv'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                " -stim_minlag 9 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                " -stim_maxlag 9 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                " -stim_file 10 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[2]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      " -stim_label 10 'PITCH'_deriv"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      " -stim_minlag 10 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      " -stim_maxlag 10 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      " -stim_file 11 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[3]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_label 11 'YAW_deriv'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_minlag 11 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_maxlag 11 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_file 12 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[4]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  " -stim_label 12 'dS_deriv'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  " -stim_minlag 12 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  " -stim_maxlag 12 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  " -stim_file 13 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[5]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        " -stim_label 13 'dL_deriv'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        " -stim_minlag 13 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        " -stim_maxlag 13 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        " -stim_file 14 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[6]'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              " -stim_label 14 'dP_deriv'"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              " -stim_minlag 14 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              " -stim_maxlag 14 0"
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              " -errts " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm +
                            " -jobs 10 "
                            " -xsave"
                            " -float"
                            " -nobucket")

            cmd.run()

        #    bandpass filtering
        utils.prettyOut("Resting State Residuals : bandpass filtering")
        if not os.path.isfile(
                subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig.HEAD"):
            cmd = utils.Cmd(" 3dBandpass"
                            " -band 0.01 0.1"
                            " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + "+orig"
                                                                                                                                                                         " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1")

            cmd.run()

        # convert residual file to tlrc at 3x3x3 resolution
        utils.prettyOut("Resting State Residuals : convert residual file to tlrc at 3x3x3 resolution")
        if not os.path.isfile(
                subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3+tlrc.HEAD"):
            # cmd = utils.Cmd(" 3dAllineate"
            # " -overwrite"
            # " -master " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.sksp+tlrc"
            # " -mast_dxyz 3.0"
            # " -1Dmatrix_apply " + subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2tlrc.aff12.1D"
            # " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig"
            # " -floatize"
            # " -final trilinear"
            # " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3" )

            nl_warp_src = subj.get_transform_path(
                study) + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc"
            t1_to_tlrc_mat = subj.get_transform_path(
                study) + study.label + "." + subj.ursi.short + ".T1_2tlrc.aff3x4.1D"
            epi_2_T1_aff_mat = subj.get_transform_path(
                study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"

            #    the reference volume output will be aligned to
            #    anat.sksp.qwarp+tlrc
            #    which is non-linear warp to TT_N27
            #    so reference volume output should look good over that
            cmd = utils.Cmd("3dNwarpApply " +
                            " -source " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig " +
                            " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
                            " -master " + subj.get_path(
                study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp+tlrc" +
                            " -interp wsinc5" +
                            " -dxyz 3.0" +
                            " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3")

            cmd.run()

    #    clean up even if nothing else was run
    utils.purge_paths(
        subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig*",
        subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.*." + run_type + "_Reg+orig* ")


def create_tissue_regressor(subj, study, run_label, tissue_classes=['wm', 'csf'], base_suffix="dist_corr"):
    """
    Generate tissue regressor files for any run type.
    Regressor file outputs into run_label directory

    """

    #    run_label indicates where processing data is
    #    run_type indicates final processing output
    run_type = run_label.rstrip("0123456789")

    #    define inputs to check

    # allow null base suffix
    if len(base_suffix):
        base_suffix = "." + base_suffix

    rest_input_file = subj.get_path(
        study) + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc" + base_suffix + "+orig.HEAD"
    rest_tx_T12epi_file = subj.get_transform_path(
        study) + study.label + "." + subj.ursi.short + "." + run_type + "_T12epi.aff12.1D"
    rest_tx_epi2T1_file = subj.get_transform_path(
        study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"

    for tissue_class in tissue_classes:

        spm_tissue_file = os.path.join(subj.get_path(study, 'anat'),
                                       subj.id_fn(study, fn="anat.spm." + tissue_class + ".mask.60+orig.HEAD"))
        if not os.path.isfile(spm_tissue_file):
            raise utils.CmdProcError("\nMissing Input " + spm_tissue_file + "\n")

            #    check preproc input
        if not os.path.isfile(rest_input_file):
            raise utils.CmdProcError("\nMissing Input " + rest_input_file + "\n")

        #   check that transform exists
        if not os.path.isfile(rest_tx_T12epi_file):
            if not os.path.isfile(rest_tx_epi2T1_file):
                raise utils.CmdProcError("\nMissing Input " + rest_tx_epi2T1_file + "\n")
            else:
                utils.Cmd("cat_matvec -ONELINE " + rest_tx_epi2T1_file + " -I > " + rest_tx_T12epi_file).run()

        #    now start creating the tissue regressor files
        o_file_1d = os.path.join(subj.get_path(study), run_label,
                                 subj.id_fn(study, fn="Tissue.Regressor." + tissue_class + ".1D"))
        o_file_reg = os.path.join(subj.get_path(study, 'anat'),
                                  subj.id_fn(study, fn="anat.spm." + tissue_class + "." + run_type + "_Reg+orig.HEAD"))

        if not os.path.isfile(o_file_1d):
            if not os.path.isfile(o_file_reg):
                utils.prettyOut("RSFC : Creating the tissue regressor files")
                utils.Cmd("3dAllineate "
                          "-master " + rest_input_file + " "
                                                         "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
                                                                                                    "-input " + spm_tissue_file + " "
                                                                                                                                  "-final NN "
                                                                                                                                  "-prefix " + o_file_reg).run()

            utils.Cmd("3dROIstats -quiet -mask " + o_file_reg + " " + rest_input_file + " > " + o_file_1d).run()

            utils.purge_paths(o_file_reg)



# def generateResiduals(subj, study, fwhm, run_label, base_suffix="dist_corr"):
#     """
#     generate residuals for resting state analyzes
#     Deprecate?
#     rsfc_create_residuals() can generate residuals with or without blurring prior
#     The only difference in this function then is that bandpass filters and spatially normalizes the residuals
#
#     """
#     #    build subject_path
#     subject_path = study.path + "/" + subj.ursi.full
#
#     #    workaround fact that regressor files contain visit
#     #    and any other differences between studies that have or have not a visit....
#     # subject_visit_path = ""
#     if subj.visit:
#         # subject_visit_path = "." + subj.visit
#         subject_path = subject_path + "/" + subj.visit
#
#     #    run_label indicates where processing data is
#     #    run_type indicates final processing outpute
#     run_type = run_label.rstrip("0123456789")
#     final_blur_fwhm = str(fwhm)
#
#     #    check for final output only before running
#     if not os.path.isfile(
#             subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3+tlrc.HEAD"):
#
#         #    deposit junk in /tmp
#         os.chdir("/tmp")
#
#         #    define inputs to check
#
#         # allow null base suffix
#         if len(base_suffix):
#             base_suffix = "." + base_suffix
#         rest_input_file = subject_path + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc" + base_suffix + "+orig"
#         rest_tx_epi2T1_file = subj.get_transform_path(
#             study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"
#         rest_tx_T12epi_file = subj.get_transform_path(
#             study) + study.label + "." + subj.ursi.short + "." + run_type + "_T12epi.aff12.1D"
#         rest_movement_file = study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D"
#         spm_wm_file = subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm.mask.60+orig"
#         spm_csf_file = subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf.mask.60+orig"
#
#         #    check inputs
#         if not os.path.isfile(rest_input_file + ".HEAD"):
#             raise utils.CmdProcError("\nMissing Input " + rest_input_file + "\n")
#
#         if not os.path.isfile(rest_tx_T12epi_file):
#
#             if not os.path.isfile(rest_tx_epi2T1_file):
#
#                 raise utils.CmdProcError("\nMissing Input " + rest_tx_epi2T1_file + "\n")
#
#             else:
#
#                 cmd = utils.Cmd("cat_matvec -ONELINE " + rest_tx_epi2T1_file + " -I > " + rest_tx_T12epi_file)
#                 cmd.run()
#
#         if not os.path.isfile(rest_movement_file):
#             raise utils.CmdProcError("\nMissing Input " + rest_movement_file + "\n")
#
#         if not os.path.isfile(spm_wm_file + ".HEAD"):
#             raise utils.CmdProcError("\nMissing Input " + spm_wm_file + "\n")
#
#         if not os.path.isfile(spm_csf_file + ".HEAD"):
#             raise utils.CmdProcError("\nMissing Input " + spm_csf_file + "\n")
#
#         #    now start creating the tissue regressor files
#         utils.prettyOut("Resting State Residuals : Creating the tissue regressor files")
#         if not os.path.isfile(
#                 subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D"):
#
#             if not os.path.isfile(
#                     subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg+orig.HEAD"):
#                 cmd = utils.Cmd("3dAllineate "
#                                 "-master " + rest_input_file + " "
#                                                                "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
#                                                                                                           "-input " + spm_wm_file + " "
#                                                                                                                                     "-final NN "
#                                                                                                                                     "-prefix " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg")
#
#                 cmd.run()
#
#             cmd = utils.Cmd("3dROIstats -quiet "
#                             "-mask " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg+orig "
#                             + rest_input_file + " > " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D")
#             cmd.run()
#
#         if not os.path.isfile(
#                 subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D"):
#
#             if not os.path.isfile(
#                     subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg+orig.HEAD"):
#                 cmd = utils.Cmd("3dAllineate "
#                                 "-master " + rest_input_file + " "
#                                                                "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
#                                                                                                           "-input " + spm_csf_file + " "
#                                                                                                                                      "-final NN "
#                                                                                                                                      "-prefix " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg")
#                 cmd.run()
#
#             cmd = utils.Cmd("3dROIstats -quiet "
#                             "-mask " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg+orig "
#                             + rest_input_file + " > " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D")
#             cmd.run()
#
#             # cmd = utils.Cmd("rm -f " +subject_path + "/anat/" + study.label + "." + subj.ursi.short+ ".anat.spm.csf." + run_type + "_Reg+orig* " )
#             # return_code = subprocess.call(cmd, shell=True)
#             # if return_code:
#             # raise utils.CmdProcError(return_code)
#
#         #    generate motion derivative regessors if they don't exist
#         if not os.path.isfile(
#                 study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D"):
#             cmd = utils.Cmd("1d_tool.py"
#                             " -infile " + rest_movement_file +
#                             " -set_nruns 1"
#                             " -derivative"
#                             " -write " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D")
#
#             cmd.run()
#
#         #    blur input to glm
#         if not os.path.isfile(
#                 subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig.HEAD"):
#             cmd = utils.Cmd("3dmerge"
#                             " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm +
#                             " -1blur_fwhm " + final_blur_fwhm +
#                             " -doall "
#                             + rest_input_file)
#
#             cmd.run()
#
#         #    run glm to generate residuals
#         utils.prettyOut("Resting State Residuals : Run glm to generate residuals")
#         if not os.path.isfile(
#                 subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + "+orig.HEAD"):
#             cmd = utils.Cmd("3dDeconvolve"
#                             " -overwrite"
#                             " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig"
#                                                                                                                                                                            " -polort 2"
#                                                                                                                                                                            " -allzero_OK"
#                                                                                                                                                                            " -GOFORIT 10"
#                                                                                                                                                                            " -num_stimts 14"
#                                                                                                                                                                            " -stim_file 1 " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D"
#                                                                                                                                                                                                                                                                            " -stim_label 1 'CSF'"
#                                                                                                                                                                                                                                                                            " -stim_minlag 1 0"
#                                                                                                                                                                                                                                                                            " -stim_maxlag 1 0"
#                                                                                                                                                                                                                                                                            " -stim_file 2 " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D"
#                                                                                                                                                                                                                                                                                                                                                                            " -stim_label 2 'WM'"
#                                                                                                                                                                                                                                                                                                                                                                            " -stim_minlag 2 0"
#                                                                                                                                                                                                                                                                                                                                                                            " -stim_maxlag 2 0"
#                                                                                                                                                                                                                                                                                                                                                                            " -stim_file 3 " + rest_movement_file + "'[1]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_label 3 'ROLL'"
#                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_minlag 3 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_maxlag 3 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_file 4 " + rest_movement_file + "'[2]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_label 4 'PITCH'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_minlag 4 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_maxlag 4 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_file 5 " + rest_movement_file + "'[3]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_label 5 'YAW'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_minlag 5 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_maxlag 5 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_file 6 " + rest_movement_file + "'[4]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_label 6 'dS'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_minlag 6 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_maxlag 6 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_file 7 " + rest_movement_file + "'[5]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_label 7 'dL'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_minlag 7 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_maxlag 7 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " -stim_file 8 " + rest_movement_file + "'[6]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_label 8 'dP'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_minlag 8 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_maxlag 8 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            " -stim_file 9 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[1]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 " -stim_label 9 'ROLL_deriv'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 " -stim_minlag 9 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 " -stim_maxlag 9 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 " -stim_file 10 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[2]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       " -stim_label 10 'PITCH'_deriv"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       " -stim_minlag 10 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       " -stim_maxlag 10 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       " -stim_file 11 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[3]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             " -stim_label 11 'YAW_deriv'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             " -stim_minlag 11 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             " -stim_maxlag 11 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             " -stim_file 12 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[4]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_label 12 'dS_deriv'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_minlag 12 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_maxlag 12 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   " -stim_file 13 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[5]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         " -stim_label 13 'dL_deriv'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         " -stim_minlag 13 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         " -stim_maxlag 13 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         " -stim_file 14 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[6]'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               " -stim_label 14 'dP_deriv'"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               " -stim_minlag 14 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               " -stim_maxlag 14 0"
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               " -errts " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm +
#                             " -jobs 10 "
#                             " -xsave"
#                             " -float"
#                             " -nobucket")
#
#             cmd.run()
#
#         #    bandpass filtering
#         utils.prettyOut("Resting State Residuals : bandpass filtering")
#         if not os.path.isfile(
#                 subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig.HEAD"):
#             cmd = utils.Cmd(" 3dBandpass"
#                             " -band 0.01 0.1"
#                             " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + "+orig"
#                                                                                                                                                                          " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1")
#
#             cmd.run()
#
#         # convert residual file to tlrc at 3x3x3 resolution
#         utils.prettyOut("Resting State Residuals : convert residual file to tlrc at 3x3x3 resolution")
#         if not os.path.isfile(
#                 subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3+tlrc.HEAD"):
#             # cmd = utils.Cmd(" 3dAllineate"
#             # " -overwrite"
#             # " -master " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.sksp+tlrc"
#             # " -mast_dxyz 3.0"
#             # " -1Dmatrix_apply " + subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2tlrc.aff12.1D"
#             # " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig"
#             # " -floatize"
#             # " -final trilinear"
#             # " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3" )
#
#             nl_warp_src = subj.get_transform_path(
#                 study) + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc"
#             t1_to_tlrc_mat = subj.get_transform_path(
#                 study) + study.label + "." + subj.ursi.short + ".T1_2tlrc.aff3x4.1D"
#             epi_2_T1_aff_mat = subj.get_transform_path(
#                 study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"
#
#             #    the reference volume output will be aligned to
#             #    anat.sksp.qwarp+tlrc
#             #    which is non-linear warp to TT_N27
#             #    so reference volume output should look good over that
#             cmd = utils.Cmd("3dNwarpApply " +
#                             " -source " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig " +
#                             " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
#                             " -master " + subj.get_path(
#                 study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp+tlrc" +
#                             " -interp wsinc5" +
#                             " -dxyz 3.0" +
#                             " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3")
#
#             cmd.run()
#
#     #    clean up even if nothing else was run
#     utils.purge_paths(
#         subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig*",
#         subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.*." + run_type + "_Reg+orig* ")


# def create_tissue_regressor(subj, study, run_label, tissue_classes=['wm', 'csf'], base_suffix="dist_corr"):
#     """
#     Generate tissue regressor files for any run type.
#     Regressor file outputs into run_label directory
#
#     """
#
#     #    run_label indicates where processing data is
#     #    run_type indicates final processing output
#     run_type = run_label.rstrip("0123456789")
#
#     #    define inputs to check
#
#     # allow null base suffix
#     if len(base_suffix):
#         base_suffix = "." + base_suffix
#
#     rest_input_file = subj.get_path(
#         study) + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc" + base_suffix + "+orig.HEAD"
#     rest_tx_T12epi_file = subj.get_transform_path(
#         study) + study.label + "." + subj.ursi.short + "." + run_type + "_T12epi.aff12.1D"
#     rest_tx_epi2T1_file = subj.get_transform_path(
#         study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"
#
#     for tissue_class in tissue_classes:
#
#         spm_tissue_file = os.path.join(subj.get_path(study, 'anat'),
#                                        subj.id_fn(study, fn="anat.spm." + tissue_class + ".mask.60+orig.HEAD"))
#         if not os.path.isfile(spm_tissue_file):
#             raise utils.CmdProcError("\nMissing Input " + spm_tissue_file + "\n")
#
#             #    check preproc input
#         if not os.path.isfile(rest_input_file):
#             raise utils.CmdProcError("\nMissing Input " + rest_input_file + "\n")
#
#         #   check that transform exists
#         if not os.path.isfile(rest_tx_T12epi_file):
#             if not os.path.isfile(rest_tx_epi2T1_file):
#                 raise utils.CmdProcError("\nMissing Input " + rest_tx_epi2T1_file + "\n")
#             else:
#                 utils.Cmd("cat_matvec -ONELINE " + rest_tx_epi2T1_file + " -I > " + rest_tx_T12epi_file).run()
#
#         #    now start creating the tissue regressor files
#         o_file_1d = os.path.join(subj.get_path(study), run_label,
#                                  subj.id_fn(study, fn="Tissue.Regressor." + tissue_class + ".1D"))
#         o_file_reg = os.path.join(subj.get_path(study, 'anat'),
#                                   subj.id_fn(study, fn="anat.spm." + tissue_class + "." + run_type + "_Reg+orig.HEAD"))
#
#         if not os.path.isfile(o_file_1d):
#             if not os.path.isfile(o_file_reg):
#                 utils.prettyOut("RSFC : Creating the tissue regressor files")
#                 utils.Cmd("3dAllineate "
#                           "-master " + rest_input_file + " "
#                                                          "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
#                                                                                                     "-input " + spm_tissue_file + " "
#                                                                                                                                   "-final NN "
#                                                                                                                                   "-prefix " + o_file_reg).run()
#
#             utils.Cmd("3dROIstats -quiet -mask " + o_file_reg + " " + rest_input_file + " > " + o_file_1d).run()
#
#             utils.purge_paths(o_file_reg)


class qc:
    
    def run_mriqc(subj):
        """
        Run MRIQC for each task, each run as needed
        """
        
        proc_map = subj.get_subj_task_map()

        for task, run_number in proc_map.items():
            for run in range(1, run_number + 1):

                if not os.path.isfile(os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'func', subj.bids_file_join('task-'+task.lower()+'_run-'+str(run)+'_bold.json'))):
        
                    # create workspace path outside of mriqc for correction group permissions
                    os.umask(0o007)
                    utils.make_path(os.path.join(subj.study.mriqcPath, 'workspace', 'func'))

                    # build path outside of container for correct permissions?
                    workdir = os.path.join(subj.study.mriqcPath, 'workspace', 'func', subj.bids_file_join("task-"+task.lower()+"_run-"+str(run)+"_bold"))
                    utils.make_path(workdir)

                    # must have this env variable set to have access to templates; does not work to include as --env 
                    os.environ['SINGULARITYENV_TEMPLATEFLOW_HOME'] = "/home/bidsapp/.cache/templateflow"
                    utils.Cmd("singularity run --cleanenv --containall --disable-cache --writable-tmpfs --bind /export:/export /export/research/analysis/human/amayer/shared/apps/containers/runtime/mayerlab_mriqc_0.16.1.sif "
                              + " --participant_label " + subj.ursi.short 
                              + " --session-id " + str(subj.visit)
                              + " --modalities bold"
                              + " --task-id " + task.lower()
                              + " --run-id " + str(run)
                              + " --work-dir " + workdir
                              + " --no-sub "
                              + " --verbose "
                              + " --nprocs 6 "
                              + " --verbose-reports "
                              + " --float32 "
                              + " --ants-nthreads 6 "
                              + " --omp-nthreads 6 "
                              + subj.study.bidsPath + " "
                              + subj.study.mriqcPath + " participant "
                              ).run()
                    
                    # have to clean up after poldrack
                    utils.purge_paths(subj.study.mriqcPath + "/workspace/func/" + subj.bids_file_join("task-"+task.lower()+"_run-"+str(run)+"_bold"))

            try:
                # fix permissions on container task output for group access
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id) )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit) )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'func') )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'func', '*.json') )
                utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_file_join('task-'+task.lower() + '*.html')) )

            except:
                pass




    def snapshot_mask(subj, task):
        
        """
        
        """
            
        if not os.path.isfile(os.path.join(subj.derivatives_path('func'), 'qc', subj.bids_file_join('task-' + task.lower(), 'mask.axi.png'))):
            
            ulay_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task, "preproc-dc.nii.gz'[0]'")) 
            pref = os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), 'mask'))
            label = subj.bids_file_join('task-' + task.lower(), 'mask')
            annot = " mask "
                
            olay_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'mask.nii.gz'))
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC Mask Snapshots" )
            
            utils.Cmd("@chauffeur_afni \
                    -ulay "+ulay_file+" \
                    -olay "+olay_file+" \
                    -prefix "+ pref +" \
                    -set_xhairs OFF \
                    -montx 10 \
                    -monty 3 \
                    -label_string "+ label +" \
                    -pbar_posonly \
                    -cbar_ncolors 1").run()
                    
                    
            for view in ['axi', 'cor', 'sag']:
                utils.Cmd("convert " + pref + '.'+view+'.png' +" -gravity North -background YellowGreen -splice 0x18 -pointsize 18 -annotate +0+2 '"+task+annot+view+"'  " + pref + '.'+view+'.png' ).run()
                    
         
                
                
                
    def snapshot_spatial_normalization(subj, task, space='MNI'):

        if not os.path.isfile(os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('func', 'task-' + task.lower(), space+'.jpg'))):

            src_file = os.path.join(subj.derivatives_path('func'), subj.bids_file_join('task-' + task.lower(), 'run-1', 'sbref-dc_masked.nii.gz'))
            spa_norm_lay = os.path.join(subj.derivatives_path('func'), 'qc', 'sbref.masked.'+space+'.nii.gz')
            apply_spatial_normalization(subj, 
                                        task, 
                                        src_file, 
                                        o_file=spa_norm_lay,
                                        space=space)
    
            edge_first_image = os.path.join(subj.derivatives_path(bidsType='func'), 'qc', 'sbref.masked.'+space+'.edge.nii.gz')
            utils.Cmd("3dedge3 -overwrite  -input "+ spa_norm_lay +" -prefix " + edge_first_image).run()
            
            if space == 'MNI':
                ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/AFNI_MNI152_T1_2009c+tlrc.HEAD"
            else:
                ulay_file = "/export/research/analysis/human/amayer/shared/apps/brains/TT_N27+tlrc.HEAD"
                
                
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC "+space+" Snapshots" )
            
            utils.Cmd("@chauffeur_afni \
                    -ulay "+ulay_file+" \
                    -olay "+edge_first_image+" \
                    -prefix "+ os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space)) +" \
                    -label_string "+ subj.bids_file_join('task-' + task.lower(), space) + " \
                    -set_xhairs OFF \
                    -delta_slices 4 4 4 \
                    -montx 6 \
                    -monty 3 \
                    -opacity 9 \
                    -pbar_posonly \
                    -cbar_ncolors 6 \
                    -cbar_topval '' \
                    -cbar '1000=blue \
                            800=cyan \
                            600=rbgyr20_10 \
                            400=rbgyr20_08 \
                            200=rbgyr20_05 \
                            100=none \
                            0=none'").run()
    
    
            for view in ['axi', 'cor', 'sag']:
                utils.Cmd("convert " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space+'.' +view+'.png')) +" -gravity North -background YellowGreen -splice 0x18 -pointsize 18 -annotate +0+2 '"+task+" "+space+" "+view+"'  " + os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('task-' + task.lower(), space+'.' +view+'.png')) ).run()
    
    
            #utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, spa_norm_lay, os.path.join(subj.derivatives_path(bidsType='func'), 'qc', subj.bids_file_join('func', 'task-' + task.lower(), space)))).run()
            
            utils.purge_paths(spa_norm_lay,edge_first_image)
        
        
    def aggregate_qc_pdf(subj):
        
        
        """
        
        """    
        import glob
        
        
        func_path = subj.derivatives_path(bidsType='func')
        
        o_file = os.path.join(func_path, 'qc', subj.bids_file_join('func-QC.pdf'))
        if not os.path.isfile(o_file):
            
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : FUNC : QC PDF" )
        
            os.chdir(os.path.join(func_path, 'qc'))
            
            files_for_pdf = ''
            jpg_list = glob.glob('*.jpg')
            if len(jpg_list):
                files_for_pdf = files_for_pdf + ' *.jpg'
            png_list = glob.glob('*.png')
            if len(png_list):
                files_for_pdf = files_for_pdf + ' *.png'
            if len(files_for_pdf):
                utils.Cmd("convert "+files_for_pdf+" -quality 100 " + o_file).run()
            
            
            
            
