function cfg = analysis_VT(i_sub,scans,cfg)
% i_sub = subject to be processed (i.e 1:nsubs). see config_subjects_VT.m
% for subject list
% scans = cell array of strings of scans to be analyzed 
% (e.g {'MVPA_OB_pre','MVPA_OB_post'}
% cfg = cfg file created with config_subjects_VT.m
%% preprocessing
cfg.do.dicom_convert = 0; % not currently recommended, use mriconvert
cfg.do.move_dummyTRs = 0;
cfg.do.slicetiming = 0;
cfg.do.realign = 0;
cfg.do.data_qc = 0;m
cfg.do.coreg = 0; 
cfg.do.check_coreg = 0;
cfg.do.segment = 0; % new segment algorithm
cfg.do.check_segment = 0;
cfg.do.normalize_mean_anat = 0;
cfg.do.check_normalize = 0;
cfg.do.normalize_functional = 0; % new normalize algorithm
cfg.do.binarize_mask1 = 0;
cfg.do.inverse_normalize = 0;
cfg.do.convert_dtype = 0;
cfg.do.binarize_mask = 0;
cfg.do.smooth = 0;

%% roi analysis for connectivity
cfg.do.regressCovariates = 0;
cfg.do.niiToMat = 0;
cfg.do.getROItimeSeries = 0;

%% first level
cfg.do.makeTimingFiles = 0;
cfg.do.firstlevel = 1;
cfg.do.firstlevel_overwrite = 0;
cfg.do.contrasts_only = 1; % if this is ticked, the firstlevel will not be executed again

%% MVPA 
cfg.do.decode = 0;
cfg.do.normalize_decode = 0;
cfg.do.smooth_decode = 0;

%% create cfg   for all subjects
if exist('cfg','var') % if cfg was passed, use it; else, create one
    cfg = config_subjects_VT(cfg);
else
    cfg = config_subjects_VT;
end

data_dir = cfg.dirs.data_dir;

%% PREPROCESSING 

%% dicom conversion
% TO DO: figure out how this conversion is different than mriconvert; voxel
% size changes when converted with spm_dicom_convert
if cfg.do.dicom_convert
    for ss=1:length(scans)
        sprintf('Importing DICOM files for subject %s %s...\n',cfg.sub(i_sub).id,scans{ss})
        input_dir = fullfile(data_dir,'DICOM',cfg.sub(i_sub).id,scans{ss}); % folder where DICOMs are located
        disp('Getting DICOM filenames...')
        P = spm_select('FPListRec',input_dir,'1*.IMA'); % recursively get DICOMs
        disp('Reading DICOM headers...')
        hdr = spm_dicom_headers(P,true); % get DICOM headers
        spm_dicom_convert(hdr,'all','series',cfg.suffix,fullfile(data_dir,'preproc',cfg.sub(i_sub).id,scans{ss}));

        clear hdr; % free memory

        disp('done')
    end
end
%% move dummy TRs
% move first n time points into new folder; set with cfg.n_dummy
if cfg.do.move_dummyTRs
    for ss=1:length(scans)
        n_runs = cfg.scans_nruns(ismember(cfg.scans,scans{ss}));
        for i_run=1:n_runs
        mkdir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,['run' num2str(i_run)],'dummyTRs'));
        cd(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,['run' num2str(i_run)]));
        !move *0001.* dummyTRs;
        !move *0002.* dummyTRs;
        !move *0003.* dummyTRs;
        !move *0004.* dummyTRs;
        end
    end
end

%% slice timing correction
if cfg.do.slicetiming
    % initializeu the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch

    n_slices = cfg.n_slices;
    TR = cfg.TR;
    prefix = 'MR'; % prefix of images to be slice time corrected
    
    for ss=1:length(scans)
        
        n_runs = cfg.scans_nruns(ismember(cfg.scans,scans{ss}));

        matlabbatch{ss}.spm.temporal.st.scans = select_files_VT(cfg,i_sub,prefix,scans{ss},n_runs);
        matlabbatch{ss}.spm.temporal.st.nslices = n_slices;
        matlabbatch{ss}.spm.temporal.st.tr = TR;
        matlabbatch{ss}.spm.temporal.st.ta = TR - (TR/n_slices);
        switch lower(cfg.sliceorder)
            case 'ascending'
                so = 1:n_slices;
            case 'descending'
                so = n_slices:-1:1 ;
            case 'interleaved'
                if mod(cfg.n_slices,2)
                    so = [1:2:n_slices 2:2:n_slices];
                else
                    warning('Using settings for Siemens, if no Siemens Scanner is used please find this message and change the script!') 
                    so = [2:2:n_slices 1:2:n_slices];
                end
            case 'interleaved descending'
                if mod(cfg.n_slices,2)
                    so = [n_slices:-2:1 n_slices-1:-2:1];
                else
                    warning('Using settings for Siemens, if no Siemens Scanner is used please find this message and change the script!')
                    so = [n_slices-1:-2:1 n_slices:-2:1];
                end
        end
        matlabbatch{ss}.spm.temporal.st.so = so;
        matlabbatch{ss}.spm.temporal.st.refslice = cfg.reference_slice;
        matlabbatch{ss}.spm.temporal.st.prefix = 'a';

    end
    spm_jobman('run',matlabbatch)
end

%% realignment 
if cfg.do.realign
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch

    prefix = 'aMR'; % prefix of images to be realigned
    preproc_dir = cfg.dirs.preproc_dir;
    for ss=1:length(scans)
        
    n_runs = cfg.scans_nruns(ismember(cfg.scans,scans{ss}));
    
    matlabbatch{ss}.spm.spatial.realign.estwrite.data = select_files_VT(cfg,i_sub,prefix,scans{ss},n_runs);
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.quality = 0.9; % default 0.9, 0.95 more accurate (but slower)
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.sep = 2; % default 4, 2 gives more accuracte (but slower) results
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{ss}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{ss}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{ss}.spm.spatial.realign.estwrite.roptions.interp = 7;
    matlabbatch{ss}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{ss}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{ss}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    
    end
    
    cd(fullfile(preproc_dir,scans{ss},cfg.sub(i_sub).id)); % cd to subject's preproc dir so movement figure is saved
    spm_figure('GetWin','Graphics'); % ensures the graphics window is open so .ps movement figure is saved
    spm_jobman('run',matlabbatch)
end

%% data qc using QC_badscanfinder (based on art toolbox)
if cfg.do.data_qc
      
    prefix = 'raMR';
    global_thresh = 1.5;
    diff_thresh = 1.0;
    mv_thresh = 0.3;
    
    
    for ss=1:length(scans)
       n_runs = cfg.scans_nruns(ismember(cfg.scans,scans{ss}));
       images = select_files_VT(cfg,i_sub,prefix,scans{ss},n_runs);
       for i_run=1:n_runs
           run_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,['run' num2str(i_run)]);
           realign_p = spm_select('fplist',run_dir,'^rp_a.*\.txt$');
           images_char = char(images{i_run});
           QC_badscanfinder(images_char,realign_p,global_thresh,diff_thresh,mv_thresh,1); 
           clear realign_p
       end 
    end
    
    
end

%% coregistration
if cfg.do.coreg
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    prefix = 'meanaMR'; % prefix of mean image
    struct_prefix = 'MR';
    
    for ss=1:length(scans)
        
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        source_path = spm_select('fplist',anat_dir,['^' struct_prefix '.*\.(nii|img)$']);
        source_path = [source_path ',1'];
        
        mean_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'run1'); % mean image is in run1 dir
        ref_path = spm_select('fplist',mean_dir,['^' prefix '.*\.(nii|img)$']);
        ref_path = [ref_path ',1'];
        
        matlabbatch{ss}.spm.spatial.coreg.estimate.ref = {ref_path};
        matlabbatch{ss}.spm.spatial.coreg.estimate.source = {source_path};
        matlabbatch{ss}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{ss}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{ss}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{ss}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{ss}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];  
    end
    
    spm_jobman('run',matlabbatch)
    
end

%% visually assess quality of coregistration 
if cfg.do.check_coreg
    for ss=1:length(scans)
        anat_prefix = 'MR';
        mean_prefix = 'meanaMR';
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        anat_path = spm_select('fplist',anat_dir,['^' anat_prefix '.*\.(img|nii)$']);
        anat_path = [anat_path ',1'];
        mean_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'run1');
        mean_path = spm_select('fplist',mean_dir,['^' mean_prefix '.*\.(img|nii)$']);
        anat_path = char(anat_path,mean_path);
        spm_check_registration(anat_path);
        
        prompt = 'Press any key when done QCing \n';
        pause;
    end
end


%% segmentation (new algorithm in spm12)
if cfg.do.segment
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    struct_prefix = 'MR';
    
    for ss=1:length(scans)
        
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        anat_path = spm_select('fplist',anat_dir,['^' struct_prefix '.*\.(nii|img)$']);
        anat_path = [anat_path ',1'];

        matlabbatch{ss}.spm.spatial.preproc.channel.vols = {anat_path};
        matlabbatch{ss}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{ss}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{ss}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(1).tpm = {[cfg.dirs.spm_dir '/tpm/TPM.nii,1']};
        matlabbatch{ss}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{ss}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(2).tpm = {[cfg.dirs.spm_dir '/tpm/TPM.nii,2']};
        matlabbatch{ss}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{ss}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(3).tpm = {[cfg.dirs.spm_dir '/tpm/TPM.nii,3']};
        matlabbatch{ss}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{ss}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(4).tpm = {[cfg.dirs.spm_dir '/tpm/TPM.nii,4']};
        matlabbatch{ss}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{ss}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(5).tpm = {[cfg.dirs.spm_dir '/tpm/TPM.nii,5']};
        matlabbatch{ss}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{ss}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(6).tpm = {[cfg.dirs.spm_dir '/tpm/TPM.nii,6']};
        matlabbatch{ss}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{ss}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{ss}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{ss}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{ss}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{ss}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{ss}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{ss}.spm.spatial.preproc.warp.samp = 2; % default is 3, 2 is more accurate 
        matlabbatch{ss}.spm.spatial.preproc.warp.write = [1 1];

        
    end
    
    spm_jobman('run',matlabbatch)
 
end

%% check visually if segmentation worked
if cfg.do.check_segment
    
    for ss=1:length(scans)
        seg_prefix = 'c';
        segment_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        segment_path = spm_select('fplist',segment_dir,['^' seg_prefix '.*\.(img|nii)$']);
        segment_path = [segment_path repmat(',1',size(segment_path,1),1)];

        spm_check_registration(segment_path)
        prompt = 'Press any key when done QCing \n';
        pause;
    end
    
end

%% normalize struct and mean to check quality of normalization
% TO DO: normalize anat outside of loop through scans, right now 
% the anat is normalized every iteration
if cfg.do.normalize_mean_anat
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    prefix = 'meanaMR'; % prefix of mean image
    cfg.suffix = '(nii|img)';
    anat_prefix = 'MR';
   
    for ss=1:length(scans)
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        normparams_path = spm_select('fplist',anat_dir,['^y_' anat_prefix '.*\.(nii|img)']);
        matlabbatch{ss}.spm.spatial.normalise.write.subj.def = {normparams_path};
        anat_path = spm_select('fplist',anat_dir,['^' anat_prefix '.*\.(nii|img)$']);
        mean_path = spm_select('fplist',fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'run1'),['^' prefix '.*\.(nii|img)$']);
        matlabbatch{ss}.spm.spatial.normalise.write.subj.resample = {[anat_path ',1']; [mean_path ',1']};
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                   78 76 85];
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.interp = 4;
    end
    
    spm_jobman('run',matlabbatch)
end

%% visually assess quality of normalization 
if cfg.do.check_normalize
    for ss=1:length(scans)
        anat_prefix = 'wMR';
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        anat_path = spm_select('fplist',anat_dir,['^' anat_prefix '.*\.(img|nii)$']);
        anat_path = [anat_path ',1'];
        ref_path = fullfile(cfg.dirs.spm_dir,'canonical','single_subj_T1.nii,1');
        anat_path = char(anat_path,ref_path);
        spm_check_registration(anat_path);

        prompt = 'Press any key when done QCing \n';
        pause;
    end
end


%% functional normalization (new algorithim in spm12)
if cfg.do.normalize_functional
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    prefix = 'raMR';
    cfg.suffix = '(nii|img)';
    anat_prefix = 'MR';

    for ss=1:length(scans)
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        normparams_path = spm_select('fplist',anat_dir,['^y_' anat_prefix '.*\.(nii|img)']);
        matlabbatch{ss}.spm.spatial.normalise.write.subj.def = {normparams_path};
        matlabbatch{ss}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPListRec',fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id),['^' prefix '.*\.(nii|img)$']));
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                   78 76 85];
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.interp = 4;
    end
    
    spm_jobman('run',matlabbatch)
    
end

%%
if cfg.do.binarize_mask1
    
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');
    
    clear matlabbatch
    
    for ss=1:length(scans)
        niiRois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen','*.nii'));
        for j=1:length(niiRois)
            clear matlabbatch
            matlabbatch{ss}.spm.util.imcalc.input = {fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen',niiRois(j).name(1:length(niiRois(j).name)));};
            matlabbatch{ss}.spm.util.imcalc.output = [niiRois(j).name(1:length(niiRois(j).name))];
            matlabbatch{ss}.spm.util.imcalc.outdir = {fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen')};
            matlabbatch{ss}.spm.util.imcalc.expression = 'i1>0.9';
            matlabbatch{ss}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{ss}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{ss}.spm.util.imcalc.options.mask = 0;
            matlabbatch{ss}.spm.util.imcalc.options.interp = 1;
            matlabbatch{ss}.spm.util.imcalc.options.dtype = 2;
            spm_jobman('run',matlabbatch)
        end

    end
end

%%
if cfg.do.convert_dtype
    
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');
    
    clear matlabbatch
    
    for ss=1:length(scans)
        niiRois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen','*.nii'));
        for j=1:length(niiRois)
            clear matlabbatch
            matlabbatch{ss}.spm.util.imcalc.input = {fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen',niiRois(j).name(1:length(niiRois(j).name)));};
            matlabbatch{ss}.spm.util.imcalc.output = ['dtypeConverted_' niiRois(j).name(1:length(niiRois(j).name))];
            matlabbatch{ss}.spm.util.imcalc.outdir = {fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen')};
            matlabbatch{ss}.spm.util.imcalc.expression = 'i1';
            matlabbatch{ss}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{ss}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{ss}.spm.util.imcalc.options.mask = 0;
            matlabbatch{ss}.spm.util.imcalc.options.interp = 1;
            matlabbatch{ss}.spm.util.imcalc.options.dtype = 2;
            spm_jobman('run',matlabbatch)
        end

    end
end

%% inverse normalization 
if cfg.do.inverse_normalize
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    prefix = 'dtypeConverted_';
    cfg.suffix = '(nii|gimg)';
    anat_prefix = 'MR';
    
    for ss=1:length(scans)
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        normparams_path = spm_select('fplist',anat_dir,['^iy_' anat_prefix '.*\.(nii|img)']);

        matlabbatch{ss}.spm.util.defs.comp{1}.def = {normparams_path};
        matlabbatch{ss}.spm.util.defs.out{1}.pull.fnames = cellstr(spm_select('FPListRec',fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen'),['^' prefix '.*\.(nii|img)$']));
        matlabbatch{ss}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
        matlabbatch{ss}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{ss}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{ss}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];

    end
    
    spm_jobman('run',matlabbatch)
    
end

%%
if cfg.do.binarize_mask
    
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');
    
    clear matlabbatch
    
    for ss=1:length(scans)
        niiRois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen','wdtypeConverted_*.nii'));
        for j=1:length(niiRois)
            clear matlabbatch
            matlabbatch{ss}.spm.util.imcalc.input = {fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen',niiRois(j).name(1:length(niiRois(j).name)));};
            matlabbatch{ss}.spm.util.imcalc.output = ['bwdtypeConverted_' niiRois(j).name(1:length(niiRois(j).name))];
            matlabbatch{ss}.spm.util.imcalc.outdir = {fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen')};
            matlabbatch{ss}.spm.util.imcalc.expression = 'i1>0.9';
            matlabbatch{ss}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{ss}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{ss}.spm.util.imcalc.options.mask = 0;
            matlabbatch{ss}.spm.util.imcalc.options.interp = 1;
            matlabbatch{ss}.spm.util.imcalc.options.dtype = 2;
            spm_jobman('run',matlabbatch)
        end

    end
end

%% smooth normalized data
if cfg.do.smooth
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    prefix = 'wraMR';
    FWHM = cfg.FWHM;
    
    for ss=1:length(scans)
        matlabbatch{ss}.spm.spatial.smooth.data = cellstr(spm_select('FPListRec',fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id),['^' prefix '.*\.(nii|img)$']));
        matlabbatch{ss}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
        matlabbatch{ss}.spm.spatial.smooth.dtype = 0 ;
        matlabbatch{ss}.spm.spatial.smooth.im = 0 ;
        matlabbatch{ss}.spm.spatial.smooth.prefix = sprintf('s%02d',FWHM);
    end
    spm_jobman('run',matlabbatch)
    
end

%% regress covariates of no interest (CSF, WM, movement parameters) using CONN toolbox 
if cfg.do.regressCovariates
    
    clear BATCH
    
    prefix = 'ra';
    anat_prefix = 'MR';
    
    
    for ss=1:length(scans)
        
        n_runs = cfg.scans_nruns(ismember(cfg.scans,scans{ss})); % get number of runs
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        
        if ~exist(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id),'dir')
            mkdir(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id));
        end
        
        BATCH{ss}.filename=fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id);
        BATCH{ss}.Setup.isnew=1;
        BATCH{ss}.Setup.nsubjects=1;
        BATCH{ss}.Setup.RT=cfg.TR;
        BATCH{ss}.Setup.voxelmask=2; %Specifies the type of mask to use for voxel-based analyses: 1 for a fixed template (explicit) mask, and 2 for a subject-specific (implicit) mask
        BATCH{ss}.Setup.voxelresolution=3; %Specifies resolution for voxel-based analyses: 1 for fixed template resolution (2mm isotropic); 2 for same resolution as structural scan; and 3 for same resolution as functional scan
        BATCH{ss}.Setup.outputfiles=[0,1,0]; % outputfiles(2) set to 1 for outputting confound corrected BOLD time series
        BATCH{ss}.Setup.covariates.names={'motion'};
        BATCH{ss}.Setup.conditions.names={'rest'}; % not really rest, but not incoporating task timing yet
        for i_run=1:n_runs
            BATCH{ss}.Setup.covariates.files{1}{1}{i_run}=cellstr(spm_select('fplist',fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,['run' num2str(i_run)]),'^rp_a.*\.txt$'));
            BATCH{ss}.Setup.conditions.onsets{1}{1}{i_run}=[0];
            BATCH{ss}.Setup.conditions.durations{1}{1}{i_run}=[inf];
            BATCH{ss}.Setup.functionals{1}{i_run}=cellstr(spm_select('FPListRec',fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,['run' num2str(i_run)]),['^' prefix '.*\.(nii|img)$']));
        end
        BATCH{ss}.Setup.structurals{1}=spm_select('fplist',anat_dir,['^' anat_prefix '.*\.(nii|img)$']);
        BATCH{ss}.Setup.masks.Grey.files{1}=spm_select('fplist',anat_dir,'^c1.*\.(nii|img)$');
        BATCH{ss}.Setup.masks.White.files{1}=spm_select('fplist',anat_dir,'^c2.*\.(nii|img)$');
        BATCH{ss}.Setup.masks.CSF.files{1}=spm_select('fplist',anat_dir,'^c3.*\.(nii|img)$');
        BATCH{ss}.Setup.rois.names={}; % initialize to empty so CONN doesn't import the default ROIs
        BATCH{ss}.Setup.rois.dimensions={};
        BATCH{ss}.Setup.rois.files='';
        BATCH{ss}.Setup.overwrite='Yes';
        BATCH{ss}.Setup.done=1; %
        BATCH{ss}.Preprocessing.confounds.names={'White','CSF','motion'}; 
        BATCH{ss}.Preprocessing.confounds.dimensions={3,3,6};
        BATCH{ss}.Preprocessing.confounds.deriv={0,0,1};
        BATCH{ss}.Preprocessing.filter=[.01,inf]; % high pass filter at 0.01 Hz
        BATCH{ss}.Preprocessing.overwrite='Yes';
        BATCH{ss}.Preprocessing.done=1;
    end
    
    conn_batch(BATCH); % run batch 
    
end

%% convert nii masks to .mat for marsbar
% shen atlas
if cfg.do.niiToMat
   for ss=1:length(scans)
       niiRois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen','b*.nii'));
       for j=1:length(niiRois)
        imgname = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen',niiRois(j).name(1:length(niiRois(j).name)));
        fname = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen',sprintf('%s_roi.mat',niiRois(j).name(1:length(niiRois(j).name)-4)));
        o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',0,'func', 'img'));
        saveroi(o, fname)
       end
   end
end

% % aal atlas
% if cfg.do.niiToMat
%    for ss=1:length(scans)
%        niiRois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_aal','wMNI_*roi.nii'));
%        for j=1:length(niiRois)
%         imgname = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_aal',sprintf('wMNI_%s_roi.nii',niiRois(j).name(6:length(niiRois(j).name)-8)));
%         fname = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_aal',sprintf('wMNI_%s_roi.mat',niiRois(j).name(6:length(niiRois(j).name)-8)));
%         o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',0,'func', 'img'));
%         saveroi(o, fname)
%        end
%    end
% end

%% extract roi time series

% shen atlas
if cfg.do.getROItimeSeries
    for ss=1:length(scans)
        %aalrois = dir(fullfile(cfg.dirs.preproc_dir,cfg.sub(i_sub).id,scans{ss},'masks_shen','wMNI_*roi.mat')); 
        shen_rois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen','b*roi.mat'));
        for j=1:length(shen_rois)
            roi_file = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_shen',sprintf('b%s_roi.mat',shen_rois(j).name(2:length(shen_rois(j).name)-8)));
            clear rois;
            rois = maroi('load_cell', roi_file);

            P = fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'results','preprocessing','niftiDATA_Subject001_Condition000.nii'); % confound corrected BOLD time course from CONN
            [timeCourse, ~, ~, ~] = getdata(rois{1}, P);
            timeCourse(:,~all(timeCourse)) = NaN; % convert any columns with all 0s into NaN so they can be ignored when averaging
            timeCourse = nanmean(timeCourse,2); % mean across all voxels, ignoring NaN

            if ~exist(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'shenfmriTimeSeries'),'dir')
                mkdir(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'shenfmriTimeSeries'));
            end
            save(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'shenfmriTimeSeries',shen_rois(j).name(34:length(shen_rois(j).name))),'timeCourse');
        end
    end
end  

% % aal atlas
% if cfg.do.getROItimeSeries
%     for ss=1:length(scans)
%         aalrois = dir(fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_aal','wMNI_*roi.mat')); 
%         for j=1:length(aalrois)
%             roi_file = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'masks_aal',sprintf('wMNI_%s_roi.mat', aalrois(j).name(6:length(aalrois(j).name)-8)));
%             clear rois;
%             rois = maroi('load_cell', roi_file);
% 
%             P = fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'results','preprocessing','niftiDATA_Subject001_Condition000.nii'); % confound corrected BOLD time course from CONN
%             [timeCourse, ~, ~, ~] = getdata(rois{1}, P);
%             timeCourse(:,~all(timeCourse)) = NaN; % convert any columns with all 0s into NaN so they can be ignored when averaging
%             timeCourse = nanmean(timeCourse,2); % mean across all voxels, ignoring NaN
% 
%             if ~exist(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'aalfmriTimeSeries'),'dir')
%                 mkdir(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'aalfmriTimeSeries'));
%             end
%             save(fullfile(cfg.dirs.results_dir,scans{ss},cfg.sub(i_sub).id,'aalfmriTimeSeries',aalrois(j).name(6:length(aalrois(j).name)-8)),'timeCourse');
%         end
%     end
% end   

%% make timing files 

if cfg.do.makeTimingFiles
   for ss=1:length(scans)
       if (strcmp('MVPA_OB_pre',scans{ss}) || strcmp('MVPA_OB_post',scans{ss}))
            makeTimingFile_VT_MVPA_OB(i_sub,scans{ss});
       elseif (strcmp('MVPA_cat_post',scans{ss}))
            makeTimingFile_VT_MVPA_cat(i_sub,scans{ss});
       end
   end
end
%% first level

if cfg.do.firstlevel
    
    spm('defaults','fmri')
    spm_jobman('initcfg');

    for ss=1:length(scans)
        prefix = 's06'; %raMR for making MVPA betas, s06 for classical univariate 
        onsname = [scans{ss} '_allMorphs']; % name of onsets file
        resname = 'NormSmooth_allMorphs'; % results path; noNormNoSmooth for MVPA betas, NormSmooth for classical univariate
        mparams = 1; % 1 = include movement parameters as regressors, 0 = do not

        if ~cfg.do.contrasts_only
            standard_firstlevel(cfg,i_sub,prefix,onsname,resname,scans{ss},mparams) % See bottom of the script
        end
        
        results_dir = fullfile(cfg.dirs.results_dir,scans{ss},'GLM',resname,cfg.sub(i_sub).id);
  
  %% contrasts - MPVA 
  
  if (strcmp('MVPA_OB_pre',scans{ss}) || strcmp('MVPA_OB_post',scans{ss}))
%       % set contrasts
        clear matlabbatch

        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(results_dir,'SPM.mat')};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'all non OB';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1 1 1 1 0];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'f_93_26';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 0 0 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'f_75_32';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 1 0 0 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'f_61_40';
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 1 0 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'f_50_50';
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 1 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'f_40_61';
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'f_32_75';
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'f_26_93';
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'OB';
        matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'skay>gark';
        matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [1 1 1 0 -1 -1 -1 0];
        matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'gark>skay';
        matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [-1 -1 -1 0 1 1 1 0];
        matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'allTrials';
        matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [1 1 1 1 1 1 1 1];
        matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'movement';
        matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = [0 0 0 0 0 0 0 0 1 1 1 1 1 1];
        matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'OB>non OB';
        matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = [(-1/7) (-1/7) (-1/7) (-1/7) (-1/7) (-1/7) (-1/7) 1];
        matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'repl';
        matlabbatch{1}.spm.stats.con.delete = 1;
        spm_jobman('run',matlabbatch)

  end 

  if (strcmp('MVPA_cat_post',scans{ss}))
      % set contrasts
      clear matlabbatch
      matlabbatch{1}.spm.stats.con.spmmat = {fullfile(results_dir,'SPM.mat')};
      matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'all morphs';
      matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1 1 1 1];
      matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'f_93_26';
      matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 0 0 0 0 0];
      matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'f_75_32';
      matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 1 0 0 0 0 0];
      matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'f_61_40';
      matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 1 0 0 0 0];
      matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'f_50_50';
      matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 1 0 0 0];
      matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'f_40_61';
      matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 1 0 0];
      matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'f_32_75';
      matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 1 0];
      matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'f_26_93';
      matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 1];
      matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'skay>gark';
      matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [1 1 1 0 -1 -1 -1];
      matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'gark>skay';
      matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [-1 -1 -1 0 1 1 1];
      matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = '50 morph > rest';
      matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [(-1/6) (-1/6) (-1/6) 1 (-1/6) (-1/6) (-1/6)];
      matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = '50 morph < rest';
      matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [(1/6) (1/6) (1/6) -1 (1/6) (1/6) (1/6)];
      matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'repl';
      matlabbatch{1}.spm.stats.con.delete = 1;
      spm_jobman('run',matlabbatch)
     
  end
  
    end
end
  %% decode 
  if cfg.do.decode
      searchlights = [4];
      decode_software = 'libsvm'; % classification algorithm to be used; can be libsvm or correlation_classifier
      resname = 'noNormNoSmooth_allMorphs'; % first level results to be decoded
      mvpa_resname = 'category_xClass2_cat';
      overwrite = 0; % overwrite results; 0=do not overwrite, 1=overwrite 
      
      for ss=1:length(scans)
        for i_searchlight=1:length(searchlights)
            results_dir = fullfile(cfg.dirs.results_dir,scans{ss},'MVPA',resname,mvpa_resname,['searchlight' num2str(searchlights(i_searchlight))],cfg.sub(i_sub).id);
            beta_dir = fullfile(cfg.dirs.results_dir,scans{ss},'GLM',resname,cfg.sub(i_sub).id);
            decode_VT_OB_MVPA(results_dir,beta_dir,scans{ss},cfg.sub(i_sub).id,decode_software,searchlights(i_searchlight),overwrite);
        end         
      end
  end
  
  %% normalize decoding 
  if cfg.do.normalize_decode
    % initialize the jobmanager
    spm('defaults','fmri')
    spm_jobman('initcfg');

    clear matlabbatch
    
    prefix = 'res';
    resname = 'noNormNoSmooth_allMorphs';
    mvpa_resname = 'category_xClass2_cat';
    searchlight = 4;
    cfg.suffix = '(nii|img)';
    anat_prefix = 'MR';

    for ss=1:length(scans)
        anat_dir = fullfile(cfg.dirs.preproc_dir,scans{ss},cfg.sub(i_sub).id,'anat');
        normparams_path = spm_select('fplist',anat_dir,['^y_' anat_prefix '.*\.(nii|img)']);
        matlabbatch{ss}.spm.spatial.normalise.write.subj.def = {normparams_path};
        matlabbatch{ss}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPListRec',fullfile(cfg.dirs.results_dir,scans{ss},'MVPA',resname,mvpa_resname,['searchlight' num2str(searchlight)],cfg.sub(i_sub).id),['^' prefix '.*\.(nii|img)$']));
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                   78 76 85];
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{ss}.spm.spatial.normalise.write.woptions.interp = 4;
    end
    
    spm_jobman('run',matlabbatch) 
          
  end
  
    %% smooth decoding data
    if cfg.do.smooth_decode
        % initialize the jobmanager
        spm('defaults','fmri')
        spm_jobman('initcfg');

        clear matlabbatch

        prefix = 'wres';
        resname = 'noNormNoSmooth_allMorphs';
        mvpa_resname = 'category_xClass2_cat';
        searchlight = 4;
        FWHM = cfg.FWHM;

        for ss=1:length(scans)
            matlabbatch{ss}.spm.spatial.smooth.data = cellstr(spm_select('FPListRec',fullfile(cfg.dirs.results_dir,scans{ss},'MVPA',resname,mvpa_resname,['searchlight' num2str(searchlight)],cfg.sub(i_sub).id),['^' prefix '.*\.(nii|img)$']));
            matlabbatch{ss}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
            matlabbatch{ss}.spm.spatial.smooth.dtype = 0 ;
            matlabbatch{ss}.spm.spatial.smooth.im = 0 ;
            matlabbatch{ss}.spm.spatial.smooth.prefix = sprintf('s%02d',FWHM);
        end
        spm_jobman('run',matlabbatch)

    end

end

%% SUBFUNCTIONS
function standard_firstlevel(cfg,i_sub,prefix,onsname,resname,scan,mparams)

sliceTiming = 1; % was slice timing performed? 1=yes,0=no
results_dir = fullfile(cfg.dirs.results_dir,scan,'GLM',resname,cfg.sub(i_sub).id);
if ~isdir(results_dir), mkdir(results_dir), end
if exist(fullfile(results_dir,'SPM.mat'),'file') && cfg.do.firstlevel_overwrite == 0 
    error('SPM.mat already exists for %s',results_dir)
else
    spm_unlink(fullfile(results_dir,'SPM.mat')); % delete SPM to overwrite without dialog
end

% Set up model
matlabbatch = firstlevel_VT(cfg,prefix,scan,results_dir,onsname,i_sub,mparams,sliceTiming); % if you need this function, let me know

spm_jobman('run',matlabbatch)
end
