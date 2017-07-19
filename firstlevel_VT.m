function matlabbatch = firstlevel_VT(cfg,prefix,scan,results_dir,onsname,i_sub,mparams,sliceTiming)

onset_dir = fullfile(cfg.dirs.timing_dir);
n_runs = cfg.scans_nruns(ismember(cfg.scans,scan)); % get number of runs 

matlabbatch{1}.spm.stats.fmri_spec.dir = {results_dir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = cfg.TR;

for i_run=1:n_runs
    % get file names
    run_dir = fullfile(cfg.dirs.preproc_dir,scan,cfg.sub(i_sub).id,['run' num2str(i_run)]);
    tmp = spm_select('FPList',run_dir,['^' prefix '.*\.nii$']);
    files = cell(size(tmp,1),1);
    for i = 1:size(tmp,1)
        files{i} = [tmp(i,:) ',1'];
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).scans = files;
    if sliceTiming
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = cfg.n_slices;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = ceil(cfg.n_slices/2);
    else
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).multi = {fullfile(onset_dir,scan,cfg.sub(i_sub).id,sprintf('sub%s_%s_timing_%d.mat',cfg.sub(i_sub).id,onsname,i_run))};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).regress = struct('name', {}, 'val', {});
    if mparams
        rp_name = spm_select('fplist',run_dir,'^rp_a.*\.txt$');
        matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).multi_reg = {rp_name};
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).multi_reg = {''};
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(i_run).hpf = 128;
end

matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
end
