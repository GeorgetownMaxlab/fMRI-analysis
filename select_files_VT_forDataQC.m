function files = select_files_VT_forDataQC(cfg,i_sub,prefix,scan,n_runs)

files = {};
preproc_dir = cfg.dirs.preproc_dir;
scan_dir = fullfile(preproc_dir,scan,cfg.sub(i_sub).id);
suffix = cfg.suffix;

for i_run=1:n_runs
    run_dir = fullfile(scan_dir,['run' num2str(i_run)]);
    tmp = spm_select('FPList',run_dir,['^' prefix '.*\.' suffix '$']);
    tmp2 = cell(size(tmp,1),1);
    for i = 1:size(tmp,1)
        tmp2{i} = [tmp(i,:) ',1'];
    end
    files = [files; tmp2];
end

end