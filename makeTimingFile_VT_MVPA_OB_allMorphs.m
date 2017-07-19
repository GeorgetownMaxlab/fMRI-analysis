function makeTimingFile_VT_MVPA_OB_allMorphs(sub,scan)

% scan arg is either 'MVPA_OB_pre' or 'MVPA_OB_post'

% names for regressors
names{1} = 'f_93_26';
names{2} = 'f_75_32';
names{3} = 'f_61_40';
names{4} = 'f_50_50';
names{5} = 'f_40_61';
names{6} = 'f_32_75';
names{7} = 'f_26_93';
names{8} = 'OB';

cfg = config_subjects_VT;

TR = cfg.TR;
n_dummyTRs = cfg.n_dummy;
% above variables used in timing calculation to account for TRs removed at
% beginning of each run

for i_run=1:cfg.scans_nruns(ismember(cfg.scans,scan))
    
    % load psychtoolbox mri timing file 
    tm_file_path = dir(fullfile(cfg.dirs.timing_dir,scan,cfg.sub(sub).id,sprintf('20*run%d.mat',i_run)));
    if size(tm_file_path,1)>1, error('More than 1 timing file for run %d found',i_run), end
    load(fullfile(cfg.dirs.timing_dir,scan,cfg.sub(sub).id,tm_file_path(1).name));
    
    % onsets in ms after beginning of scan
    onsets = cell(1,8);

for i=1:length(trialOutput)
    if trialOutput(i).oddball == 1, k=8; % oddball trial
    elseif trialOutput(i).sResp == 1, k=8; % code non-oddball trials where subject responded (false alarms) as OB trial
    elseif round(trialOutput(i).stimuli{1}(1)) == 93, k=1;
    elseif round(trialOutput(i).stimuli{1}(1)) == 76, k=2;
    elseif round(trialOutput(i).stimuli{1}(1)) == 62, k=3;
    elseif round(trialOutput(i).stimuli{1}(1)) == 50, k=4;
    elseif round(trialOutput(i).stimuli{1}(1)) == 41, k=5;
    elseif round(trialOutput(i).stimuli{1}(1)) == 33, k=6;
    elseif round(trialOutput(i).stimuli{1}(1)) == 27, k=7;
    end
    onsets{1,k} = [onsets{1,k} (single((trialOutput(i).FixationOnsetTime(1) - exptdesign.scanStart))) - (TR*n_dummyTRs)];
end
    
    durations{1} = repmat(6,1,length(onsets{1}))';
    durations{2} = repmat(6,1,length(onsets{2}))';
    durations{3} = repmat(6,1,length(onsets{3}))';
    durations{4} = repmat(6,1,length(onsets{4}))';
    durations{5} = repmat(6,1,length(onsets{5}))';
    durations{6} = repmat(6,1,length(onsets{6}))';
    durations{7} = repmat(6,1,length(onsets{7}))';
    durations{8} = repmat(6,1,length(onsets{8}))';
    
    save(fullfile(cfg.dirs.timing_dir,scan,cfg.sub(sub).id,['sub' cfg.sub(sub).id '_' scan '_allMorphs_timing_' num2str(i_run)]),'names','onsets','durations');
    clear durations onsets
end


end