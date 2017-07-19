function cfg = config_subjects_VT(cfg)

%% set analysis specifc variables 
cfg.name = 'vibrotactile'; % project name
cfg.scans = {'MVPA_OB_pre', 'MVPA_OB_post', 'MVPA_cat_post'};
cfg.scans_nruns = [6 6 6]; % corresponds to cfg.scans 
cfg.suffix = 'img'; % img or nii

%% set paths 

[~,hostname] = system('hostname');
if strcmp(cellstr(hostname),'hmaxinator') 
    data_dir = '/data/Data/Data_MBP/vibrotactile';
    spm_dir = '/home/malone/Documents/MATLAB/spm12';
else 
    data_dir = 'C:\Users\Patrick Malone\Documents\Data\vibrotactile';
    spm_dir = 'C:\Users\Patrick Malone\Documents\MATLAB\spm12';
end

cfg.dirs.spm_dir = spm_dir;
cfg.dirs.data_dir = data_dir;
cfg.dirs.preproc_dir = fullfile(data_dir,'preproc');
cfg.dirs.PIRE_preproc_dir = fullfile(data_dir,'connectivity_PIRE'); % for PIRE connectivity analyis 
cfg.dirs.results_dir = fullfile(data_dir,'results');
cfg.dirs.behav_dir = fullfile(data_dir,'behavData');
cfg.dirs.timing_dir = fullfile(cfg.dirs.preproc_dir,'timing');
cfg.dirs.mask_dir = fullfile(data_dir,'masks');

%% set scanner-specific and basic analysis variables

cfg.TE = 0.03; % in seconds
cfg.TR = 2;
cfg.n_slices = 33;
cfg.voxdims = [3.2 3.2 4.2];
cfg.sliceorder = 'interleaved'; % ascending, descending, interleaved, or interleaved descending
cfg.reference_slice = 33; % set to middle slice in acquisition sequence
cfg.FWHM = 6; % smoothing kernel (mm)
cfg.n_dummy = 4; % how many TRs to remove at beginning of run

%% set subjects and runs to use in analysis
% TO DO: this is currently not implemented
excluded_subjects = [];
excluded_runs = {};

%% set subject variables 

% Subject 1
cfg.sub(1).id = 'MR1034'; % subject id anonymized
cfg.sub(1).age = 27;
cfg.sub(1).gender = 0;  % 0=male, 1=female
cfg.sub(1).handedness = 1; % 0=left, 1=right
cfg.sub(1).preproc.corrupted_slices = [];

% Subject 2
cfg.sub(2).id = 'MR1035'; % subject id anonymized
cfg.sub(2).age = 19;
cfg.sub(2).gender = 0;  % 0=male, 1=female
cfg.sub(2).handedness = 1; % 0=left, 1=right
cfg.sub(2).preproc.corrupted_slices = [];

% Subject 3
cfg.sub(3).id = 'MR1036'; % subject id anonymized
cfg.sub(3).age = 20;
cfg.sub(3).gender = 0; % 0=male, 1=female
cfg.sub(3).handedness = 1; % 0=left, 1=right
cfg.sub(3).preproc.corrupted_slices = [];

% Subject 4
cfg.sub(4).id = 'MR940'; % subject id anonymized
cfg.sub(4).age = 21;
cfg.sub(4).gender = 0; % 0=male, 1=female
cfg.sub(4).handedness = 1; % 0=left, 1=right
cfg.sub(4).preproc.corrupted_slices = [];

% Subject 5
cfg.sub(5).id = 'MR1038'; % subject id anonymized
cfg.sub(5).age = 20;
cfg.sub(5).gender = 0; % 0=male, 1=female
cfg.sub(5).handedness = 1; % 0=left, 1=right
cfg.sub(5).preproc.corrupted_slices = [];

% Subject 6
cfg.sub(6).id = 'MR1043'; % subject id anonymized
cfg.sub(6).age = 22;
cfg.sub(6).gender = 0; % 0=male, 1=female
cfg.sub(6).handedness = 1; % 0=left, 1=right
cfg.sub(6).preproc.corrupted_slices = [];

% Subject 7
cfg.sub(7).id = 'MR1040'; % subject id anonymized
cfg.sub(7).age = 23;
cfg.sub(7).gender = 0; % 0=male, 1=female
cfg.sub(7).handedness = 1; % 0=left, 1=right
cfg.sub(7).preproc.corrupted_slices = [];

% Subject 8
cfg.sub(8).id = 'MR1058'; % subject id anonymized
cfg.sub(8).age = 24;
cfg.sub(8).gender = 1; % 0=male, 1=female
cfg.sub(8).handedness = 1; % 0=left, 1=right
cfg.sub(8).preproc.corrupted_slices = [];

% Subject 9
cfg.sub(9).id = 'MR1064'; % subject id anonymized
cfg.sub(9).age = 23;
cfg.sub(9).gender = 0; % 0=male, 1=female
cfg.sub(9).handedness = 1; % 0=left, 1=right
cfg.sub(9).preproc.corrupted_slices = [];

% Subject 10
cfg.sub(10).id = 'MR1028'; % subject id anonymized
cfg.sub(10).age = 18;
cfg.sub(10).gender = 0; % 0=male, 1=female
cfg.sub(10).handedness = 1; % 0=left, 1=right
cfg.sub(10).preproc.corrupted_slices = [];

% Subject 11
cfg.sub(11).id = 'MR1062'; % subject id anonymized
cfg.sub(11).age = 28;
cfg.sub(11).gender = 1; % 0=male, 1=female
cfg.sub(11).handedness = 1; % 0=left, 1=right
cfg.sub(11).preproc.corrupted_slices = [];

% Subject 12
cfg.sub(12).id = 'MR1072'; % subject id anonymized
cfg.sub(12).age = 28;
cfg.sub(12).gender = 1; % 0=male, 1=female
cfg.sub(12).handedness = 1; % 0=left, 1=right
cfg.sub(12).preproc.corrupted_slices = [];

% Subject 13
cfg.sub(13).id = 'MR1073'; % subject id anonymized
cfg.sub(13).age = 28;
cfg.sub(13).gender = 1; % 0=male, 1=female
cfg.sub(13).handedness = 1; % 0=left, 1=right
cfg.sub(13).preproc.corrupted_slices = [];

end

