function cfg = eeg_etParams_275(varargin)

if ~isstruct(varargin{1})
    % default parameter for experiment EEG features
    cfg.expname                 = 'E283';                         % used for eeg population stats and meta-analysis across experiments, keep the same name for all analysis done with the same experimental data

    % paths
%     cfg.expfolder               = '/Users/jossando/trabajo/E275/';
    cfg.expfolder               = '/net/store/nbp/projects/EEG/E283/'    
    cfg.edfreadpath             = '/net/store/users/jossando/edfread/build/linux64/';
   
    % triggers and trial definition
    cfg.trial_trig_eeg          = {'S  1','S  2','S  5','S  6','S  9','S 10','S 13','S 14','S 96'};                          % 'S  1','S  2','S  3', vm task eeg triggers 
    cfg.trial_trig_et           = {'1','2','5','6','9','10','13','14','96'};                                % vm task EDF triggers
    cfg.trial_time              = [0 10000];                                    % experiment trial window

    % channel that need to be remove or switch
    cfg.correct_chan            = [];                                           % [neworder] if channels need to be re-arrenged (for example if electrodes were connected wrongly to the amps)
    cfg.elim_chan               = [];
    
    % cleaning data
    cfg.clean_ica_correct       = 'no'; 
    cfg.clean_name              = 'pre';
    cfg.clean_movwin_length     = 0.256;
    cfg.clean_mov_step          = 0.005;
    cfg.clean_exclude_eye       = 0;
    cfg.clean_minclean_interval = 200;
    
    cfg.clean_foi               = 30:5:120;
    cfg.clean_freq_threshold    = 30;
    cfg.clean_range_threshold   = 350;
    cfg.clean_trend_threshold   = .2;
    
    cfg.clean_bad_channel_criteria  = .5;
%     cfg.artifact_reject         = 'continous';                                  % trial or continous data 
%     cfg.artifact_reject_type    = 'automatic';                                  % automatic or visual
%     cfg.artifact_chunks_length  = 512;                                          % in ms, defines size of non-overlapping chunks for continous automatic artifact rejection
%     cfg.artifact_auto_stat      = 'yes';
%     cfg.thresholds              = 'subj';                                        % 'pop' - bases in overall population statistics 'subj' = based only in subject statistics 'file' = base only in current file statistics
%     cfg.thresholds_otherexp     = 'no';                                         % 'yes' - use also data from other eeg-eyetracking experiment to determine therhsolds, 'no', only data from current exp
%     cfg.datastats               = 'all';  
%     
    % ICA
%     cfg.ica_it                  = [];                                           % How expica.m do ICA: [] - once,  # - number of std to define outlier trial, if there are outlier trials after one ICA it perform a new ICA iteration
%     cfg.ica_type                = 'exp';                                        % ICA weigths obtained with data from: 'preexp' - preexperiment , 'exp' - experiment data, 'both' 
    cfg.ica_data                = 'all';
    cfg.ica_chunks_length       = 512;
    % use ICA and overall eeg stats for finding artifact and artifact component
    % in other experiments
%     cfg.keepforotherexp         = 'no' ;                                        % saves ICA weigths for artifact correction of current experiment in a common folder so they can be used for component selection by topography clustering (not implemented yet)  
%     cfg.trim_percent            = 1;

    % eyedata
    cfg.imagefield              = 'image';                                      % name of the image field in EDF data
    cfg.eyedata                 = 'yes';                                        %there is an EDF file associated
    cfg.conditionfield          = {'ETtrigger',1};                              % field that define trial condition. 'ETtrigger' correspond to the time and data send with the eye-tracker to eeg trigger 
    cfg.resolution              = 41;                               % image resolution in pixels per degree
    cfg.recalculate_eye         = 'no';                                         % 'yes' to recalculate eye movements with Engbert algorithm (or with a fixed velocity threshold)
    cfg.eyes                    = 'monocular';


    % analysis
    cfg.analysisname            = 'cleaning';                                % name of the analysis, results are saved in corresponding analysis/* folder. preanalysis folder is used for the files that contain the artifact segments that are going to be elimined

    cfg.sujid                   = '1';                                        % id of the subject as how is in EdF and EEG files names
    % cfg.subjects                = [9:10,12:14,16:24,26:28, 30, 32:39];                          % subjects in the experiment

    % cfg.session                 = 1;
   % cfg.task_id                 = 'vm';
   % cfg.filename                 = 'jo01vm01';
    % cfg.task_num                = 1;
    for nv = 1:2:length(varargin)                                           % redefintion of eeg_etParams
        cfg.(varargin{nv}) = varargin{nv+1};
    end
else
    cfg = varargin{1};
    for nv = 2:2:length(varargin)                                           % redefintion of eeg_etParams
        cfg.(varargin{nv}) = varargin{nv+1};
    end
end
cfg.datapath            = [cfg.expfolder 'data/'];
cfg.chanloc             = [cfg.expfolder 'channels/easycapM1E275.txt'];              % text file with positions
cfg.chanlocs            = [cfg.expfolder 'channels/chanlocseasycapM1E275'];             % eeglab ready format
cfg.chanfile            = [cfg.expfolder 'channels/eleceasycapM1E275'];
cfg.analysisfolder      = [cfg.expfolder 'analysis/'];

%cfg.masterfile           = [cfg.expfolder 'subjects_master_files/' upper(cfg.sujid) 'wc'];

cfg.EDFfolder           = [cfg.datapath cfg.sujid  '/'];
%cfg.xensor              = [cfg.datapath upper(cfg.sujid) '/xensor/'];
cfg.eegfolder           = [cfg.datapath cfg.sujid '/'];
cfg.matfolder           = [cfg.datapath cfg.sujid '/'];
cfg.eyeanalysisfolder   = [cfg.analysisfolder 'eyedata/' cfg.sujid '/'];
cfg.eegstats            = [cfg.analysisfolder 'eeg_stats/' cfg.sujid '/'];
%cfg.channelcorfolder    = [cfg.expfolder 'subjects_master_files/'];
% cfg.tasks_info          = [    cfg.expfolder 'subjects_master_files/tasks_triggers.mat']; 

% create the folders if they don't exist
if ~isdir(cfg.eyeanalysisfolder), mkdir(cfg.eyeanalysisfolder),end
% if ~isdir(cfg.eeganalysisfolder), mkdir(cfg.eeganalysisfolder),end
if ~isdir(cfg.analysisfolder), mkdir(cfg.analysisfolder),end
%     if ~isdir([cfg.analysisfolder 'expstats']), mkdir([cfg.analysisfolder 'expstats']),end
if ~isdir([cfg.analysisfolder cfg.analysisname]), mkdir([cfg.analysisfolder cfg.analysisname]),end
if ~isdir([cfg.analysisfolder 'ICAm/' cfg.sujid]), mkdir([cfg.analysisfolder 'ICAm/' cfg.sujid]),end
if ~isdir([cfg.analysisfolder 'cleaning/' cfg.sujid]), mkdir([cfg.analysisfolder 'cleaning/' cfg.sujid]),end



