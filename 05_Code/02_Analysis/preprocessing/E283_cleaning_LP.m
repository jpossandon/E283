%clean_path
%%
<<<<<<< HEAD:code/analysis/preprocessing/E283_cleaning_LP.m
MACpath = '/Volumes/nibaldo/trabajo/E283/';
%suj = str2num(getenv('SGE_TASK_ID'));
for suj             = [43,44 ,59];
=======
MACpath = '/Users/jossando/trabajo/E283/';
% datafileStruct

for suj             = [68,69,71,72];
>>>>>>> ebb813b44ef59b5af415b489a9aba4243dc36d94:05_Code/02_Analysis/preprocessing/E283_cleaning_LP.m
eegfilename     = sprintf('s%02dvs',suj);
suj             = sprintf('s%02dvs',suj);
if ismac    
    cfg             = eeg_etParams_E283('sujid',suj,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
else
      cfg             = eeg_etParams_E283('sujid',suj,'expfolder','C:\Users\jpo\trabajo\E283\');
end
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean eye correction from this specific file
% pause(rand(1)*30)   % to get rid of that random error that seem to be cause
% by many computers accesing and saving the same file at the same moment 

clean_channel_corrections(cfg,eegfilename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% correct channels (inversions) (this is inside all code)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a first sweep to remove absurdely bad segments independently of eye
% movements
cfg             = eeg_etParams_E283(cfg,...
                                   'task_id','fv_touch',...
                                   'filename',eegfilename,...
                                   'event',[eegfilename '.vmrk'],...
                                   'analysisname','cleaning',...
                                   'clean_exclude_eye',0,...
                                   'clean_foi',20:2:42,...
                                   'clean_freq_threshold',600,...
                                   'clean_range_threshold',[5 600],...
                                   'clean_ica_correct',[],...
                                   'clean_trend_threshold',200,...
                                   'clean_minclean_interval',100,...
                                   'clean_movwin_length',.256,...
                                   'clean_mov_step',.006,...
                                   'lpfreq',42,...
                                   'ratio20thresh',1,...
                                    'clean_name','pre');
                                
% bad because of gamma
[value_g,tot_sample_g]              = freq_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_g,badchans_g,channelbad_g]     = badsegments(cfg,value_g,tot_sample_g,cfg.clean_freq_threshold);
% 
% % bad because of amplitude
[value_a,tot_sample_a]              = range_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_a,badchans_a,channelbad_a]     = badsegments(cfg,value_a,tot_sample_a,cfg.clean_range_threshold);
% bad = bad_a; badchans = badchans_a; channelbad = channelbad_a;
[bad,badchans]                      = combine_bad({bad_a;bad_g},{badchans_a,badchans_g},cfg.clean_minclean_interval);
channelbad                          = combine_bad([channelbad_a;channelbad_g],[],cfg.clean_minclean_interval);

cfg_clean                           = cfg;

% save([cfg.analysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','channelbad','cfg_clean','value_g','tot_sample_g','value_a','tot_sample_a') % we save all this for now until we now the good seetings
save([cfg.preprocanalysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','channelbad','cfg_clean') % TODO: info about the cleaning parameters

% here we check if there is a channel that is continuously bad even for the
% lax criteria were using
cfg                                 = eeg_etParams_E283(cfg,'clean_bad_channel_criteria',.15,...
                                                                'clean_bad_channel_denominator',{'S 96',[5000 17000]});
check_session(cfg)

% re-check bad segments now whitout taking in accound bad channels (done in badsegments and check_session)
 cfg                                = eeg_etParams_E283(cfg,'clean_movwin_length',.256,'clean_mov_step',.006);
[bad_g,badchans_g]                  = badsegments(cfg,value_g,tot_sample_g,cfg.clean_freq_threshold);
[bad_a,badchans_a]                  = badsegments(cfg,value_a,tot_sample_a,cfg.clean_range_threshold);
% bad = bad_a; badchans = badchans_a;
[bad,badchans]                      = combine_bad({bad_a;bad_g},{badchans_a,badchans_g},cfg.clean_minclean_interval);
save([cfg.preprocanalysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','bad_a','badchans_a','-append') % TODO: info about the cleaning parameters

% run first ICA
expica(cfg)

%%
% the second sweep is over data clean from eye-movement compontent and muscle artifact components, we can
% use narrower thresholds here, we check visually that everything is ok and then we run ICA again
clear cfg
if ismac    
    cfg             = eeg_etParams_E283('sujid',suj,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
else
   cfg             = eeg_etParams_E283('sujid',suj,'expfolder','C:\Users\jpo\trabajo\E283\');
end

cfg             = eeg_etParams_E283(cfg,...
                                   'task_id','fv_touch',...
                                   'filename',eegfilename,...
                                   'event',[eegfilename '.vmrk'],...
                                   'analysisname','cleaning',...
                                   'clean_exclude_eye',0,...
                                   'clean_foi',20:2:42,...
                                   'clean_freq_threshold',300,...
                                   'clean_range_threshold',[1 125],...
                                   'clean_trend_threshold',70,...
                                   'clean_minclean_interval',500,...
                                   'clean_ica_correct','yes',...
                                   'clean_movwin_length',.256,...
                                   'clean_mov_step',.006,...
                                   'lpfreq',42,...
                                   'ratio20thresh',1,...
                                   'clean_name','general');

% bad because of gamma
[value_g,tot_sample_g]              = freq_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_g,badchans_g,channelbad_g]     = badsegments(cfg,value_g,tot_sample_g,cfg.clean_freq_threshold);
% 
% % bad because of amplitude
[value_a,tot_sample_a]              = range_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_a,badchans_a,channelbad_a]     = badsegments(cfg,value_a,tot_sample_a,cfg.clean_range_threshold);
% 
% % bad because of trend, longer
 cfg                                = eeg_etParams_E283(cfg,'clean_movwin_length',1,...
                                          'clean_mov_step',.06);
[value,tot_sample]                  = trend_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_t,badchans_t,channelbad_t]     = badsegments(cfg,value,tot_sample,cfg.clean_trend_threshold); %TODO: borders need to be readjusted to hthe length of gthe moving window
%                    
% % combine info and save,[channelbad_a;channelbad_g;channelbad_t]
[bad,badchans]                      = combine_bad({bad_a;bad_g;bad_t},{badchans_a,badchans_g,badchans_t},cfg.clean_minclean_interval);
channelbad                          = combine_bad([channelbad_a;channelbad_g;channelbad_t],[],cfg.clean_minclean_interval);

cfg_clean                           = cfg;

% save([cfg.analysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','channelbad','cfg_clean','value_g','tot_sample_g','value_a','tot_sample_a') % we save all this for now until we now the good seetings
save([cfg.preprocanalysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','channelbad','cfg_clean') % TODO: info about the cleaning parameters

cfg                                 = eeg_etParams_E283(cfg,'clean_bad_channel_criteria',.10,...
                                                                'clean_bad_channel_denominator',{'S 96',[5000 17000]}); % to take in account only relevant peaces, 5 second before and after the max length of a trial (12s) 
check_session(cfg)

% re-check bad segments now whitout taking in accound bad channels (done in badsegments)
cfg                                = eeg_etParams_E283(cfg,'clean_movwin_length',.256,'clean_mov_step',.006);
[bad_g,badchans_g]                 = badsegments(cfg,value_g,tot_sample_g,cfg.clean_freq_threshold);
[bad_a2,badchans_a2]                 = badsegments(cfg,value_a,tot_sample_a,cfg.clean_range_threshold);
cfg                                = eeg_etParams_E283(cfg,'clean_movwin_length',1,'clean_mov_step',.06);
[bad_t,badchans_t]                 = badsegments(cfg,value,tot_sample,cfg.clean_trend_threshold); %TODO: borders need to be readjusted to hthe length of gthe moving window

% we reuse the very bad segments becuase of range of step1 to get rid of
% problems with some dataset where there is a very bad segment that gets
% celan by ICA
load([cfg.preprocanalysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename 'pre'],'bad_a','badchans_a') % TODO: info about the cleaning parameters

[bad,badchans]                     = combine_bad({bad_a;bad_a2;bad_g;bad_t},{badchans_a,badchans_a2,badchans_g,badchans_t},cfg.clean_minclean_interval);
save([cfg.preprocanalysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','-append') % TODO: info about the cleaning parameters

%%
%        cfgvis             = eeg_etParams_E283(cfg,...
%                                              'remove_eye',0,...
%                                              'remove_m',0,'raw',1); 
%       visual_clean(cfgvis)
%% 
% run second definitive ICA
 expica(cfg)

%%
% and we do the cleaning one more time with the final ICA weights
clear cfg
if ismac    
    cfg             = eeg_etParams_E283('sujid',suj,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
else
   cfg             = eeg_etParams_E283('sujid',suj,'expfolder','C:\Users\jpo\trabajo\E283\');
end

cfg             = eeg_etParams_E283(cfg,...
                                   'task_id','fv_touch',...
                                   'filename',eegfilename,...
                                   'event',[eegfilename '.vmrk'],...
                                   'analysisname','cleaning',...
                                   'clean_exclude_eye',0,...
                                   'clean_foi',20:2:42,...
                                   'clean_freq_threshold',300,...
                                   'clean_range_threshold',[1 125],...
                                   'clean_trend_threshold',70,...
                                   'clean_minclean_interval',500,...
                                   'clean_ica_correct','yes',...
                                   'clean_movwin_length',.256,...
                                   'clean_mov_step',.006,...
                                   'lpfreq',42,...
                                   'ratio20thresh',1,...
                                   'clean_name','final');

% bad because of gamma
[value_g,tot_sample_g]              = freq_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_g,badchans_g,channelbad_g]     = badsegments(cfg,value_g,tot_sample_g,cfg.clean_freq_threshold);
% 
% % bad because of amplitude
[value_a,tot_sample_a]              = range_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_a,badchans_a,channelbad_a]     = badsegments(cfg,value_a,tot_sample_a,cfg.clean_range_threshold);
% 
% % bad because of trend, longer
 cfg                                = eeg_etParams_E283(cfg,'clean_movwin_length',1,...
                                          'clean_mov_step',.06);
[value,tot_sample]                  = trend_artifact(cfg,'lpfilter','yes','lpfreq',cfg.lpfreq);
[bad_t,badchans_t,channelbad_t]     = badsegments(cfg,value,tot_sample,cfg.clean_trend_threshold); %TODO: borders need to be readjusted to hthe length of gthe moving window
%                    
% % combine info and save,[channelbad_a;channelbad_g;channelbad_t]
[bad,badchans]                      = combine_bad({bad_a;bad_g;bad_t},{badchans_a,badchans_g,badchans_t},cfg.clean_minclean_interval);
channelbad                          = combine_bad([channelbad_a;channelbad_g;channelbad_t],[],cfg.clean_minclean_interval);

cfg_clean                           = cfg;
save([cfg.preprocanalysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','channelbad','cfg_clean')  % TODO: info about the cleaning parameters
% save([cfg.analysisfolder cfg.analysisname '/' cfg.sujid '/' cfg.filename cfg.clean_name],'bad','badchans','channelbad','cfg_clean','value_g','tot_sample_g','value_a','tot_sample_a','value','tot_sample')  % TODO: info about the cleaning parameters
end
% 

% so there is the problem of very bad periods during eye-movement 
% (possible solutions: embracing periods or ahigher non eye movement
% threshold)
% Can eegplot show whithin periods different causes?

% tag bad channels (bad>25% time within experimental period)

% save preinfo

% 1st ICA with corrected channels (dead_crazy) 
% removal of eye movements components
% bad segment with eye movements
% visual inspection
% 2nd ICA
% removal of eye movements components nad other artifactual components
% bad segmetns with criteria for ERP and TFR
% need to decide other components to remove or mark like heart beats
