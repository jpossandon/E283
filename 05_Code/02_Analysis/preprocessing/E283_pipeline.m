% E283 pipeline

%%
% Read the eyetracker file
clear
MACpath = '/Volumes/nibaldo/trabajo/E283/';
MACpath = '/Users/jossando/trabajo/E283/';
oldsubjects             = [2,4,5,6,7,8,10,11,12,13,14,15,17,18,19,22,20,21,23,24,25,26,27,28,29,31,41:47,50,51,52,53,54,55,56,57,58,59,60,62,63,64,66];
     
addsubjects             = [68,69,71,72];
% eeg                     = [0,0,1,1];
for s=addsubjects
    eegfilename     = sprintf('s%02dvs',s);
    suj             = sprintf('s%02dvs',s);

    if ismac    
        cfg             = eeg_etParams_E283('sujid',suj,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
    else
        cfg             = eeg_etParams_E283('sujid',suj,'expfolder','C:\Users\jpo\trabajo\E283\');
    end
    cfg             = eeg_etParams_E283(cfg,'sujid',suj,...%'expfolder','/net/store/nbp/projects/EEG/E275/',...      % to run things in different environments
                                    'task_id','fv_touch',...
                                    'filename',eegfilename,...
                                    'event',[eegfilename '.vmrk'],...
                                    'trial_trig_eeg',{'S 96'},...
                                    'trial_trig_et',{'96'});      % experiment parameters 
    if ~isempty(oldsubjects)
        load([cfg.analysisfolder '03_Eye/eyedata/alleyedata'],'data','stim')
        load([cfg.analysisfolder '03_Eye/eyedata/allsampledata'],'sample')
    end                            
    if ~ismember(s,oldsubjects)
        load([cfg.EDFfolder suj 'eye.mat'])
                
        auxstim             = eyedata.marks.value(strcmp(eyedata.marks.type,'ETtrigger'));
        auxstimtime         = eyedata.marks.time(strcmp(eyedata.marks.type,'ETtrigger'));
        auxstimtrial        = eyedata.marks.trial(strcmp(eyedata.marks.type,'ETtrigger'));
        stimdata.value      = auxstim(auxstim<15 & auxstim>0);  % take in account only stimulation start and only stimulation during the trial and not at the start (values 1,2,3 - left, right and bilateral respectively; value 10 is for stim stop and there is NaNs when there was nono stimulation during the complete trial)(initial stimulation is always at time 150 or 151)
        stimdata.time       = auxstimtime(auxstim<15 & auxstim>0); %
        stimdata.trial      = auxstimtrial(auxstim<15 & auxstim>0);
        stimdata.subject    = s*ones(1,length(stimdata.time));
        sampledata          = eyedata.samples;
        sampledata.subject  = s*ones(1,length(sampledata.time));
        data                = struct_up('data',eyedata.events,2);
        stim                = struct_up('stim',stimdata,2);
        sample              = struct_up('sample',sampledata,2);
        
        eyedata         = synchronEYEz(cfg, eyedata);   
        save(sprintf('%s%seye',cfg.eyeanalysisfolder,suj),'eyedata')
    end
      clear stimdata sampledata eyedata
      save([cfg.analysisfolder '03_Eye/eyedata/alleyedata'],'data','stim')
save([cfg.analysisfolder '03_Eye/eyedata/allsampledata'],'sample')
end

%%
% to fix eyedata files 
% eyedata.marks = struct_elim(eyedata.marks,find(eyedata.marks.trial<27),2,0)
% eyedata.samples = struct_elim(eyedata.samples,find(eyedata.samples.trial<27),2,0)
% eyedata.events = struct_elim(eyedata.events,find(eyedata.events.trial<27),2,0)
