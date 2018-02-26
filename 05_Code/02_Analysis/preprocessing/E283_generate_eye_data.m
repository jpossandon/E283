function E283_generate_eye_data(subjects)
% generate data experiment touch
cfg.edfreadpath     = '/home/th/code/edfread/build/linux64/';
cfg.eyes            = 'monoocular';
for s = subjects
    cfg.filename        = sprintf('s%02dvs',s);
%     cfg.EDFfolder       = sprintf('/home/th/Experiments/E275/data/s%02dvs/',s);
%     cfg.EDFfolder       = sprintf('C:\\Users\\jpo\\trabajo\\E283\\data\\s%02dvs\\',s);
     cfg.EDFfolder       = sprintf('/Users/jossando/trabajo/E283/06_RawData/s%02dvs/',s);
    eyedata             = eyeread2(cfg); 
    auxdata             = eyedata.events;
    
    auxdata.block       = zeros(1,length(auxdata.start));
    auxdata.cue         = zeros(1,length(auxdata.start));
    auxdata.target      = zeros(1,length(auxdata.start));
    auxdata.blockstart  = zeros(1,length(auxdata.start));
    auxdata.tpos        = zeros(1,length(auxdata.start));
    auxdata.stim        = zeros(1,length(auxdata.start));
    auxdata.latestim    = zeros(1,length(auxdata.start));
    
    % actual stimulations
    auxstim             = eyedata.marks.value(strcmp(eyedata.marks.type,'ETtrigger'));
    auxstimtime         = eyedata.marks.time(strcmp(eyedata.marks.type,'ETtrigger'));
    auxstimtrial        = eyedata.marks.trial(strcmp(eyedata.marks.type,'ETtrigger'));
%     stimdata.value      = auxstim(auxstim<5 & auxstim>0);  % take in account only stimulation start and only stimulation during the trial and not at the start (values 1,2,3 - left, right and bilateral respectively; value 10 is for stim stop and there is NaNs when there was nono stimulation during the complete trial)(initial stimulation is always at time 150 or 151)
%     stimdata.time       = auxstimtime(auxstim<5 & auxstim>0); %
     stimdata.trial      = auxstimtrial(auxstim<15 & auxstim>0);
%     stimdata.subject    = s*ones(1,length(stimdata.trial));
%     sampledata          = eyedata.samples;
%     sampledata.subject  = s*ones(1,length(sampledata.time));
    for t = unique(eyedata.events.trial)
        % which block type (0 - uncrossed ; 1 - crossed)
        auxvalue        = eyedata.marks.value(strcmp(eyedata.marks.type,'block') & eyedata.marks.trial == t);
        auxdata.block(auxdata.trial == t) = auxvalue;
        % was block start (0 - no ; 1 - yes)
        auxvalue        = eyedata.marks.value(strcmp(eyedata.marks.type,'block_start') & eyedata.marks.trial == t);
        auxdata.blockstart(auxdata.trial == t) = auxvalue;
        % which cue (0 - uniformative; 1 - informative)
        auxvalue        = eyedata.marks.value(strcmp(eyedata.marks.type,'cue') & eyedata.marks.trial == t);
        auxdata.cue(auxdata.trial == t) = auxvalue;
        % which tpos (64 different numbers)
        auxvalue        = eyedata.marks.value(strcmp(eyedata.marks.type,'tpos') & eyedata.marks.trial == t);
        auxdata.tpos(auxdata.trial == t) = auxvalue;
        % if there is stimulation during the trial (0 - no; 1 - yes) , as there can be more than one stimulation that information will be in the stimdata structure       
         if ~isempty(find(stimdata.trial==t))
             auxdata.latestim(auxdata.trial == t) = 1;
         end
    end
    auxdata.subject     = s*ones(1,length(auxdata.trial));
    eyedata.events      = auxdata;
   % data                = struct_up('data',auxdata,2);
   % stim                = struct_up('stim',stimdata,2);
   % sample              = struct_up('sample',sampledata,2);
       save(sprintf('%ss%02dvseye',cfg.EDFfolder,s),'eyedata') 
end
     
% save('/net/store/nbp/touch/data/alleyedata','data','stim')
% save('/net/store/nbp/touch/data/allsampledata','sample')
    