% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'deconvTS';
%%
% subject configuration and data
 
if ismac 
run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')        
else
run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')   
end    
%p.subj              = [7,11,12,15];
for tk = p.subj
    tk
%  for tk = p.subj;
% tk = str2num(getenv('SGE_TASK_ID'));
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jpo/trabajo/E283/');
    end
    
    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname',p.analysisname);    % single experiment/session parameters 
   
    % get relevant epochevents
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata               
    [trl,events]           = define_event(cfg_eeg,eyedata,2,{'&origstart','>0';'&origstart','<7000'},...
                                [800 100],{-1,1,'origstart','>0'}); 
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.type        = cell(1,length(events.start));
    epochevents.type(events.type==1) = repmat({'fix'},1,sum(events.type==1));
    epochevents.type(events.type==2) = repmat({'sac'},1,sum(events.type==2));

    epochevents.pxini       = (events.posinix-960)/45;            
    epochevents.pyini       = (events.posiniy-540)/45;     
    epochevents.pxend       = (events.posendx-960)/45;            
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;  
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini; 
    epochevents.side        = nan(1,length(events.start));    
    epochevents.cross       = nan(1,length(events.start));    
    epochevents.inst        = nan(1,length(events.start)); 
    
    [trl,events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','>0'},[1500 1000]);
    events                  = struct_elim(events,find(~ismember(events.value,[1:6,9,10,13,14,96])),2,0);
    epochevents.latency     = [epochevents.latency,events.time];
    ETttype                 = cell(1,length(events.value));
    ETttype(events.value==96)   = repmat({'image'},1,sum(events.value==96));
    ETttype(events.value<96)    = repmat({'stim'},1,sum(events.value<96));
    epochevents.type            = [epochevents.type,ETttype];
    ETside                  = nan(1,length(events.value));
    ETside(ismember(events.value,[1 3 5 9 13])) = -1;    % left
    ETside(ismember(events.value,[2 4 6 10 14])) = 1;
    epochevents.side        = [epochevents.side,ETside];
    ETcross                 = nan(1,length(events.value));
    ETcross(ismember(events.value,[1 2 9 10])) = -1;    % uncross
    ETcross(ismember(events.value,[3:6 13 14])) = 1;
    epochevents.cross       = [epochevents.cross,ETcross];
    ETinst                 = nan(1,length(events.value));
    ETinst(ismember(events.value,[1 2 5 6])) = 1;  % instructive
    ETinst(ismember(events.value,[9 10 13 14])) = -1;
    epochevents.inst       = [epochevents.inst,ETinst];
    
    epochevents.pxini       = [epochevents.pxini,nan(1,length(events.value))];    
    epochevents.pyini       = [epochevents.pyini,nan(1,length(events.value))];       
    epochevents.pxend       = [epochevents.pxend,nan(1,length(events.value))];          
    epochevents.pyend       = [epochevents.pyend,nan(1,length(events.value))];  
    epochevents.pxdiff      = [epochevents.pxdiff,nan(1,length(events.value))];  
    epochevents.pydiff      = [epochevents.pydiff,nan(1,length(events.value))];  
  
    % getting the data in EEGlab format
    [EEG,winrej] = getDataDeconv(cfg_eeg,epochevents,100);  

    % deconvolution design
    cfgDesign           = [];
    cfgDesign.eventtype = {'fix','sac','image','stim'};
    cfgDesign.formula   = {'y ~pxini+pyini','y~pxdiff+pydiff','y~1','y~ side*cross*inst'};
    model               = 'Fxy_Sxdyd_IM_STsc';
    EEG                 = dc_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-.75,.75];tic
    EEG                 = dc_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = dc_continuousArtifactExclude(EEG,struct('winrej',winrej,'zerodata',0));
    EEG                 = dc_glmfit(EEG);
 
    unfold              = dc_beta2unfold(EEG);
    for nep = 1:length(unfold.epoch)
        if iscell(unfold.epoch(nep).event)
            unfold.epoch(nep).event = cell2mat(unfold.epoch(nep).event);
        end
    end
    
  % ploting beta averages
     mkdir(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,model,'figures_subjects',cfg_eeg.sujid))
     
     B = unfold.beta;
     collim = [-.2 .2];
     p.coeff = strrep({unfold.epoch.name},':','_');
     p.coeff = strrep(p.coeff,'(','');
     p.coeff = strrep(p.coeff,')','');
     etype   = {unfold.epoch.event};
     for b = 1:size(B,3);
        betas.dof   = 1;
        betas.n     = 1;
        betas.avg   = permute(B(:,:,b),[1,3,2]);
        collim      = [-6*std(betas.avg(:)) 6*std(betas.avg(:))]; 
        
        betas.time      = unfold(1).times; 
        auxresult.time  =  unfold(1).times;
        fh = plot_stat(cfg_eeg,auxresult,betas,[],[-.64 .64 .02],collim,.05,sprintf('Beta: %s %s',strrep(p.coeff{b},'_',' | '),etype{b}),1);
           doimage(fh,fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,model,'figures_subjects',cfg_eeg.sujid),'png',...
                [datestr(now,'ddmmyy') cfg_eeg.sujid '_'  etype{b} '_' p.coeff{b}],1)
     end

      mkdir(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,model,'glm'))
      save(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
      clear unfold
end 
%%
% %2nd level analysis
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'deconvTSmirr';
 if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', 'deconvTS'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname', 'deconvTS');
 end
    
stimB = [];
model               = 'Fxy_Sxdyd_IM_STsc';
for tk = p.subj
     cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    load(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    if any(strfind(p.analysisname,'mirr'))
        mirindx         = mirrindex({unfold.chanlocs.labels},[cfg_eeg.expfolder '/channels/mirror_chans']); 
        stimB = cat(4,stimB,permute(unfold.beta-unfold.beta(mirindx ,:,:),[1,3,2]));
    else
        stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
    end
end
load(cfg_eeg.chanfile)
result      = regmodel2ndstat(stimB,unfold.times,elec,1000,'signpermT','cluster');

mkdir(fullfile(cfg_eeg.analysisfolder,p.analysisname ,model,'glm'))
save(fullfile(cfg_eeg.analysisfolder,p.analysisname ,model,'glm','glmALL'),'result')
%%
interval = [-.64 .64 .02];

pathfig = fullfile(cfg_eeg.analysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
coeffs  = strrep({unfold.epoch.name},':','XX');
coeffs  = strrep(coeffs,'(','');
coeffs  = strrep(coeffs,')','');
coeffs = strcat({unfold(1).epoch.event}','_',coeffs');

glm_betaplots(cfg_eeg,stimB,result,interval,pathfig,coeffs)

