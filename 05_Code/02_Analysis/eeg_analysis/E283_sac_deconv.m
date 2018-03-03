% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'TL_dc_Stim_500';
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
try
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
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'&origstart','<100';'&latposStim','>-1000'},...
                                [800 100],{-1,2,'origstart','>0'}); 
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.type        = cell(1,length(events.start));
    epochevents.type(events.type==1) = repmat({'fix'},1,sum(events.type==1));
    epochevents.type(events.type==2) = repmat({'sac'},1,sum(events.type==2));

    epochevents.orderposStim= events.orderposStim;
    epochevents.pxini       = (events.posinix-960)/45;            
    epochevents.pyini       = (events.posiniy-540)/45;     
    epochevents.pxend       = (events.posendx-960)/45;            
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;  
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini; 
    epochevents.pxdiff(2:2:end) = epochevents.pxdiff(1:2:end);         % fixation vector is the same sa previous saccade
    epochevents.pydiff(2:2:end) = epochevents.pydiff(1:2:end);
    epochevents.side        = nan(1,length(events.start));    
    epochevents.cross       = nan(1,length(events.start));    
    epochevents.inst        = nan(1,length(events.start)); 
    
    [trl,events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','>0'},[1500 900]);
    events                  = struct_elim(events,find(~ismember(events.value,[1:6,9,10,13,14,96])),2,0);
    epochevents.latency     = [epochevents.latency,events.time];
    ETttype                 = cell(1,length(events.value));
    ETttype(events.value==96)   = repmat({'image'},1,sum(events.value==96));
    ETttype(events.value<96)    = repmat({'stim'},1,sum(events.value<96));
    epochevents.type            = [epochevents.type,ETttype];
    ETside                      = nan(1,length(events.value));
    ETside(ismember(events.value,[1 3 5 9 13])) = -1;    % left
    ETside(ismember(events.value,[2 4 6 10 14])) = 1;
    epochevents.side            = [epochevents.side,ETside];
    ETcross                     = nan(1,length(events.value));
    ETcross(ismember(events.value,[1 2 9 10])) = -1;    % uncross
    ETcross(ismember(events.value,[3:6 13 14])) = 1;
    epochevents.cross           = [epochevents.cross,ETcross];
    ETinst                      = nan(1,length(events.value));
    ETinst(ismember(events.value,[1 2 5 6])) = 1;  % instructive
    ETinst(ismember(events.value,[9 10 13 14])) = -1;
    epochevents.inst            = [epochevents.inst,ETinst];
    
    epochevents.pxini           = [epochevents.pxini,nan(1,length(events.value))];    
    epochevents.pyini           = [epochevents.pyini,nan(1,length(events.value))];       
    epochevents.pxend           = [epochevents.pxend,nan(1,length(events.value))];          
    epochevents.pyend           = [epochevents.pyend,nan(1,length(events.value))];  
    epochevents.pxdiff          = [epochevents.pxdiff,nan(1,length(events.value))];  
    epochevents.pydiff          = [epochevents.pydiff,nan(1,length(events.value))];  
    epochevents.orderposStim    = [epochevents.orderposStim,nan(1,length(events.value))];  
    % getting the data in EEGlab format
    [EEG,winrej]            = getDataDeconv(cfg_eeg,epochevents,500); 
    mirindx                 = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
             
    if any(strfind(p.analysisname,'CI'))
        LstimTimes = epochevents.latency(strcmp(epochevents.type,'stim') & epochevents.side == -1);
        for tt  = 1:length(LstimTimes)
            EEGst = find(EEG.times>LstimTimes(tt),1);
            mirsamples = EEGst-floor(EEG.srate*.4):EEGst+ceil(EEG.srate*1.6);
            EEG.data(:,mirsamples) = EEG.data(mirindx,mirsamples);
            ETst  = find(epochevents.latency>LstimTimes(tt)-400 & epochevents.latency<LstimTimes(tt)+900 & ...
                (strcmp(epochevents.type,'sac') | strcmp(epochevents.type,'fix')));
            if ~isempty(ETst)
                [epochevents.pxini(ETst) epochevents.pxend(ETst) epochevents.pxdiff(ETst)] = deal(-epochevents.pxini(ETst),-epochevents.pxend(ETst),-epochevents.pxdiff(ETst));
            end
        end
    end
    if any(strfind(p.analysisname,'mirr'))
        EEG.data = EEG.data-EEG.data(mirindx,:);   
    end
    % deconvolution design
    cfgDesign           = [];
%      cfgDesign.eventtypes= {'sac','image','stim'};
%      cfgDesign.formula   = {'y~pxdiff+pydiff+pxend+pyend','y~1','y~cross*inst*side'};
%      model               = 'Sxydxyend_IM_STsci';
    cfgDesign.eventtypes= {'image','stim'};
    cfgDesign.formula   = {'y~1','y~cross*inst*side'};
    model               = 'IM_STsci';
    
    EEG                 = dc_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-.5,.5];tic
    EEG                 = dc_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = dc_continuousArtifactExclude(EEG,struct('winrej',winrej,'zerodata',0));
    EEG                 = dc_glmfit(EEG);
 
    unfold              = dc_beta2unfold(EEG);

    mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
    save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    clear unfold
catch
    tk
end
end 
%%
% %2nd level analysis
% clear
% E283_params                                 % basic experimental parameters               % 
% p.analysisname  = 'deconvTSCImirr';
 if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname', p.analysisname);
 end
    
stimB = [];
% model               = 'Sxydxyend_IM_STsci';
 model               = 'IM_STsci';
bslcor    = [];
for tk = p.subj
     cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    try
     load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    if ~isempty(bslcor) 
        stimB = cat(4,stimB,permute(unfold.beta-...
            repmat(mean(unfold.beta(:,find(unfold.times>bslcor(1) & unfold.times<bslcor(2)),:),2),1,size(unfold.beta,2),1),[1,3,2]));
    else
        stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
    end
    catch
        tk
      end
end
%
load(cfg_eeg.chanfile)
result      = regmodel2ndstat(stimB,unfold.times,elec,1000,'signpermT','cluster');
coeffs  = strrep({unfold.param.name},':','XX');
coeffs  = strrep(coeffs,'(','');
coeffs  = strrep(coeffs,')','');
coeffs = strcat([unfold.param.event]','_',coeffs');
result.coeffs = coeffs;
mkdir(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm'))
if ~isempty(bslcor) 
 save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALLbslcorr'),'result')
else
    save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALL'),'result')
end%
%%
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr']);
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
end
mkdir(pathfig)
plotinterval = [-.25  .0 .01;0 .25 .01];
setAbsoluteFigureSize
for b=1:size(result.B,2)
    betas.dof   = 1;
    betas.n     = size(stimB,4);
    betas.avg   = squeeze(mean(stimB(:,b,:,:),4));
    betas.time  = result.clusters(1).time;
    
    % topoplot across time according to interval with significant
    % clusters
    collim      =[-6*std(betas.avg(:)) 6*std(betas.avg(:))]; 
%     collim      =[-6 6]; 
    for pint = 1:size(plotinterval,1)
        fh       = topomitlines(cfg_eeg,result.clusters(b),betas,plotinterval(pint,:),collim);
        figsize  = [17.6 17.6*fh.Position(4)/fh.Position(3)];
        doimage(gcf,pathfig,'pdf',[result.coeffs{b} '_' strjoin('_',{num2str(plotinterval(pint,1)),num2str(plotinterval(pint,2))})],figsize,1)
   
    end
end

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';
% cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
   
load(cfg_eeg.chanlocs)
timeLim         = [-.3 .7];
% chnstoPlot      = {{'FC1','FCz','FC2','C1','Cz','C2','CP1','CPz','CP2'},{'O1','Oz','O2'},...
 %    {'P5','P7','PO7'},{'P6','P8','PO8'}};
 chnstoPlot      = {{'C3','C5'},{'C4','C6'},{'P5','P7','PO7'},{'P6','P8','PO8'}};
lineColors      = cbrewer('qual','Set1',9);%[134 16 9;22 79 134;11 93 24]/255;
Bnames        = {'stim_2_Intercept'};
xaxis           = result.clusters(1).time;
axLim           = [-.2 .4 -2 2];


for ch = 1:length(chnstoPlot)
    auxChns         = find(ismember({chanlocs.labels},chnstoPlot{ch}));
    datatoPlot  = [];
    for b = 1:length(Bnames)
        ixB = strmatch(Bnames{b},result.coeffs,'exact');
        datatoPlot  = cat(3,datatoPlot,permute(squeeze(mean(result.B(auxChns,ixB,:,:))),[2 1]));
    end
    lineNames   = strtok(Bnames,'_');
%         result.clusters(ixB).mask
%         signf        = squeeze(any(any(statUCIci.mask(auxChns,freqs,:),1),2))';
    fh = fillPlot(datatoPlot,[],xaxis,axLim,'mean',lineColors,lineNames);
%     fh.Name = Bnames{b};
    xlabel('Time (s)','FontSize',12)
    ylabel('\muV','Interpreter','tex','FontSize',12)
    figsize     = [8 8*fh.Position(4)/fh.Position(3)];
%     doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_Inf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
end