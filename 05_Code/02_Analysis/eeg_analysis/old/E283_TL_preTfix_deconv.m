% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'TL_dc_preT';
if ismac 
run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')        
else
run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')   
end  
%%
% subject configuration and data
 
  
%p.subj              = [7,11,12,15];
for tk =p.subj
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
    eyedata.events.nextToTarg = any([eyedata.events.nextToTargH;eyedata.events.nextToTargV;eyedata.events.nextToTargD]);
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'orderPreT','<6'},...
                                [800 100]);
    
   
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.dur         = events.dur; 
    epochevents.type        = cell(1,length(events.start));
%     epochevents.prefixdur   = events.prefixdur; 
  
   %ADD NEXT SACCADE DURATION
      epochevents.type(events.type==1 & events.orderPreT==0) = repmat({'fixpreT0'},1,sum(events.type==1 & events.orderPreT==0));
      epochevents.type(events.type==1 & events.orderPreT==1 & events.onTarg)    = repmat({'fixpreT1onTarg'},1,sum(events.type==1  & events.orderPreT==1 & events.onTarg));
      epochevents.type(events.type==1 & events.orderPreT==1 & ~events.onTarg)   = repmat({'fixpreT1'},1,sum(events.type==1  & events.orderPreT==1 & ~events.onTarg));
      epochevents.type(events.type==1 & events.orderPreT==2)   = repmat({'fixpreT2'},1,sum(events.type==1  & events.orderPreT==2));
      epochevents.type(events.type==1 & events.orderPreT==3)   = repmat({'fixpreT3'},1,sum(events.type==1  & events.orderPreT==3 ));
     epochevents.type(events.type==1 & events.orderPreT>3)   = repmat({'fixpre'},1,sum(events.type==1  & events.orderPreT>3));
    

    epochevents.orderposStim= events.orderposStim;
    epochevents.pxini       = (events.posinix-960)/45;            
    epochevents.pyini       = (events.posiniy-540)/45;     
    epochevents.pxend       = (events.posendx-960)/45;            
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;  
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini;
    events.angle(events.angle<0) = 360+events.angle(events.angle<0);           % transform angles to 0-360
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
    
    epochevents.dur       = [epochevents.dur,nan(1,length(events.value))];
    epochevents.pxini       = [epochevents.pxini,nan(1,length(events.value))];    
    epochevents.pyini       = [epochevents.pyini,nan(1,length(events.value))];       
    epochevents.pxend       = [epochevents.pxend,nan(1,length(events.value))];          
    epochevents.pyend       = [epochevents.pyend,nan(1,length(events.value))];  
    epochevents.pxdiff      = [epochevents.pxdiff,nan(1,length(events.value))];  
    epochevents.pydiff      = [epochevents.pydiff,nan(1,length(events.value))];  
     epochevents.orderposStim= [epochevents.orderposStim,nan(1,length(events.value))];  
    % getting the data in EEGlab format
    [EEG,winrej] = getDataDeconv(cfg_eeg,epochevents,200); 
    EEG = pop_eegfiltnew(EEG, [],45,300);
    mirindx         = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
       
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

    cfgDesign.eventtypes = {'fixpreT0','fixpreT1onTarg',...
                              'fixpreT1','fixpreT2','fixpreT3',...
                              'fixpre',...
                               'image','stim'};

    cfgDesign.formula   = {'y ~pxini','y ~pxini+spl(dur,10)',...
                              'y ~pxini+spl(dur,10)','y ~pxini+spl(dur,10)','y ~pxini+spl(dur,10)',...
                              'y ~pxini+spl(dur,10)',...
                              'y~1','y~cross*inst*side'};

    model               = 'F_preT_onT_T0_T1_T2_T3_xini_durspl_IM_STsc';
% model               = 'test';
    EEG                 = dc_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-1,1];tic
    EEG                 = dc_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = dc_continuousArtifactExclude(EEG,struct('winrej',winrej,'zerodata',0));
    EEG                 = dc_glmfit(EEG);
 
    unfold              = dc_beta2unfold(EEG);

      mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
      save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
      clear unfold
end 
%%
% %2nd level analysis

 if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname', p.analysisname);
 end
    
stimB = [];

model               = 'F_preT_onT_T0_T1_T2_T3_xini_durspl_IM_STsc';
bslcor    = [];
stimSpl = [];
splmes =   linspace(120,260,8);
for tk = p.subj
     cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
%     if any(strfind(p.analysisname,'mirr'))
%         mirindx         = mirrindex({unfold.chanlocs.labels},[cfg_eeg.expfolder '/channels/mirror_chans']); 
%         stimB = cat(4,stimB,permute(unfold.beta-unfold.beta(mirindx ,:,:),[1,3,2]));
%     else
    if ~isempty(bslcor) 
        stimB = cat(4,stimB,permute(unfold.beta-...
            repmat(mean(unfold.beta(:,find(unfold.times>bslcor(1) & unfold.times<bslcor(2)),:),2),1,size(unfold.beta,2),1),[1,3,2]));
    else
        stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
    end
      for us = 2:5
         unfold.deconv.splines{us}.name = [num2str(us+1) '_dur'];
      end
      
     splnames = {'dur','3_dur','4_dur','5_dur','6_dur'};
     
     auxstimSpl = [];
      for ss = 1:length(splnames)
        cfgu.plotParam  = {splnames{ss}};
        cfgu.pred_value = {{splnames{ss},splmes}};
        unf             = dc_getParam(unfold,cfgu);
        unf             = dc_addmarginal(unf);
         auxstimSpl = cat(4, auxstimSpl,unf.beta(:,:,strmatch(splnames{ss},{unf.param.name})));
      end
     stimSpl = cat(5,stimSpl,auxstimSpl);
end

%%
% SPLINES COMPARISONS

% splavg = mean(stimSpl,5);
% diffspl = mean(squeeze(stimSpl(:,:,:,3,:)-stimSpl(:,:,:,4,:)),4);
% cmap = cbrewer('qual','Paired',length(splmes));
% 
% for ch=[7]
% figure,hold on
%     for ss = 1:length(splmes)
%          h = plot(unf.times,squeeze(splavg(ch,:,ss,3)),'Color',cmap(ss,:));,hold on
%          plot(unf.times,squeeze(splavg(ch,:,ss,4)),'--','Color',cmap(ss,:))
%         plot(unf.times,squeeze(diffspl(ch,:,ss)),':','Color',cmap(ss,:))
%     end
%     title(unfold.chanlocs(ch).labels)
% % set(h, {'color'}, num2cell(cmap, 2));
% vline(splmes/1000,':k')
% end
% 
splavg = mean(stimSpl,5);
% diffspl = mean(squeeze(stimSpl(:,:,:,3,:)-stimSpl(:,:,:,4,:)),4);
cmap = cbrewer('qual','Paired',length(splmes));
cmap = cbrewer('qual','Set1',9);
for ch=[4]
figure,hold on
    for spl =1:5
        for ss = 1:length(splmes)
             h = plot(unf.times,squeeze(splavg(ch,:,ss,spl)),'Color',cmap(spl,:));,hold on
%              plot(unf.times,squeeze(splavg(ch,:,ss,4)),'--','Color',cmap(ss,:))
%             plot(unf.times,squeeze(diffspl(ch,:,ss)),':','Color',cmap(ss,:))
        end
    end
    title(unfold.chanlocs(ch).labels)
% set(h, {'color'}, num2cell(cmap, 2));
vline(splmes/1000,':k')
end
%  ylim([-6 6])
%%
load(cfg_eeg.chanfile)
result      = regmodel2ndstat(stimB,unfold.times,elec,100,'signpermT','cluster');

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
end
%
%%
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr']);
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
end
mkdir(pathfig)
plotinterval = [-.3  .2 .02;.2 .7 .02];
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
setAbsoluteFigureSize
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';
cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
   
load(cfg_eeg.chanlocs)
timeLim         = [-.3 .7];
chnstoPlot      = {{'Cz','CPz'},{'Pz','POz'},{'O1','Oz','O2'},...
    {'P5','P7','PO7'},{'P6','P8','PO8'},{'C3','CP3'}};
lineColors      = cbrewer('qual','Set1',9);%[134 16 9;22 79 134;11 93 24]/255;
  Bnames        = {'sacpreT0onTarg_Intercept','sacpreT0_2_Intercept',...
                              'sacpreT1_3_Intercept','sacpreT2_4_Intercept','sacpreT3_5_Intercept','sacpre_6_Intercept','sacpre_6_pxdiff','image_13_Intercept','stim_14_Intercept'}
%   Bnames        = {'fixpreT0_7_Intercept','fixpreT1onTarg_8_Intercept','fixpreT1_9_Intercept',...
%        'fixpreT2_10_Intercept','fixpreT3_11_Intercept','fixpre_12_Intercept','fixpre_12_pxini','image_13_Intercept'};
Bnames        = {'sacpreT0_2_Intercept','fixpreT0_7_Intercept','sacpreT1_3_Intercept','fixpreT1_9_Intercept'}
xaxis           = result.clusters(1).time;
axLim           = [-.7 .9 -15 15];

for ch = 1:length(chnstoPlot)
    auxChns         = find(ismember({chanlocs.labels},chnstoPlot{ch}));
    datatoPlot  = [];
    for b = 1:length(Bnames)
        ixB = strmatch(Bnames{b},result.coeffs);
        datatoPlot  = cat(3,datatoPlot,permute(squeeze(mean(result.B(auxChns,ixB,:,:))),[2 1]));
    end
    lineNames   = strtok(Bnames,'_');
%         result.clusters(ixB).mask
%         signf        = squeeze(any(any(statUCIci.mask(auxChns,freqs,:),1),2))';
    fh = fillPlot(datatoPlot,[],xaxis,axLim,'mean',lineColors,lineNames);
    fh.Name = strjoin(' ',chnstoPlot{ch});
    xlabel('Time (s)','FontSize',12)
    ylabel('\muV','Interpreter','tex','FontSize',12)
    figsize     = [8 8*fh.Position(4)/fh.Position(3)];
%     doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_Inf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
end
%%
for ch = 1:length(chnstoPlot)
    figure
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
%     doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',['Channs_' strjoin('_',chnstoPlot{ch})],[2 2],1)
end
