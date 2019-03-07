% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                  % basic experimental parameters               %
p.analysisname  = 'TL_dc_preT';
if ismac
    run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')
else
    run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')
end


 model = 'F_2dsplvec_IM_STsc_eff';
% model = 'F_2dsplvec_2dsplend_orderPrecat_IM_STsc_eff';
% model = 'F_2Dsplvec_orderPrecat_IM_STsc_eff';
% model = 'F_orderPrecat_spldiffxy_IM_STsc_eff';
% model = 'F_orderPrecat_DToTargspl_DToTargspl2_spldiffxy_IM_STsc_eff';
% model = 'F_T0_Tx_DToTargspl_spldiffxy_IM_STsc_eff';
% model = 'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff';
% model = 'F_T0_T1_T2_Tx_DToTargspl_spldiffxy_IM_STsc_eff';
%   model = 'F_T0_T1_Tx_DToTargspl_2Dsplvec_IM_STsc_eff';
%%
% subject configuration and data
%%{
if ismac
    load('/Users/jossando/trabajo/E283/07_Analysis/03_Eye/eyedata/alleyedataFULLevidence')
else
    load('/Users/jpo/trabajo/E283/07_Analysis/03_Eye/eyedata/alleyedataFULLevidence')
end
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
    % check eyedad is the same as full data
    normev = std(data.evidence(1,data.subject==tk & data.type==1));
    if length(eyedata.events.start)~=sum(data.subject==tk)
        error('data individual subject file not equal to the one from all')
    else
        eyedata.events.evidence = data.evidence(1,data.subject==tk);%./normev;
        eyedata.events.evidence2 = data.evidence(2,data.subject==tk);%./normev;
        eyedata.events.evidence3 = data.evidence(3,data.subject==tk);%./normev;
        eyedata.events.evidence4 = data.evidence(4,data.subject==tk);%./normev;
    end
    nmovs = 40;
    eyedata.events.nextToTarg = any([eyedata.events.nextToTargH;eyedata.events.nextToTargV;eyedata.events.nextToTargD]);
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'orderPreT',['<' num2str(nmovs)]},...
        [800 100],{-1,2,'origstart','>0';-2,1,'origstart','>0'});
%     events.DToTarg = events.DToTarg./45; 
%     events.DToTarg2 = zeros(1,length(events.DToTarg)); 
%     events.DToTarg2(3:3:end) = events.DToTarg(1:3:end);
    events.onTarg(2:3:end)  = events.onTarg(1:3:end);
%     events.DToTarg(2:3:end) = events.DToTarg(1:3:end);  
%     events.DToTarg2(2:3:end) = events.DToTarg(1:3:end); this is not ready for saccades 
     

    
    events = struct_elim(events,[1:3:length(events.type)],2);
%     events.amp(2:2:end)     = events.amp(1:2:end);
    
    epochevents             = [];
     epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
      epochevents.dur         = events.dur;
    epochevents.type        = cell(1,length(events.start));
     epochevents.trial      = events.trial;
    epochevents.orderPreT   =  events.orderPreT;
%     epochevents.DToTarg     = events.DToTarg;
%     epochevents.DToTarg2     = events.DToTarg2;
%    
   if strcmp(model,'F_2dsplvec_2dsplend_orderPrecat_IM_STsc_eff') || strcmp(model,'F_orderPrecat_spldiffxy_IM_STsc_eff') ||...
           strcmp(model,'F_2dsplvec_orderPrecat_IM_STsc_eff') || strcmp(model,'F_orderPrecat_DToTargspl_DToTargspl2_spldiffxy_IM_STsc_eff') || ...
           strcmp(model,'F_2dsplvec_IM_STsc_eff')
         epochevents.type(events.type==1 & events.orderPreT>=0)   = repmat({'fixpre'},1,sum(events.type==1 & events.orderPreT>=0));
   elseif strcmp(model,'F_T0_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
         epochevents.type(events.type==1 & events.orderPreT==0) = repmat({'fixpreT0'},1,sum(events.type==1 & events.orderPreT==0));
%         epochevents.type(events.type==1 & events.orderPreT==1 & events.onTarg)    = repmat({'fixpreT1onTarg'},1,sum(events.type==1  & events.orderPreT==1 & events.onTarg));
        epochevents.type(events.type==1 & events.orderPreT>0)   = repmat({'fixpre'},1,sum(events.type==1 & events.orderPreT>0));
    elseif strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff') || strcmp(model,'F_T0_T1_Tx_DToTargspl_2Dsplvec_IM_STsc_eff')
         epochevents.type(events.type==1 & events.orderPreT==0) = repmat({'fixpreT0'},1,sum(events.type==1 & events.orderPreT==0));
        epochevents.type(events.type==1 & events.orderPreT==1)   = repmat({'fixpreT1'},1,sum(events.type==1 & events.orderPreT==1));
        epochevents.type(events.type==1 & events.orderPreT>1)   = repmat({'fixpre'},1,sum(events.type==1 & events.orderPreT>1));
   elseif strcmp(model,'F_T0_T1_T2_Tx_DToTargspl_spldiffxy_IM_STsc_eff') 
         epochevents.type(events.type==1 & events.orderPreT==0) = repmat({'fixpreT0'},1,sum(events.type==1 & events.orderPreT==0));
        epochevents.type(events.type==1 & events.orderPreT==1)   = repmat({'fixpreT1'},1,sum(events.type==1 & events.orderPreT==1));
        epochevents.type(events.type==1 & events.orderPreT==2)   = repmat({'fixpreT2'},1,sum(events.type==1 & events.orderPreT==2));
        epochevents.type(events.type==1 & events.orderPreT>2)   = repmat({'fixpre'},1,sum(events.type==1 & events.orderPreT>2));
     
   end
    
    
%     epochevents.orderposStim= events.orderposStim;
    epochevents.pxini       = (events.posinix-960)/45;
    epochevents.pyini       = (events.posiniy-540)/45;
    epochevents.pxend       = (events.posendx-960)/45;
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini;
    epochevents.pxdiff(2:2:end)  = epochevents.pxdiff(1:2:end);
    epochevents.pydiff(2:2:end)  = epochevents.pydiff(1:2:end);
%     epochevents.amp       = events.amp;
    epochevents.side        = nan(1,length(events.start));
    epochevents.cross       = nan(1,length(events.start));
    epochevents.inst        = nan(1,length(events.start));
    
    
    if strncmp(model,'F_',2)
        epochevents = struct_elim(epochevents,[1:2:length(epochevents.type)],2);
    end
    
    if any(epochevents.orderPreT==nmovs)
        epochevents = struct_elim(epochevents,find(epochevents.orderPreT==nmovs),2);
    end
        
    
    [trl,events]            = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','>0'},[1500 900]);
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
    epochevents.trial      = [epochevents.trial,events.trial];
        
     epochevents.dur       = [epochevents.dur,nan(1,length(events.value))];
    epochevents.pxini       = [epochevents.pxini,nan(1,length(events.value))];
    epochevents.pyini       = [epochevents.pyini,nan(1,length(events.value))];
    epochevents.pxend       = [epochevents.pxend,nan(1,length(events.value))];
    epochevents.pyend       = [epochevents.pyend,nan(1,length(events.value))];
    epochevents.pxdiff      = [epochevents.pxdiff,nan(1,length(events.value))];
    epochevents.pydiff      = [epochevents.pydiff,nan(1,length(events.value))];
    epochevents.orderPreT       = [epochevents.orderPreT nan(1,length(events.value))];
%     epochevents.amp             = [epochevents.amp,nan(1,length(events.value))];
%     epochevents.DToTarg        = [epochevents.DToTarg,nan(1,length(events.value))];
%     epochevents.DToTarg2        = [epochevents.DToTarg2,nan(1,length(events.value))];
    % getting the data in EEGlab format
    % without trials that have a prefixation on target
%     if any(strfind(model,'S_'))
%         epochevents = struct_elim(epochevents,find(ismember(epochevents.trial,epochevents.trial(strcmp(epochevents.type,'sacpreT0onTarg')))),2);
%     elseif any(strfind(model,'F_'))
%         epochevents = struct_elim(epochevents,find(ismember(epochevents.trial,epochevents.trial(strcmp(epochevents.type,'fixpreT1onTarg')))),2);
%     end
    [EEG,winrej]    = getDataDeconv(cfg_eeg,epochevents,500,1);
    mirindx         = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']);
    
    
    % deconvolution design
    cfgDesign           = [];
    if strcmp(model,'F_2dsplvec_2dsplend_orderPrecat_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+2dspl(pxdiff,pydiff,8)+2dspl(pxend,pyend,8)+cat(orderPreT)','y~1','y~cross*inst*side'};
    elseif strcmp(model,'F_2dsplvec_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+2dspl(pxdiff,pydiff,8)','y~1','y~cross*inst*side'};
   
    elseif  strcmp(model,'F_2dsplvec_orderPrecat_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+2dspl(pxdiff,pydiff,10)+cat(orderPreT)','y~1','y~cross*inst*side'}; 
    elseif strcmp(model,'F_orderPrecat_spldiffxy_IM_STsc_eff') 
         cfgDesign.eventtypes = {'fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)+cat(orderPreT)','y~1','y~cross*inst*side'};
    elseif strcmp(model,'F_orderPrecat_DToTargspl_DToTargspl2_spldiffxy_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)+spl(DToTarg2,10)+cat(orderPreT)','y~1','y~cross*inst*side'};
    elseif strcmp(model,'F_T0_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y~1','y~cross*inst*side'};
    elseif strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpreT1','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y~1','y~cross*inst*side'};
    elseif  strcmp(model,'F_T0_T1_Tx_DToTargspl_2Dsplvec_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpreT1','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+2dspl(pxdiff,pydiff,10)','y ~1+2dspl(pxdiff,pydiff,10)+spl(DToTarg,10)','y ~1+2dspl(pxdiff,pydiff,10)+spl(DToTarg,10)','y~1','y~cross*inst*side'};
    elseif  strcmp(model,'F_T0_T1_T2_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpreT1','fixpreT2','fixpre','image','stim'};
         cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y~1','y~cross*inst*side'};
    end
    
   
    cfgDesign.codingschema = 'effects';
    
    EEG                 = uf_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-.8,.8];tic
    EEG                 = uf_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));
    EEG                 = uf_glmfit(EEG);
    
    unfold              = uf_condense(EEG);
    
    mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
    save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    clear unfold
end
%}
%%
% %2nd level analysis

if ismac
    cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
else
    cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname', p.analysisname);
end

stimB = [];

%  model               = 'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc';
bslcor      = [-.8 -.4];
%  bslcor      = [];
stimSpl     = [];
stimBspl    = [];
ufpredict   = [];

for tk = p.subj
    cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
   
    
    if  strcmp(model,'F_2dsplvec_2dsplend_orderPrecat_IM_STsc_eff') || strcmp(model,'F_orderPrecat_spldiffxy_IM_STsc_eff') || ...
            strcmp(model,'F_orderPrecat_DToTargspl_DToTargspl2_spldiffxy_IM_STsc_eff') || strcmp(model,'F_T0_Tx_DToTargspl_spldiffxy_IM_STsc_eff') || ...
            strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff') || strcmp(model,'F_T0_T1_T2_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
          pxd_pred    = -6:1.5:6;   pyd_pred    = -6:1.5:6;
          pxe_pred    = -12:3:12;   pye_pred    = -12:3:12;
          pDToTatg = [.5:.5:8];
        if  strcmp(model,'F_2dsplvec_2dsplend_orderPrecat_IM_STsc_eff')  
            ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{'pxdiffpydiff',pxd_pred,pyd_pred},{'pxendpyend',pxe_pred,pye_pred}});
        end
        if  strcmp(model,'F_orderPrecat_spldiffxy_IM_STsc_eff') || strcmp(model,'F_orderPrecat_DToTargspl_DToTargspl2_spldiffxy_IM_STsc_eff') 
            ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{'pxdiff',pxd_pred,},{'pydiff',pyd_pred},{'DToTarg',pDToTatg},{'DToTarg2',pDToTatg}});
        end
       if strcmp(model,'F_T0_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{'pxdiff',pxd_pred,},{'pydiff',pyd_pred},{'2_pxdiff',pxd_pred,},{'2_pydiff',pyd_pred},{'DToTarg',pDToTatg}});
        
       end
       
       if strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{'pxdiff',pxd_pred,},{'pydiff',pyd_pred},{'2_pxdiff',pxd_pred,},{'2_pydiff',pyd_pred},...
                {'3_pxdiff',pxd_pred,},{'3_pydiff',pyd_pred},{'DToTarg',pDToTatg},{'3_DToTarg',pDToTatg}});
        
       end
       if strcmp(model,'F_T0_T1_T2_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{'pxdiff',pxd_pred,},{'pydiff',pyd_pred},{'2_pxdiff',pxd_pred,},{'2_pydiff',pyd_pred},...
                {'3_pxdiff',pxd_pred,},{'3_pydiff',pyd_pred},{'4_pxdiff',pxd_pred,},{'4_pydiff',pyd_pred},{'DToTarg',pDToTatg},{'3_DToTarg',pDToTatg},{'4_DToTarg',pDToTatg}});
        
       end
    end
    
  
    % remove splines from coefficients estimates
    if any(strfind(model,'spl'))
        rmspl                   = strmatch('spline',{unfold.param.type});
        unfold.beta(:,:,rmspl)  = [];
        unfold.param(rmspl)     = [];
        
        rmnospl                   = setdiff(1:length(ufpredict.param),strmatch('spline',{ufpredict.param.type}));
        ufpredict.beta(:,:,rmnospl) = [];
        ufpredict.param(rmnospl)    = [];
    end
    if ~isempty(bslcor)
        stimB       = cat(4,stimB,permute(unfold.beta-...
            repmat(mean(unfold.beta(:,find(unfold.times>bslcor(1) & unfold.times<bslcor(2)),:),2),1,size(unfold.beta,2),1),[1,3,2]));
        if ~isempty(ufpredict)
            stimBspl    = cat(4,stimBspl,permute(ufpredict.beta-...
                repmat(mean(ufpredict.beta(:,find(ufpredict.times>bslcor(1) & ufpredict.times<bslcor(2)),:),2),1,size(ufpredict.beta,2),1),[1,3,2]));
        end
    else
        stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
        if ~isempty(ufpredict)
            stimBspl    = cat(4,stimBspl,permute(ufpredict.beta(:,:,:),[1,3,2]));
        end
    end
    
end

%
% run 2nd level stat and save
stattype = 'signpermT';
mc = 'cluster';
np= 1000
 E283_run2nd_save

%%
% plot betas
% bslcor = [];
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr']);
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
end
plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.3  .6 .05],[],0)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';
cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop

chnstoPlot      = {{'O1','Oz','O2'},{'P1','Pz','P2','PO3','POz','PO4'},{'P5','P7','PO7'},{'P6','P8','PO8'}};
 Bnames          = {'fixpreT0_Intercept','fixpreT1_2_Intercept','fixpreT2_3_Intercept','fixpre_4_Intercept'}
 Blabels         = {'order  0','order 1','order 2','order >1'}
% Bnames          = {'fixpre_orderPreT_0','fixpre_orderPreT_1','fixpre_orderPreT_2','fixpre_orderPreT_3','fixpre_orderPreT_4','fixpre_orderPreT_5','fixpre_orderPreT_6'}
%  Blabels         = {'order  0','order -1','order -2','order -3','order -4','order -5','order -6'}

%  Blabels         = {[0:1:6]};
axLim           = [-.8 .8 -8 8];
filled = ones(1,7);
% lineColors      = cbrewer('qual','Set1',length(Bnames)+1);
lineColors      = cbrewer('div','RdYlGn',length(Bnames));
lineColors = [1 0 0;0 0 1;0 1 0;0 1 1]
plotBetasChannels(cfg_eeg,'mean',result,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'fixpreTarget')

axLim           = [-.8 .8 -6 8];
 Bnames = {result.coeffs(1:8),result.coeffs(1:8),result.coeffs(1:8),result.coeffs(1:8),result.coeffs(1:8),result.coeffs(1:8),result.coeffs(1:8)}
 signs1 = [1  1 0 0 0 0 0 0;1 0 1 0 0 0 0 0;1 0 0 1 0 0 0 0;1 0 0 0 1 0 0 0;1 0 0 0 0 1 0 0;1 0 0 0 0 0 1 0;1 0 0 0 0 0 0 1]; 
 lineColors      = cbrewer('div','RdYlGn',length(Bnames));
Blabels         = {[0:1:6]};
 whstat          = [0 2 3 4 5 6 7];
 ploteffectChannels(cfg_eeg,result,whstat,chnstoPlot,signs1,Bnames,Blabels,ones(1,7),lineColors,pathfig,axLim,'fixpreTargetInt')
 
 %%
 axLim           = [-.6 .9 -4 8];
 Bnames = {result.coeffs(1:7),result.coeffs(1:7),result.coeffs(1:7),result.coeffs(1:7),result.coeffs(1:7),result.coeffs(1:7)}
 signs1 = [1 0 0 0 0 0 0;0 1 1 0 0 0 0;0 1 0 1 0 0 0;0 1 0 0 1 0 0;0 1 0 0 0 1 0;0 1 0 0 0 0 1]; 
 lineColors      = cbrewer('div','RdYlGn',length(Bnames));
Blabels         = {[0:1:5]};
 ploteffectChannels(cfg_eeg,result,chnstoPlot,signs1,Bnames,Blabels,ones(1,6),lineColors,pathfig,axLim,'sacpreTargetInt')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],'splines');
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],'splines');
end

if any(strfind(p.analysisname,'mirr'))
    plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3  .6 .05],[],1,0)
else
    plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3  .6 .05],[],0,0)

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],'splines');
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],'splines');
end

if any(strfind(p.analysisname,'mirr'))
    plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3  .6 .05],[],1,0)
else
    plotSplinesTopos2D(cfg_eeg,resultSplines,'mean',pathfig,[-.2  .4 .005],[-7 7],0,0)

end

% combine 2 splines
% comb2spl
% plotSplinesTopos2D(cfg_eeg,resultSplinescomb,'mean',pathfig,[.06  .16 .02],[-7 7],0,0)

%%


% %resultSplines.B = resultSplines.B; 
% chnstoPlot      = {{'P1','Pz','P2','PO3','POz','PO4'},{'O1','Oz','O2'}};
% axLim           = [-.3 .6 -2.5 2.5];
% Bnames = {resultSplines.coeffs{strmatch('sacpre_DToTarg',resultSplines.coeffs)}};
% Blabels = {[1:10]};
% filled = zeros(1,10);
% lineColors      = cbrewer('div','RdYlGn',length(Bnames));
% plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'splineDToTarg')
% 
% Bnames = {resultSplines.coeffs{strmatch('sacpre_amp',resultSplines.coeffs)}};
% Blabels = {[1:10]};
% filled = zeros(1,10);
% lineColors      = cbrewer('div','RdYlGn',length(Bnames));
% plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'sacamp')
% 
% axLim           = [-.3 .6 -5 5];
%  Bnames = {resultSplines.coeffs{strmatch('sacpre_pydiff',resultSplines.coeffs)}};
%  Blabels = {[-10:1:10]};
%  filled = zeros(1,21);
%  lineColors      = cbrewer('div','RdYlGn',length(Bnames));
%  plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'pydif')
% 
% axLim           = [-.3 .6 -2 4];
% chnstoPlot      = {{'P1','Pz','P2','PO3','POz','PO4'},{'O1','PO7','PO3'},{'O2','PO8','PO4'}};
%  Bnames = {resultSplines.coeffs{strmatch('sacpre_pxdiff',resultSplines.coeffs)}};
%  Blabels = {[-10:1:10]};
%  filled = zeros(1,21);
%  lineColors      = cbrewer('div','RdYlGn',length(Bnames));
%  plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'pxdif')
% 
%  axLim           = [-.3 .6 -2 4];
% chnstoPlot      = {{'P1','Pz','P2','PO3','POz','PO4'},{'O1','PO7','PO3'},{'O2','PO8','PO4'}};
%  Bnames = {resultSplines.coeffs{strmatch('fixpre_DToTarg',resultSplines.coeffs)}};
%  Blabels = {[1:1:20]};
%  filled = zeros(1,20);
%  lineColors      = cbrewer('div','RdYlGn',length(Bnames));
%  plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'fixpre_DToTarg')

 %%
 % FIXATION
 spltoplot = {'fixpreT1_DToTarg_','fixpreT2_3_DToTarg_','fixpre_4_DToTarg_'};
 for sp = 1:length(spltoplot)
    chnstoPlot      = {{'P1','Pz','P2','PO3','POz','PO4'},{'O1','Oz','O2'}};
    axLim           = [-.8 .8 -5 5];
    Bnames = {resultSplines.coeffs{strmatch(spltoplot{sp},resultSplines.coeffs)}};
    Blabels = {[.5:.5:8]};
    filled = zeros(1,16);
    % lineColors      = cbrewer('div','RdYlGn',length(Bnames));
    lineColors     = cmocean('thermal',length(Bnames)); 
    plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,spltoplot{sp})
 end



 %%
for ch = 1:length(chnstoPlot)
    fh=figure;
    fh.Position = fh.Position/2
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
    doimage(gcf,['/Users/jossando/trabajo/E283/02_Planning/report'],'tiff',['Channs_' strjoin('_',chnstoPlot{ch})],'opengl',[2 2],1)
end
    