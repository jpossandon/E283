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
 % model               = 'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc';
%model = 'S_orderPrecat_xdifspl_IM_STsc_eff'; 
 model = 'F_orderPrecat_xdifspl_IM_STsc_eff'; 
 % model = 'F_nextto_xdifspl_IM_STsc_eff'; 
%%
% subject configuration and data
%  p.subj = p.subj(4:end)
  
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
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'orderPreT','<8'},...
                                [800 100],{-1,2,'origstart','>0';-2,1,'origstart','>0'});
    
    events.prefixdur        = nan(1,length(events.start));  
    events.prefixdur(2:3:end)  = events.dur(1:3:end); 
    events.onTarg(2:3:end)  = events.onTarg(1:3:end);
    events = struct_elim(events,[1:3:length(events.type)],2);
  
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.dur         = events.dur; 
    epochevents.type        = cell(1,length(events.start));
    epochevents.prefixdur   = events.prefixdur; 
    epochevents.orderPreT   =  events.orderPreT;
     epochevents.type(events.type==2 & events.orderPreT==0 & events.onTarg) = repmat({'sacpreT0onTarg'},1,sum(events.type==2 & events.orderPreT==0 & events.onTarg));
     epochevents.type(events.type==2 & events.orderPreT==0 & ~events.onTarg) = repmat({'sacpre'},1,sum(events.type==2 & events.orderPreT==0 & ~events.onTarg));
     epochevents.type(events.type==2 & events.orderPreT>0) = repmat({'sacpre'},1,sum(events.type==2 & events.orderPreT>0));
    
      epochevents.type(events.type==1 & (events.orderPreT==0 | (events.orderPreT==1 & events.onTarg)) ) = repmat({'fixpreT0T1onTarg'},1,sum(events.type==1 & (events.orderPreT==0 | (events.orderPreT==1 & events.onTarg)) ));
   %   epochevents.type(events.type==1 & events.orderPreT==1 & events.onTarg)    = repmat({'fixpreT1onTarg'},1,sum(events.type==1  & events.orderPreT==1 & events.onTarg));
      epochevents.type(events.type==1 & events.orderPreT==1 & ~events.onTarg)    = repmat({'fixpreT1'},1,sum(events.type==1 & events.orderPreT==1 & ~events.onTarg));
      
      epochevents.type(events.type==1 & (events.orderPreT>1))   = repmat({'fixpre'},1,sum(events.type==1 & (events.orderPreT>1)));
    
      

    epochevents.orderposStim= events.orderposStim;
    epochevents.pxini       = (events.posinix-960)/45;            
    epochevents.pyini       = (events.posiniy-540)/45;     
    epochevents.pxend       = (events.posendx-960)/45;            
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;  
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini; 
    epochevents.pxdiff(2:2:end)  = epochevents.pxdiff(1:2:end); 
    epochevents.pydiff(2:2:end)  = epochevents.pydiff(1:2:end); 
    epochevents.dirtotarg(epochevents.pxdiff>0) = 1;
    epochevents.dirtotarg(epochevents.pxdiff<0) = 2;
     epochevents.dirtotarg(2:2:end-2)  = epochevents.dirtotarg(3:2:end); 
     
    events.angle(events.angle<0) = 360+events.angle(events.angle<0);           % transform angles to 0-360
    epochevents.angle       = events.angle;
    epochevents.amp       = events.amp;
    epochevents.side        = nan(1,length(events.start));    
    epochevents.cross       = nan(1,length(events.start));    
    epochevents.inst        = nan(1,length(events.start)); 
 
    if strcmp(model,'S_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff') || strcmp(model,'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc')
         epochevents = struct_elim(epochevents,[2:2:length(epochevents.type)],2);
         if any(epochevents.orderPreT==8)
            epochevents = struct_elim(epochevents,find(epochevents.orderPreT==8),2);
         end
    end
    if strcmp(model,'F_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff')
        epochevents = struct_elim(epochevents,[1:2:length(epochevents.type)],2);
         if any(epochevents.orderPreT==8)
            epochevents = struct_elim(epochevents,find(epochevents.orderPreT==8),2);
         end
    end
    
   
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
    epochevents.prefixdur      = [epochevents.prefixdur ,nan(1,length(events.value))];  
    epochevents.orderposStim= [epochevents.orderposStim,nan(1,length(events.value))];
    epochevents.angle       = [epochevents.angle,nan(1,length(events.value))];
         epochevents.orderPreT        = [epochevents.orderPreT nan(1,length(events.value))];
     epochevents.amp       = [epochevents.amp,nan(1,length(events.value))];
   epochevents.dirtotarg       = [epochevents.dirtotarg,nan(1,length(events.value))];
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

    if strcmp(model,'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc')
        cfgDesign.eventtypes = {'sacpreT0onTarg','sacpreT0',...
                            'sacpreT1','sacpreT2','sacpreT3',...
                            'sacpre','image','stim'};

        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)',...
                            'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)',...
                            'y ~1+spl(pxdiff,10)+spl(pydiff,10)', 'y~1','y~cross*inst*side'};
    end

     if strcmp(model,'F_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0T1onTarg','fixpreT1','fixpre','image','stim'};

        cfgDesign.formula   = {'y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(dirtotarg)+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+cat(dirtotarg)+spl(pxdiff,10)+spl(pydiff,10)','y~1','y~cross*inst*side'};
        if strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff')
            cfgDesign.codingschema = 'effects';
        end
     end
  
     if strcmp(model,'S_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'sacpreT0onTarg','sacpre','image','stim'};

        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)','y~1','y~cross*inst*side'};
        if strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
            cfgDesign.codingschema = 'effects';
        end
     end
     
% model               = 'test';
    EEG                 = uf_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-1,1];tic
    EEG                 = uf_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));
    EEG                 = uf_glmfit(EEG);
 
    unfold              = uf_condense(EEG);

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

%  model               = 'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc';
bslcor      = [];
stimSpl     = [];
stimBspl    = [];
ufpredict   = [];

for tk = p.subj
     cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
%     if any(strfind(p.analysisname,'mirr'))
%         mirindx         = mirrindex({unfold.chanlocs.labels},[cfg_eeg.expfolder '/channels/mirror_chans']); 
%         stimB = cat(4,stimB,permute(unfold.beta-unfold.beta(mirindx ,:,:),[1,3,2]));
%     else
    if  strcmp(model,'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc') 
        pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
        pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         
        ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred},...
            {'2_pxdiff',pxd_pred},{'2_pydiff',pyd_pred},{'3_pxdiff',pxd_pred},{'3_pydiff',pyd_pred},{'4_pxdiff',pxd_pred},{'4_pydiff',pyd_pred},...
            {'5_pxdiff',pxd_pred},{'5_pydiff',pyd_pred},{'6_pxdiff',pxd_pred},{'6_pydiff',pyd_pred}});
     
    end
    if  strcmp(model,'S_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
        pxdM        = 2;   pxd_pred0    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
        pydM        = 2;   pyd_pred0    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
        pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
        pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         
        ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred0},{'pydiff',pyd_pred0},...
            {'2_pxdiff',pxd_pred},{'2_pydiff',pyd_pred}});
     
    end
    if  strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff') 
        pxdM        = 2;   pxd_pred0    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
        pydM        = 2;   pyd_pred0    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
        pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
        pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         
        ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred0},{'pydiff',pyd_pred0},...
            {'2_pxdiff',pxd_pred},{'2_pydiff',pyd_pred},{'3_pxdiff',pxd_pred},{'3_pydiff',pyd_pred}});
     
    end
       % remove splines from coefficients estimates
     if strcmp(model,'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc') ||  strcmp(model,'S_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
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
E283_run2nd_save

%%
% plot betas
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr']);
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
end
plotBetasTopos(cfg_eeg,result,pathfig,[-.1  .6 .025],[])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';
cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop

%chnstoPlot      = {{'Cz','CPz'},{'Pz','POz'},{'O1','Oz','O2'},...
 %   {'P5','P7','PO7'},{'P6','P8','PO8'},{'C3','CP3'}};
chnstoPlot      = {{'FCz','Cz','CPz'},{'O1','Oz','O2'},{'P1','Pz','P2','PO3','POz','PO4'},{'P5','P7','PO7'}};
Bnames        = {'sacpre_2_Intercept','sacpre_orderPreT_0','sacpre_orderPreT_1','sacpre_orderPreT_2','sacpre_orderPreT_3','sacpre_orderPreT_4','sacpre_orderPreT_5','sacpreT0onTarg_Intercept'}
%  Bnames        = {'fixpreT0_Intercept','fixpre_3_Intercept','fixpre_orderPreT_1','fixpre_orderPreT_2','fixpre_orderPreT_3','fixpre_orderPreT_4','fixpre_orderPreT_5','fixpre_orderPreT_6'}
axLim           = [-1 .9 -2 5];
plotBetasChannels(cfg_eeg,result,chnstoPlot,Bnames,pathfig,axLim,'preTarget')


%%
for ch = 1:length(chnstoPlot)
    figure
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
%     doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',['Channs_' strjoin('_',chnstoPlot{ch})],[2 2],1)
end
