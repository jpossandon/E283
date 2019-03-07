
% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'deconvTFCI';
%model           = 'IM_STcsi';
model           = 'IM_STcsi_SdxyCI';
%
% subject configuration and data
 
if ismac 
run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')        
else
run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')   
end    
%%
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
    [trl,events]           = define_event(cfg_eeg,eyedata,2,{'&origstart','>0';'&origstart','<800'},...
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
    fprintf('\nMinimum intertrial latency subject %d:%d',tk,min(diff(events.time(events.value==96))))

    epochevents.pxini       = [epochevents.pxini,nan(1,length(events.value))];    
    epochevents.pyini       = [epochevents.pyini,nan(1,length(events.value))];       
    epochevents.pxend       = [epochevents.pxend,nan(1,length(events.value))];          
    epochevents.pyend       = [epochevents.pyend,nan(1,length(events.value))];  
    epochevents.pxdiff      = [epochevents.pxdiff,nan(1,length(events.value))];  
    epochevents.pydiff      = [epochevents.pydiff,nan(1,length(events.value))];  
    % getting the data in EEGlab format
    [EEG,winrej] = getDataDeconv(cfg_eeg,epochevents,250,0);  
     mirindx         = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
       
    if any(strfind(p.analysisname,'CI'))
        LstimTimes = epochevents.latency(find(strcmp(epochevents.type,'stim') & epochevents.side == -1)-1);
        for tt  = 1:length(LstimTimes)
            EEGst = find(EEG.times>LstimTimes(tt),1);
            mirsamples = EEGst-floor(EEG.srate*.45):EEGst+ceil(EEG.srate*1.6);
            EEG.data(:,mirsamples) = EEG.data(mirindx,mirsamples);
            ETst  = find(epochevents.latency>LstimTimes(tt) & epochevents.latency<LstimTimes(tt)+850 & ...
                (strcmp(epochevents.type,'sac') | strcmp(epochevents.type,'fix')));
            if ~isempty(ETst)
                [epochevents.pxini(ETst) epochevents.pxend(ETst) epochevents.pxdiff(ETst)] = deal(-epochevents.pxini(ETst),-epochevents.pxend(ETst),-epochevents.pxdiff(ETst));
            end
        end
    end
   
    cfgDesign           = [];
    if strcmp(model,'IM_STcsi')
        cfgDesign.eventtypes= {'image','stim'};
        cfgDesign.formula   = {'y~1','y~cat(cross)*cat(inst)*cat(side)'};
    elseif strcmp(model,'IM_STcsi_Sdxy')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)',...
                                'y~1','y~cat(cross)*cat(inst)*cat(side)'};
    elseif strcmp(model,'IM_STcsi_SdxyCI')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)',...
                                'y~1','y~cat(cross)*cat(inst)'};
  
    end
    cfgDesign.codingschema = 'effects';
    %freqbands       = {'alfa','beta'};
    %bplim           = [9 15;16 26];
    
    freqbands       = {'alfa','beta'};
    bplim           = [9 15;16 25];
    filtPnts        = [368,254]; % check this
    trSt            = strmatch('S 96',{EEG.event.type});
    evLat           = [EEG.event.latency];
    lat_trSt        = round(evLat(trSt));
    bsl             = [-.450 0];
    winsize         = .3/(1/EEG.srate);
%     for tt = winsize/2+1:1:EEG.pnts-winsize/2-1
%         auxd = EEG.data(:,tt-winsize/2:tt+winsize/2-1);

    for fb = 1:length(freqbands)
      
      
        % this two are quite equivalent but holber path is faster
        EEGaux              = pop_eegfiltnew(EEG, bplim(fb,1), bplim(fb,2), [], 0, [], 0);
        EEGaux.data         = abs(hilbert(EEGaux.data')').^2; %'power'
          %foi = bplim(fb,1):bplim(fb,2)
         % EEGaux= EEG;
%         for ch = 1:EEG.nbchan
%             [spectrum_mtmconvol,ntaper,foi,toi] = ft_specest_mtmconvol(EEG.data(ch,:),EEG.times/1000,...
%                 'timeoi',EEG.times/1000,'timwin',3./foi ,'taper','hanning',...
%                 'pad',[],'freqoi',foi,'polyorder',0,'dimord','chan_time_freqtap');
%             EEGaux.data(ch,:) =  sum(abs(spectrum_mtmconvol).^2,3);
%         end
        bslsample           = floor(bsl.*EEGaux.srate);
%         sampltobsl            = [-.4 1.6].*EEGaux.srate;
        for ll = 1:length(lat_trSt)
            trlbslpow(:,ll)       = mean(EEGaux.data(:,bslsample(1)+lat_trSt(ll):bslsample(2)+lat_trSt(ll)),2);
        end
        EEGaux.data = 10*log10(EEGaux.data./repmat(mean(trlbslpow,2),1,EEGaux.pnts));
    %        EEGaux.data(:,sampltobsl(1)+lat_trSt(ll):sampltobsl(2)+lat_trSt(ll)) = 10*log10(EEGaux.data(:,sampltobsl(1)+lat_trSt(ll):sampltobsl(2)+lat_trSt(ll))./...
     %           repmat(trlbslpow,1,diff(sampltobsl)+1)); 
    
        %EEGaux.data         = 10*log10(EEGaux.data./repmat(nanmean(EEGaux.data,2),1,size(EEGaux.data,2))); 
        if any(strfind(p.analysisname,'mirr'))
            EEGaux.data = EEGaux.data-EEGaux.data(mirindx,:); 
        end
        EEGaux              = uf_designmat(EEGaux,cfgDesign);
        cfgTexp             = [];
        cfgTexp.timelimits  = [-.5,.8];tic
        EEGaux              = uf_timeexpandDesignmat(EEGaux,cfgTexp);toc
        EEGaux              = uf_continuousArtifactExclude(EEGaux,struct('winrej',winrej));
        EEGaux              = uf_glmfit(EEGaux);
 
        unfold.(freqbands{fb})          = uf_condense(EEGaux);
        
    end
  for fb = 1:length(freqbands)
      unfold.(freqbands{fb}).bplim = bplim(fb,:);
  end
      mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
      save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
      clear unfold
end 
%%
% %2nd level analysis
 if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname',p.analysisname);
 end
    
freqbands       = {'alfa','beta'};
%freqbands       = {'alfa'};
for fb = 1%:2%:length(freqbands)
    stimB       = [];
    bslcor      = [-.4 -.1];
    bslcorspl      = [-.5 -.3];
    stimBspl    = [];
    ufpredict   = [];
    
    for tk = p.subj
         cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/' model '/glm/' cfg_eeg.sujid '_' model],'unfold')
        unfold = unfold.(freqbands{fb});
         if  strcmp(model,'IM_STcsi_Sdxy') | strcmp(model,'IM_STcsi_SdxyCI')
            pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
            pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
            pxd_pred    = [[-8.7:2:.7],0,[.7:2:8.7]];
            pyd_pred    = [[-4.38:1:.38],0,[.38:1:4.38]];
            ufpredict       = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred}});
         end
    
        if strcmp(model,'IM_STcsi_Sdxy')| strcmp(model,'IM_STcsi_SdxyCI')
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
        %  stimB       = cat(4,stimB,permute(unfold.(freqbands{fb}).beta./...
         %   repmat(mean(unfold.(freqbands{fb}).beta(:,find(unfold.(freqbands{fb}).times>bslcor(1) & ...
          %      unfold.(freqbands{fb}).times<bslcor(2)),:),2),1,size(unfold.(freqbands{fb}).beta,2),1),[1,3,2]));
        % FALTA aQUI      
            if ~isempty(ufpredict)
                stimBspl    = cat(4,stimBspl,permute(ufpredict.beta-...
                    repmat(mean(ufpredict.beta(:,find(ufpredict.times>bslcorspl(1) & ufpredict.times<bslcorspl(2)),:),2),1,size(ufpredict.beta,2),1),[1,3,2]));
            end
        else
            stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
            if ~isempty(ufpredict)
                stimBspl    = cat(4,stimBspl,permute(ufpredict.beta(:,:,:),[1,3,2]));
            end
        end
    end
%   stattype = 'boottrimet';
     stattype = 'signpermT';
    mc       = 'cluster';
    E283_run2nd_save
end


%%
% resultALL = result;
for fb = 1
    
    if isempty(bslcor)
        load(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',['glmALL_' freqbands{fb} stattype '_' mc]),'result')
    else
        load(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',['glmALLbslcorr_' freqbands{fb} stattype '_' mc]),'result')
    end
    if ~isempty(bslcor)
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],(freqbands{fb}));
    else
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],(freqbands{fb}));
    end

%     auxresult = resultALL.(freqbands{fb});
result.B(:,4:end,:,:) = -1*result.B(:,4:end,:,:);
    if any(strfind(p.analysisname,'mirr'))
        plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.3  .6 .0375],[],1);
    else
        plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.3  .6 .075],[],0);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';

lineColors      = cbrewer('qual','Set1',9);
cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
if strcmp(p.analysisname,'deconvTFCI')
   chnstoPlot      = {{'P6','P8','PO8'},{'C4','C6'},{'P5','P7','PO7'},{'C3','C5'}};
%    axLim           = [-.2 .8 -8 1];
%    Bnames          = {result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end)}
%    signs           = [1 1 1;1 -1 1;1 1 -1;1 -1 -1]; signs1 = [signs,signs(:,2).*signs(:,3)];
%    ploteffectChannels(cfg_eeg,result,chnstoPlot,signs1,Bnames,{'UunI','CunI','UI','CI'},lineColors,pathfig,axLim,'condition')
%    
   cb1      = cbrewer('qual','Set1',5)
   cb2      = cbrewer('qual','Pastel1',3)
   lineColors = [cb1(1,:);cb2(1,:);cb1(2,:);cb2(2,:);cb1(1,:);cb1(2,:)]; 
   Bnames          = {result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end)}
   signs           = [1 1 1;1 -1 1;1 1 -1;1 -1 -1;1 1 1;1 1 -1]; signs1 = [signs,signs(:,2).*signs(:,3)];
   signs           = [0 0 0;0 0 0;0 0 0;0 0 0;1 -1 1;1 -1 -1]; signs2 = [signs,signs(:,2).*signs(:,3)];
   filled          = [0 0 0 0 1 1];
   whstat          = [0 0 0 0 6 0];
   axLim           = [-.15 .75 -7 2];
   
%    ploteffectChannels(cfg_eeg,result,whstat,chnstoPlot,cat(3,signs1,signs2),Bnames,{'UunI','CunI','UI','CI','U-CunI','U-CI'},filled,lineColors,pathfig,axLim,'conditionandUCcontrast')
interval = [.3 .75 .45];
   Bnames          = {result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end)};
signs           = [1 1 1;1 -1 1;1 1 -1;1 -1 -1]; signs1 = [signs,signs(:,2).*signs(:,3)];
plotsingleTopo(cfg_eeg,result,signs1,interval,Bnames,{'UunI','CunI','UI','CI'},[-5 1],pathfig,'topo')
     Bnames          = {result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end)}
   signs           = [1 1 1;1 1 -1;1 1 -1;1 -1 -1]; signs1 = [signs,signs(:,2).*signs(:,3)];
   signs           = [1 -1 1;1 -1 -1;1 1 1;1 -1 1]; signs2 = [signs,signs(:,2).*signs(:,3)];
 plotsingleTopo(cfg_eeg,result,cat(3,signs1,signs2),interval,Bnames,{'U-CunI','U-CI','Uinf-uni','Cinf-uni'},[-1.2 1.2],pathfig,'topodiff')
   
else 
    chnstoPlot      = {{'O1','Oz','O2'},{'P1','Pz','P2','PO3','POz','PO4'},{'P5','P7','PO7'},{'P6','P8','PO8'}};
    chnstoPlot      = {{'P5','P7','PO7'},{'C3','C5'}};
    axLim           = [-.4 .8 -2 2];
%     Bnames        = {'image_2_Intercept','stim_3_Intercept','stim_side_-1'};
%     Blabels = Bnames
     cb1      = cbrewer('qual','Set1',5)
     cb2      = cbrewer('qual','Pastel1',3)
     lineColors = [cb1(1,:);cb2(1,:);cb1(2,:);cb2(2,:)]; 
%     plotBetasChannels(cfg_eeg,'mean',result,chnstoPlot,Bnames,Blabels,lineColors,pathfig,axLim,'stimLockTF_intercepts')
   
    Bnames = {result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end),result.coeffs(3:end)}
    signs = [1  1 1  1;1  1 1 -1;1 -1 1  1;1 -1 1 -1]; signs1 = [signs,signs(:,2).*signs(:,3),signs(:,2).*signs(:,4),signs(:,3).*signs(:,4),signs(:,2).*signs(:,3).*signs(:,4)]
    ploteffectChannels(cfg_eeg,result,chnstoPlot,signs1,Bnames,{'ULuni','URuni','CLuni','CRuni'},lineColors,pathfig,axLim,'mirr_uni')
    signs = [1  1 -1  1;1  1 -1 -1;1 -1 -1  1;1 -1 -1 -1]; signs2 = [signs,signs(:,2).*signs(:,3),signs(:,2).*signs(:,4),signs(:,3).*signs(:,4),signs(:,2).*signs(:,3).*signs(:,4)]
    ploteffectChannels(cfg_eeg,result,chnstoPlot,signs2,Bnames,{'ULi','URi','CLi','CRi'},lineColors,pathfig,axLim,'mirr_inf')
    signs = [1  1 1  1;1  -1 1 1;1  1 -1  1;1  -1 -1 1]; signs1 = [signs,signs(:,2).*signs(:,3),signs(:,2).*signs(:,4),signs(:,3).*signs(:,4),signs(:,2).*signs(:,3).*signs(:,4)]
    signs = [1  1 1  -1;1  -1 1 -1;1  1 -1  -1;1  -1 -1 -1]; signs2 = [signs,signs(:,2).*signs(:,3),signs(:,2).*signs(:,4),signs(:,3).*signs(:,4),signs(:,2).*signs(:,3).*signs(:,4)]

    ploteffectChannels(cfg_eeg,result,chnstoPlot,cat(3,signs1,signs2),Bnames,{'UL-URuni','CL-CRuni','UL-URi','CL-CRi'},lineColors,pathfig,axLim,'mirr_LR')

    signs = [1  1 1  1;1  1 1 -1;1  1 -1  1;1  1 -1 -1]; signs1 = [signs,signs(:,2).*signs(:,3),signs(:,2).*signs(:,4),signs(:,3).*signs(:,4),signs(:,2).*signs(:,3).*signs(:,4)]
    signs = [1  -1 1  1;1  -1 1 -1;1  -1 -1  1;1  -1 -1 -1]; signs2 = [signs,signs(:,2).*signs(:,3),signs(:,2).*signs(:,4),signs(:,3).*signs(:,4),signs(:,2).*signs(:,3).*signs(:,4)]

    ploteffectChannels(cfg_eeg,result,chnstoPlot,cat(3,signs1,signs2),Bnames,{'UL-CLuni','UR-CRuni','UL-CLi','UR-CRi'},lineColors,pathfig,axLim,'mirr_UC')
end

%%
 chnstoPlot      = {{'P6','P8','PO8','C4','C6','P5','P7','PO7','C3','C5'}};
for ch = 1:length(chnstoPlot)
    fh=figure;
    fh.Position = fh.Position/2;
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
    doimage(gcf,['/Users/jossando/trabajo/E283/08_Publications/manuscript1/figures/'],'pdf',['Channs_' strjoin('_',chnstoPlot{ch})],'1200','opengl',[1 1],1)
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for fb = 1

    if isempty(bslcor)
        load(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',['glmALL_' freqbands{fb} stattype '_' mc]),(freqbands{fb}),'resultSplines')
    else
%         load(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',['glmALLbslcorr_' freqbands{fb} stattype '_' mc]),(freqbands{fb}),'resultSplines')
         load(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',['glmALLbslcorr_' freqbands{fb} stattype '_' mc '_spldiff']),(freqbands{fb}),'resultSplines','resultSpldiff')
    end
    if ~isempty(bslcor)
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],(freqbands{fb}),'splines');
    else
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],(freqbands{fb}),'splines');
    end
    if any(strfind(p.analysisname,'mirr'))
        plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3  .6 .0375],[],1,1)
   resultSpldiff.times = resultSpldiff.clusters(1).time
        plotBetasTopos(cfg_eeg,resultSpldiff,'mean',pathfig,[-.3  .6 .0375],[],1);

    else
        plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3  .6 .075],[],0,1)
        resultSpldiff.times = resultSpldiff.clusters(1).time
        plotBetasTopos(cfg_eeg,resultSpldiff,'mean',pathfig,[-.3  .6 .075],[],0);
%         plotSplinesTopos(cfg_eeg,resultSpldiff,'mean',pathfig,[-.3  .6 .075],[],0,1)

    end
end
%%
chnstoPlot      = {{'P6','P8','PO8'},{'C3','C5'}};
axLim           = [-.3 .6 -1 1];
Bnames = {resultSplines.coeffs{strmatch('sac_pxdiff',resultSplines.coeffs)}};
Blabels = {[-8.7:2:-.7,0,.7:2:8.7]};
filled = zeros(1,11);
lineColors      = cbrewer('div','RdYlGn',length(Bnames));
plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'splinepxdiff')

chnstoPlot      = {{'P6','P8','PO8'},{'C3','C5'}};
axLim           = [-.3 .6 -1 1];
Bnames = {resultSplines.coeffs{strmatch('sac_pydiff',resultSplines.coeffs)}};
Blabels = {[-4.38:1:-.38,0,.38:2:4.38]};
filled = zeros(1,11);
lineColors      = cbrewer('div','RdYlGn',length(Bnames));
plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,filled,lineColors,pathfig,axLim,'splineydiff')

%%

%chnstoPlot      = {{'Cz','CPz'},{'Pz','POz'},{'O1','Oz','O2'},...
 %   {'P5','P7','PO7'},{'P6','P8','PO8'},{'C3','CP3'}};
chnstoPlot      = {{'FCz','Cz','CPz'},{'O1','Oz','O2'}};
Bnames        = {'image_2_Intercept','stim_3_Intercept','stim_side_-1'};
axLim           = [-1 .9 -10 10];
for fb = 1:2
    result = resultALL.(freqbands{fb});
    plotBetasChannels(cfg_eeg,result,chnstoPlot,Bnames,axLim)
end