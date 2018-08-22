% E283
% Basic TF analysis aligned to saccade pre tartge
% - Simple TF charts
% - Selection fo peak frequencies at theta, alpha and beta bands
% - GlM analysis per frequency

%%
% TFR
% sge               = str2num(getenv('SGE_TASK_ID'));
clear
E283_params
MACpath = '/Users/jossando/trabajo/E283/';
p.analysisname = 'TFR_SacPreT_hann';
%%
s=1
%p.analysisname = 'saclockTFRdps';

for tk = p.subj(s:end) % subject number

    % Analysis parameters
    p.times_tflock              = [700 700];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-.550 -.150]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'yes';
    p.cfgTFR.toi                = (-p.times_tflock(1):10:p.times_tflock(2))/1000;	
    p.cfgTFR.method             = 'mtmconvol';%'wavelet';%%
    
    % wavelet parameters       
    %time of a cycle in ms times number of cicles divide by two for baseline limit :1/fo*width/2 
    %bandwith sf = f0/width   st = 1/(2*pi*f0/width)
%     p.cfgTFR.width              = 4; 
    
    % single tapers
    if any(strfind(p.analysisname,'hann'))
        p.cfgTFR.foi                = 4:1:30;
        p.cfgTFR.taper              = 'hanning';%'dpss';
        p.cfgTFR.pad                = 4;
        p.cfgTFR.t_ftimwin          = .3*ones(1,length(p.cfgTFR.foi));      %    single taper
    end
    if any(strfind(p.analysisname,'dps'))
        p.cfgTFR.foi                = 2.^(3:.125:6);%(2.^([3:.25:5.25]));% %6:1:40	
        p.cfgTFR.taper              = 'dpss';
        p.cfgTFR.t_ftimwin          = 4./p.cfgTFR.foi;
        p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
        %     plottp(p.cfgTFR)
    end


    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk));
    end

    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname',p.analysisname);       % single experiment/session parameters 
  
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])                         
    eyedata.events.diffx = eyedata.events.posendx-eyedata.events.posinix;
    
    if strcmp(p.analysisname,'TF_Sac_dps')
        fieldstoav      = {'Left','Right'};%,'IM'
        [trls.Left,events]           = define_event(cfg_eeg,eyedata,2,{'origstart','>0';...
                                'origstart','<1000';...                            
                                'diffx','<0'},...
                                p.times_tflock); 
        [trls.Right,events]           = define_event(cfg_eeg,eyedata,2,{'origstart','>0';...
                                'origstart','<1000';...                            
                                'diffx','>0'},...
                                p.times_tflock); 
    elseif strcmp(p.analysisname,'TFR_SacPreT_hann') || strcmp(p.analysisname,'TFR_SacPreT_dps')
        % NEXT BY DIFFX SIDE?
        [trl,events]           = define_event(cfg_eeg,eyedata,2,{'&origstart','>0';'orderPreT','<5';'diffx','<0'},...
                                p.times_tflock,{-1,1,'origstart','>0'}); 
        events.onTarg(2:2:end)  = events.onTarg(1:2:end); 
        trls.LsacpreT0onTarg = trl(events.orderPreT(2:2:end)==0 & events.onTarg(2:2:end),:);
        trls.LsacpreT0       = trl(events.orderPreT(2:2:end)==0 & ~events.onTarg(2:2:end),:);
        trls.LsacpreT1       = trl(events.orderPreT(2:2:end)==1,:);
        trls.Lsacpre         = trl(events.orderPreT(2:2:end)>1,:);
        [trl,events]           = define_event(cfg_eeg,eyedata,2,{'&origstart','>0';'orderPreT','<5';'diffx','>0'},...
                                p.times_tflock,{-1,1,'origstart','>0'}); 
        events.onTarg(2:2:end)  = events.onTarg(1:2:end); 
        trls.RsacpreT0onTarg = trl(events.orderPreT(2:2:end)==0 & events.onTarg(2:2:end),:);
        trls.RsacpreT0       = trl(events.orderPreT(2:2:end)==0 & ~events.onTarg(2:2:end),:);
        trls.RsacpreT1       = trl(events.orderPreT(2:2:end)==1,:);
        trls.Rsacpre         = trl(events.orderPreT(2:2:end)>1,:);
    end
    
    fieldstoav = fields(trls);
    at  = 1;
    for f = 1:length(fieldstoav)
         [auxTFR toelim]    = getTFRsfromtrl({cfg_eeg},{trls.(fieldstoav{f})},...
             p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
         auxTFR.N_tot       = [size(trls.(fieldstoav{f}),1)-length(toelim{1}), size(trls.(fieldstoav{f}),1)];
         N_allSubj(s,f,:)   = [size(trls.(fieldstoav{f}),1)-length(toelim{1}), size(trls.(fieldstoav{f}),1)];
         TFR.(fieldstoav{f}) = auxTFR;
         clear auxTFR
     end

%      mirroring left stimulation data to make contra/ipsi plots
    mirindx         = mirrindex(TFR.(fieldstoav{1}).(p.analysis_type{1}).label,[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
%     Left_mirr   = TFR.Left;
%     Left_mirr.(p.analysis_type{at}).powspctrm = TFR.Left.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
%     
%     cfgs            = [];
%     cfgs.parameter  = 'powspctrm';
%     TFR.CI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.Right.(p.analysis_type{1}),Left_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
%     clear Left_mirr
%     fieldstoav = fields(TFR);
%     
     for f = 1:length(fieldstoav)
         TFRav.(fieldstoav{f}).(p.analysis_type{1}) = ft_freqdescriptives([], TFR.(fieldstoav{f}).(p.analysis_type{1}));
     end
    clear TFR
    mkdir([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/'])
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_' p.analysis_type{at}],'TFRav','cfg_eeg','p')
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/Nall'],'N_allSubj','fieldstoav')
    %     save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/' cfg_eeg.sujid '_imcoh_stim_' p.analysis_type{at}],'imcoh','cfg_eeg','p')
%     save(['/Users/jossando/trabajo/E283/analysis/' cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')
s=s+1;
end

%%
% grand averages
at                  = 1;
for b = 1%:2
    p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    cfgr                = [];
    p.bsl               = [-.550 -.300];  
    if b==1
        cfgr.baseline       = p.bsl;
        cfgr.baselinetype   = 'db';
    else
        cfgr.baseline       = 'no';
    end
    cfgr.keepindividual = 'yes';
    s=1;
    MACpath = '/Users/jossando/trabajo/E283/';
    % MACpath = '/Volumes/nibaldo/trabajo/E283/';
    for tk = p.subj; % subject number
        if ismac    
            cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname',p.analysisname,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
        else
            cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname',p.analysisname);
        end
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_' p.analysis_type{at}],'TFRav')
        fTFR    = fields(TFRav);
        for ff=1:length(fTFR)
            faux(s,ff) = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
        end
        s=s+1;
    end
    cfgga.keepindividual = 'yes';
    for ff=1:length(fTFR)
        str_GA = 'GA.(fTFR{ff}) = ft_freqgrandaverage(cfgga';
        for ss = 1:size(faux(:,ff),1)
            str_GA = [str_GA, ',faux(', num2str(ss), ',ff)'];
        end
        str_GA = [str_GA,');'];
        eval(str_GA)
    end

%     mirindx         = mirrindex(GA.U_I.label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
%   mirindx         = mirrindex(GA.Left.label,[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
  
%     GA.Leftci             = GA.Left;
%     GA.Leftci.powspctrm   = GA.Left.powspctrm-GA.Left.powspctrm(:,mirindx,:,:);
%     GA.Rightci             = GA.Right;
%     GA.Rightci.powspctrm   = GA.Right.powspctrm-GA.Right.powspctrm(:,mirindx,:,:);
    
%     GA.CIci             = GA.CI;
%     GA.CIci.powspctrm   = GA.CI.powspctrm-GA.CI.powspctrm(:,mirindx,:,:);
    
    if b==1
        GAbsl = GA;
    else
        GAnobsl = GA;
    end
    clear faux
end
%%
cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';

%  GAbsl.LeftvsRight       = ft_math(cfgs,GAbsl.Left,GAbsl.Right);
%  GAbsl.LeftvsRightci       = ft_math(cfgs,GAbsl.Leftci,GAbsl.Rightci);
GAbsl.LvsRpreT1      = ft_math(cfgs,GAbsl.LsacpreT1,GAbsl.RsacpreT1); 
GAbsl.LvsRpreT0      = ft_math(cfgs,GAbsl.LsacpreT0,GAbsl.RsacpreT0);
GAbsl.LvsRpreT0onTarg      = ft_math(cfgs,GAbsl.LsacpreT0onTarg,GAbsl.RsacpreT0onTarg);
  GAbsl.LvsRpre      = ft_math(cfgs,GAbsl.Lsacpre,GAbsl.Rsacpre);
GAbsl.LpreT0vsLpre      = ft_math(cfgs,GAbsl.LsacpreT0,GAbsl.Lsacpre);
GAbsl.RpreT0vsRpre      = ft_math(cfgs,GAbsl.RsacpreT0,GAbsl.Rsacpre);

%%
load(cfg_eeg.chanfile)
cfgp                = [];
cfgp.showlabels     = 'no'; 
cfgp.fontsize       = 12; 
cfgp.elec           = elec;
cfgp.interactive    = 'yes';
% cfgp.channel        = mirindx(1:38);
% cfgp.trials         = 4
% cfgp.baseline       = p.bsl;
% cfgp.baselinetype   = 'db';
% cfgp.ylim           = [0 40];
%   cfgp.xlim           = [-.75 .1];
%  cfgp.zlim           = [-.5 .5];
%   cfgp.maskparameter  = 'mask';
%   cfgp.maskalpha      = 1;
% cfgp.parameter      = 'stat';

   data = GAbsl.LsacpreT0
% data = TFRav.LsacpreT1.ICAem
%      data =GAbsl.C_Ici;
% %  data.mask = statUCIci.mask;
figure,ft_multiplotTFR(cfgp,data)
% 2.^([3:.5:5]+-.25)

%%
% load('cmapjp','cmap') 
%   cmap = cmocean('thermal');
 cmap = flipud(cbrewer('div','RdBu',128)); 
cfgp                = [];
cfgp.xlim           = [-.3 .9];
cfgp.zlim           = [-1.5 1.5];
% cfgp.ylim           = [8 30];
cfgp.colorbar       = 'no';
datatoplot      = {'U_Ici','C_Ici','diffUCIci','U_unIci','C_unIci','diffUCunIci','diffUC_IvsUnici'};
chnstoPlot      = {{'C3','CP3'},{'P5','P7','PO7'}};

for ch = 1:length(chnstoPlot)
    for ff = 1:length(datatoplot)
        cfgp.channel        = chnstoPlot{ch};
        figure
        ft_singleplotTFR(cfgp,GAbsl.(datatoplot{ff}))
        colormap(cmap),
        vline(0,'k')
        figsize     = [8 4];
        doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',...
            [datestr(now,'ddmmyy') '_TF_' datatoplot{ff} '_' strjoin('_',cfgp.channel)],figsize,1)
    end
end
cfgp.channel        = {'C3','CP3'};
figure,ft_singleplotTFR(cfgp,GAbsl.U_Ici),colormap(cmap),vline(0,'k')
figsize     = [8 4];
doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_TF_U_Ici_' strjoin('_',cfgp.channel)],figsize,1)
figure,ft_singleplotTFR(cfgp,GAbsl.C_Ici),colormap(cmap),vline(0,'k')
doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_TF_C_Ici_' strjoin('_',cfgp.channel)],figsize,1)
figure,ft_singleplotTFR(cfgp,GAbsl.diffUCIci),colormap(cmap),vline(0,'k')
doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_TF_diffUC_Ici_' strjoin('_',cfgp.channel)],figsize,1)

cfgp.channel        = {'P5','P7','PO7'};
figure,ft_singleplotTFR(cfgp,GAbsl.U_Ici),colormap(cmap),vline(0,'k')
figure,ft_singleplotTFR(cfgp,GAbsl.C_Ici),colormap(cmap),vline(0,'k')
figure,ft_singleplotTFR(cfgp,GAbsl.diffUCIci),colormap(cmap),vline(0,'k')

cfgp.channel        = {'C3','CP3'};
figure,ft_singleplotTFR(cfgp,GAbsl.U_unIci),colormap(cmap),vline(0,'k')
figure,ft_singleplotTFR(cfgp,GAbsl.C_unIci),colormap(cmap),vline(0,'k')
figure,ft_singleplotTFR(cfgp,GAbsl.diffUCunIci),colormap(cmap),vline(0,'k')
cfgp.channel        = {'P5','P7','PO7'};
figure,ft_singleplotTFR(cfgp,GAbsl.U_unIci),colormap(cmap),vline(0,'k')
figure,ft_singleplotTFR(cfgp,GAbsl.C_unIci),colormap(cmap),vline(0,'k')
figure,ft_singleplotTFR(cfgp,GAbsl.diffUCunIci),colormap(cmap),vline(0,'k')
%%
% Comparisons stat
% chtoplot = find(ismember(difffreq2.label,{'P5', 'P7', 'PO7'}));
load(cfg_eeg.chanfile)
statUCIci       = freqpermWS(GAbsl.U_Ici,GAbsl.C_Ici,elec,1000);
statUCunIci     = freqpermWS(GAbsl.U_unIci,GAbsl.C_unIci,elec,1000);
statUCIvsUnici  = freqpermWS(GAbsl.diffUCIci,GAbsl.diffUCunIci,elec,1000);
save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/stat_hann_' datestr(now,'ddmmyy')],'statUCIci','statUCunIci','statUCIvsUnici')
% statUC = freqpermWS(GAbsl.U_I,GAbsl.C_I,elec,1000)
% statBsl = freqpermBSL(GAnobsl.U_I,p.bsl,elec,500)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOPOPLOT SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('stat_hann_220218.mat')
clim            = [-2.5 2.5];
plotinterval    = [-.1 .8 .05];
bands           = [8 16;15 27;20 27;36 40];
bandnames       = {'alpha','beta','hbeta','lgamma'};
datatoplot      = {'U_Ici','C_Ici','diffUCIci','U_unIci','C_unIci','diffUCunIci','diffUC_IvsUnici'};
stattoplot      = {[],[],'statUCIci',[],[],'statUCunIci','statUCIvsUnici'};
setAbsoluteFigureSize
for b = 1:length(bandnames)
    for ff = 1:length(datatoplot)
        if isempty(stattoplot{ff})
            stat = [];
            stat.time   = GAbsl.(datatoplot{ff}).time;
        else
            eval(['stat=' stattoplot{ff} ';'])
        end
       fh          = topomitlinesTFR(cfg_eeg,stat,GAbsl.(datatoplot{ff}),plotinterval,clim,bands(b,:));
         figsize     = [17.6 17.6*fh.Position(4)/fh.Position(3)];
        doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' datatoplot{ff} '_' bandnames{b}],figsize,1)
    end
end

%%
setAbsoluteFigureSize
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('stat_hann_160118.mat')
load(cfg_eeg.chanlocs)
timeLim         = [-.1 .8];
bands           = [8 16;15 27];
chnstoPlot      = {{'C3','CP3'},{'P5','P7','PO7'}};
chnsLabel       = {'central','parietal'};
bandnames       = {'alpha','beta'};
% datatoplot      = {'U_Ici','C_Ici','diffUCIci','U_unIci','C_unIci','diffUCunIci','diffUC_IvsUnici'};
% stattoplot      = {[],[],'statUCIci',[],[],'statUCunIci','statUCIvsUnici'};
lineColors      = [134 16 9;22 79 134;11 93 24]/255;

xaxis           = GAbsl.U_Ici.time;
axLim           = [-.2 1 -2 2];


for b = 1:length(bandnames)
    for ch = 1:length(chnstoPlot)
        freqs           = find(GAbsl.U_Ici.freq>bands(b,1) & GAbsl.U_Ici.freq<bands(b,2));
        auxChns         = find(ismember({chanlocs.labels},chnstoPlot{ch}));
        
        datatoPlot  = cat(3,squeeze(nanmean(nanmean(GAbsl.U_Ici.powspctrm(:,auxChns,freqs,:),2),3)),...
            squeeze(nanmean(nanmean(GAbsl.C_Ici.powspctrm(:,auxChns,freqs,:),2),3)),...
            squeeze(nanmean(nanmean(GAbsl.diffUCIci.powspctrm(:,auxChns,freqs,:),2),3)));
        lineNames   = {'|| Inf','X Inf','|| - X Inf'};
        signf        = squeeze(any(any(statUCIci.mask(auxChns,freqs,:),1),2))';
        fh = fillPlot(datatoPlot,signf,xaxis,axLim,'mean',lineColors,lineNames);
        xlabel('Time (s)','FontSize',12)
        ylabel('Chan-mirrChan Pow','Interpreter','tex','FontSize',12)
        figsize     = [8 8*fh.Position(4)/fh.Position(3)];
        doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_Inf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
    
         datatoPlot  = cat(3,squeeze(nanmean(nanmean(GAbsl.U_unIci.powspctrm(:,auxChns,freqs,:),2),3)),...
            squeeze(nanmean(nanmean(GAbsl.C_unIci.powspctrm(:,auxChns,freqs,:),2),3)),...
            squeeze(nanmean(nanmean(GAbsl.diffUCunIci.powspctrm(:,auxChns,freqs,:),2),3)));
        lineNames   = {'|| unInf','X unInf','|| - X nuInf'};
        signf        = squeeze(any(any(statUCunIci.mask(auxChns,freqs,:),1),2))';
        fh = fillPlot(datatoPlot,signf,xaxis,axLim,'mean',lineColors,lineNames);
        xlabel('Time (s)','FontSize',12)
        ylabel('Chan-mirrChan Pow','Interpreter','tex','FontSize',12)
        figsize     = [8 8*fh.Position(4)/fh.Position(3)];
        doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_unInf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
    
	end
end
for ch = 1:length(chnstoPlot)
    figure
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
    doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',['Channs_' strjoin('_',chnstoPlot{ch})],[2 2],1)
end
    
  