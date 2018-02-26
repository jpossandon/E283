% E283
% Basic TF analysis aligned to stimulus
% - Simple TF charts
% - Selection fo peak frequencies at theta, alpha and beta bands
% - GlM analysis per frequency

%%
% TFR
% sge               = str2num(getenv('SGE_TASK_ID'));
clear
E283_params
MACpath = '/Users/jossando/trabajo/E283/';

%%
s=1
for tk = p.subj % subject number

    % Analysis parameters
    p.times_tflock              = [2000 1500];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-1.5 -1]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'no';
    p.cfgTFR.toi                = (-p.times_tflock(1):10:p.times_tflock(2))/1000;	
    p.cfgTFR.foi                = 6:1:40;%6:1:80;%2.^(3:.125:5);%(2.^([3:.25:5.25]));% %6:1:40	
    p.cfgTFR.method             = 'wavelet';%%'mtmconvol';%
    p.cfgTFR.output             = 'pow';
    % wavelet parameters       
    %time of a cycle in ms times number of cicles divide by two for baseline limit :1/fo*width/2 
    %bandwith sf = f0/width   st = 1/(2*pi*f0/width)
%     p.cfgTFR.width              = 4; 
    
    % single tapers
%     p.cfgTFR.taper              = 'hanning';%'dpss';
%     p.cfgTFR.pad                = 4;
%     p.cfgTFR.t_ftimwin          = .250*ones(1,length(p.cfgTFR.foi));      %    single taper

%    single tapers

    p.cfgTFR.taper              = 'dpss';
    p.cfgTFR.t_ftimwin          = 4./p.cfgTFR.foi;
    p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
%    plottp(p.cfgTFR)

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
                                            'analysisname','preT_TFRwavlong');       % single experiment/session parameters 
  
%     mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])                         
    eyedata.events.nextToTarg =any(eyedata.events.nextToTarg);
    fieldstoav = {'fixpreT0',...
                'fixpreT1OnTarg','fixpreT1nextT','fixpreT1',...
                'fixpreT2nextT','fixpreT2',...
                'fixpreT3nextT','fixpreT3',...
                'fixpreT4nextT','fixpreT4',...
                'fixpreT5nextT','fixpreT5'};
    oPT      = [0 1 1 1 2 2 3 3 4 4 5 5];
    onT      = [1 1 0 0 0 0 0 0 0 0 0 0];
    ntT      = [0 0 1 0 1 0 1 0 1 0 1 0];
    
    at=1;
    for tt = [1,8]%:length(fieldstoav)
        [trls.(fieldstoav{tt}),events]           = define_event(cfg_eeg,eyedata,1,{'origstart','>0';...
                                'orderPreT',['==' num2str(oPT(tt))];...
                                'onTarg',['==' num2str(onT(tt))];...
                                'nextToTarg',['==' num2str(ntT(tt))]},...
                                p.times_tflock); 
    end
    
   for tt = [1,8]%:length(fieldstoav)
        [auxTFR toelim]    = getTFRsfromtrl({cfg_eeg},{trls.(fieldstoav{tt})},...
                            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
         auxTFR.N_tot       = [size(trls.(fieldstoav{tt}),1)-length(toelim{1}), size(trls.(fieldstoav{tt}),1)];
        N_allSubj(s,tt,:)   = [size(trls.(fieldstoav{tt}),1)-length(toelim{1}), size(trls.(fieldstoav{tt}),1)];
        TFR.(fieldstoav{tt}) = auxTFR;
   end

%      mirroring left stimulation data to make contra/ipsi plots
%     mirindx         = mirrindex(TFR.LU_I.(p.analysis_type{1}).label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
%    
%     fieldstoav = fields(TFR);
%     for f = 1:length(fieldstoav)
%         TFRav.(fieldstoav{f}).(p.analysis_type{1}) = ft_freqdescriptives([], TFR.(fieldstoav{f}).(p.analysis_type{1}));
%     end
TFRav = TFR;
    if ~isdir([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/']);mkdir([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/']);end
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/Nall'],'N_allSubj','fieldstoav')
   s=s+1;
end

%%
% grand averages
clear
E283_params
at                  = 1;

for b = 1%:2
    p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    cfgr                = [];
    p.bsl               = [-2 -1.5];  
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
            cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','preT_TFRwavlong','expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
        else
            cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','preT_TFRwavlong');
        end
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav')
        fTFR    = fields(TFRav);
        for ff=1:length(fTFR)
            auxf = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
            if ~isfield(auxf,'elec')
                auxf.elec = faux(s-1,ff).elec;
            end
            faux(s,ff) =orderfields(auxf);
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
% 
%     GA.U_Ici             = GA.U_I;
%     GA.U_Ici.powspctrm   = GA.U_I.powspctrm-GA.U_I.powspctrm(:,mirindx,:,:);
%     GA.U_unIci           = GA.U_unI;
%     GA.U_unIci.powspctrm = GA.U_unI.powspctrm-GA.U_unI.powspctrm(:,mirindx,:,:);
%     GA.C_Ici             = GA.C_I;
%     GA.C_Ici.powspctrm   = GA.C_I.powspctrm-GA.C_I.powspctrm(:,mirindx,:,:);
%     GA.C_unIci           = GA.C_unI;
%     GA.C_unIci.powspctrm = GA.C_unI.powspctrm-GA.C_unI.powspctrm(:,mirindx,:,:);
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

 GAbsl.diff       = ft_math(cfgs,GAbsl.fixpreT0,GAbsl.fixpreT3);
%   difffreq2.mask   = statUC.mask;
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
%  cfgp.xlim           = [-.2 .4];
%  cfgp.zlim           = [-.5 .5];
%   cfgp.maskparameter  = 'mask';
%   cfgp.maskalpha      = 1;
% cfgp.parameter      = 'stat';

  data = GAbsl.diff;
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
        doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',...
            [datestr(now,'ddmmyy') '_TF_' datatoplot{ff} '_' strjoin('_',cfgp.channel)],figsize,1)
    end
end
cfgp.channel        = {'C3','CP3'};
figure,ft_singleplotTFR(cfgp,GAbsl.U_Ici),colormap(cmap),vline(0,'k')
figsize     = [8 4];
doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_TF_U_Ici_' strjoin('_',cfgp.channel)],figsize,1)
figure,ft_singleplotTFR(cfgp,GAbsl.C_Ici),colormap(cmap),vline(0,'k')
doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_TF_C_Ici_' strjoin('_',cfgp.channel)],figsize,1)
figure,ft_singleplotTFR(cfgp,GAbsl.diffUCIci),colormap(cmap),vline(0,'k')
doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_TF_diffUC_Ici_' strjoin('_',cfgp.channel)],figsize,1)

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
chtoplot = find(ismember(difffreq2.label,{'P5', 'P7', 'PO7'}));
load(cfg_eeg.chanfile)
statUCIci       = freqpermWS(GAbsl.U_Ici,GAbsl.C_Ici,elec,1000);
statUCunIci     = freqpermWS(GAbsl.U_unIci,GAbsl.C_unIci,elec,1000);
statUCIvsUnici  = freqpermWS(difffreq1,difffreq2,elec,1000);
save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/stat_hann_' datestr(now,'ddmmyy')],'statUCIci','statUCunIci','statUCIvsUnici')
% statUC = freqpermWS(GAbsl.U_I,GAbsl.C_I,elec,1000)
% statBsl = freqpermBSL(GAnobsl.U_I,p.bsl,elec,500)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOPOPLOT SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('stat_hann_160118.mat')
clim            = [-2.5 2.5];
plotinterval    = [-.1 .8 .05];
bands           = [8 16;15 27];
bandnames       = {'alpha','beta'};
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
        doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' datatoplot{ff} '_' bandnames{b}],figsize,1)
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
        doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_Inf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
    
         datatoPlot  = cat(3,squeeze(nanmean(nanmean(GAbsl.U_unIci.powspctrm(:,auxChns,freqs,:),2),3)),...
            squeeze(nanmean(nanmean(GAbsl.C_unIci.powspctrm(:,auxChns,freqs,:),2),3)),...
            squeeze(nanmean(nanmean(GAbsl.diffUCunIci.powspctrm(:,auxChns,freqs,:),2),3)));
        lineNames   = {'|| unInf','X unInf','|| - X nuInf'};
        signf        = squeeze(any(any(statUCunIci.mask(auxChns,freqs,:),1),2))';
        fh = fillPlot(datatoPlot,signf,xaxis,axLim,'mean',lineColors,lineNames);
        xlabel('Time (s)','FontSize',12)
        ylabel('Chan-mirrChan Pow','Interpreter','tex','FontSize',12)
        figsize     = [8 8*fh.Position(4)/fh.Position(3)];
        doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_unInf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
    
	end
end
for ch = 1:length(chnstoPlot)
    figure
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
    doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',['Channs_' strjoin('_',chnstoPlot{ch})],[2 2],1)
end
    
  