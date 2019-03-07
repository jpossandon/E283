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
p.analysisname = 'TF_Stim_hann';
load(fullfile(MACpath,'07_Analysis','03_Eye','eyedata','tgtPos.mat')) 
%%
s=1
for tk = p.subj % subject number

    % Analysis parameters
    p.times_tflock              = [1100 1100];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-.650 -.350]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'yes';
    p.cfgTFR.toi                = (-p.times_tflock(1):15:p.times_tflock(2))/1000;	
    %p.cfgTFR.foi                = 8:1:40;%2.^(3:.125:7);%(2.^([3:.25:5.25]));% %6:1:40	
    p.cfgTFR.method             = 'mtmconvol';%'wavelet';%%
    
    % wavelet parameters       
    %time of a cycle in ms times number of cicles divide by two for baseline limit :1/fo*width/2 
    %bandwith sf = f0/width   st = 1/(2*pi*f0/width)
%     p.cfgTFR.width              = 4; 
    
    % single tapers
     if any(strfind(p.analysisname,'hann'))
          p.cfgTFR.foi                = 2:1:40;
        p.cfgTFR.taper              = 'hanning';%'dpss';
        p.cfgTFR.pad                = 4;
        p.cfgTFR.t_ftimwin          = .5.*ones(1,length(p.cfgTFR.foi));      %    single taper
     elseif any(strfind(p.analysisname,'dpss'))
        %p.cfgTFR.foi                = 2.^(5:.125:7);
          p.cfgTFR.foi                = 25:5:150;
        p.cfgTFR.taper              = 'dpss';
        p.cfgTFR.t_ftimwin          = 5./p.cfgTFR.foi;
        p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
        %plottp(p.cfgTFR)
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
  
%     mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])                         
    
    fieldstoav      = {'LU_I','RU_I','LC_I','RC_I','LU_unI','RU_unI','LC_unI','RC_unI'};%,'IM'
    trigs           = [1,2,5,6,9,10,13,14,96];
    for f = 1:8
        [trls.(fieldstoav{f}),events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value',['==' num2str(trigs(f))]},p.times_tflock);            
        if strfind(fieldstoav{f},'unI')
            trls.([fieldstoav{f} '_L']) = trls.(fieldstoav{f});
            trls.([fieldstoav{f} '_R']) = trls.(fieldstoav{f});
            if strfind(fieldstoav{f},'LU') | strfind(fieldstoav{f},'RU')
                addlat = 0;
            else
                addlat =200;
            end
            RR = []; LL = [];
            for ev = 1:length(events.time)
            auxsac = find(eyedata.events.type==2 & eyedata.events.start>events.time(ev)+addlat & eyedata.events.trial==events.trial(ev),1,'first');
                if eyedata.events.posendx(auxsac)>posVec.stimRes(1)/2
                    RR = [RR,ev];
                else
                    LL = [LL,ev];
                end
            end
            trls.([fieldstoav{f} '_L']) = trls.([fieldstoav{f} '_L'])(LL,:);
            trls.([fieldstoav{f} '_R']) = trls.([fieldstoav{f} '_R'])(RR,:);
        end
        
    end
    
    % Locked to stimulus
    at  = 1;
    fieldstoav = fields(trls);
%     mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/TFR_' p.analysis_type{at} '/'])
    for f = 1:length(fieldstoav)
        p.cfgTFR.output        = 'fourier';
%          p.cfgTFR.foi           = 8:1:40;
        [auxTFR toelim]         = getTFRsfromtrl({cfg_eeg},{trls.(fieldstoav{f})},...
             p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
        auxTFR.ICAem.powspctrm  = abs(auxTFR.ICAem.fourierspctrm).^2;
        auxTFR.ICAem.itpc       = squeeze(abs(sum(auxTFR.ICAem.fourierspctrm./abs(auxTFR.ICAem.fourierspctrm),1))/size(auxTFR.ICAem.fourierspctrm,1)); % remove the first singleton dimension
        auxTFR.ICAem            = rmfield(auxTFR.ICAem,'fourierspctrm');
        auxTFR.ICAem.dimord     = 'rpt_chan_freq_time';
        auxTFR.N_tot       = [size(trls.(fieldstoav{f}),1)-length(toelim{1}), size(trls.(fieldstoav{f}),1)];
        N_allSubj(s,f,:)   = [size(trls.(fieldstoav{f}),1)-length(toelim{1}), size(trls.(fieldstoav{f}),1)];
        TFR.(fieldstoav{f}) = auxTFR;
    end

%      mirroring left stimulation data to make contra/ipsi plots
    mirindx         = mirrindex(TFR.LU_I.(p.analysis_type{1}).label,[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
    for f = 1:length(fieldstoav)
        TFR.([fieldstoav{f} '_mirr']) = TFR.([fieldstoav{f}]);
        TFR.([fieldstoav{f} '_mirr']).(p.analysis_type{at}).powspctrm = TFR.([fieldstoav{f}]).(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
        TFR.([fieldstoav{f} 'ci'])   = TFR.([fieldstoav{f}]);
        TFR.([fieldstoav{f} 'ci']).(p.analysis_type{at}).powspctrm  = TFR.([fieldstoav{f}]).(p.analysis_type{at}).powspctrm-TFR.([fieldstoav{f} '_mirr']).(p.analysis_type{at}).powspctrm; 
    
    end
    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR.U_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RU_I.(p.analysis_type{1}),TFR.LU_I_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR.C_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RC_I.(p.analysis_type{1}),TFR.LC_I_mirr.(p.analysis_type{1})); %SAME
    TFR.U_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RU_unI.(p.analysis_type{1}),TFR.LU_unI_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR.U_unI_ipsi.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RU_unI_R.(p.analysis_type{1}),TFR.LU_unI_L_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR.U_unI_contra.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RU_unI_L.(p.analysis_type{1}),TFR.LU_unI_R_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    
    
    TFR.C_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RC_unI.(p.analysis_type{1}),TFR.LC_unI_mirr.(p.analysis_type{1})); %SAME
    TFR.C_unI_contra.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RC_unI_L.(p.analysis_type{1}),TFR.LC_unI_R_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR.C_unI_ipsi.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RC_unI_R.(p.analysis_type{1}),TFR.LC_unI_L_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
   
    TFR.ALLUC.(p.analysis_type{1})   = ft_appendfreq(cfgs,TFR.U_I.(p.analysis_type{1}),TFR.C_I.(p.analysis_type{1}),TFR.U_unI.(p.analysis_type{1}),TFR.C_unI.(p.analysis_type{1}));
    bands               = {'alpha','beta'};
    for bb = 1:2
        if strcmp(bands{bb},'alpha')
            findx = find(TFR.U_I.ICAem.freq>8 &TFR.U_I.ICAem.freq<16);
        elseif strcmp(bands{bb},'beta')
            findx = find(TFR.U_I.ICAem.freq>15 &TFR.U_I.ICAem.freq<27);
        end
        cfgr                = [];
        cfgr.baseline       =  [-.650 -.150];
        cfgr.baselinetype   = 'db';
        aux1 = ft_freqbaseline(cfgr,TFR.U_I.ICAem);
        aux2 = ft_freqbaseline(cfgr,TFR.C_I.ICAem);
        aux3 = ft_freqbaseline(cfgr,TFR.U_unI.ICAem);
        aux4 = ft_freqbaseline(cfgr,TFR.C_unI.ICAem);
        Y = [squeeze(mean(aux1.powspctrm(:,:,findx,:),3));
        squeeze(mean(aux2.powspctrm(:,:,findx,:),3));
        squeeze(mean(aux3.powspctrm(:,:,findx,:),3));
        squeeze(mean(aux4.powspctrm(:,:,findx,:),3))];
        
        XY = [[ones(size(TFR.U_I.ICAem.powspctrm,1),1);-1*ones(size(TFR.C_I.ICAem.powspctrm,1),1);
            ones(size(TFR.U_unI.ICAem.powspctrm,1),1);-1*ones(size(TFR.C_unI.ICAem.powspctrm,1),1)],...
            [ones(size(TFR.U_I.ICAem.powspctrm,1)+size(TFR.C_I.ICAem.powspctrm,1),1);
            -1.*ones(size(TFR.U_unI.ICAem.powspctrm,1)+size(TFR.C_unI.ICAem.powspctrm,1),1)]];
        XY(:,3) = XY(:,1).*XY(:,2); 
        [sujreg.(bands{bb}).B,sujreg.(bands{bb}).Bt,sujreg.(bands{bb}).STATS] = regntcfe(Y,XY,1,'effect',[],0);    
    end   
    sujreg.times = TFR.U_I.ICAem.time;
    fieldstoav = fields(TFR);
    for f = 1:length(fieldstoav)
        TFRav.(fieldstoav{f}).(p.analysis_type{1}) = ft_freqdescriptives([], TFR.(fieldstoav{f}).(p.analysis_type{1}));
        try
            TFRav.([fieldstoav{f}]).(p.analysis_type{1}).erppow = squeeze(TFR.(fieldstoav{f}).(p.analysis_type{1}).erppow); 
            TFRav.([fieldstoav{f}]).(p.analysis_type{1}).itpc   = squeeze(TFR.(fieldstoav{f}).(p.analysis_type{1}).itpc); 
        catch
         f
        end
    end
    mkdir([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/'])
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','sujreg','cfg_eeg','p')
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/Nall'],'N_allSubj','fieldstoav')
    %     save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/' cfg_eeg.sujid '_imcoh_stim_' p.analysis_type{at}],'imcoh','cfg_eeg','p')
%     save(['/Users/jossando/trabajo/E283/analysis/' cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')
s=s+1;
end

%%
% grand averages
clear
E283_params
at                  = 1;
p.analysisname = 'TF_Stim_hann';
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
for b = 1:2
    p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    cfgr                = [];
    p.bsl               = [-.600 -.400];  
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
    clear faux
    for tk = p.subj; % subject number
%          try
            if ismac    
                cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname',p.analysisname,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
            else
                cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname',p.analysisname);
            end
            load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','sujreg')
            fTFR    = fields(TFRav);
            mirrones = find(~cellfun(@isempty,strfind(fTFR,'_mirr')));
            fTFR(mirrones)= []
            for ff=1:length(fTFR)
                if any(strfind(fTFR{ff},'mirr'))
                    continue
                elseif any(strfind(fTFR{ff},'ci'))
                   bla = TFRav.(fTFR{ff}).(p.analysis_type{1});
                   bla = rmfield(bla,'cumtapcnt');
                   bla.erppow = [];
                   bla.itpc = [];
                   bla = orderfields(bla,fields(faux(1,1)));
                else
                    cfgr.parameter = 'powspctrm';
                    bla = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
                    if isfield(bla,'cumtapcnt'),bla = rmfield(bla,'cumtapcnt');end
                    if isfield(TFRav.(fTFR{ff}).(p.analysis_type{1}),'erppow')
                        cfgr.parameter = 'erppow';
                        bla2 = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
                        bla.erppow = bla2.erppow;
                    else
                        bla.erppow = [];
                    end
                    if isfield(TFRav.(fTFR{ff}).(p.analysis_type{1}),'itpc')
                        bla.itpc = TFRav.(fTFR{ff}).(p.analysis_type{1}).itpc;
                    else
                        bla.itpc = [];
                    end
                    
                end
                faux(s,ff) = bla;
            end
            s=s+1;
%          catch
%              sprintf('no data for suj %d',tk)
%          end
    end
    cfgga.keepindividual = 'yes';
    for ff=1:length(fTFR)
        cfgga.parameter = 'powspctrm';
        str_GA = ' = ft_freqgrandaverage(cfgga';
        for ss = 1:size(faux(:,ff),1)
            str_GA = [str_GA, ',faux(', num2str(ss), ',ff)'];
        end
        str_GA1 = ['GA.(fTFR{ff})' str_GA,');'];
        eval(str_GA1)
        if ~isempty(faux(1,ff).erppow)
            cfgga.parameter = 'erppow';
            str_GA2 = ['auxGA' str_GA,');'];
            eval(str_GA2)
            GA.(fTFR{ff}).erppow = auxGA.erppow;
        else
            GA.(fTFR{ff}).erppow = [];
        end
        if ~isempty(faux(1,ff).itpc)
            cfgga.parameter = 'itpc';
            str_GA3 = ['auxGA' str_GA,');'];
            eval(str_GA3)
            GA.(fTFR{ff}).itpc = auxGA.itpc;
        else
            GA.(fTFR{ff}).itpc = [];
        end
    end

%     mirindx         = mirrindex(GA.U_I.label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
  mirindx         = mirrindex(GA.U_I.label,[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
  doci = {'LU_I','LC_I','LU_unI','LC_unI','RU_I','RC_I','RU_unI','RC_unI','U_I','U_unI','C_I','C_unI','ALLUC'};
  
  for dc = 1:length(doci)
      GA.([doci{dc} 'ci']) = GA.(doci{dc});
      GA.([doci{dc} 'ci']).powspctrm = GA.(doci{dc}).powspctrm-GA.(doci{dc}).powspctrm(:,mirindx,:,:);
  end
    if b==1
        GAbsl = GA;
    else
        GAnobsl = GA;
    end
    clear faux GA
end
%%
cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';

 GAbsl.diffUIvsunI       = ft_math(cfgs,GAbsl.U_I,GAbsl.U_unI);
 GAbsl.diffCIvsunI       = ft_math(cfgs,GAbsl.C_I,GAbsl.C_unI);
    
 GAbsl.diffUIvsunIci       = ft_math(cfgs,GAbsl.U_Ici,GAbsl.U_unIci);
 GAbsl.diffCIvsunIci       = ft_math(cfgs,GAbsl.C_Ici,GAbsl.C_unIci);

 GAbsl.diffUCI       = ft_math(cfgs,GAbsl.U_I,GAbsl.C_I);
 GAbsl.diffUCunI       = ft_math(cfgs,GAbsl.U_unI,GAbsl.C_unI);
 
 GAbsl.diffLUCI       = ft_math(cfgs,GAbsl.LU_I,GAbsl.LC_I);
 GAbsl.diffLUCIci       = ft_math(cfgs,GAbsl.LU_Ici,GAbsl.LC_Ici);
 GAbsl.diffLUCunI       = ft_math(cfgs,GAbsl.LU_unI,GAbsl.LC_unI);
 GAbsl.diffLUCunIci       = ft_math(cfgs,GAbsl.LU_unIci,GAbsl.LC_unIci);
 
 GAbsl.diffLRUI       = ft_math(cfgs,GAbsl.LU_I,GAbsl.RU_I);
 GAbsl.diffLRCI       = ft_math(cfgs,GAbsl.LC_I,GAbsl.RC_I);
  GAbsl.diffLRUunI       = ft_math(cfgs,GAbsl.LU_unI,GAbsl.RU_unI);
 GAbsl.diffLRCunI       = ft_math(cfgs,GAbsl.LC_unI,GAbsl.RC_unI);
 
 GAbsl.diffRUCI       = ft_math(cfgs,GAbsl.RU_I,GAbsl.RC_I);
 GAbsl.diffRUCIci       = ft_math(cfgs,GAbsl.RU_Ici,GAbsl.RC_Ici);
 
 GAbsl.diffUCIci       = ft_math(cfgs,GAbsl.U_Ici,GAbsl.C_Ici);
 GAbsl.diffUCunIci       = ft_math(cfgs,GAbsl.U_unIci,GAbsl.C_unIci);
 GAbsl.diffUC_IvsUni      = ft_math(cfgs,GAbsl.diffUCI,GAbsl.diffUCunI);
 GAbsl.diffUC_IvsUnici      = ft_math(cfgs,GAbsl.diffUCIci,GAbsl.diffUCunIci);
% difffreq2.mask   = statUC.mask;
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
 % cfgp.xlim           = [-.4 .6];
%  cfgp.zlim           = [-1 1 ];
 cfgp.maskparameter  = 'mask';
   cfgp.maskalpha      = 1;
 %cfgp.parameter      = 'erppow';

  data =  GAbsl.ALLUC
%      data =GAbsl.C_Ici;`
    data.mask = statBslALLUC.mask;

figure,ft_multiplotTFR(cfgp,data)


%%
 load('cmapjp','cmap') 
%   cmap = cmocean('thermal');
% cmap = flipud(cbrewer('div','RdBu',128)); 
cfgp                = [];
cfgp.xlim           = [-.3 .75];
cfgp.zlim           = [-4 4];
% cfgp.ylim           = [8 30];
cfgp.colorbar       = 'no';
datatoplot      = {'U_Ici','C_Ici','diffUCIci','U_unIci','C_unIci','diffUCunIci','diffUC_IvsUnici'};
chnstoPlot      = {{'C3','CP3'},{'P5','P7','PO7'}};
datatoplot      = {'ALLUC'}
chnstoPlot      = {GAbsl.U_I.label};
for ch = 1:length(chnstoPlot)
    for ff = 1:length(datatoplot)
        cfgp.channel        = chnstoPlot{ch};
        fh = figure;
        ft_singleplotTFR(cfgp,GAbsl.(datatoplot{ff}))
        colormap(cmap),
        line([0 0],[2.1 30],'LineWidth',.25,'Color',[0 0 0])
%         vline(0,'k')
        figsize     = [8 4];
        axis([-.3 .75 2 30])
        axis square
        box off
        set(gca,'Title',[],'XTick',-0.3:.150:750,'XTickLabel',{'300','','0','','300','','600',''},...
            'YTick',0:5:30,'YTickLabel',{'0','','10','','20','','30'},'FontSize',6,'LineWidth',.25)
        ylabel('Frequency (Hz)','FontSize',8)
        xlabel('Time (ms)','FontSize',8)
%         tightfig
        hc = colorbar;
        hc.Position = [.85 .25 .05 .6];
        hc.LineWidth = .25;
        hc.Ticks = -4:2:4;
        hc.TickLabels = {'-4','','0','','4'};
          figsize  = [17.6/4 17.6/4*fh.Position(4)/fh.Position(3)];
        doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',...
            [datestr(now,'ddmmyy') '_TF_' datatoplot{ff} '_allchannels'],'600','opengl',figsize,1)
          %doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',...
        %    [datestr(now,'ddmmyy') '_TF_' datatoplot{ff} '_' strjoin('_',cfgp.channel)],'600','painters',figsize,1)
    end
end
%%
% simple spect5ra 500 ms and 500 ms after stim
t1 = [57];
t2 = [91];
ff = find(GAnobsl.ALLUC.freq>0 & GAnobsl.ALLUC.freq<31);
cmap1 = cbrewer('qual','Paired',12);
fh=figure;
   ,plot(GAnobsl.ALLUC.freq(ff),squeeze(nanmean(nanmean(GAnobsl.ALLUC.powspctrm(:,:,ff,t1),2),4)),'Color',cmap1(3,:),'LineWidth',.1)
hold on,plot(GAnobsl.ALLUC.freq(ff),squeeze(nanmean(nanmean(GAnobsl.ALLUC.powspctrm(:,:,ff,t2),2),4)),'Color',cmap1(7,:),'LineWidth',.1)
% hold on,plot(GAnobsl.ALLUC.freq(ff),squeeze(nanmean(nanmean(GAnobsl.ALLUC.powspctrm(:,:,:,t1),2),4))-squeeze(nanmean(nanmean(GAnobsl.ALLUC.powspctrm(:,:,:,t2),2),4)),'Color',[.7 .7 .7])

hold on,plot(GAnobsl.ALLUC.freq(ff),squeeze(nanmean(nanmean(nanmean(GAnobsl.ALLUC.powspctrm(:,:,ff,t1),1),2),4)),'Color',cmap1(4,:),'LineWidth',1)
hold on,plot(GAnobsl.ALLUC.freq(ff),squeeze(nanmean(nanmean(nanmean(GAnobsl.ALLUC.powspctrm(:,:,ff,t2),1),2),4)),'Color',cmap1(8,:),'LineWidth',1)
 axis square
box off
axis([0 30 0 20])
set(gca,'XTick',0:2:30,'XTickLabel',{'0','','','','','10','','','','','20','','','','','30'},...
    'YTick',0:5:20,'FontSize',6, 'Layer','top','LineWidth',.25)
xlabel('Frequency (Hz)','FontSize',8)
ylabel('Absolute Power','FontSize',8)
% tightfig
     figsize  = [17.6/4 17.6/4*fh.Position(4)/fh.Position(3)];
     
  doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',...
            [datestr(now,'ddmmyy') '_spectra_ALLcond_allchannels'],'600','opengl',figsize,1)
%%
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
statLRUI       = freqpermWS(GAbsl.LU_I,GAbsl.RU_I,elec,1000);
statLRCI       = freqpermWS(GAbsl.LC_I,GAbsl.RC_I,elec,1000);
statUCI       = freqpermWS(GAbsl.U_I,GAbsl.C_I,elec,1000);
statUCIci       = freqpermWS(GAbsl.U_Ici,GAbsl.C_Ici,elec,1000);
statUCunIci     = freqpermWS(GAbsl.U_unIci,GAbsl.C_unIci,elec,1000);
statUCIvsUnici  = freqpermWS(GAbsl.diffUCIci,GAbsl.diffUCunIci,elec,1000);
save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/stat_hann_' datestr(now,'ddmmyy')],'statUCIci','statUCunIci','statUCIvsUnici')
% statUC = freqpermWS(GAbsl.U_I,GAbsl.C_I,elec,1000)
 statBslALLUC = freqpermBSL(GAbsl.ALLUC,p.bsl,elec,1000)
 statBslALLUCci = freqpermBSL(GAbsl.ALLUCci,p.bsl,elec,500)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOPOPLOT SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('stat_hann_220218.mat')
clim            = [-2.5 2.5];
plotinterval    = [-.1 .8 .05];
bands           = [8 16;15 27;20 27;36 40];
bandnames       = {'alpha'}%,'beta','hbeta','lgamma'};
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
chnstoPlot      = {{'C3'},{'c4'}};
chnsLabel       = {'cL','cR'};
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
%         doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_unInf_' chnsLabel{ch} '_' bandnames{b}],figsize,1)
    
	end
end
for ch = 1:length(chnstoPlot)
    figure
    tp = topo_markCh(cfg_eeg,chnstoPlot{ch});
    tightfig
    doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',['Channs_' strjoin('_',chnstoPlot{ch})],[2 2],1)
end
    
  