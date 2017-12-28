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
for tk = p.subj; % subject number

    % Analysis parameters
    p.times_tflock              = [500 1000];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-.5 -.25]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'yes';	                
    p.cfgTFR.method             = 'mtmconvol';
    p.cfgTFR.taper              = 'hanning';
    % p.cfgTFR.width              = 5; 
    p.cfgTFR.output             = 'pow';	
    % p.cfgTFR.foi                = 4:2:35;	
    % p.cfgTFR.t_ftimwin          = .512*ones(1,length(p.cfgTFR.foi));
    p.cfgTFR.toi                = (-p.times_tflock(1):10:p.times_tflock(2))/1000;	

    p.cfgTFR.foi               = 4:1:40;	
    p.cfgTFR.pad                = 4;
    p.cfgTFR.t_ftimwin         = .250.*ones(1,length(p.cfgTFR.foi));
%     p.cfgTFR.t_ftimwin          = 3./p.cfgTFR.foi;
%     p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
%     plottp(p.cfgTFR)


    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Volumes/nibaldo/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk));
    end

    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname','stimlockTFR');       % single experiment/session parameters 
  
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])                         
    
    fieldstoav      = {'LU_I','RU_I','LC_I','RC_I','LU_unI','RU_unI','LC_unI','RC_unI','IM','U_I','C_I','U_unI','C_unI'};
    trigs           = [1,2,5,6,9,10,13,14,96];
    for f = 1:9
        [trls.(fieldstoav{f}),events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value',['==' num2str(trigs(f))]},p.times_tflock);            
    end
%     [trls.trlLU_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==1'},p.times_tflock);            
%     [trls.trlRU_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==2'},p.times_tflock);            
%     [trls.trlLC_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==5'},p.times_tflock);            
%     [trls.trlRC_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==6'},p.times_tflock);            
%     [trls.image,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==96'},p.times_tflock);  
    
    % Locked to stimulus
    at  = 1;
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/TFR_' p.analysis_type{at} '/'])

     for f = 1:9
         TFR.(fieldstoav{f})  = getTFRsfromtrl({cfg_eeg},{trls.(fieldstoav{f})},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%          eval(['[TFR_' fieldstoav{f} '] = getTFRsfromtrl({cfg_eeg},{trls.' fieldstoav{f} '},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);'])
     end

%      mirroring left stimulation data to make contra/ipsi plots
    mirindx         = mirrindex(TFR.LU_I.(p.analysis_type{1}).label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
    TFR.LU_I_mirr      = TFR.LU_I;
    TFR.LC_I_mirr      = TFR.LC_I;
    TFR.LU_I_mirr.(p.analysis_type{at}).powspctrm = TFR.LU_I.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
    TFR.LC_I_mirr.(p.analysis_type{at}).powspctrm = TFR.LC_I.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);

    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR.U_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RU_I.(p.analysis_type{1}),TFR.LU_I_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR.C_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RC_I.(p.analysis_type{1}),TFR.LC_I_mirr.(p.analysis_type{1})); %SAME

    TFR.LU_unI_mirr      = TFR.LU_unI;
    TFR.LC_unI_mirr      = TFR.LC_unI;
    TFR.LU_unI_mirr.(p.analysis_type{at}).powspctrm = TFR.LU_unI.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
    TFR.LC_unI_mirr.(p.analysis_type{at}).powspctrm = TFR.LC_unI.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);

    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR.U_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RU_unI.(p.analysis_type{1}),TFR.LU_unI_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR.C_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR.RC_unI.(p.analysis_type{1}),TFR.LC_unI_mirr.(p.analysis_type{1})); %SAME

    
    for f = 1:length(fieldstoav)
        TFRav.(fieldstoav{f}).(p.analysis_type{1}) = ft_freqdescriptives([], TFR.(fieldstoav{f}).(p.analysis_type{1}));
    end
    save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')

end

%%
% grand averages
at                  = 1;
p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
cfgr                = [];
p.bsl               = [-.5 -.25]; 
cfgr.baseline       = p.bsl;
cfgr.baselinetype   = 'db';
cfgr.keepindividual = 'yes';
s=1
for tk = p.subj; % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlockTFR','expfolder','/Volumes/nibaldo/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlockTFR');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg')
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
%%
% [freq1]             = ft_freqbaseline(cfgr,GA.LU_I);
% [freq2]             = ft_freqbaseline(cfgr,GA.LC_I);
% 

% freq1av   = ft_freqdescriptives([], GA.U_unI);
% freq2av   = ft_freqdescriptives([], GA.C_unI);
% freq3av   = ft_freqdescriptives([], GA.U_I);
% freq4av   = ft_freqdescriptives([], GA.C_I);
mirindx         = mirrindex(GA.U_I.label,[cfg_eeg.expfolder '/channels/mirror_chans']); 

GA.U_Ici             = GA.U_I;
GA.U_Ici.powspctrm   = GA.U_I.powspctrm-GA.U_I.powspctrm(:,mirindx,:,:);
GA.U_unIci             = GA.U_unI;
GA.U_unIci.powspctrm = GA.U_unI.powspctrm-GA.U_unI.powspctrm(:,mirindx,:,:);
GA.C_Ici             = GA.C_I;
GA.C_Ici.powspctrm   = GA.C_I.powspctrm-GA.C_I.powspctrm(:,mirindx,:,:);
GA.C_unIci             = GA.C_unI;
GA.C_unIci.powspctrm = GA.C_unI.powspctrm-GA.C_unI.powspctrm(:,mirindx,:,:);

%%
cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';
difffreq1       = ft_math(cfgs,GA.U_unI,GA.C_unI);
difffreq2       = ft_math(cfgs,GA.U_I,GA.C_I);
difffreq1       = ft_math(cfgs,GA.U_unIci,GA.C_unIci);
difffreq2       = ft_math(cfgs,GA.U_Ici,GA.C_Ici);
difffreq3      = ft_math(cfgs,difffreq2,difffreq1);
% difffreq1   = ft_freqdescriptives([],difffreq1);
% difffreq2   = ft_freqdescriptives([], difffreq2);
% difffreq3   = ft_freqdescriptives([], difffreq3);
% difffreq.mask   = statUC.mask;
%%
load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
% cfgp.trials     =4
%                  cfgp.baseline       = p.bsl;
%                   cfgp.baselinetype   = 'db';
% cfgp.ylim           = [0 40];
%  cfgp.xlim           = [-.75 1.25];
%         cfgp.zlim           = [-3 3];
%     cfgp.maskparameter = 'mask';
%       cfgp.maskalpha = .3
%     data = TFRav.C_I.ICAem;
%        data = ft_freqdescriptives([],GA.C_unIci);
        data =difffreq3;

figure,ft_multiplotTFR(cfgp,data)

%%
% Comparisons stat
% cfgr                = [];
% cfgr.baseline       = p.bsl;
% cfgr.baselinetype   = 'db';
% [freq1]             = ft_freqbaseline(cfgr, TFR.RU.(p.analysis_type{at}));
% [freq2]             = ft_freqbaseline(cfgr, TFR.RC.(p.analysis_type{at}));
%     find(ismember(difffreq2.label,{'P5', 'P7', 'PO7'}))
load(cfg_eeg.chanfile)
statUC = freqpermWS(GA.U_Ici,GA.C_Ici,elec,2000)
plot_topos_TFR(cfg_eeg,ft_freqdescriptives([],difffreq2),[-.3 .8 .025],[11 15],[],[-3 3],'alpha I U-C')
plot_topos_TFR(cfg_eeg,ft_freqdescriptives([],difffreq1),[-.3 .8 .025],[11 15],[],[-3 3],'alpha unI U-C')
plot_topos_TFR(cfg_eeg,ft_freqdescriptives([],difffreq2),[-.3 .8 .025],[18 24],[],[-2 2],'beta I U-C')
plot_topos_TFR(cfg_eeg,ft_freqdescriptives([],difffreq1),[-.3 .8 .025],[15 25],[],[-3 3],'beta unI U-C')
%%
figure
for ff=4:2:25
    fffs = ff:ff+2;
%  figure,plot(difffreq2.time,squeeze(mean(mean(difffreq2.powspctrm(:,[58,65,66],fffs,:),2),3)),'Color',[0.9 0.9 0.9])
% hold on,plot(difffreq2.time,squeeze(mean(mean(mean(difffreq2.powspctrm(:,[58,65,66],fffs,:),2),3))),'LineWidth',3,'Color',[0 0 0])
figure
Mval = squeeze(mean(mean(mean(difffreq2.powspctrm(:,[58,65,66],fffs,:),2),3)));
SE   = squeeze(std(mean(mean(difffreq2.powspctrm(:,[58,65,66],fffs,:),2),3)))./sqrt(size(difffreq2.powspctrm,1));
jbfill(difffreq2.time,Mval'+SE',Mval'-SE',[1 .7 .7],[1 .9 .9],1,.8)
hold on
plot(difffreq2.time,Mval,'LineWidth',2,'Color',[1 0 0])
% plot(difffreq2.time,squeeze(mean(mean(difffreq1.powspctrm(:,[58,65,66],fffs,:),2),3)),'Color',[1 0.9 0.9])
% hold on,plot(difffreq2.time,squeeze(mean(mean(m`ean(difffreq1.powspctrm(:,[58,65,66],fffs,:),2),3))),'LineWidth',3,'Color',[1 0 0])
hline(0)
title(sprintf('freq %d',difffreq2.freq(ff+1)))
ylim([-2.5 1])
end