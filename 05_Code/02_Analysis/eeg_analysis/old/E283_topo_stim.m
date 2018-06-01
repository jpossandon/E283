% Time-series analysis locked to tactile stimulation
%%
clear
E283_params                                 % basic experimental parameters               % 

for tk = p.subj
% tk = str2num(getenv('SGE_TASK_ID'));
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk));
    end
    
    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname','stimlock');    % single experiment/session parameters 
   
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata               
    
    % ERP topoplots for all four conditions 
    p.times                 = [600 1200];
    E283_base_trl_event_def_stim  
    % p.data                  = 'stim';
    p.rref                  = 'yes';
    p.plot                  = 1;
    p.analysis_type         = {'ICAem'};
    p.bsl                   = [-.4 -.25];
    p.interval              = [-.15 .6 .025]; % [start end step]
    p.colorlim              = [-10 10];
    p.analysisname          = 'stimLock';   
    p.colorlim              = [-10 10];
    
    p.keep                  = 'no';
    for e = 1:8%length(p.trls_stim)
        [ERP.(p.trls_stim{e})] = getERPsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,p.analysis_type{1},42,p.keep);
    end

    % difference plot U vs C pooled hand by mirroring
    % clear ERPs
    p.keep          = 'yes';
    for e = 1:8 %length(p.trls_stim)
        [ERPaux(e)] = getERPsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,p.analysis_type{1},42,p.keep);
        mirindx         = mirrindex(ERPaux(1).(p.analysis_type{1}).label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
        if ismember(e,[1,2,5,6])
            ERPaux(e).(p.analysis_type{1}).trial = ERPaux(e).(p.analysis_type{1}).trial(:,mirindx,:);
        end
    end

    difflabels      = {'U_I','C_I','U_unI','C_unI'};
    compidx         = [3 1;4 2;7 5;8 6];
   
    for f=1:size(compidx,1)
        ERP.(difflabels{f}).(p.analysis_type{1})   = ft_timelockanalysis([],ft_appenddata([],ERPaux(compidx(f,1)).(p.analysis_type{1}),ERPaux(compidx(f,2)).(p.analysis_type{1})));
    end

    difflabels      = {'U_Ici','C_Ici','U_unIci','C_unIci'};
    for e = 1:8
%         if ismember(e,[3,4,7,8])
%              ERPaux(e).(p.analysis_type{1}).trial = ERPaux(e).(p.analysis_type{1}).trial(:,mirindx,:);
%         end
        ERPaux(e).(p.analysis_type{1}).trial = ERPaux(e).(p.analysis_type{1}).trial-ERPaux(e).(p.analysis_type{1}).trial(:,mirindx,:);
    end
    for f=1:size(compidx,1)
        ERP.(difflabels{f}).(p.analysis_type{1})   = ft_timelockanalysis([],ft_appenddata([],ERPaux(compidx(f,1)).(p.analysis_type{1}),ERPaux(compidx(f,2)).(p.analysis_type{1})));
    end

    
    save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/ERP_' cfg_eeg.sujid '_stimlock'],'ERP','p')
end

%%
% Grand averages
clear
E283_params  
at                  = 1;
p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
cfgr                = [];
s=1
for tk = p.subj; % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlock','expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlock');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/ERP_' cfg_eeg.sujid '_stimlock'],'ERP')
    
    ERPall(s) = ERP;
    s=s+1;
end

fERP    = fields(ERPall);
cfgga.keepindividual = 'yes';
for ff=1:length(fERP)
    str_GA = 'GA.(fERP{ff}) = ft_timelockgrandaverage(cfgga'
    for ss = 1:length(ERPall)
        str_GA = [str_GA, ',ERPall(' num2str(ss) ').' fERP{ff} '.' p.analysis_type{1} ''];
    end
    str_GA = [str_GA,');'];
    eval(str_GA)
end
%%
mirindx         = mirrindex(GA.U_I.label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
        
cfgs            = [];
cfgs.parameter  = 'individual';
cfgs.operation  = 'subtract';
GA.UvsC_I       = ft_math(cfgs,GA.U_I,GA.C_I);
GA.UvsC_unI     = ft_math(cfgs,GA.U_unI,GA.C_unI);
GA.UvsC_Ici     = ft_math(cfgs,GA.U_Ici,GA.C_Ici);
GA.UvsC_unIci   = ft_math(cfgs,GA.U_unIci,GA.C_unIci);
GA.LUvsRU_I     = ft_math(cfgs,GA.LU_I,GA.RU_I);
GA.LCvsRC_I     = ft_math(cfgs,GA.LC_I,GA.RC_I);
GA.LUvsRU_unI   = ft_math(cfgs,GA.LU_unI,GA.RU_unI);
GA.LCvsRC_unI   = ft_math(cfgs,GA.LC_unI,GA.RC_unI);
GA.LUvsLC_I     = ft_math(cfgs,GA.LU_I,GA.LC_I);
GA.RUvsRC_I     = ft_math(cfgs,GA.RU_I,GA.RC_I);
GA.UvsC_IcivcUvsC_unIci     = ft_math(cfgs,GA.UvsC_Ici,GA.UvsC_unIci);

%%
% GA FIGURE
fERP    = fields(GA);
p.interval = [-.1 1 .02]
for ff=[27]%length(fERP)
    fh            = plot_topos(cfg_eeg,GA.(fERP{ff}),p.interval,p.bsl,[-3 3],['all ' fERP{ff} ' '  p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],...
             'tiff',['all_imlock_' fERP{ff}],1)
end

%%
% Multiplot figure

load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
cfgp.graphcolor  = 'rbgkmc';% cbrewer('qual','Set1',6);
figure
  ft_multiplotER(cfgp,GA.U_I,GA.C_I,GA.UvsC_I)
% ft_multiplotER(cfgp,GA.U_Ici,GA.C_Ici,GA.UvsC_Ici,GA.U_unIci,GA.C_unIci,GA.UvsC_unIci)
% ft_multiplotER(cfgp,ERPall(2).U_Ici.ICAem,ERPall(2).C_Ici.ICAem)

%  ft_multiplotER(cfgp,GA.LU_I,GA.LC_I,GA.LUvsLC_I,GA.RU_I,GA.RC_I,GA.RUvsRC_I)