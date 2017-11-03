% Time-series analysis locked to tactile stimulation
%%
clear
E283_params                                 % basic experimental parameters               % 
fmodel                      = 1;            % wich glm model
E283_models                                 % the models

%%
% subject configuration and data
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

    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
     mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/' cfg_eeg.sujid '/'])
      E283_base_trl_event_def_stim                                        % trial configuration  
%%
% ERP topoplots for all four conditions 

p.colorlim      = [-10 10];
for e = 1:8%length(p.trls_stim)
    [ERP.(p.trls_stim{e})] = getERPsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,p.analysis_type{1},p.keep);
    if p.plot
        fh = plot_topos(cfg_eeg,ERP.(p.trls_stim{e}).(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' imlock ' p.trls_stim{e} ' / ' p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))]);
         saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' cfg_eeg.sujid '_imlock_' p.analysis_type{1} '_' p.trls_stim{e}],'fig')
        close(fh)
    end
end

%%    
% difference plots
cfgs            = [];
cfgs.parameter  = 'avg';
cfgs.operation  = 'subtract';

difflabels      = {'LU_IvsLC_I','RU_IvsRC_I','LU_IvsRU_I','LC_IvsRC_I',...
    'LU_unIvsLC_unI','RU_unIvsRC_unI','LU_unIvsRU_unI','LC_unIvsRC_unI'};
compidx         = [1 2;3 4;1 3;2 4;5 6;7 8;5 7;6 8];
p.colorlim      = [-5 5];
erpfield        = fields(ERP);
for e = 1:length(difflabels)
    ERP.(difflabels{e}).(p.analysis_type{1})         = ft_math(cfgs,ERP.(erpfield{compidx(e,1)}).(p.analysis_type{1}),ERP.(erpfield{compidx(e,2)}).(p.analysis_type{1}));
    fh                  = plot_topos(cfg_eeg,ERP.(difflabels{e}).(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' imlock ' difflabels{e} ' / ' p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))]);
    saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' cfg_eeg.sujid '_imlock_' p.analysis_type{1} '_' difflabels{e}],'fig')
close(fh)
end
% [diffdata.name]           = difflabels{:};

%%
% difference plot U vs C pooled hand by mirroring
% clear ERPs
p.keep          = 'yes';
for e = 1:8 %length(p.trls_stim)
    [ERPaux(e)] = getERPsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,p.analysis_type{1},p.keep);
    mirindx         = mirrindex(ERPaux(1).(p.analysis_type{1}).label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
    if ismember(e,[1,2,5,6])
        ERPaux(e).(p.analysis_type{1}).trial = ERPaux(e).(p.analysis_type{1}).trial(:,mirindx,:);
    end
end

difflabels      = {'UncI','CrossI','UncunI','CrossunI'};
compidx         = [3 1;4 2;7 5;8 6];
p.colorlim      = [-10 10];

for f=1:size(compidx,1)
    ERP.(difflabels{f}).(p.analysis_type{1})   = ft_timelockanalysis([],ft_appenddata([],ERPaux(compidx(f,1)).(p.analysis_type{1}),ERPaux(compidx(f,2)).(p.analysis_type{1})));
    fh              = plot_topos(cfg_eeg,ERP.(difflabels{f}).(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' imlock ' (difflabels{f}) ' / ' p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))]);
    saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' cfg_eeg.sujid '_imlock_' p.analysis_type{1} '_' (difflabels{f})],'fig')
close(fh)
end

%continue here
% U vs C
ERP.UvsCI.(p.analysis_type{1}) = ft_math(cfgs,ERP.UncI.(p.analysis_type{1}),ERP.CrossI.(p.analysis_type{1}));
p.colorlim      = [-5 5];
fh            = plot_topos(cfg_eeg,ERP.UvsCI.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' imlock U minus C Inf / ' p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))]);
saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' cfg_eeg.sujid '_imlock_' p.analysis_type{1} '_UvsCI'],'fig')
ERP.UvsCunI.(p.analysis_type{1}) = ft_math(cfgs,ERP.UncunI.(p.analysis_type{1}),ERP.CrossunI.(p.analysis_type{1}));
p.colorlim      = [-5 5];
fh            = plot_topos(cfg_eeg,ERP.UvsCunI.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' imlock U minus C unInf / ' p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))]);
saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' cfg_eeg.sujid '_imlock_' p.analysis_type{1} '_UvsCunI'],'fig')

close all
save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/' cfg_eeg.sujid '/ERP_' cfg_eeg.sujid '_stimlock'],'ERP','p')
end

%%
% Grand averages
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
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/' cfg_eeg.sujid '/ERP_' cfg_eeg.sujid '_stimlock'],'ERP','p')
    
    ERPall(s) = ERP;
    s=s+1;
end

fERP    = fields(ERPall);
for ff=1:length(fERP)
    str_GA = 'GA.(fERP{ff}) = ft_timelockgrandaverage([]'
    for ss = 1:length(ERPall)
        str_GA = [str_GA, ',ERPall(' num2str(ss) ').' fERP{ff} '.' p.analysis_type{1} ''];
    end
    str_GA = [str_GA,');'];
    eval(str_GA)
end

%%
% GA FIGURE
fERP    = fields(GA);
for ff=1:length(fERP)
fh            = plot_topos(cfg_eeg,GA.(fERP{ff}),p.interval,p.bsl,p.colorlim,['all ' fERP{ff} ' '  p.analysis_type{1} ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))]);
saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/figures/GA/all_imlock_' p.analysis_type{1} '_' fERP{ff}],'fig')
close(fh)
end

%%
% Multiplot figure

load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';

ft_multiplotER(cfgp,GA.LU,GA.LC,GA.LUvsLC)