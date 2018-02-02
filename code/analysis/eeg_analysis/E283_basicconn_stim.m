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
for tk = p.subj % subject number

      % Analysis parameters
    p.times_tflock              = [600 1100];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-.400 -.150]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.output             = 'fourier';
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'yes';
    p.cfgTFR.toi                = (-p.times_tflock(1):20:p.times_tflock(2))/1000;	
    p.cfgTFR.foi                = 8:1:40;%2.^(3:.125:5);%(2.^([3:.25:5.25]));% %6:1:40	
    p.cfgTFR.method             = 'mtmconvol';%'wavelet';%%
    
    % wavelet parameters       
    %time of a cycle in ms times number of cicles divide by two for baseline limit :1/fo*width/2 
    %bandwith sf = f0/width   st = 1/(2*pi*f0/width)
%     p.cfgTFR.width              = 4; 
    
    % single tapers
    p.cfgTFR.taper              = 'hanning';%'dpss';
    p.cfgTFR.pad                = 4;
   
%    single tapers

%    p.cfgTFR.taper              = 'dpss';
%    p.cfgTFR.t_ftimwin          = 3./p.cfgTFR.foi;
%    p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
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
                                            'analysisname','stimlockcoh');       % single experiment/session parameters 
  
     load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])                         
    
    fieldstoav      = {'LU_I','RU_I','LC_I','RC_I','LU_unI','RU_unI','LC_unI','RC_unI','U_I','C_I','U_unI','C_unI'};%,'IM'
    trigs           = [1,2,5,6,9,10,13,14,96];
    for f = 1:8
        [trls.(fieldstoav{f}),events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value',['==' num2str(trigs(f))]},p.times_tflock);            
    end
    
    % Locked to stimulus
    at  = 1;
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/TFR_' p.analysis_type{at} '/'])
    bandFreqs   = [9 15;16 26];
    cfgcoh=[];     
    cfgcoh.method       = 'coh';
    cfgcoh.complex      = 'complex';
    bands               = {'alpha','beta'};
    for b = 1:2
        p.cfgTFR.foi        	= bandFreqs(b,1):1:bandFreqs(b,2);
        p.cfgTFR.t_ftimwin      = .250*ones(1,length(p.cfgTFR.foi));      %    single taper

        for f = 1:8
            auxfour                = getTFRsfromtrl({cfg_eeg},{trls.(fieldstoav{f})},...
                                p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
            imcoh.(bands{b}).(fieldstoav{f}) = ft_connectivityanalysis(cfgcoh,auxfour.ICAem);
            imcoh.(bands{b}).(fieldstoav{f}).cohspctrm = mean(imcoh.(bands{b}).(fieldstoav{f}).cohspctrm,3);
%             imcoh.(bands{b}).(fieldstoav{f}).Zcohspctrm = imcoh.(bands{b}).(fieldstoav{f}).cohspctrm./...
%                 abs(imcoh.(bands{b}).(fieldstoav{f}).cohspctrm).*atanh(abs(imcoh.(bands{b}).(fieldstoav{f}).cohspctrm));
%             coherency./abs(coherency).*atanh(abs(coherency));
            imcoh.(bands{b}).(fieldstoav{f}).freq = mean(imcoh.(bands{b}).(fieldstoav{f}).freq) ;
        end
    end
    
    mirindx         = mirrindex(imcoh.alpha.LU_I.label,[cfg_eeg.expfolder '/channels/mirror_chans']);
    for b = 1:2
        LU_I_mirr   = imcoh.(bands{b}).LU_I;
        LC_I_mirr   = imcoh.(bands{b}).LC_I;
        LU_I_mirr.cohspctrm = imcoh.(bands{b}).LU_I.cohspctrm(mirindx,mirindx,:,:);
        LC_I_mirr.cohspctrm = imcoh.(bands{b}).LC_I.cohspctrm(mirindx,mirindx,:,:);

        imcoh.(bands{b}).U_I           = imcoh.(bands{b}).RU_I;
        imcoh.(bands{b}).C_I           = imcoh.(bands{b}).RC_I;
        imcoh.(bands{b}).U_I.cohspctrm  = mean(cat(3,imcoh.(bands{b}).RU_I.cohspctrm,LU_I_mirr.cohspctrm),3);
        imcoh.(bands{b}).C_I.cohspctrm  = mean(cat(3,imcoh.(bands{b}).RC_I.cohspctrm,LC_I_mirr.cohspctrm),3);

        LU_unI_mirr   = imcoh.(bands{b}).LU_unI;
        LC_unI_mirr   = imcoh.(bands{b}).LC_unI;
        LU_unI_mirr.cohspctrm = imcoh.(bands{b}).LU_unI.cohspctrm(mirindx,mirindx,:,:);
        LC_unI_mirr.cohspctrm = imcoh.(bands{b}).LC_unI.cohspctrm(mirindx,mirindx,:,:);

        imcoh.(bands{b}).U_unI           = imcoh.(bands{b}).RU_unI;
        imcoh.(bands{b}).C_unI           = imcoh.(bands{b}).RC_I;
        imcoh.(bands{b}).U_unI.cohspctrm  = mean(cat(3,imcoh.(bands{b}).RU_unI.cohspctrm,LU_unI_mirr.cohspctrm),3);
        imcoh.(bands{b}).C_unI.cohspctrm  = mean(cat(3,imcoh.(bands{b}).RC_unI.cohspctrm,LC_unI_mirr.cohspctrm),3);
    end
    save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/' cfg_eeg.sujid '_imcoh_stim_' p.analysis_type{at}],'imcoh','cfg_eeg','p')
%     save(['/Users/jossando/trabajo/E283/analysis/' cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')

end

%%
% grand averages
clear
E283_params
at                  = 1;
p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 

MACpath = '/Users/jossando/trabajo/E283/';
clear faux
s=1
% MACpath = '/Volumes/nibaldo/trabajo/E283/';
for tk = p.subj; % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlockcoh','expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlockcoh');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/' cfg_eeg.sujid '_imcoh_stim_' p.analysis_type{at}],'imcoh')
    bands               = {'alpha','beta'};
    for b = 1:2
        fTFR    = fields(imcoh.(bands{b}));
            for ff=1:length(fTFR)
%         faux(s,ff) = ft_freqbaseline(cfgr,imcoh.(fTFR{ff}).(p.analysis_type{1}));
%         if ~isfield(imcoh.(fTFR{ff}),'elec')
%             imcoh.(fTFR{ff}).elec = faux(s-1,ff).elec;
%             imcoh.(fTFR{ff})      = orderfields(imcoh.(fTFR{ff}),faux(s-1,ff));
%         end
%         faux(s,ff) =imcoh.(fTFR{ff});

             faux(:,:,:,:,s,ff,b) = imcoh.(bands{b}).(fTFR{ff}).cohspctrm;
            end
        
    end
    s=s+1;
end
%%
%     faux = squeeze(nanmean(faux,5));
    for b = 1:2
        for ff=1:length(fTFR)
            GA.(bands{b}).(fTFR{ff}) = imcoh.(bands{b}).(fTFR{ff});
            GA.(bands{b}).(fTFR{ff}) = rmfield(GA.(bands{b}).(fTFR{ff}),'cohspctrm');
            coherency                = squeeze(faux(:,:,:,:,:,ff,b));
            GA.(bands{b}).(fTFR{ff}).coherency = coherency;
            GA.(bands{b}).(fTFR{ff}).avgImagCoherency = mean(imag(coherency),4);
            
%             anlstdCohy =(1-abs(mean(coherency,4)).^2).*atanh(abs(mean(coherency,4))).^2./abs(mean(coherency,4)).^2;
%             stdZcoh    = sqrt(1./2./size(coherency,4)THIS IS WRONG THIS N CORRESPOND TO THE NUMBER OF OBSERVATION TO CALCULATE THE COHERENCY.*(anlstdCohy.*cos(angle(mean(coherency,4))).^2+sin(angle(mean(coherency,4))).^2));
%             avgZcoherency = mean(coherency./abs(coherency).*atanh(abs(coherency)),4);
            
            GA.(bands{b}).(fTFR{ff}).avgZcoherency = mean(coherency./...
                abs(coherency).*atanh(abs(coherency)),4);
            GA.(bands{b}).(fTFR{ff}).dimord = 'chan_chan_time_subj';
        end
    end
%     clear faux imcoh
%%
diffFields = {'LU_I','RU_I';'LU_I','LC_I';'RU_I','RC_I';'LC_I','RC_I';...
    'LU_unI','RU_unI';'LU_unI','LC_unI';'RU_unI','RC_unI';'LC_unI','RC_unI';...
    'U_I','C_I';'U_unI','C_unI'};
diffName   = {'LvsR_U_I','LUvsC_I','RUvsC_I','LvsR_C_I',...
    'LvsR_U_unI','LUvsC_unI','RUvsC_unI','LvsR_C_unI',...
    'UvsC_I','UvsC_unI'};
for b = 1:2
    for ff = 1:length(diffName)
        GAdiff.(bands{b}).(diffName{ff})   = GA.(bands{b}).(diffFields{ff,1});
%         [H,P,~,STATS] = ttest(permute(GA.(bands{b}).(diffFields{ff,1}).cohspctrm,[4,1,2,3]),permute(GA.(bands{b}).(diffFields{ff,2}).cohspctrm,[4,1,2,3]));
%         [GAdiff.(bands{b}).(diffName{ff}).H,GAdiff.(bands{b}).(diffName{ff}).p,GAdiff.(bands{b}).(diffName{ff}).t] = ...
%             deal(squeeze(H),squeeze(P),squeeze(STATS.tstat));
%            GAdiff.(bands{b}).(diffName{ff}).avgZcoherency = mean(GA.(bands{b}).(diffFields{ff,1}).coherency./...
%                abs(GA.(bands{b}).(diffFields{ff,1}).coherency).*...
%                atanh(abs(GA.(bands{b}).(diffFields{ff,1}).coherency))-...
%                GA.(bands{b}).(diffFields{ff,2}).coherency./...
%                abs(GA.(bands{b}).(diffFields{ff,2}).coherency).*...
%                atanh(abs(GA.(bands{b}).(diffFields{ff,2}).coherency)),4)

%         GAdiff.(bands{b}).(diffName {ff}).cohspctrm = squeeze(mean(GA.(bands{b}).(diffFields{ff,1}).cohspctrm-GA.(bands{b}).(diffFields{ff,2}).cohspctrm,4));
        GAdiff.(bands{b}).(diffName{ff}).dimord = 'chan_chan_time';
        fprintf('\n%s %s p<0.001 = %2.4f',bands{b},diffName{ff},sum(GAdiff.(bands{b}).(diffName {ff}).p(:)<.001)./sum(~isnan(GAdiff.(bands{b}).(diffName {ff}).p(:))))

    end
end
% save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/' 'GA_' datestr(now,'ddmmyy') ],'GA','GAdiff')
%%

    
%%
alfa     = .001;
collim   = [-.15 .15];
times    = [-.1 .9 .05];

for b = 1
    for ff = 1:length(diffName)
        auxGA   = GAdiff.(bands{b}).(diffName{ff});
        fh = seq_cohplots(cfg_eeg,auxGA,times,collim,alfa);
        fh.Name = diffName{ff};
    end
end
%HACER MODELOS SEPARADOS POR SACADAS 1,2,3 AND TARGET FIXATOIN
%%
alfa     = .0005;
collim   = [-.3 .3];
times    = [-.2 .7 .05];
for b = 1
    fTFR    = fields(GA.(bands{b}));
    for ff = 9:10%1:length(fTFR)
        auxGA   = GA.(bands{b}).(fTFR{ff});
        auxGA.p = nan(size(auxGA.cohspctrm));
        auxGA.p(abs(auxGA.cohspctrm)>.1) = .0001;
        fh = seq_cohplots(cfg_eeg,auxGA,times,collim,alfa);
        fh.Name = fTFR{ff};
    end
end
%%
% stat
diffFields = {'LU_I','RU_I';'LU_I','LC_I';'RU_I','RC_I';'LC_I','RC_I';...
    'LU_unI','RU_unI';'LU_unI','LC_unI';'RU_unI','RC_unI';'LC_unI','RC_unI';...
    'U_I','C_I';'U_unI','C_unI'};
diffName   = {'LvsR_U_I','LUvsC_I','RUvsC_I','LvsR_C_I',...
    'LvsR_U_unI','LUvsC_unI','RUvsC_unI','LvsR_C_unI',...
    'UvsC_I','UvsC_unI'};
alfa     = .05;
collim   = [-.7 .7];
times    = [-.2 .9 .05];
for b = 1:2
    for ff =1:length(diffName)
%         auxpermT = coh_cluster(cfg_eeg,GA.(bands{b}).(diffFields{ff,1}).cohspctrm,...
%             GA.(bands{b}).(diffFields{ff,2}).cohspctrm,1000); 
%         auxpermT.time = GA.(bands{b}).(diffFields{ff,1}).time;
%         result.(bands{b}).(diffName{ff}) = auxpermT;
        fh      = seq_cohplots(cfg_eeg,result.(bands{b}).(diffName{ff}),times,collim,alfa);
        fh.Name = diffName{ff};
        figsize     = [17.6 17.6*fh.Position(4)/fh.Position(3)];

        doimage(gcf,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' bands{b} '_' diffName{ff}],figsize,1)
    
    end
%     save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/diff_stat'],'result')
end