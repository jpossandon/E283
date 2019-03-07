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
load(fullfile(MACpath,'07_Analysis','03_Eye','eyedata','tgtPos.mat')) 
for tk = p.subj % subject number

     % Analysis parameters
    p.times_tflock              = [800 1200];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-.650 -.350]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'yes';
    p.cfgTFR.output             = 'fourier';
    p.cfgTFR.toi                = (-p.times_tflock(1):15:p.times_tflock(2))/1000;	
    %p.cfgTFR.foi                = 8:1:40;%2.^(3:.125:7);%(2.^([3:.25:5.25]));% %6:1:40	
    p.cfgTFR.method             = 'mtmconvol';%'wavelet';%%
    
    % wavelet parameters       
    %time of a cycle in ms times number of cicles divide by two for baseline limit :1/fo*width/2 
    %bandwith sf = f0/width   st = 1/(2*pi*f0/width)
%     p.cfgTFR.width              = 4; 
    
    % single tapers
%      if any(strfind(p.analysisname,'hann'))
%           p.cfgTFR.foi                = 8:1:30;
        p.cfgTFR.taper              = 'hanning';%'dpss';
        p.cfgTFR.pad                = 4;
%         p.cfgTFR.t_ftimwin          = 3./p.cfgTFR.foi;      %    single taper
%      elseif any(strfind(p.analysisname,'dpss'))
%         %p.cfgTFR.foi                = 2.^(5:.125:7);
%           p.cfgTFR.foi                = 25:5:150;
%         p.cfgTFR.taper              = 'dpss';
%         p.cfgTFR.t_ftimwin          = 5./p.cfgTFR.foi;
%         p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
%         %plottp(p.cfgTFR)
%      end
     

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
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/TFR_' p.analysis_type{at} '/'])
    bandFreqs   = [9 15;16 26];
    cfgcoh=[];     
    cfgcoh.method       = 'coh';
    cfgcoh.complex      = 'complex';
    bands               = {'alpha'};%,'beta'};
    conn_meas           = {'imcoh'};%{'imcoh','wpli'};
    fieldstoav = fields(trls);
    
    for b = 1:length(bands)
        p.cfgTFR.foi        	= bandFreqs(b,1):1:bandFreqs(b,2);
        p.cfgTFR.t_ftimwin      = 3./p.cfgTFR.foi;      %    single taper
        for f = 1:length(fieldstoav)
            auxfour                = getTFRsfromtrl({cfg_eeg},{trls.(fieldstoav{f})},...
                                p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
            for connm = conn_meas
                if strcmp(connm,'imcoh')
                    cfgcoh.method       = 'coh';
                    cfgcoh.complex      = 'complex';
                elseif strcmp(connm,'wpli')
                    cfgcoh              = [];
                    cfgcoh.method       = 'wpli';
                end
                conn.(connm{:}).(bands{b}).(fieldstoav{f}) = ft_connectivityanalysis(cfgcoh,auxfour.ICAem);
                conn.(connm{:}).(bands{b}).(fieldstoav{f}).([cfgcoh.method 'spctrm']) = mean(conn.(connm{:}).(bands{b}).(fieldstoav{f}).([cfgcoh.method 'spctrm']),3);  % mean across frequencies
                conn.(connm{:}).(bands{b}).(fieldstoav{f}).freq = mean(conn.(connm{:}).(bands{b}).(fieldstoav{f}).freq) ;
            end
        end
    end
    
    mirindx         = mirrindex(conn.imcoh.alpha.LU_I.label,[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']);
    for connm = conn_meas
        if strcmp(connm,'imcoh')
            cfgcoh.method       = 'coh';
        elseif strcmp(connm,'wpli')
            cfgcoh.method       = 'wpli';
        end
        for b = 1:length(bands)
            LU_I_mirr   = conn.(connm{:}).(bands{b}).LU_I;
            LC_I_mirr   = conn.(connm{:}).(bands{b}).LC_I;
            LU_I_mirr.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LU_I.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);
            LC_I_mirr.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LC_I.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);

            conn.(connm{:}).(bands{b}).U_I           = conn.(connm{:}).(bands{b}).RU_I;
            conn.(connm{:}).(bands{b}).C_I           = conn.(connm{:}).(bands{b}).RC_I;
            conn.(connm{:}).(bands{b}).U_I.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RU_I.([cfgcoh.method 'spctrm']),LU_I_mirr.([cfgcoh.method 'spctrm'])),3);
            conn.(connm{:}).(bands{b}).C_I.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RC_I.([cfgcoh.method 'spctrm']),LC_I_mirr.([cfgcoh.method 'spctrm'])),3);

            LU_unI_mirr   = conn.(connm{:}).(bands{b}).LU_unI;
            LC_unI_mirr   = conn.(connm{:}).(bands{b}).LC_unI;
            LU_unI_mirr.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LU_unI.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);
            LC_unI_mirr.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LC_unI.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);

            conn.(connm{:}).(bands{b}).U_unI           = conn.(connm{:}).(bands{b}).RU_unI;
            conn.(connm{:}).(bands{b}).C_unI           = conn.(connm{:}).(bands{b}).RC_unI;
            conn.(connm{:}).(bands{b}).U_unI.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RU_unI.([cfgcoh.method 'spctrm']),LU_unI_mirr.([cfgcoh.method 'spctrm'])),3);
            conn.(connm{:}).(bands{b}).C_unI.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RC_unI.([cfgcoh.method 'spctrm']),LC_unI_mirr.([cfgcoh.method 'spctrm'])),3);
            
            LU_unI_mirr_L   = conn.(connm{:}).(bands{b}).LU_unI_L;
            LC_unI_mirr_L   = conn.(connm{:}).(bands{b}).LC_unI_L;
            LU_unI_mirr_L.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LU_unI_L.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);
            LC_unI_mirr_L.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LC_unI_L.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);

            
            LU_unI_mirr_R   = conn.(connm{:}).(bands{b}).LU_unI_R;
            LC_unI_mirr_R   = conn.(connm{:}).(bands{b}).LC_unI_R;
            LU_unI_mirr_R.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LU_unI_R.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);
            LC_unI_mirr_R.([cfgcoh.method 'spctrm']) = conn.(connm{:}).(bands{b}).LC_unI_R.([cfgcoh.method 'spctrm'])(mirindx,mirindx,:,:);

            conn.(connm{:}).(bands{b}).U_unI_ipsi        = conn.(connm{:}).(bands{b}).RU_unI_R;
            conn.(connm{:}).(bands{b}).U_unI_contra      = conn.(connm{:}).(bands{b}).RU_unI_L;
            conn.(connm{:}).(bands{b}).C_unI_contra      = conn.(connm{:}).(bands{b}).RC_unI_L;
            conn.(connm{:}).(bands{b}).C_unI_ipsi        = conn.(connm{:}).(bands{b}).RC_unI_R;
           
            conn.(connm{:}).(bands{b}).U_unI_ipsi.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RU_unI_R.([cfgcoh.method 'spctrm']),LU_unI_mirr_L.([cfgcoh.method 'spctrm'])),3);
            conn.(connm{:}).(bands{b}).U_unI_contra.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RU_unI_L.([cfgcoh.method 'spctrm']),LU_unI_mirr_R.([cfgcoh.method 'spctrm'])),3);
            conn.(connm{:}).(bands{b}).C_unI_contra.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RC_unI_L.([cfgcoh.method 'spctrm']),LC_unI_mirr_R.([cfgcoh.method 'spctrm'])),3);
            conn.(connm{:}).(bands{b}).C_unI_ipsi.([cfgcoh.method 'spctrm'])  = mean(cat(3,conn.(connm{:}).(bands{b}).RC_unI_R.([cfgcoh.method 'spctrm']),LC_unI_mirr_L.([cfgcoh.method 'spctrm'])),3);
        end
    end
    save([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/imcoh/' cfg_eeg.sujid '_imcoh_stim_' p.analysis_type{at}],'conn','cfg_eeg','p')
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
    load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/imcoh/' cfg_eeg.sujid '_imcoh_stim_' p.analysis_type{at}],'conn')
    bands               = {'alpha','beta'};
    for b = 1
        fTFR    = fields( conn.imcoh.(bands{b}));
            for ff=1:length(fTFR)
%         faux(s,ff) = ft_freqbaseline(cfgr,imcoh.(fTFR{ff}).(p.analysis_type{1}));
%         if ~isfield(imcoh.(fTFR{ff}),'elec')
%             imcoh.(fTFR{ff}).elec = faux(s-1,ff).elec;
%             imcoh.(fTFR{ff})      = orderfields(imcoh.(fTFR{ff}),faux(s-1,ff));
%         end
%         faux(s,ff) =imcoh.(fTFR{ff});

             faux(:,:,:,:,s,ff,b) =  imag(conn.imcoh.(bands{b}).(fTFR{ff}).cohspctrm);
            % fauxwpli(:,:,:,:,s,ff,b) =  conn.wpli.(bands{b}).(fTFR{ff}).wplispctrm;
            end
        
    end
    s=s+1
end
% clear conn
%
%     faux = squeeze(nanmean(faux,5));
    for b = 1%:2
        for ff=1:length(fTFR)
            GA.(bands{b}).(fTFR{ff}) = conn.imcoh.(bands{b}).(fTFR{ff});
          %  GA.(bands{b}).(fTFR{ff}) = conn.wpli.(bands{b}).(fTFR{ff});
            GA.(bands{b}).(fTFR{ff}).imagcohspctrm = squeeze(faux(:,:,:,:,:,ff,b));
           % GA.(bands{b}).(fTFR{ff}).wplispctrm = squeeze(fauxwpli(:,:,:,:,:,ff,b));
             coherency                = squeeze(faux(:,:,:,:,:,ff,b));
    %        GAimcoh.(bands{b}).(fTFR{ff}).coherency = coherency;
             GA.(bands{b}).(fTFR{ff}).imagcoherency = coherency;
     %       GAimcoh.(bands{b}).(fTFR{ff}).avgImagCoherency = mean(imag(coherency),4);
            
%             anlstdCohy =(1-abs(mean(coherency,4)).^2).*atanh(abs(mean(coherency,4))).^2./abs(mean(coherency,4)).^2;
%             stdZcoh    = sqrt(1./2./size(coherency,4)THIS IS WRONG THIS N CORRESPOND TO THE NUMBER OF OBSERVATION TO CALCULATE THE COHERENCY.*(anlstdCohy.*cos(angle(mean(coherency,4))).^2+sin(angle(mean(coherency,4))).^2));
%             avgZcoherency = mean(coherency./abs(coherency).*atanh(abs(coherency)),4);
            
   %         GAimcoh.(bands{b}).(fTFR{ff}).avgZcoherency = mean(coherency./...
    %            abs(coherency).*atanh(abs(coherency)),4);
            GAimcoh.(bands{b}).(fTFR{ff}).dimord = 'chan_chan_time_subj';
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
         [H,P,~,STATS] = ttest(permute(GA.(bands{b}).(diffFields{ff,1}).imagcoherency,[4,1,2,3]),permute(GA.(bands{b}).(diffFields{ff,2}).imagcoherency,[4,1,2,3]));
         [GAdiff.(bands{b}).(diffName{ff}).H,GAdiff.(bands{b}).(diffName{ff}).p,GAdiff.(bands{b}).(diffName{ff}).t] = ...
             deal(squeeze(H),squeeze(P),squeeze(STATS.tstat));
%            GAdiff.(bands{b}).(diffName{ff}).avgZcoherency = mean(GA.(bands{b}).(diffFields{ff,1}).coherency./...
%                abs(GA.(bands{b}).(diffFields{ff,1}).coherency).*...
%                atanh(abs(GA.(bands{b}).(diffFields{ff,1}).coherency))-...
%                GA.(bands{b}).(diffFields{ff,2}).coherency./...
%                abs(GA.(bands{b}).(diffFields{ff,2}).coherency).*...
%                atanh(abs(GA.(bands{b}).(diffFields{ff,2}).coherency)),4)

         GAdiff.(bands{b}).(diffName {ff}).cohspctrm = squeeze(mean(GA.(bands{b}).(diffFields{ff,1}).imagcoherency-GA.(bands{b}).(diffFields{ff,2}).imagcoherency,4));
        GAdiff.(bands{b}).(diffName{ff}).dimord = 'chan_chan_time';
        fprintf('\n%s %s p<0.001 = %2.4f',bands{b},diffName{ff},sum(GAdiff.(bands{b}).(diffName {ff}).p(:)<.001)./sum(~isnan(GAdiff.(bands{b}).(diffName {ff}).p(:))))

    end
end
clear faux
% save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/' 'GA_' datestr(now,'ddmmyy') ],'GA','GAdiff')
%%

    
%%
alfa     = .0005;
collim   = [-.1 .1];
times    = [-.1 .9 .05];

for b = 1
    for ff = [9]
        auxGA   = GAdiff.(bands{b}).(diffName{ff});
        fh = seq_cohplots(cfg_eeg,auxGA,times,collim,alfa);
        fh.Name = diffName{ff};
    end
end
%HACER MODELOS SEPARADOS POR SACADAS 1,2,3 AND TARGET FIXATOIN
%%
alfa     = .0005;
collim   = [-.15 .15];
times    = [-.2 .7 .05];
for b = 1
    fTFR    = fields(GA.(bands{b}));
    for ff = 1:4%1:length(fTFR)
        auxGA   = GA.(bands{b}).(fTFR{ff});
        auxGA.cohspctrm   = mean(GA.(bands{b}).(fTFR{ff}).imagcoherency,4 );
        auxGA.p = nan(size(auxGA.cohspctrm));
        auxGA.p(abs(auxGA.cohspctrm)>.085) = .0001;
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

for b = 1
    for ff =[9 10] %1:length(diffName)
         auxpermT = coh_cluster(cfg_eeg,GA.(bands{b}).(diffFields{ff,1}).imagcoherency,...
             GA.(bands{b}).(diffFields{ff,2}).imagcoherency,10); 
         auxpermT.time = GA.(bands{b}).(diffFields{ff,1}).time;
         result.(bands{b}).(diffName{ff}) = auxpermT;
        fh      = seq_cohplots(cfg_eeg,result.(bands{b}).(diffName{ff}),times,collim,alfa);
        fh.Name = diffName{ff};
        figsize     = [17.6 17.6*fh.Position(4)/fh.Position(3)];

%         doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' bands{b} '_' diffName{ff}],figsize,1)
    
    end
%     save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/imcoh/diff_stat'],'result')
end











%%
%%
% connectivity
dattocomp = {'U_unI','C_unI';'U_I','C_I';'U_unI_ipsi','C_unI_contra';'U_unI_contra','C_unI_ipsi'};%'U_Ici','C_Ici';
             %'LU_I','LC_I';'RU_I','RC_I';
             %'LU_unI','LC_unI';'RU_unI','RC_unI';}
             %'LU_unI','RU_unI';'LU_unIci','RU_unIci';'LC_unI','RC_unI';'LC_unIci','RC_unIci';
             %'LU_unI','LC_unI';'LU_unIci','LC_unIci';'RU_unI','RC_unI';'RU_unIci','RC_unIci';
             %'U_unI','C_unI';'U_unIci','C_unIci'};
bands               = {'alpha','beta'};
conn_meas           = 'imagcohspctrm';%'wplispctrm';%'imagcohspctrm';%,};
tiempo = GA.(bands{1}).(dattocomp{1,1}).time;
chans = {[50],[18],[58,65,66],[26,33,34]}
chanlab = {'C3','C4','P5-P7-PO7' ,'P6-P8-PO8'}
cluslab = {'pos','neg'};
combs = nchoosek(1:length(chans),2);
% combs([4,6],:) = []
N = size(GA.(bands{1}).(dattocomp{1,1}).imagcohspctrm,4);
setAbsoluteFigureSize
sppos=[3,6,7,9,10,13];
sppos=1:6

cb1      = cbrewer('qual','Set1',5)
 cb2      = cbrewer('qual','Pastel1',3)
lineColors = [cb1(1,:);cb2(1,:);cb1(2,:);cb2(2,:);cb1(1,:);cb2(1,:);cb1(1,:);cb2(1,:)]; 
%    
lc=1
for bb = 1:length(bands)
    for pp = 1:size(dattocomp,1)
        fh=figure;
        fh.Position = [1 1 800 100];

        for sp=1:size(combs,1)
            aux1 = squeeze(mean(mean(GA.(bands{bb}).(dattocomp{pp,1}).(conn_meas)(chans{combs(sp,1)},chans{combs(sp,2)},:,:),1),2));
            aux2 = squeeze(mean(mean(GA.(bands{bb}).(dattocomp{pp,2}).(conn_meas)(chans{combs(sp,1)},chans{combs(sp,2)},:,:),1),2));
          
             subplot(1,6,sppos(sp));,hold on
%             figure, hold on
          %  plot(tiempo,aux1,'Color',[1 .7 .7])
          %  plot(tiempo,aux2,'Color',[.7 .7 1])
            jbfill(tiempo,[mean(aux1,2)+std(aux1,1,2)/sqrt(N)]',[mean(aux1,2)-std(aux1,1,2)/sqrt(N)]',lineColors(lc,:),lineColors(lc,:),1,.7,.7);,hold on
            plot(tiempo,mean(aux1,2),'Color',lineColors(lc,:),'LineWidth',1);
            jbfill(tiempo,[mean(aux2,2)+std(aux2,1,2)/sqrt(N)]',[mean(aux2,2)-std(aux2,1,2)/sqrt(N)]',lineColors(lc+1,:),lineColors(lc+1,:),1,.7,.7);,hold on
            plot(tiempo,mean(aux2,2),'Color',lineColors(lc+1,:),'LineWidth',2);
            hline(0,'k-')
            axis([-.2 .8 -.1 .15])
            %axis square
            box off
%             title([chanlab{combs(sp,1)} ' - ' chanlab{combs(sp,2)}])
            
%             subplot(2,size(combs,1),sp+size(combs,1)),hold on
            hline(0,'r-')
            %plot(tiempo,aux1-aux2,'Color',[.7 .7 .7])
            jbfill(tiempo,[mean(aux1-aux2,2)+std(aux1-aux2,1,2)/sqrt(N)]',[mean(aux1-aux2,2)-std(aux1-aux2,1,2)/sqrt(N)]',[.3 .3 0],[.15 .15 0],1,.7,.7),hold on
            [h,p] = ttest(aux1',aux2');
         %   plot(tiempo,p)
            elec.channeighbstructmat = 0;
            [result(bb,pp,sp)] = regmodel2ndstat(reshape(aux1,[1 1 size(aux1)])-reshape(aux2,[1 1 size(aux1)]),tiempo,elec,2000,'signpermT','cluster'); %need check this with clusters and singnperm to
%             plot(tiempo,result(bb,pp,sp).pvalnc,'k')
%             plot(tiempo,median(aux1-aux2,2),':k')
            plot(tiempo,mean(aux1-aux2,2),'k')
            for pn = 1:2
                if ~isempty(result(bb,pp,sp).clusters.([cluslab{pn} 'clusters']))
                    for cc = 1:length(result(bb,pp,sp).clusters.([cluslab{pn} 'clusters']))
                        if result(bb,pp,sp).clusters.([cluslab{pn} 'clusters'])(cc).prob < 0.05
                            auxclus = find(result(bb,pp,sp).clusters.([cluslab{pn} 'clusterslabelmat'])==cc);
                        auxdiff = mean(aux1(auxclus,:)-aux2(auxclus,:),2);
                            if result(bb,pp,sp).clusters.([cluslab{pn} 'clusters'])(cc).prob < 0.05
                                plot(tiempo(auxclus),auxdiff,'r','LineWidth',1)
                            end
                            if strcmp(cluslab{pn},'neg')
                                [~,di] = min(mean(aux1(auxclus,:)-aux2(auxclus,:),2));
                            elseif strcmp(cluslab{pn},'pos')
                                [~,di] = max(mean(aux1(auxclus,:)-aux2(auxclus,:),2));
                            end
                             text(tiempo(auxclus(di)),auxdiff(di)+.02*sign(auxdiff(di)),sprintf('p = %1.3f',result(bb,pp,sp).clusters.([cluslab{pn} 'clusters'])(cc).prob))
                        end
                    end
                end
            end
             hline(0,'k:')
            vline(0,'k:')
            set(gca,'XTick',-.2:.2:.8,'XTickLabel',{'-.2','0','','.4','','.8'},'YTick',-.1:.05:0.15,'YTickLabels',{'-.1','','0','','','.15'})
           % axis([-.6 1 -.1 .2])
%         figsize     = [17.6 17.6*fh.Position(4)/fh.Position(3)];
%         doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' bands{bb} '_' dattocomp{pp,1} '_' dattocomp{pp,2}],figsize,1)
        
        end,
          figsize     = [17.6 17.6*fh.Position(4)/fh.Position(3)];
          doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' bands{bb} '_' dattocomp{pp,1} '_' dattocomp{pp,2}],'600','painters',figsize,1)
lc = lc+2;
    end
end

%%
% the channels
for ch = 1:size(combs,1)
    chantoplot = [chans{combs(ch,1)} chans{combs(ch,2)}];
    fh = figure;
     fh.Position = fh.Position/2;
    tp = topo_markCh(cfg_eeg,chantoplot);
    tightfig
    doimage(fh,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'tiff',['Channs ' strjoin(' X ',chanlab(combs(ch,:)))],'1200','opengl',[1 1],1)
end    
  
%%
% linear model the same
for bb = 1:length(bands)
    for sp=1:size(combs,1)
         auxT = [squeeze(mean(mean(GA.(bands{bb}).U_I.(conn_meas)(chans{combs(sp,1)},chans{combs(sp,2)},:,:),1),2))';...
        squeeze(mean(mean(GA.(bands{bb}).C_I.(conn_meas)(chans{combs(sp,1)},chans{combs(sp,2)},:,:),1),2))';...
        squeeze(mean(mean(GA.(bands{bb}).U_unI.(conn_meas)(chans{combs(sp,1)},chans{combs(sp,2)},:,:),1),2))';...
        squeeze(mean(mean(GA.(bands{bb}).C_unI.(conn_meas)(chans{combs(sp,1)},chans{combs(sp,2)},:,:),1),2))'];
   auxT = reshape(auxT,[size(auxT,1) 1 size(auxT,2)])
    XY = [[ones(31,1);-1*ones(31,1);ones(31,1);-1*ones(31,1)],... % hand crossing
       [ones(62,1);-1*ones(62,1)]];
   XY(:,3) = XY(:,1).*XY(:,2);
    [B,Bt,STATS,T,residuals] = regntcfe(auxT,XY,1,'effect',[],0)
    end
end