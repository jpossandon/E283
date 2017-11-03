clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'spectra';
%%
% subject configuration and data
 s=1;
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
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata               
   
    % get relevant epochevents
    load([cfg_eeg.analysisfolder 'cleaning/' cfg_eeg.sujid '/' cfg_eeg.filename cfg_eeg.clean_name],'bad');
    [trl,events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==96'},[0 8000]);
    epochevents             = [];
    epochevents.latency     = [events.time;events.time+4000];                       % fixation start, here the important thing is the ini pos
    epochevents.latency     = epochevents.latency(:)';
    epochevents.type        = repmat({'seg'},1,length(epochevents.latency));
    blockst                 = [1:length(events.time)];
    epochevents.latency     = [epochevents.latency,events.time(blockst)-4000];
    epochevents.type        = [epochevents.type,repmat({'bstart'},1,length(blockst))];
 
    rem = [];
    for e = 1:length(epochevents.latency)
        if any((epochevents.latency(e)>bad(:,1) & epochevents.latency(e)<bad(:,2)) |...
                (epochevents.latency(e)+4000>bad(:,1) & epochevents.latency(e)+4000<bad(:,2)))
            rem = [rem,e];
        end
        if epochevents.latency(e)<1
            rem = [rem,e];
        end
    end
        epochevents.latency(rem)             = [];
    epochevents.type(rem)             = [];
    [EEG,winrej] = getDataDeconv(cfg_eeg,epochevents,200);  
     EEGepoch.seg = pop_epoch( EEG, {  'seg'  }, [0  4], 'newname', ' repochs', 'epochinfo', 'yes');
    EEGepoch.bst = pop_epoch( EEG, {  'bstart'  }, [0  4], 'newname', ' repochs', 'epochinfo', 'yes');
    epTypes = {'seg','bst'};
    for eT = 1:length(epTypes)
        for ch = 1:EEGepoch.(epTypes{eT}).nbchan
            [Pxx,F] = periodogram(squeeze(EEGepoch.(epTypes{eT}).data(ch,:,:)),[],400,EEGepoch.(epTypes{eT}).srate,'power');
            spctM(ch,:) = mean(Pxx,2);
            spctSTD(ch,:) = std(Pxx,1,2);
            [~,locs] = findpeaks(double(spctM(ch,find(F>2 & F<40))),double(F(F>2 & F<40)),'MinPeakProminence',.05,'MinPeakDistance',3,'Annotate','extents');
            chpeaks{ch} = locs' ;
        end
        Pfreqs = unique(cell2mat(chpeaks));
        countsPfreqs = histc(cell2mat(chpeaks),Pfreqs);
    
        fh = figure,
        set(gcf,'Position',[10 350, 1300 600])
         subplot(2,3,1)
        plot(F,spctM'), hold on
        plot(F,mean(spctM),'r','LineWidth',3)
        ylim([0 10])
    
        subplot(2,3,2)
        loglog(F,spctM'), hold on
        loglog(F,mean(spctM),'r','LineWidth',3)
        ylim([.01 20])
        xlim([0 100])
    
        subplot(2,3,3)
        bar(Pfreqs,countsPfreqs)
        ylim([0 76])
        set(gca,'XTick',Pfreqs)
    
        load(cfg_eeg.chanlocs)
        load('cmapjp','cmap') 
        sbd = [3 8.5;7.5 16;15.5 30];
        sbdNames    = {'theta','alpha','beta'};
        for sb = 1:size(sbd,1)
            ix = find(Pfreqs>sbd(sb,1) & Pfreqs<sbd(sb,2) );

            if ~isempty(ix)
                [a,b] = max(countsPfreqs(ix));
                frtoplot = Pfreqs(ix(b));
                fpstr.(sbdNames{sb}) = frtoplot;
                [Y,I] = sort(spctM(:,find(F == frtoplot)),1,'descend');
                fpstr.([sbdNames{sb} '_chans']) = I(1:4);
                subplot(2,3,3+sb)    
                topoplot(spctM(:,F==frtoplot),chanlocs,'colormap',cmap,'headrad','rim','electrodes','on');
                colorbar
                title(sprintf('%2.1f Hz | %d peaks',frtoplot,a))
             else
                subplot(2,3,3+sb)    
                topoplot(mean(spctM(:,find(F>sbd(sb,1) & F<sbd(sb,2))),2),chanlocs,'colormap',cmap,'headrad','rim','electrodes','on');
                colorbar
                title(sprintf('%2.1f-%2.1f Hz | no peaks',sbd(sb,1),sbd(sb,2)))
                fpstr.(sbdNames{sb}) = [];
                fpstr.([sbdNames{sb} '_chans']) = [];
                
            end
        end
        mkdir(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'figures_subjects'))
        doimage(fh,fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'figures_subjects'),'png',...
                [datestr(now,'ddmmyy') cfg_eeg.sujid '_' epTypes{eT}],1)    
        spectra.(epTypes{eT})(s).id            = tk;
        spectra.(epTypes{eT})(s).spctM         = spctM;
        spectra.(epTypes{eT})(s).spctSTD       = spctSTD;
        spectra.(epTypes{eT})(s).Pfreqs        = Pfreqs;
        spectra.(epTypes{eT})(s).countsPfreqs  = countsPfreqs;
        spectra.(epTypes{eT})(s).F             = F;
        spectra.(epTypes{eT})(s).fpstr         = fpstr;
    end
    s = s+1
    save(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'allspectra'),'spectra')
end

%%

%%
%v
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'spectra';
if ismac    
    cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname',p.analysisname); % this is just to being able to do analysis at work and with my laptop
else
    cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname',p.analysisname);
end
load(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'allspectra'),'spectra')

 load(cfg_eeg.chanlocs)
 load('cmapjp','cmap') 
 sbd = [3 8;7.5 16;15.5 30];
%% 
sbd         = [3 8;7.5 16;15.5 30];
sbdNames    = {'theta','alpha','beta'};
  epTypes = {'seg','bst'};
%%   
fh1 = figure; 
set(gcf,'Position',[300 300, 1000 500])

for eT = 1:length(epTypes)
    F = spectra.(epTypes{eT})(1).F;

    for s = 1:length(spectra.(epTypes{eT}))
        allChS(s,:) = mean( spectra.(epTypes{eT})(s).spctM );
    end
  

     figure(fh1)
     subplot(2,2,2*eT-1),hold on
     plot(F,allChS,'Color',[.7 .7 .7])
     plot(F,mean(allChS),'k.-','LineWidth',2,'MarkerSize',6)
     axis([0 40 0 10])
     title('All Channels mean Spectra per Subj')
     subplot(2,2,2*eT, 'XScale', 'log', 'YScale', 'log'),hold on
     loglog(F,allChS,'Color',[.7 .7 .7])
     loglog(F,mean(allChS),'k.-','LineWidth',2,'MarkerSize',6)
     title('All Channels mean loglog Spectra per Subj')
     axis([0 40 0.01 10])
     mkdir(fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'figures'))
     
     auxsb = [spectra.(epTypes{eT}).fpstr]
     fh2 = figure;
     for sb = 1:size(sbd,1)
        auxpeak   = [auxsb.(sbdNames{sb})];
        auxchpeak =[auxsb.([sbdNames{sb} '_chans'])];
        auxfreqs.(sbdNames{sb})(eT,:)    = {auxsb.(sbdNames{sb})};
  
        set(gcf,'Position',[300 300, 1000 500])
        subplot(3,3,3*sb-2),hold on
        bar(unique(auxpeak),histc(auxpeak,unique(auxpeak)))
        title(sprintf('%s %s Subjects Peaks',epTypes{eT},sbdNames{sb}))
        xlim([sbd(sb,1)-1 sbd(sb,2)+1])
        
         subplot(3,3,3*sb-1),hold on
         for s = 1:length(spectra.(epTypes{eT}))
            plot(F,mean(spectra.(epTypes{eT})(s).spctM(spectra.(epTypes{eT})(s).fpstr.([sbdNames{sb} '_chans']),:)))
         end
         xlim([sbd(sb,1)-1 sbd(sb,2)+1])
     
        subplot(3,3,3*sb)
        topoplot(accumarray(auxchpeak(:),ones(numel(auxchpeak,1)),[length(chanlocs),1]),chanlocs,'colormap',cmap,'headrad','rim','electrodes','on');
        caxis([-15 15])
        title(sprintf('%s %s 4 peak channels',epTypes{eT},sbdNames{sb}))
        colorbar
     end
     tightfig
     doimage(fh2,fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'figures'),'png',...
                    [datestr(now,'ddmmyy') '_averageSpectra' epTypes{eT}],1)    
       
end
figure(fh1)
tightfig
doimage(fh1,fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'figures'),'png',...
                    [datestr(now,'ddmmyy') '_averageSpectra' epTypes{eT}],1)    

figure
set(gcf,'Position',[300 300, 1200 400])
for sb = 1:size(sbd,1)
subplot(1,3,sb)
auxfreqs.(sbdNames{sb})(cellfun('isempty',auxfreqs.(sbdNames{sb})))={NaN};
plot(cell2mat(auxfreqs.(sbdNames{sb})(1,:)),cell2mat(auxfreqs.(sbdNames{sb})(2,:)),'.k','MarkerSize',16),hold on
axis([sbd(sb,1)-1 sbd(sb,2)+1 sbd(sb,1)-1 sbd(sb,2)+1])
line([sbd(sb,1)-1 sbd(sb,2)+1],[sbd(sb,1)-1 sbd(sb,2)+1])
title(sprintf('%s vs %s %s Peaks',epTypes{1},epTypes{2},sbdNames{sb}))
end
tightfig
doimage(gcf,fullfile(cfg_eeg.analysisfolder,cfg_eeg.analysisname,'figures'),'png',...
                    [datestr(now,'ddmmyy') '_intrafreqCor'],1)    
