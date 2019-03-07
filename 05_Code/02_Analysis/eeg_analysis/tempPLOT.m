%%

dattocomp = {'U_unI_ipsi','C_unI_contra';'U_unI_contra','C_unI_ipsi';'U_I','C_I';'U_unI','C_unI'}%;'U_Ici','C_Ici';
            % 'LU_I','RU_I';'LU_Ici','RU_Ici';'LC_I','RC_I';'LC_Ici','RC_Ici';
           %  'LU_I','LC_I';'LU_Ici','LC_Ici';'RU_I','RC_I';'RU_Ici','RC_Ici';
             %'LU_unI','RU_unI';'LU_unIci','RU_unIci';'LC_unI','RC_unI';'LC_unIci','RC_unIci';
             %'LU_unI','LC_unI';'LU_unIci','LC_unIci';'RU_unI','RC_unI';'RU_unIci','RC_unIci';
        %     'U_Ici','C_Ici';'U_unIci','C_unIci'};
bands               = {'alpha','beta'};

chans = {[50],[18],[58,65,66],[26,33,34]}
chansLabels = {'C3','C4','P5 P7 P07','P6 P8 P08'}
cluslab = {'pos','neg'};

tiempo =   GAbsl.(dattocomp{1,1}).time;       
findx = find(GAbsl.LU_I.freq>16 &GAbsl.LU_I.freq<26);
for bb = 1:2
    if strcmp(bands{bb},'alpha')
        findx = find(GAbsl.LU_I.freq>8 &GAbsl.LU_I.freq<16);
    elseif strcmp(bands{bb},'beta')
         findx = find(GAbsl.LU_I.freq>15 &GAbsl.LU_I.freq<27);
    end
    for pp = 1:size(dattocomp,1)
        fh=figure;
       fh.Position = [1 1 1440 805];
sp =0;
        for ch = 1:length(chans)
            aux1 = squeeze(mean(mean(GAbsl.(dattocomp{pp,1}).powspctrm(:,chans{ch},findx,:),2),3));
            aux2 = squeeze(mean(mean(GAbsl.(dattocomp{pp,2}).powspctrm(:,chans{ch},findx,:),2),3));
            N = size(aux1,1);
            subplot(2,2,1+sp), hold on
            %plot(GAbsl.(dattocomp{pp,1}).time,aux1,'Color',[1 .7 .7])
            %plot(GAbsl.(dattocomp{pp,1}).time,aux2,'Color',[.7 .7 1])
            jbfill(tiempo,mean(aux1)+std(aux1)/sqrt(N),mean(aux1)-std(aux1)/sqrt(N),[1 .3 .3],[1 .15 .15],1,.7);hold on
            plot(tiempo,mean(aux1),'r','LineWidth',2),hold on
            jbfill(tiempo,mean(aux2)+std(aux2)/sqrt(N),mean(aux2)-std(aux2)/sqrt(N),[.3 .3 1],[.15 .15 1],1,.7);hold on
            plot(tiempo,mean(aux2),'b','LineWidth',2)
            axis([-.6 .8 -6 2])
            title(sprintf('%s %s Ch %s', dattocomp{pp,1},dattocomp{pp,2},chansLabels{ch}), 'Interpreter', 'none')
           
             %plot(tiempo,aux1-aux2,'Color',[.7 .7 .7])
            jbfill(tiempo,[mean(aux1-aux2)+std(aux1-aux2)/sqrt(N)],[mean(aux1-aux2)-std(aux1-aux2)/sqrt(N)],[.3 .3 0],[.15 .15 0],1,.7),hold on
            plot(tiempo,median(aux1-aux2),'k--','LineWidth',2)
            plot(tiempo,mean(aux1-aux2),'k:','LineWidth',2)
            [h,p] = ttest(aux1,aux2);
           % plot(tiempo,p)
            elec.channeighbstructmat = 0;
            [result] = regmodel2ndstat(reshape(aux1',[1 1 size(aux1')])-reshape(aux2',[1 1 size(aux2')]),tiempo,elec,2000,'signpermT','cluster'); %need check this with clusters and singnperm to
            plot(tiempo,result.pvalnc,'k')
            for pn = 1:2
                    if ~isempty(result.clusters.([cluslab{pn} 'clusters']))
                        for cc = 1:length(result.clusters.([cluslab{pn} 'clusters']))
                            if result.clusters.([cluslab{pn} 'clusters'])(cc).prob < 0.1
                                auxclus = find(result.clusters.([cluslab{pn} 'clusterslabelmat'])==cc);
                            auxdiff = mean(aux1(:,auxclus)-aux2(:,auxclus));
                                if result.clusters.([cluslab{pn} 'clusters'])(cc).prob < 0.05
                                    plot(tiempo(auxclus),auxdiff,'r','LineWidth',2)
                                end
                                if strcmp(cluslab{pn},'neg')
                                    [~,di] = min(mean(aux1(:,auxclus)-aux2(:,auxclus)));
                                elseif strcmp(cluslab{pn},'pos')
                                    [~,di] = max(mean(aux1(:,auxclus)-aux2(:,auxclus)));
                              
                                end
                                text(tiempo(auxclus(di)),auxdiff(di)+.2*sign(auxdiff(di)),sprintf('p = %1.4f',result.clusters.([cluslab{pn} 'clusters'])(cc).prob))
                            end
                        end
                    end
                end   
            %axis([-.8 1 -1 1])
            %title(sprintf('%s-%s Ch %s', dattocomp{pp,1},dattocomp{pp,2},chansLabels{ch}), 'Interpreter', 'none')
            hline(0,'k:')
            vline(0,'k:')
            sp=sp+1;
        end
         figsize     = [17.6 17.6*fh.Position(4)/fh.Position(3)];
            doimage(fh,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_' bands{bb} '_' dattocomp{pp,1} '_' dattocomp{pp,2}],figsize,1)

    end
end

%%
% the channels
for ch = 1:length(chans)
    figure
    tp = topo_markCh(cfg_eeg,chans {ch});
    tightfig
    doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',['Channs_' strjoin('_',chansLabels(ch))],[2 2],1)
end

%%
GAfields = { 'LU_I','LC_I','RU_I','RC_I','LU_unI','LC_unI','RU_unI','RC_unI'};

findx = find(GAbsl.LU_I.freq>9 &GAbsl.LU_I.freq<15);
freqs = GAbsl.LU_I.freq;
auxtBSL  = find(GAbsl.LU_I.time>-.65 &GAbsl.LU_I.time<-.15);
auxt  = find(GAbsl.LU_I.time>.15);
chns  = [18,50];
ch=18
ff=1
clear bslaft befbsl
for ch = 1:length(chns)
 for ff = 1:length(GAfields)
   % aux = squeeze(nanmean(GAbsl.(GAfields{ff}).powspctrm(:,ch,:,auxt),4))
   % figure,plot(freqs,aux','.-')
   % hold on,plot(freqs,mean(aux),'.-k','LineWidth',3)
   % aux = squeeze(nanmean(GAnobsl.(GAfields{ff}).powspctrm(:,ch,:,auxtBSL),4))
    %   figure,plot(freqs,aux','.-')
   % hold on,plot(freqs,mean(aux),'.-k','LineWidth',3)
    
    bslaft(:,ff,ch) = squeeze(nanmean(nanmean(GAbsl.(GAfields{ff}).powspctrm(:,chns(ch),findx,auxt),4),3));
    nobslbef(:,ff,ch) = squeeze(nanmean(nanmean(GAnobsl.(GAfields{ff}).powspctrm(:,chns(ch),findx,auxtBSL),4),3));
    nobslaf(:,ff,ch) = squeeze(nanmean(nanmean(GAnobsl.(GAfields{ff}).powspctrm(:,chns(ch),findx,auxt),4),3));
   % figure,plot(befbsl,bslaft,'.')
   % xlabel('baseline power')
   % ylabel('stim red db')
 end
end

%%
figure,
ch = [2,1,2,1];
chLb = {'C4','C3'}
comp = [1 2;1 2;5 6;5 6];
for sp = 1:4
    subplot(2,2,sp)
    plot(bslaft(:,comp(sp,1),ch(sp)),bslaft(:,comp(sp,2),ch(sp)),'.')
    line([-8 4],[-8 4],'Color',[0 0 0],'LineStyle',':')
    xlabel(GAfields{comp(sp,1)},'interpreter','none'),ylabel(GAfields{comp(sp,2)},'interpreter','none')
%axis([-8 4 -8 4])
    title(chLb{ch(sp)})
end

figure,
ch = [2,1,2,1];
chLb = {'C4','C3'}
comp = [3 4;3 4;7 8;7 8];
for sp = 1:4
    subplot(2,2,sp)
    plot(bslaft(:,comp(sp,1),ch(sp)),bslaft(:,comp(sp,2),ch(sp)),'.')
    line([-8 4],[-8 4],'Color',[0 0 0],'LineStyle',':')
    xlabel(GAfields{comp(sp,1)},'interpreter','none'),ylabel(GAfields{comp(sp,2)},'interpreter','none')
    axis([-10 5 -10 5])
    title(chLb{ch(sp)})
end

%%
% join ipsi/contra data
setAbsoluteFigureSize
uc   = cat(3,(bslaft(:,[2,4,7,8],2)+bslaft(:,[1,2,5,6],1))/2,(bslaft(:,[1,2,5,6],2)+bslaft(:,[2,4,7,8],1))/2);
ucdiff = [uc(:,1,:)-uc(:,2,:) uc(:,3,:)-uc(:,4,:)]; % row - subjects, col - 1 inf 2 uninf, 3rd - 1 ipsi, 2 contra

splabels = {'Informative','Uninformative'};
fh = figure;
fh.Position = [9 479 465 212];
lims = [-3.5 3.5];
for sp = 1:2
    
M = squeeze(mean(ucdiff(:,sp,:)));
Med = squeeze(median(ucdiff(:,sp,:)));
SE = squeeze(std(ucdiff(:,sp,:)))./sqrt(size(ucdiff,1));
subplot(1,2,sp),hold on
line([lims],[lims],'Color',[.8 .8 .8],'LineStyle',':')
line([lims],[0 0],'Color',[.8 .8 .8],'LineStyle',':')
line([0 0],[lims],'Color',[.8 .8 .8],'LineStyle',':')
plot(ucdiff(:,sp,1),ucdiff(:,sp,2),'ok','MarkerSize',6)
line([M(1)-SE(1) M(1)+SE(1)],[M(2) M(2)],'LineWidth',1,'Color',[1 0 0])
line([M(1) M(1)],[M(2)-SE(2) M(2)+SE(2)],'LineWidth',1,'Color',[1 0 0])
set(gca,'XTick',-3:1.5:3,'YTick',-3:1.5:3)
axis([lims lims])
axis square
xlabel('U - C ipsi (dB)')
ylabel('U - C contra (dB)')
box off
title(splabels{sp})
end
 figsize     = [8.8 8.8*fh.Position(4)/fh.Position(3)];
%         doimage(gcf,[cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/figures/GA/'],'pdf',[datestr(now,'ddmmyy') '_UvsC_alpha'],figsize,1)

%%
% model with constant and interaction, with the means of Subjects
'U_I','C_I';'U_unI','C_unI'
ch=1
chans = {[50],[18],[58,65,66],[26,33,34]}
chansLabels = {'C3','C4','P5 P7 P07','P6 P8 P08'}
baseXY = [1 1;-1 1;1 -1;-1 -1];
N = size(GAbsl.LU_I.powspctrm,1);
baseXYS = [eye(N-1); -1*ones(1,N-1)];

for bb = 1
    if strcmp(bands{bb},'alpha')
        findx = find(GAbsl.LU_I.freq>8 &GAbsl.LU_I.freq<16);
    elseif strcmp(bands{bb},'beta')
         findx = find(GAbsl.LU_I.freq>15 &GAbsl.LU_I.freq<27);
    end
    for ch = 1:length(chans)
        for pp = 1:501
            XY = []; Y = [];
            for suj= 1:size(GAbsl.LU_I.powspctrm,1)
                sujdata = [squeeze(mean(mean(GAbsl.U_I.powspctrm(suj,chans{ch},findx,:),2),3))';
                            squeeze(mean(mean(GAbsl.C_I.powspctrm(suj,chans{ch},findx,:),2),3))';
                            squeeze(mean(mean(GAbsl.U_unI.powspctrm(suj,chans{ch},findx,:),2),3))';
                            squeeze(mean(mean(GAbsl.C_unI.powspctrm(suj,chans{ch},findx,:),2),3))'];
                if pp == 1
                    order = 1:4;
                else
                    order = randsample(4,4);
                end
                Y = [Y;sujdata];
                auxXY = baseXY(order,:);
                auxXY(:,3) = auxXY(:,1).*auxXY(:,2);
                auxXY = [auxXY,repmat(baseXYS(suj,:),size(auxXY,1),1)];
                auxXY = [auxXY,repmat(auxXY(:,1),1,size(baseXYS,2)).*repmat(baseXYS(suj,:),size(auxXY,1),1),...
                    repmat(auxXY(:,2),1,size(baseXYS,2)).*repmat(baseXYS(suj,:),size(auxXY,1),1)];
                XY = [XY;auxXY];
            end
             [B,Bt,STATS,T,residuals] = regntcfe(Y,XY,1,'effect',[],0)    
        end
    end
end
