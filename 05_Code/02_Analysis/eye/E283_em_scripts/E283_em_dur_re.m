s=1;
SUMMARY.lagBINS = [1.5:1:40,inf];
for ss = subjects

    auxdat    = struct_select(data,{'subject'},{['==' num2str(ss)]},2);
    lastfix    = [diff(auxdat.trial) 1];
    nextToTarg = any([auxdat.nextToTargH;auxdat.nextToTargV;auxdat.nextToTargD]);
    
    SUMMARY.refixpTr(s) =sum(auxdat.refix & ~lastfix)./length(unique(auxdat.trial));
    SUMMARY.revisitpTr(s) =sum(auxdat.revisit>0 & ~lastfix)./length(unique(auxdat.trial));
    
    SUMMARY.newpTr(s) =sum(~auxdat.refix & ~lastfix & auxdat.start>0 & auxdat.revisit==0)./length(unique(auxdat.trial));
   
    onlyrefix                           =  find(auxdat.refix & ~lastfix &auxdat.type==1);           % only refixation that are njot hte last fixation
    preonlyrefix                        =  find(auxdat.torefix & ~auxdat.refix);                                        % the fixation previous to that
                                                           
    SUMMARY.refixMdur(s)                = mean(auxdat.dur(onlyrefix));
    SUMMARY.prerefixMdur(s)             = mean(auxdat.dur(preonlyrefix));
    
    revisit                             = find(auxdat.revisit>0 & auxdat.type==1 & ~lastfix);
    
    prerevisit                          = revisit-2;
    prerevisitLag                       = auxdat.revisit(revisit);
    prerevisitLag(auxdat.type(prerevisit)~=1) = [];
    prerevisit(auxdat.type(prerevisit)~=1)    = [];
    
    
    SUMMARY.revisitfixMdur(s)           = mean(auxdat.dur(revisit));
    SUMMARY.revisitfixMdur(s)           = mean(auxdat.dur(prerevisit));
    SUMMARY.torevisitfixMdur(s)         = mean(auxdat.dur(auxdat.torevisit>0));
    auxprerevisit                       = zeros(1,length(auxdat.dur));
    auxprerevisit(prerevisit)           = 1;
    SUMMARY.fixdurnoreMdur(s)           = mean(auxdat.dur(~auxdat.torefix & ...
                                            ~auxdat.refix & ~lastfix & auxdat.type==1 & ...
                                            ~auxdat.revisit & auxdat.start>0 & ~auxprerevisit));
    %revisit by lag                                    
    SUMMARY.revisitpLag(s,:)            = histc(auxdat.revisit(revisit),SUMMARY.lagBINS);                 
    for lag = 2:5
       SUMMARY.revisitpLagMdur(s,lag-1)    = mean(auxdat.dur(auxdat.revisit == lag & ~lastfix & auxdat.type==1));
       SUMMARY.prerevisitpLagMdur(s,lag-1) = mean(auxdat.dur(prerevisit(prerevisitLag== lag)));
       SUMMARY.torevisitpLagMdur(s,lag-1) = mean(auxdat.dur(auxdat.torevisit == lag  & auxdat.type==1));
    end
       SUMMARY.revisitpLagMdur(s,lag)    = mean(auxdat.dur(auxdat.revisit > lag & ~lastfix & auxdat.type==1));
       SUMMARY.prerevisitpLagMdur(s,lag) = mean(auxdat.dur(prerevisit(prerevisitLag > lag)));
       SUMMARY.torevisitpLagMdur(s,lag) = mean(auxdat.dur(auxdat.torevisit > lag  & auxdat.type==1));
    
    % ditances refix to element                                    
    SUMMARY.preonlyrefixMdist(s)        = mean(sqrt((auxdat.posx(preonlyrefix)-posVec.scr(1,auxdat.elfix(preonlyrefix))).^2 + ...
        (auxdat.posy(preonlyrefix)-posVec.scr(2,auxdat.elfix(preonlyrefix))).^2));
    SUMMARY.onlyrefixMdist(s)           = mean(sqrt((auxdat.posx(onlyrefix)-posVec.scr(1,auxdat.elfix(onlyrefix))).^2 + ...
        (auxdat.posy(onlyrefix)-posVec.scr(2,auxdat.elfix(onlyrefix))).^2));
    SUMMARY.allfixMdist(s)              =  mean(sqrt((auxdat.posx(~isnan(auxdat.elfix))-posVec.scr(1,auxdat.elfix(~isnan(auxdat.elfix)))).^2 + ...
        (auxdat.posy(~isnan(auxdat.elfix))-posVec.scr(2,auxdat.elfix(~isnan(auxdat.elfix)))).^2));
    
    % fixation duration by order
    for eo = 1:40
        SUMMARY.fixperEOMdur(s,eo) = mean(auxdat.dur(auxdat.event_order==eo & auxdat.type==1 & ~lastfix));
    end
        SUMMARY.fixperEOMdur(s,eo+1) = mean(auxdat.dur(auxdat.event_order>eo & auxdat.type==1 & ~lastfix));
    
    % fixation duration by order pre target
    for eo = 1:20
        SUMMARY.fixperOPTdur(s,eo) = mean(auxdat.dur(auxdat.orderPreT==eo-1 & auxdat.type==1));
        if eo>1 & eo<8
            SUMMARY.fixperOPTonTdur(s,eo-1) = nanmean(auxdat.dur(auxdat.orderPreT==eo-1 & auxdat.type==1 & auxdat.onTarg));
            SUMMARY.fixperOPnextTdur(s,eo-1) = nanmean(auxdat.dur(nextToTarg & auxdat.orderPreT==eo-1 & auxdat.type==1)); 
            SUMMARY.fixperOPnotnextTdur(s,eo-1) = nanmean(auxdat.dur(~nextToTarg & auxdat.orderPreT==eo-1 & auxdat.type==1 & ~auxdat.onTarg)); 
        end
    end
    
    % saccade amplitude by order pre target
    for eo = 1:20
        SUMMARY.sacperOPTamp(s,eo) = mean(auxdat.amp(auxdat.orderPreT==eo-1 & auxdat.type==2 & auxdat.amp<50));
        
    end
    s = s+1;
end
    
[r,p] = corr(SUMMARY.rT',SUMMARY.refixpTr');
SUMMARY.corrRT_refix = [r, p];
[r,p] = corr(SUMMARY.rT',SUMMARY.revisitpTr');
SUMMARY.corrRT_revisit = [r, p];
[r,p] = corr(SUMMARY.refixpTr',SUMMARY.revisitpTr');
SUMMARY.corrrefix_revisit = [r, p];
[r,p] = corr(SUMMARY.newpTr',SUMMARY.revisitpTr');
SUMMARY.corrnew_revisit = [r, p];
[r,p] = corr(SUMMARY.refixpTr',SUMMARY.newpTr');
SUMMARY.corrrefix_new = [r, p];

%%
% figure correlation between reaction times and number of refix and
% revisits per trial
figure, hold on
plot(SUMMARY.rT,SUMMARY.refixpTr,'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[1 .4 .4])
plot(SUMMARY.rT,SUMMARY.revisitpTr,'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.4 .4 1])
plot(1.25,2.3,'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[1 .4 .4])
text(1.31,2.3,sprintf('refix, corr: %2.1f  p : %1.3f',SUMMARY.corrRT_refix(1),SUMMARY.corrRT_refix(2)),'FontSize',10)
plot(1.25,2.1,'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.4 .4 1])
text(1.31,2.1,sprintf('revisit, corr: %2.1f  p : %1.3f',SUMMARY.corrRT_revisit(1),SUMMARY.corrRT_revisit(2)),'FontSize',10)

axis([1 4 0 2.5])
set(gca,'FontSize',10)
xlabel('Reaction Time (s)','FontSize',12)
ylabel('# refix/revisit per trial','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['corrRT_refix_revisit_' namegr],[],1)

%%
% figure correlation between refix nd revisit

figure, hold on
line([0 3],[0 3],'Color',[0.5 0.5 0.5],'LineStyle',':')
plot(SUMMARY.revisitpTr,SUMMARY.refixpTr,'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.4 .4 .4])
text(0.25,2.25,sprintf('corr: %2.1f  p : %1.3f',SUMMARY.corrrefix_revisit(1),SUMMARY.corrrefix_revisit(2)),'FontSize',10)

axis([0 2.5 0 2.5])
set(gca,'FontSize',10)
xlabel('Reaction Time (s)','FontSize',12)
ylabel('Reaction Time (s)','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
doimage(gcf,fullfile(patheye,'figures'),...
              'tiff',['corrrefix_revisit_' namegr],[],1)

%%
% revisit per lag
figure, hold on

plot(SUMMARY.lagBINS+unique(diff(SUMMARY.lagBINS(1:end-1)))/2,SUMMARY.revisitpLag./repmat(SUMMARY.nT',1,length(SUMMARY.lagBINS)),...
    'Color',[.7 .7 .7])
plot(SUMMARY.lagBINS+unique(diff(SUMMARY.lagBINS(1:end-1)))/2,mean(SUMMARY.revisitpLag./repmat(SUMMARY.nT',1,length(SUMMARY.lagBINS))),'-',...
    'Color',[0 0 0],'LineWidth',3,'MarkerSize',12)
plot(SUMMARY.lagBINS+unique(diff(SUMMARY.lagBINS(1:end-1)))/2,mean(SUMMARY.revisitpLag./repmat(SUMMARY.nT',1,length(SUMMARY.lagBINS))),'.',...
    'Color',[1 0 0],'LineWidth',3,'MarkerSize',8)
set(gca,'XTick',[[2:5],10:10:40],...
     'FontSize',10)
axis([0 40 0 .25])
xlabel('Revisit lag','FontSize',12)
ylabel('revisit pre trial','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
% doimage(gcf,fullfile(patheye,'figures'),...
%              'tiff',['revisitlag_' namegr],[],1)
        

%%
% figure fixation duration different type of fixation
labelsST = {'refixMdur','prerefixMdur','revisitfixMdur','torevisitfixMdur','fixdurnoreMdur'};
labels   = {'refix','prerefix','revisit','torevisit','fixdurOther'};
figure
hold on
jitter = (rand(1,Ns)-.5)/10;
for ll = 1:length(labelsST)
    plot(ll+jitter,SUMMARY.(labelsST{ll}),'ok','MarkerSize',4,'LineWidth',1,'MarkerFaceColor',[.7 .7 .7])
     errorbar(ll,mean(SUMMARY.(labelsST{ll})),std(SUMMARY.(labelsST{ll}))./sqrt(Ns),'LineWidth',2,'Color',[1 0 0])
    plot(ll,mean(SUMMARY.(labelsST{ll})),'sk','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[1 0 0])
   
end
axis([.5 5.5 0 400])
set(gca,'FontSize',10,'XTick',1:length(labels),'XTickLabel',labels)
ylabel('Fixation duration (ms)','FontSize',12)
tightfig
doimage(gcf,fullfile(patheye,'figures'),...
            'tiff',['fixdurs_' namegr],[],1)

%%
% figure fixation duration at different event orders
figure, hold on
plot(1:size(SUMMARY.fixperEOMdur,2)-1,SUMMARY.fixperEOMdur(:,2:end),...
    'Color',[.7 .7 .7])
plot(1:size(SUMMARY.fixperEOMdur,2)-1,mean(SUMMARY.fixperEOMdur(:,2:end)),'-',...
    'Color',[0 0 0],'LineWidth',3,'MarkerSize',12)
errorbar(1:size(SUMMARY.fixperEOMdur,2)-1,mean(SUMMARY.fixperEOMdur(:,2:end)),std(SUMMARY.fixperEOMdur(:,2:end),1,1),...
    'Color',[0 0 0],'LineWidth',2)
plot(1:size(SUMMARY.fixperEOMdur,2)-1,mean(SUMMARY.fixperEOMdur(:,2:end)),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[1 .4 .4])


set(gca,'FontSize',10)
axis([0 41 100 300])
xlabel('Event order','FontSize',12)
ylabel('Fixation duration (ms)','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
doimage(gcf,fullfile(patheye,'figures'),...
              'tiff',['fixdurperEO_' namegr],[],1)

%%
% figure fixation duration at different event orders pre T
fh = figure, hold on
plot(1:10,[SUMMARY.fixperOPTdur(:,2:10) mean(SUMMARY.fixperOPTdur(:,11:end),2)],...
    'Color',[.7 .7 .7],'LineWidth',.5)
plot(1:10,[mean(SUMMARY.fixperOPTdur(:,2:10)) mean(mean(SUMMARY.fixperOPTdur(:,11:end),2))] ,'-',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',12)
errorbar(1:10,[mean(SUMMARY.fixperOPTdur(:,2:10)) mean(mean(SUMMARY.fixperOPTdur(:,11:end),2))],...
    [std(SUMMARY.fixperOPTdur(:,2:10)) std(mean(SUMMARY.fixperOPTdur(:,11:end),2))],...
    'Color',[0 0 0],'LineWidth',.5)
plot(1:10,[mean(SUMMARY.fixperOPTdur(:,2:10)) mean(mean(SUMMARY.fixperOPTdur(:,11:end),2))],'s',...
    'Color',[0 0 0],'LineWidth',.5,'MarkerSize',4,'MarkerFaceColor',[1 .4 .4])


set(gca,'FontSize',6,'YTick',100:50:300,'XTick',[1:2:9, 10],'XTickLabels',{'1', '3', '5', '7', '9', '>9'})
axis([0 10.5 100 300])
xlabel('Event order pre Target','FontSize',8)
ylabel('Fixation duration (ms)','FontSize',8)
% title(sprintf('N = %d',Ns),'FontSize',12)
% tightfig
figsize     = [4.6 4.6*fh.Position(4)/fh.Position(3)];
doimage(fh,fullfile(patheye,'figures'),...
                'pdf',['fixdurperOPT_' namegr],'painters',figsize,1)
           
%%
% saccade amplitude at different event orders pre T
fh = figure, hold on
plot(1:10,[SUMMARY.sacperOPTamp(:,2:10) mean(SUMMARY.sacperOPTamp(:,11:end),2)],...
    'Color',[.7 .7 .7],'LineWidth',.5)
plot(1:10,[mean(SUMMARY.sacperOPTamp(:,2:10)) mean(mean(SUMMARY.sacperOPTamp(:,11:end),2))] ,'-',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',12)
errorbar(1:10,[mean(SUMMARY.sacperOPTamp(:,2:10)) mean(mean(SUMMARY.sacperOPTamp(:,11:end),2))],...
    [std(SUMMARY.sacperOPTamp(:,2:10)) std(mean(SUMMARY.sacperOPTamp(:,11:end),2))],...
    'Color',[0 0 0],'LineWidth',.5)
plot(1:10,[mean(SUMMARY.sacperOPTamp(:,2:10)) mean(mean(SUMMARY.sacperOPTamp(:,11:end),2))],'s',...
    'Color',[0 0 0],'LineWidth',.5,'MarkerSize',4,'MarkerFaceColor',[1 .4 .4])


set(gca,'FontSize',6,'YTick',0:2:8,'XTick',[1:2:9, 10],'XTickLabels',{'1', '3', '5', '7', '9', '>9'})
axis([0 10.5 2 8])
xlabel('Event order pre Target','FontSize',8)
ylabel('Saccade Amplitude (\circ)','FontSize',8)
% title(sprintf('N = %d',Ns),'FontSize',12)
% tightfig
figsize     = [4.6 4.6*fh.Position(4)/fh.Position(3)];
doimage(fh,fullfile(patheye,'figures'),...
                'pdf',['sacampperOPT_' namegr],'painters',figsize,1)
                       
%%
%% 
% figure fixation duration on target (but not last), next to target and not
% next to target
figure, hold on
plot([1:size(SUMMARY.fixperOPTonTdur,2)]-.15,nanmean(SUMMARY.fixperOPTonTdur(:,1:end)),'-',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',12)
errorbar([1:size(SUMMARY.fixperOPTonTdur,2)]-.15,nanmean(SUMMARY.fixperOPTonTdur(:,1:end)),nanstd(SUMMARY.fixperOPTonTdur(:,1:end),1,1),...
    'Color',[0 0 0],'LineWidth',1)
hp(1) = plot([1:size(SUMMARY.fixperOPTonTdur,2)]-.15,nanmean(SUMMARY.fixperOPTonTdur(:,1:end)),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[1 .4 .4])

plot([1:size(SUMMARY.fixperOPnextTdur,2)],nanmean(SUMMARY.fixperOPnextTdur(:,1:end)),'-',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',12)
errorbar([1:size(SUMMARY.fixperOPnextTdur,2)],nanmean(SUMMARY.fixperOPnextTdur(:,1:end)),nanstd(SUMMARY.fixperOPnextTdur(:,1:end),1,1),...
    'Color',[0 0 0],'LineWidth',1)
hp(2) =plot([1:size(SUMMARY.fixperOPnextTdur,2)],nanmean(SUMMARY.fixperOPnextTdur(:,1:end)),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[.4 1 .4])

plot([1:size(SUMMARY.fixperOPnotnextTdur,2)]+.15,nanmean(SUMMARY.fixperOPnotnextTdur(:,1:end)),'-',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',12)
errorbar([1:size(SUMMARY.fixperOPnotnextTdur,2)]+.15,nanmean(SUMMARY.fixperOPnotnextTdur(:,1:end)),nanstd(SUMMARY.fixperOPnotnextTdur(:,1:end),1,1),...
    'Color',[0 0 0],'LineWidth',1)
hp(3) =plot([1:size(SUMMARY.fixperOPnotnextTdur,2)]+.15,nanmean(SUMMARY.fixperOPnotnextTdur(:,1:end)),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[.4 .4 1])

legend(hp,{'onTarget','nextTarget','notNext'})
set(gca,'FontSize',10)
axis([0 7 000 400])
xlabel('Event order','FontSize',12)
ylabel('Fixation duration (ms)','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['fixdurperOPTpertype_' namegr],[],1)

%%
% figure fixation duration at different revisit lags and the duration of the fixation before 
figure
hold on
jitter = (rand(1,Ns)-.5)/10;
for ll = 1:5
    plot(ll+jitter,SUMMARY.revisitpLagMdur(:,ll),'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.7 .7 .7])
    errorbar(ll,mean(SUMMARY.revisitpLagMdur(:,ll)),std(SUMMARY.revisitpLagMdur(:,ll))./sqrt(Ns),'LineWidth',2,'Color',[0 0 0])
    plot(ll,mean(SUMMARY.revisitpLagMdur(:,ll)),'sk','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 0 0])
  
    plot(ll+6+jitter,SUMMARY.prerevisitpLagMdur(:,ll),'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.7 .7 .7])
    errorbar(ll+6,mean(SUMMARY.prerevisitpLagMdur(:,ll)),std(SUMMARY.prerevisitpLagMdur(:,ll))./sqrt(Ns),'LineWidth',2,'Color',[0 0 0])
    plot(ll+6,mean(SUMMARY.prerevisitpLagMdur(:,ll)),'sk','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 0 0])
  
    plot(ll+12+jitter,SUMMARY.torevisitpLagMdur(:,ll),'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.7 .7 .7])
    errorbar(ll+12,mean(SUMMARY.torevisitpLagMdur(:,ll)),std(SUMMARY.torevisitpLagMdur(:,ll))./sqrt(Ns),'LineWidth',2,'Color',[0 0 0])
    plot(ll+12,mean(SUMMARY.torevisitpLagMdur(:,ll)),'sk','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 0 0])
  
end

plot(ll+14+jitter,SUMMARY.fixdurnoreMdur,'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.7 .7 .7])
    errorbar(ll+14,mean(SUMMARY.fixdurnoreMdur),std(SUMMARY.fixdurnoreMdur)./sqrt(Ns),'LineWidth',2,'Color',[0 0 0])
    plot(ll+14,mean(SUMMARY.fixdurnoreMdur),'sk','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 0 0])
  hline(mean(SUMMARY.fixdurnoreMdur))

axis([0 20 0 300])
set(gca,'FontSize',10,'XTick',[1:5,7:11,13:17,19],'XTickLabel',{'2','3','4','5','>6','2','3','4','5','>6','2','3','4','5','>6','nore'})
ylabel('Fixation duration (ms)','FontSize',12)
xlabel('Duration revisit at lag          Duration pre-revisit at lag         Duration to-revisit at lag')
tightfig
  doimage(gcf,fullfile(patheye,'figures'),...
               'tiff',['fixdursrevisit_' namegr],[],1)