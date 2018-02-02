%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUAL SEARCH E283 ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% GET DATA
clear
pathexp = '/Users/jossando/trabajo/E283/';
group   = 2;    % 1-  subject 6 to 27, 2, 28 to 59, 3 all
load(fullfile(pathexp,'analysis','eyedata','alleyedataFULL.mat'),'data')    % info from subjects EDF files, augments with augmentinfALLdata.m 
load(fullfile(pathexp,'analysis','eyedata','allRTFULL.mat'))                    % info from subjects matlab files 
load(fullfile(pathexp,'analysis','eyedata','tgtPos.mat'))                   % position of targets on the screen

if group==1
    data    = struct_elim(data,data.subject>27,2,1);
    result  = struct_elim(result,result.subject>27,2,1);
    namegr  = 'randStim';
elseif group==2
    data    = struct_elim(data,data.subject<28,2,1);
    result  = struct_elim(result,result.subject<28,2,1);
    namegr  = 'earlyStim';
else
    namegr  = 'all';
end
subjects    = unique(data.subject);
Ns          = length(subjects);
subjectswin = unique(result.subject);
if ~all(subjects==subjectswin)
    error('subjects in data and result files do not match')
end
clear subjectswin



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT AND PERFORMANCE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SUMMARY

s=1;
SUMMARY.rTBINS   = 0:.5:12;
for ss = subjects
    SUMMARY.subjIndx(s) = ss;
    auxres  = struct_select(result,{'subject'},{['==' num2str(ss)]},2);
     if length(auxres.perf) ~= 592
        display(sprintf('\nSubject %d has only %d trials',ss,length(auxres.perf))) %  subject 20 and 22 the EEG recording was started late
     end
    SUMMARY.nT(s)         = length(auxres.perf);
    SUMMARY.perf(s)       = sum(auxres.perf)./length(auxres.perf);
    SUMMARY.rT(s)         = mean(auxres.rT(auxres.perf==1));
    SUMMARY.rTpbin(s,:)   = histc(auxres.rT(auxres.perf==1),SUMMARY.rTBINS);
    
    % rt per position
    for pos = 1:48
        SUMMARY.rtpPos(s,pos) = mean(auxres.rT(auxres.perf==1 & auxres.tpos==pos));
    end
    %nfixs
   
    if group==2
        % value: 1 - LUI, 2 - RUI, 5 - LCI, 6 - RCI, 9 - LUunI, 10 - RUuni
        %  13 - LCunI, 14 - RCuni
        % cue: 0/1 - informative, 1/2 - uninformative ,2 test trials
        SUMMARY.rTpCond(s,:) = [mean(auxres.rT(auxres.perf==1 & auxres.value==1)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==2)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==5)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==6)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==9)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==10)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==13)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==14))];
    end
     s = s+1;
end
SUMMARY.condLabels  = {'LUI','RUI','LCI','RCI','LUunI','RUunI','LCunI','RCunI'};
SUMMARY.condValues = [1,2,5,6,9,10,13,14];
SUMMARY.infoLabels  = {'Info','Info','Info','Info','unInfo','unInfo','unInfo','unInfo'};
SUMMARY.meanRT      = mean(SUMMARY.rT);
SUMMARY.medianRT    = median(SUMMARY.rT);
SUMMARY.sdRT        = std(SUMMARY.rT);
SUMMARY.meanRTse    = SUMMARY.sdRT./sqrt(Ns);

SUMMARY.meanRTpCond     = mean(SUMMARY.rTpCond);
SUMMARY.medianRTpCond   = median(SUMMARY.rTpCond);
SUMMARY.sdRTpCond       = std(SUMMARY.rTpCond);
SUMMARY.seRTCond        = SUMMARY.sdRTpCond./sqrt(Ns);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple check that eyedata correspond with target positions and border
% between elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, hold on
set(gcf,'Position',[0 0 1280 700])
axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
plot(data.posx(data.type==1 & data.start>0),data.posy(data.type==1 & data.start>0),'.k','MarkerSize',1)
plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',16)
axis ij
vline(posVec.hLims)
hline(posVec.vLims,'r:')
tightfig
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
            'tiff',['allFixs_' namegr],1)
        
% this is for the different types of cue
cuevalues = [1 2 5 6;9 10 13 14];
infLabel  = {'Info','unInfo'};
for ii = 1:2

    figure, hold on
    set(gcf,'Position',[0 0 1280 700])
    for sp = [1,3,2,4]
        subplot(2,2,sp), hold on
        axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
        plot(data.posx(data.type==1 & data.value==cuevalues(ii,sp) & data.start>0),...
            data.posy(data.type==1 & data.value==cuevalues(ii,sp) & data.start>0),'.k','MarkerSize',1)
        plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',16)
        axis ij
        vline(posVec.hLims)
        hline(posVec.vLims,'r:')
        title(SUMMARY.condLabels{sp+(ii-1)*4})
        
    end
    tightfig
    doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
                    'tiff',['allFixs_' namegr '_' infLabel{ii}],1)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT AND PERFORMANCE FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, hold on
plot(SUMMARY.rTBINS+unique(diff(SUMMARY.rTBINS))/2,SUMMARY.rTpbin./repmat(sum(SUMMARY.rTpbin,2),1,size(SUMMARY.rTpbin,2)),...
    'Color',[.7 .7 .7])
plot(SUMMARY.rTBINS+unique(diff(SUMMARY.rTBINS))/2,sum(SUMMARY.rTpbin./repmat(sum(SUMMARY.rTpbin,2),1,size(SUMMARY.rTpbin,2)))./Ns,...
    'Color',[0 0 0],'LineWidth',3)
plot(SUMMARY.rT,0.225+(rand(1,Ns)-.5)/150,'ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[1 .9 .9])
line([SUMMARY.meanRT-SUMMARY.meanRTse SUMMARY.meanRT+SUMMARY.meanRTse],[.225 .225],'LineWidth',3,'Color',[0 0 0])
plot(SUMMARY.meanRT,0.225,'sk','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[1 0 0])
text(max(SUMMARY.rT)+range(SUMMARY.rT)*.2,0.225,sprintf('mean : %2.2f + - %2.2f',SUMMARY.meanRT,SUMMARY.meanRTse),'FontSize',10)
plot(13+(rand(1,Ns)-.5)/5,1-SUMMARY.perf,'ok','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[1 .9 .9])
set(gca,'XTick',[SUMMARY.rTBINS(1:4:end) 13],'XTickLabel',{'0','2','4','6','8','10','12', 'misses'},...
    'FontSize',10)
axis([0 14 0 .25])
xlabel('Reaction Time (s)','FontSize',12)
ylabel('Frequency','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
 doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
             'tiff',['performance_' namegr],1)

% per condition  
figure
hold on
jitter = (rand(1,Ns)-.5)/10
for ll = 1:length(SUMMARY.condLabels)
    plot(ll+jitter,SUMMARY.rTpCond(:,ll),'ok','MarkerSize',6,'LineWidth',.5,...
        'MarkerFaceColor',[.9 .9 .9],'MarkerEdgeColor',[.5 .5 .5])
end 
errorbar(1:4,SUMMARY.meanRTpCond(1:4),SUMMARY.seRTCond(1:4),'LineWidth',2,'Color',[.1 .1 .1])
plot(1:4,SUMMARY.meanRTpCond(1:4),'sk','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[.3 .3 1])
errorbar(5:8,SUMMARY.meanRTpCond(5:8),SUMMARY.seRTCond(5:8),'LineWidth',2,'Color',[.1 .1 .1])
plot(5:8,SUMMARY.meanRTpCond(5:8),'sk','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[1 .3 .3])
    
axis([.5 8.5 0 6])
set(gca,'FontSize',10,'XTick',1:length(SUMMARY.condLabels),'XTickLabel',SUMMARY.condLabels)
ylabel('Reaction Time (s)','FontSize',12)

tightfig
 doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
             'tiff',['RTpercond_' namegr],1)

%%
% figure reaction time per position
cmap = cmocean('haline');
figure
set(gcf,'Position',[50 50 1100 1100])
subplot(7,9,[1:9:46]), hold on
errorbar(1:6,mean(mean(reshape(SUMMARY.rtpPos',[6 8 Ns]),2),3),std(mean(reshape(SUMMARY.rtpPos',[6 8 Ns]),2),1,3)./sqrt(Ns),...
    'Color',[0 0 0])
plot(1:6,mean(mean(reshape(SUMMARY.rtpPos',[6 8 Ns]),2),3),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[1 .4 .4])
axis([.5 6.5 2 4])
view([90 90])
set(gca,'FontSize',10)

subplot(7,9,repmat(2:9,1,6)+reshape(repmat([0 9:9:45],8,1),[1 48]))
imagesc(reshape(mean(SUMMARY.rtpPos),[6,8]))
axis off
caxis([1.5 3.5])
 colorbar
 title('Reaction Time per position (s)','FontSize',12)
 
subplot(7,9,56:63)
hold on
errorbar(1:8,mean(mean(reshape(SUMMARY.rtpPos',[6 8 Ns])),3),std(mean(reshape(SUMMARY.rtpPos',[6 8 Ns])),1,3)./sqrt(Ns),...
    'Color',[0 0 0])
plot(1:8,mean(mean(reshape(SUMMARY.rtpPos',[6 8 Ns])),3),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[1 .4 .4])
axis([.5 8.5 2 4])
 colorbar
set(gca,'FontSize',10)

colormap(cmap)

doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
             'tiff',['RTperPos_' namegr],1)        

%%
% figure reaction time per distance to the target
lgscale = 0;
ixfix   = find(data.start>0 & data.type==1 & data.tToTarg<-100 & data.tToTarg>-10000);
ixfix   = ixfix(1:5:end);
figure,hold on
cmap    = flipud(cmocean('thermal'));
for ff = ixfix
%     ff
if lgscale
    plot(data.xToTarg(ff),data.yToTarg(ff),'.','Color',cmap(round((log10(-data.tToTarg(ff)/1000)+1)/2*256),:),'MarkerSize',4)
else
    plot(data.xToTarg(ff),data.yToTarg(ff),'.','Color',cmap(round(-data.tToTarg(ff)/10100*256),:),'MarkerSize',4)
end
 end
axis([-1250 1250 -750 750])
vline(0,':k')
hline(0,':k')
if lgscale
    hc = colorbar('colormap',cmap,'XTick',(log10([.1:.1:.9,1:10])+1)/2,'XTickLabel',[.1:.1:.9,1:10]);
    doimage(gcf,fullfile(pathexp,'figures','behaviour'),'tiff',['RTdistTlog_' namegr],1)  
else
    hc = colorbar('colormap',cmap,'XTick',0:.1:1,'XTickLabel',[0:1:10]);
        doimage(gcf,fullfile(pathexp,'figures','behaviour'),'tiff',['RTdistT_' namegr],1)  
end   

%%
% as a surface
clear Z
ixfix   = find(data.start>0 & data.type==1 & data.tToTarg<-100 & data.tToTarg>-10000);
xxpos = data.xToTarg(ixfix);
yypos = data.yToTarg(ixfix);
tt    = data.tToTarg(ixfix);
xlin = linspace(-1250,1250,30);
ylin = linspace(-750,750,30);
for xx = 1:length(xlin)-1
    for yy = 1:length(ylin)-1
        auxix = find(xxpos>xlin(xx) & xxpos<xlin(xx+1) & yypos>ylin(yy) & yypos<ylin(yy+1));
        Z(xx,yy) = median(-tt(auxix)/1000);
    end
end
% [X,Y] = meshgrid(xlin(1:end-1)+diff(xlin)/2,ylin(1:end-1)+diff(ylin)/2);
% figure,surf(X,Y,Z)
figure,
imagesc(xlin(1:end-1)+diff(xlin)/2,ylin(1:end-1)+diff(ylin)/2,log10(Z),[-1 1])
   hc = colorbar('colormap',cmap,'XTick',log10([.1:.1:.9,1:10]),'XTickLabel',[.1:.1:.9,1:10]);
  
colormap(cmap)
% f = scatteredInterpolant(data.xToTarg(ixfix)',data.yToTarg(ixfix)',log10(-data.tToTarg(ixfix)'/1000));
% Z = f(X,Y)
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXATION DURATIONS, REFIXATIONS AND REVISITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=1;
SUMMARY.lagBINS = [1.5:1:40,inf];
for ss = subjects

    auxdat    = struct_select(data,{'subject'},{['==' num2str(ss)]},2);
    lastfix                             = [diff(auxdat.trial) 1];
    
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
% figure correlation between reaction times and refix nd revisit
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
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
             'tiff',['corrRT_refix_revisit_' namegr],1)

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
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
              'tiff',['corrrefix_revisit_' namegr],1)

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
axis([0 40 0 .6])
xlabel('Revisit lag','FontSize',12)
ylabel('revisit pre trial','FontSize',12)
title(sprintf('N = %d',Ns),'FontSize',12)
tightfig
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
             'tiff',['revisitlag_' namegr],1)
        

%%
% figure fixation duration different type of fixation
labelsST = {'refixMdur','prerefixMdur','revisitfixMdur','torevisitfixMdur','fixdurnoreMdur'};
labels   = {'refix','prerefix','revisit','torevisit','fixdurOther'};
figure
hold on
jitter = (rand(1,Ns)-.5)/10;
for ll = 1:length(labelsST)
    plot(ll+jitter,SUMMARY.(labelsST{ll}),'ok','MarkerSize',8,'LineWidth',1,'MarkerFaceColor',[.7 .7 .7])
     errorbar(ll,mean(SUMMARY.(labelsST{ll})),std(SUMMARY.(labelsST{ll}))./sqrt(Ns),'LineWidth',2,'Color',[0 0 0])
    plot(ll,mean(SUMMARY.(labelsST{ll})),'sk','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 0 0])
   
end
axis([.5 5.5 0 300])
set(gca,'FontSize',10,'XTick',1:length(labels),'XTickLabel',labels)
ylabel('Fixation duration (ms)','FontSize',12)
tightfig
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
            'tiff',['fixdurs_' namegr],1)

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
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
              'tiff',['fixdurperEO_' namegr],1)
        

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
  doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
               'tiff',['fixdursrevisit_' namegr],1)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATENCY TO MOVE AFTER IMAGE APPEARANCE AND STIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins            = [0:25:1500,inf];
binspostStim    = [0:25:1000,inf];
binspostStimT    = [0:50:1000,inf];
s=1;
nsacs = 5;
clear binstarts binstPostStim
for ss = subjects
     auxdat    = struct_select(data,{'subject'},{['==' num2str(ss)]},2);
%      auxdat.latposStim = auxdat.start-auxdat.stimtime;
     starts     = [];   orders       = [];  condition       = [];        
     stimtime   = [];   startPosStim = [];  orderPosStim    = [];
     posPosStim = [];   condPosStim  = [];  posOrdPosStim   = [];
     for tt = unique(auxdat.trial)
         auxstarts      = auxdat.start(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         auxstartpost   = auxdat.latposStim(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         cond           = auxdat.value(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         pos            = auxdat.posendx(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         stimtime       = [stimtime,auxdat.stimtime(find(auxdat.trial==tt & auxdat.start>0 & auxdat.type==2,nsacs))];
         starts         = [starts,auxstarts];
         auxorder       = 1:length(auxstarts);
         orders         = [orders,1:length(auxstarts)];
         condition      = [condition,cond];
         
         afterSt        = find(auxstartpost>0);
         startPosStim   = [startPosStim auxstartpost(afterSt)];
         posOrdPosStim  = [posOrdPosStim 1:length(afterSt)];
         orderPosStim   = [orderPosStim auxorder(afterSt)];
         posPosStim     = [posPosStim pos(afterSt)];
         condPosStim    = [condPosStim cond(afterSt)];
     end
     % proportion of first saccades after stimulation that are the first,
     % second or third after trial start
     porderPosStim(s,:)     = [sum(posOrdPosStim==1 & orderPosStim ==1) sum(posOrdPosStim==1 & orderPosStim ==2) ...
                                    sum(posOrdPosStim==1 & orderPosStim ==3)]./sum(posOrdPosStim==1); 
     % proportion of first saccades after stimulation that start before 100 ms after stimulation and that they are the frist or second after trial start                           
     porderPosStimEarly(s,:) = [sum(startPosStim(posOrdPosStim==1 & orderPosStim ==1)<100)./sum(posOrdPosStim==1 & orderPosStim ==1) ...
         sum(startPosStim(posOrdPosStim==1 & orderPosStim ==2)<100)./sum(posOrdPosStim==1 & orderPosStim ==2)];
     %proportion of 1st,2nd,3rd saccades after stimulation that are to the
     %right and the ones tat start after 100 ms
     for pops = 1:3   
        pRightpCond(s,:,pops)        = accumarray(condPosStim(posOrdPosStim==pops)',posPosStim(posOrdPosStim==pops)>posVec.stimRes(1)/2)'./accumarray(condPosStim(posOrdPosStim==pops)',ones(sum(posOrdPosStim==pops),1))';
        pRightpCondLate(s,:,pops)    = accumarray(condPosStim(startPosStim>100 & posOrdPosStim==pops)',posPosStim(startPosStim>100 & posOrdPosStim==pops)>posVec.stimRes(1)/2)'./...
        accumarray(condPosStim(startPosStim>100 & posOrdPosStim==pops)',ones(sum(startPosStim>100 & posOrdPosStim==pops),1))';
     end
     % start time after trial start
     for cV = SUMMARY.condValues
         for sacor = 1:nsacs
            auxindx = orders==sacor & condition == cV;
            binstarts(s,:,sacor,cV) = histc(starts(auxindx),bins)./sum(auxindx);
         end
     end
     % start time after stimulation and p(right after stimulation)
     for cV = SUMMARY.condValues
         for sacor = 1:3
             auxindx = posOrdPosStim==sacor & condPosStim == cV;
            binstPostStim(s,:,sacor,cV) = histc(startPosStim(auxindx),binspostStim)./sum(auxindx);
         end
         auxindx = condPosStim == cV;
         [N,BIN] = histc(startPosStim(auxindx),binspostStimT);
         auxpos  = posPosStim(auxindx);
         pRightpCondTime(s,:,cV) = [accumarray(BIN',auxpos>posVec.stimRes(1)/2)'./N(1:max(BIN)) zeros(1,length(N)-1-max(BIN))];
       
     end
     s = s+1;
end

%%
% Distribution latency to move plot
cmap1 = cbrewer('qual','Set1',nsacs);
cmap2 = cbrewer('qual','Pastel1',nsacs);
figure,hold on 
for sacor = 1:nsacs
    plot(bins+12.5,squeeze(nanmean(binstarts(:,:,sacor,:),4)),'Color',cmap2(sacor,:))
end
for sacor = 1:nsacs
    h(sacor) = plot(bins+12.5,squeeze(mean(nanmean(binstarts(:,:,sacor,:),4))),'Color',cmap1(sacor,:),...
        'LineWidth',2);
end
xlabel('Latency to # movement','FontSize',12)
ylabel('Relative Frequency','FontSize',12)
legend(h,{'#1 Mov','#2 Mov','#3 Mov','#4 Mov','#5 Mov'})

clear h
  doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
                'tiff',['latToMove_' namegr],1)

% and per condition
figure, hold on
for cV = SUMMARY.condValues
     for sacor = 1:5
         if cV<9
             h(sacor) = plot(bins+12.5,squeeze(mean(binstarts(:,:,sacor,cV))),'Color',cmap1(sacor,:),...
        'LineWidth',2);
         else
             h2(sacor) = plot(bins+12.5,squeeze(mean(binstarts(:,:,sacor,cV))),'Color',cmap2(sacor,:),...
        'LineWidth',2);
         end
     end
end
xlabel('Latency to # movement per condition','FontSize',12)
ylabel('Relative Frequency','FontSize',12)
legend(h,{'#1 Mov Info','#2 Mov Info','#3 Mov Info','#4 Mov Info','#5 Mov Info'})
  doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
                'tiff',['latToMovePerCond_' namegr],1)
clear h
% latency to move after stimulatio
figure, hold on
for sacor = 1:3
    plot(binspostStim+12.5,nanmean(binstPostStim(:,:,sacor,:),4),'Color',cmap2(sacor,:))
end
for sacor = 1:3
    h(sacor) = plot(binspostStim+12.5,squeeze(mean(nanmean(binstPostStim(:,:,sacor,:),4))),'Color',cmap1(sacor,:),...
        'LineWidth',2);
end
xlabel('Latency to # movement after stimulus','FontSize',12)
ylabel('Relative Frequency','FontSize',12)
legend(h,{'#1 Mov','#2 Mov','#3 Mov'})
 doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
                 'tiff',['latToMoveAfterStim_' namegr],1)
clear h

% and per condition
figure, hold on
for cV = SUMMARY.condValues
     for sacor = 1:3
         if cV<9
            h(sacor) = plot(binspostStim+12.5,squeeze(mean(binstPostStim(:,:,sacor,cV))),'Color',cmap1(sacor,:),...
            'LineWidth',2);
         else
            h(sacor) = plot(binspostStim+12.5,squeeze(mean(binstPostStim(:,:,sacor,cV))),'Color',cmap2(sacor,:),...
            'LineWidth',2);  
         end
     end
end
xlabel('Latency to # movement after stimulus','FontSize',12)
ylabel('Relative Frequency','FontSize',12)
legend(h,{'#1 Mov Info','#2 Mov Info','#3 Mov Info'})
 doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
                 'tiff',['latToMoveAfterStimPerCond_' namegr],1)
             
%%
% proportion of movements

% order of movement after stimulation
paftStim = [porderPosStim(:,1), porderPosStim(:,1).*porderPosStimEarly(:,1), porderPosStim(:,1).*(1-porderPosStimEarly(:,1)),...
    porderPosStim(:,2), porderPosStim(:,2).*porderPosStimEarly(:,2), porderPosStim(:,2).*(1-porderPosStimEarly(:,2))];
jitter = (rand(1,Ns)-.5)/10;
figure,hold on
xpos = [1,2,3,5,6,7];
labels = {'1st','1st <100','1st >100','2nd','2nd <100','2nd >100'};
for e=1:6
    plot(xpos(e)+jitter,paftStim(:,e),'ok','MarkerSize',6,'LineWidth',.5,...
        'MarkerFaceColor',[.9 .9 .9],'MarkerEdgeColor',[.5 .5 .5])
end
    errorbar(xpos,mean(paftStim),std(paftStim)./sqrt(Ns),'LineWidth',2,'Color',[.1 .1 .1])

plot(xpos,mean(paftStim),'sk','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[1 .3 .3])

set(gca,'XTick',xpos,'XTickLabels',labels,'FontSize',10)
ylabel('p()','FontSize',12)
xlabel('Order from trial start of first movement after stim','FontSize',12)
 doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
                'tiff',['OrderMovAfterStim_' namegr],1)

%%
% movements to the right by order after stim
xpos    = [1,2,4,5,8,9,11,12];
cV      = SUMMARY.condValues;
labels  = SUMMARY.condLabels;
cmap    = cbrewer('div','RdYlBu',11);
figure, hold on
bR = bar(squeeze(mean(pRightpCond(:,cV,:))));
bL = bar(-1+squeeze(mean(pRightpCond(:,cV,:))),'FaceColor',[0 0 1]);
for e = 1:3, 
    bR(e).XData = xpos;
    bR(e).FaceColor = cmap(6-e,:);
    bL(e).XData = xpos;
    bL(e).FaceColor = cmap(6+e,:);
end
view([90 90])
leg = legend(bR,{'#1','#2','#3'},'Location','SouthEast');
legend(gca,'boxoff')

ylabel('p(left/right field)','FontSize',12)
set(gca,'XTick',xpos,'XTickLabel',labels)

  doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
               'tiff',['probRight_' namegr],1)

%%
% movement to the right by time after stimulation

cmap1 = cbrewer('qual','Set1',4);
cmap2 = cbrewer('qual','Pastel1',4);
figure, hold on
px = 1; clear h
for cV = SUMMARY.condValues(1:4)
    Mval    = nanmean(pRightpCondTime(:,:,cV));
    SE      = nanstd(pRightpCondTime(:,:,cV))./sqrt(Ns);
%     plot(binspostStimT(1:end-1)+25,pRightpCondTime(:,:,cV),'-','LineWidth',.5,'Color',cmap2(px,:));
    jbfill(binspostStimT(1:end-1)+25,Mval+SE,Mval-SE,cmap2(px,:),cmap2(px,:),1,.5),hold on
    h(px)   = plot(binspostStimT(1:end-1)+25,Mval,'.-','LineWidth',2,'Color',cmap1(px,:));

    px = px+1;
end
hline(0.5)
view([90 90])
axis([0 1000 0 1])
ylabel('p(right)','FontSize',12)
xlabel('time (ms)','FontSize',12)
legend(h,SUMMARY.condLabels(1:4),'Location','Best')
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
               'tiff',['probRightAfterStimInf_' namegr],1)
figure, hold on
px = 1; clear h
for cV = SUMMARY.condValues(5:8)
    Mval    = nanmean(pRightpCondTime(:,:,cV));
    SE      = nanstd(pRightpCondTime(:,:,cV))./sqrt(Ns);
%     plot(binspostStimT(1:end-1)+25,pRightpCondTime(:,:,cV),'-','LineWidth',.5,'Color',cmap2(px,:));
    jbfill(binspostStimT(1:end-1)+25,Mval+SE,Mval-SE,cmap2(px,:),cmap2(px,:),1,.5),hold on
    h(px)   = plot(binspostStimT(1:end-1)+25,Mval,'.-','LineWidth',2,'Color',cmap1(px,:));

    px = px+1;
end
hline(0.5)
view([90 90])
axis([0 1000 0 1])
legend(h,SUMMARY.condLabels(5:8),'Location','Best')
ylabel('p(right)','FontSize',12)
xlabel('time (ms)','FontSize',12)
doimage(gcf,fullfile(pathexp,'figures','behaviour'),...
               'tiff',['probRightAfterStimUni_' namegr],1)
%%
% saccades gaze after stimulation per condition, moveemnt order and latency
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATENCY TO MOVE AFTER IMAGE APPEARANCE AND STIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=1;
nsacs = 5;
load(fullfile(pathexp,'analysis','eyedata','allsampledata.mat')) 
clear allgazes
for ss = subjects
    auxdat      = struct_select(data,{'subject'},{['==' num2str(ss)]},2);
    auxsam      = struct_select(sample,{'subject'},{['==' num2str(ss)]},2);

    indxSac     = auxdat.type==2 & (auxdat.orderposStim ==1 |auxdat.orderposStim ==2);
    latposStim  = auxdat.latposStim(indxSac);
    sttimes     = auxdat.start(indxSac);
    endtimes    = auxdat.end(indxSac);
    trial       = auxdat.trial(indxSac);
    cond        = auxdat.value(indxSac);
    order       = auxdat.orderposStim(indxSac);
    clear gaze gazesac
    for ev = 1:sum(indxSac)
        indxev   = find(auxsam.trial == trial(ev) & auxsam.time>sttimes(ev) & auxsam.time<endtimes(ev));
        gaze{ev} = auxsam.pos(:,indxev);
        indxev   = find(auxsam.trial == trial(ev) & auxsam.time>sttimes(ev),1);
        gazesac{ev} = auxsam.pos(:,indxev:indxev+50);
    end
    allgazes(s).gaze        = gaze;
    allgazes(s).gazesac     = gazesac;
    allgazes(s).latposStim  = latposStim;
    allgazes(s).cond        = cond;
    allgazes(s).order       = order;
    s=s+1;
end
%%
subpltmap = [1 2 5 6 3 4 7 8];
cV        = SUMMARY.condValues;
for ss = 1:length(allgazes)
    for gOr = 1:2
        figure
        set(gcf,'Position',[30 53 1251 652])
        
        for ev = 1:length(allgazes(ss).gaze)
            if allgazes(ss).order(ev)==gOr
                wsP = find(cV==allgazes(ss).cond(ev));
                subplot(2,4,subpltmap(wsP)),hold on
                gazetoplot      = allgazes(ss).gazesac{ev};
                if ~isempty(gazetoplot)
                    gazetoplot      = gazetoplot-repmat(gazetoplot(:,1),1,size(gazetoplot,2));
                    gazetoplot(2,:) = gazetoplot(2,:)+allgazes(ss).latposStim(ev);
                    if any(gazetoplot(1,:)<0) & any(gazetoplot(1,:)>0)
                        plot(gazetoplot(1,:),gazetoplot(2,:),'.','Color',[1 0 0],'MarkerSize',2)
                    else
                        plot(gazetoplot(1,:),gazetoplot(2,:),'.','Color',[0 0 1],'MarkerSize',2)
                    end
                end
            end
        end
    
        for wsP=1:8
            subplot(2,4,subpltmap(wsP))
            axis([-400 400 -300 800])
            vline(0)
            hline(0)
            title(SUMMARY.condLabels(wsP))
        end
         doimage(gcf,fullfile(pathexp,'figures','behaviour','sacpersubj'),...
                   'tiff',sprintf('s%d_postStim%d_%s',subjects(ss),gOr,namegr),1)
    end
end