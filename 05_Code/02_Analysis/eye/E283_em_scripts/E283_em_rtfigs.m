% DISTRIBUTION OF RT, for each subject and pooled, plus mean RT and
% performance per subject and total

fh = figure, hold on
fh.Units = 'centimeters'
fh.Position = [100 100 4.6 4.6*fh.Position(4)/fh.Position(3)];
plot(SUMMARY.rTBINS+unique(diff(SUMMARY.rTBINS))/2,[SUMMARY.rTpbin./repmat(sum(SUMMARY.rTpbin,2),1,size(SUMMARY.rTpbin,2))]',...
    'Color',[.7 .7 .7],'LineWidth',.5)
plot(SUMMARY.rTBINS+unique(diff(SUMMARY.rTBINS))/2,sum(SUMMARY.rTpbin./repmat(sum(SUMMARY.rTpbin,2),1,size(SUMMARY.rTpbin,2)))./Ns,...
    'Color',[0 0 0],'LineWidth',2)
plot(SUMMARY.rT,0.225+(rand(1,Ns)-.5)/100,'ok','MarkerSize',2,'LineWidth',.5,'MarkerFaceColor',[.9 .9 .9])
line([SUMMARY.meanRT-SUMMARY.meanRTse SUMMARY.meanRT+SUMMARY.meanRTse],[.225 .225],'LineWidth',.5,'Color',[1 0 0])
line([SUMMARY.meanRT SUMMARY.meanRT],[.225-.005 .225+.005],'LineWidth',.5,'Color',[1 0 0])
%plot(SUMMARY.meanRT,0.225,'+k','MarkerSize',6,'LineWidth',1,'MarkerFaceColor',[1 0 0])
text(max(SUMMARY.rT)+range(SUMMARY.rT)*.3,0.225,sprintf('mean: %2.2f +- %2.2f',SUMMARY.meanRT,SUMMARY.meanRTse),'FontSize',6)
plot(13+(rand(1,Ns)-.5)/3,1-SUMMARY.perf,'ok','MarkerSize',2,'LineWidth',.5,'MarkerFaceColor',[.9 .9 .9])
set(gca,'XTick',[SUMMARY.rTBINS(1:4:end) 13],'XTickLabel',{'0','2','4','6','8','10','', 'misses'},...
    'FontSize',6)
axis([0 14 0 .25])
%axis square
xlabel('Reaction Time (s)','FontSize',8)
ylabel('Frequency','FontSize',8)
% title(sprintf('N = %d',Ns),'FontSize',8)
%set(gca,'Position',[0.05 0.1 0.95 0.6])
%tightfig
figsize     = [4.6 4.6*fh.Position(4)/fh.Position(3)];
%  doimage(gcf,fullfile(patheye,'figures'),...
%              'pdf',['performance_' namegr],figsize  ,1)

         %%
% RT per condition (stimulation side, limb crossing, and informativeness) per subject meand and SE  
fh = figure
hold on
jitter = (rand(1,Ns)-.5)/10
xpos = [1:4,6:9];
cmap2 = cbrewer('qual','Pastel1',2);
for ll = 1:length(SUMMARY.condLabels)
    if ll<5
    col = cmap2(2,:)
    else
        col = cmap2(1,:)
    end
    plot1 = plot(xpos(ll)+jitter,SUMMARY.rTpCond(:,ll),'ok','MarkerSize',4,'LineWidth',.5,...
        'MarkerFaceColor',col,'MarkerEdgeColor',[.5 .5 .5])
end 
   cmap1 = cbrewer('qual','Set1',2);
errorbar(1:4,SUMMARY.meanRTpCond(1:4),SUMMARY.seRTCond(1:4),'LineWidth',1,'Color',[0 0 0])
%plot(1:4,SUMMARY.meanRTpCond(1:4),'sk','MarkerSize',4,'LineWidth',.5,'MarkerFaceColor',[.3 .3 1])
errorbar(6:9,SUMMARY.meanRTpCond(5:8),SUMMARY.seRTCond(5:8),'LineWidth',1,'Color',[0 0 0])
%plot(6:9,SUMMARY.meanRTpCond(5:8),'sk','MarkerSize',4,'LineWidth',.5,'MarkerFaceColor',[1 .3 .3])
    
    axis([0 10 0 5.2])
   % axis square
set(gca,'FontSize',6,'XTick',xpos,'XTickLabel',{'LU','RU','LC','RC','LU','RU','LC','RC'},'YTick',0:2.5:5)
ylabel('Reaction Time (s)','FontSize',8)
%xlabel('Informative    Uninformative','FontSize',10)
%tightfig
figsize     = [4.6 4.6*fh.Position(4)/fh.Position(3)];
%  doimage(gcf,fullfile(patheye,'figures'),...
%              'pdf',['RTpercond_' namegr],figsize,1)
RTtable =array2table(SUMMARY.rTpCond,'VariableNames',SUMMARY.condLabels)
Within = table({'L';'R';'L';'R';'L';'R';'L';'R'},{'U';'U';'C';'C';'U';'U';'C';'C'},{'I';'I';'I';'I';'unI';'unI';'unI';'unI'},'VariableNames',{'Side','Crossing','Info'}) 
RTrm = fitrm(RTtable,'LUI-RCunI~1','WithinDesign',Within)
[h,p,ci,stats]=ttest(mean(SUMMARY.rTpCond(:,1:4),2),mean(SUMMARY.rTpCond(:,5:8),2))
RTrm.multcompare('Info')
RTanova = ranova(RTrm,'WithinModel','Side*Crossing*Info')

RTrmI = fitrm(RTtable(:,1:4),'LUI-RCI~1','WithinDesign',Within(1:4,:))
RTanovaI = ranova(RTrmI,'WithinModel','Side*Crossing')
[h,p,ci,stats]=ttest(mean(SUMMARY.rTpCond(:,3:4),2),mean(SUMMARY.rTpCond(:,1:2),2))

RTrmI.multcompare('Side','By','Crossing')

RTrmunI = fitrm(RTtable(:,5:8),'LUunI-RCunI~1','WithinDesign',Within(5:8,:))
RTanovaunI = ranova(RTrmunI,'WithinModel','Side*Crossing')

%%
% figure reaction time per target position with horizontal and vertical
% marginals

% color plot
cmap = cmocean('haline');
fh = figure;
set(gcf,'Position',[50 50 600 400])
subplot(7,9,[1:9:46]), hold on
errorbar(1:6,mean(mean(reshape(SUMMARY.rtpPos',[6 8 Ns]),2),3),std(mean(reshape(SUMMARY.rtpPos',[6 8 Ns]),2),1,3)./sqrt(Ns),...
    'Color',[0 0 0])
plot(1:6,mean(mean(reshape(SUMMARY.rtpPos',[6 8 Ns]),2),3),'s',...
    'Color',[0 0 0],'LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[1 .4 .4])
axis([.5 6.5 2 4])
view([90 90])
set(gca,'FontSize',10)

% marginals
subplot(7,9,repmat(2:9,1,6)+reshape(repmat([0 9:9:45],8,1),[1 48]))
imagesc(reshape(mean(SUMMARY.rtpPos),[6,8]))
axis off
caxis([1.5 3.5])
 colorbar
 title('RT per target position (s)','FontSize',12)
 
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
figsize     = [8 8*fh.Position(4)/fh.Position(3)];
doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['RTperPos_' namegr],[],1)        

%%
% figure reaction time at any fixation per distance to the target
lgscale = 1;
trimming = 5;
ixfix   = find(data.start>0 & data.type==1 & data.tToTarg<-100 & data.tToTarg>-10000);
ixfix   = ixfix(1:trimming:end);
figure,hold on
cmap    = flipud(cmocean('thermal'));
for ff = ixfix % the loop is used to plot each point with a different color
    if lgscale
        plot(data.xToTarg(ff),data.yToTarg(ff),'.','Color',cmap(round((log10(-data.tToTarg(ff)/1000)+1)/2*256),:),'MarkerSize',4)
    else
        plot(data.xToTarg(ff),data.yToTarg(ff),'.','Color',cmap(round(-data.tToTarg(ff)/10100*256),:),'MarkerSize',4)
    end
end
axis([-1250 1250 -750 750])
vline(-7*posVec.dsX+posVec.dsX/2:posVec.dsX:7*posVec.dsX-posVec.dsX/2,'k')
hline(-6*posVec.dsY+posVec.dsY/2:posVec.dsY:6*posVec.dsY-posVec.dsY/2,'k')
title('RT from fix at x dist')
if lgscale
    hc = colorbar('colormap',cmap,'XTick',(log10([.1:.1:.9,1:10])+1)/2,'XTickLabel',[.1:.1:.9,1:10]);
    doimage(gcf,fullfile(patheye,'figures'),'tiff',['RTdistTlog_' namegr],[],1)  
else
    hc = colorbar('colormap',cmap,'XTick',0:.1:1,'XTickLabel',[0:1:10]);
        doimage(gcf,fullfile(patheye,'figures'),'pdf',['RTdistT_' namegr],[],1)  
end   

% same but binned
clear Z
ixfix   = find(data.start>0 & data.type==1 & data.tToTarg<-100 & data.tToTarg>-10000);
xxpos = data.xToTarg(ixfix);
yypos = data.yToTarg(ixfix);
tt    = data.tToTarg(ixfix);
xlin = -7*posVec.dsX+posVec.dsX/2:posVec.dsX/2:7*posVec.dsX-posVec.dsX/2;
ylin = -6*posVec.dsY+posVec.dsY/2:posVec.dsY/2:6*posVec.dsY-posVec.dsY/2;
for xx = 1:length(xlin)-1
    for yy = 1:length(ylin)-1
        auxix = find(xxpos>xlin(xx) & xxpos<xlin(xx+1) & yypos>ylin(yy) & yypos<ylin(yy+1));
        Z(xx,yy) = median(-tt(auxix)/1000);
    end
end
% [X,Y] = meshgrid(xlin(1:end-1)+diff(xlin)/2,ylin(1:end-1)+diff(ylin)/2);
% figure,surf(X,Y,Z)
figure,
imagesc(xlin(1:end-1)+diff(xlin)/2,ylin(1:end-1)+diff(ylin)/2,log10(Z)',[-1 1])
vline(xlin(1:2:end),':k')
hline(ylin(1:2:end),':k')
title('RT from fix at x dist')
colormap(cmap)
if lgscale
    hc = colorbar('colormap',cmap,'XTick',(log10([.1:.1:.9,1:10])+1)/2,'XTickLabel',[.1:.1:.9,1:10]);
    doimage(gcf,fullfile(patheye,'figures'),'tiff',['RTdistTbinlog_' namegr],[],1)  
else
    hc = colorbar('colormap',cmap,'XTick',0:.1:1,'XTickLabel',[0:1:10]);
        doimage(gcf,fullfile(patheye,'figures'),'pdf',['RTdistbinT_' namegr],[],1)  
end  
% f = scatteredInterpolant(data.xToTarg(ixfix)',data.yToTarg(ixfix)',log10(-data.tToTarg(ixfix)'/1000));
% Z = f(X,Y)