% DISTRIBUTION OF RT, for each subject and pooled, plus mean RT and
% performance per subject and total

figure, hold on
plot(SUMMARY.rTBINS+unique(diff(SUMMARY.rTBINS))/2,[SUMMARY.rTpbin./repmat(sum(SUMMARY.rTpbin,2),1,size(SUMMARY.rTpbin,2))]',...
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
 doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['performance_' namegr],[],1)

% RT per condition (stimulation side, limb crossing, and informativeness) per subject meand and SE  
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
 doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['RTpercond_' namegr],[],1)

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