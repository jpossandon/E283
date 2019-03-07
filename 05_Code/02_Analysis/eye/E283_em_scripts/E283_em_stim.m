mcoef      = [0 1 0 1 0 1 0 1;... %side stim
                   0 0 1 1 0 0 1 1;...  %crossing 
                   1 1 1 1 0 0 0 0]';      % info
binsize         = 50;
bins            = [0:binsize:1500,inf];
binspostStim    = [0:binsize:1000,inf];
binspostStimT    = [0:binsize:1000,inf];
RTrightpostStimCond = [];
 RTleftpostStimCond = [];
 pRightpCondTime = [];
s=1;
nsacs = 5;
clear binstarts binstPostStim
 YY = []; XY = []; 
for ss = subjects
     auxdat    = struct_select(data,{'subject'},{['==' num2str(ss)]},2);
%      auxdat.latposStim = auxdat.start-auxdat.stimtime;
     starts     = [];   orders       = [];  condition       = [];        
     stimtime   = [];   startPosStim = [];  orderPosStim    = [];
     posPosStim = [];   condPosStim  = [];  posOrdPosStim   = [];
     dirPosStim = [];
     for tt = unique(auxdat.trial)
         auxstarts      = auxdat.start(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         auxstartpost   = auxdat.latposStim(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         cond           = auxdat.value(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         pos            = auxdat.posendx(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
         dirpos         = auxdat.posendx(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs))-auxdat.posinix(find(auxdat.trial==tt  & auxdat.start>0 & auxdat.type==2,nsacs));
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
         dirPosStim     = [dirPosStim dirpos(afterSt)];
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
     
    cVV = 1;
     for cV = SUMMARY.condValues
         for sacor = 1:3
             auxindx = posOrdPosStim==sacor & condPosStim == cV;
            binstPostStim(s,:,sacor,cV) = histc(startPosStim(auxindx),binspostStim)./sum(auxindx);
            auxindxR = posOrdPosStim==sacor & condPosStim == cV & posPosStim>posVec.stimRes(1)/2;
            auxindxL = posOrdPosStim==sacor & condPosStim == cV & posPosStim<posVec.stimRes(1)/2;
            %auxindxR = posOrdPosStim==sacor & condPosStim == cV & dirPosStim>0;
            %auxindxL = posOrdPosStim==sacor & condPosStim == cV & dirPosStim<0;
            [NR,BINR]  = histc(startPosStim(auxindxR),binspostStimT);
            [NL,BINL]  = histc(startPosStim(auxindxL),binspostStimT);
            RTrightpostStimCond(s,:,sacor,cV) = NR./(sum(auxindxR)+sum(auxindxL));
            RTleftpostStimCond(s,:,sacor,cV) = NL./(sum(auxindxR)+sum(auxindxL));
       %     if sacor==1
                 YY = [YY;[BINR';BINL'],[ones(sum(NR),1);zeros(sum(NL),1)]];
                 XY = [XY;[repmat(mcoef(cVV,:),length(BINR)+length(BINL),1) s*ones(length(BINR)+length(BINL),1)]];   
        %    end
         end
         auxindx = condPosStim == cV;
         [N,BIN] = histc(startPosStim(auxindx),binspostStimT);
         auxpos  = posPosStim(auxindx);
         pRightpCondTime(s,:,cV) = [accumarray(BIN',auxpos>posVec.stimRes(1)/2)'./N(1:max(BIN)) zeros(1,length(N)-1-max(BIN))];
       cVV = cVV+1;
     end
     
    
     s = s+1;
end

%%
% Distribution latency to move for different orders of movement after trials start
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
%   doimage(gcf,fullfile(patheye,'figures'),...
%                 'tiff',['latToMove_' namegr],[],1)

%%
% and per info/noInfo condition
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
%   doimage(gcf,fullfile(patheye,'figures'),...
%                 'tiff',['latToMovePerCond_' namegr],[],1)
clear h
%%
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
%  doimage(gcf,fullfile(patheye,'figures'),...
%                  'tiff',['latToMoveAfterStim_' namegr],[],1)
clear h

%%
% and per info/noinfo condition
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
%   doimage(gcf,fullfile(patheye,'figures'),...
%                   'tiff',['latToMoveAfterStimPerCond_' namegr],[],1)
%%
% latency to move per condition and side of the screen landing
fh = figure;
barxpos    = [-120 -90 -60];

fh.Position = [198 354 1400 600];
px = 1;
% cmap1 = cbrewer('qual','Set2',4);
cmap1 = flipud(cbrewer('seq','Blues',9));
cmap2 = cmap1(2:2:6,:);
cmap1 = cmap1(1:2:5,:);
% cmap2 = cbrewer('qual','Pastel2',4);
% cmap3 = cbrewer('qual','Set1',4);
plim = .25;
for cV = SUMMARY.condValues
    if cV<9, colname = 'Blues';else
        colname = 'Reds';
    end
    cmap1 = flipud(cbrewer('seq',colname,9));
    cmap2 = cmap1(2:2:6,:);
    cmap1 = cmap1(1:2:5,:);
    subplot(2,4,px), hold on
    
    bR = bar(plim.*squeeze(mean(pRightpCond(:,cV,:)))');
    bL = bar(-plim+plim.*squeeze(mean(pRightpCond(:,cV,:)))');
   
        bR.XData = barxpos;
        bR.FaceColor = cmap1(1,:);
        bL.XData = barxpos;
        bL.FaceColor = cmap1(1,:);
    bR.BarWidth = .9;
     bL.BarWidth = .9;
     bR.LineWidth = .1;
     bL.LineWidth = .1;
   for sacor = 1:3
        Mval    = nanmean(RTrightpostStimCond(:,(1:end-1),sacor,cV));
        SE      = nanstd(RTrightpostStimCond(:,(1:end-1),sacor,cV))./sqrt(Ns);
        %jbfill(binspostStimT(1:end-1)+25,Mval+SE,Mval-SE,cmap2(sacor,:),cmap2(sacor,:),1,.5);,hold on
        h   = plot(binspostStimT(1:end-1)+binsize/2,Mval,'.:','LineWidth',.5,'Color',cmap1(sacor,:));
        Mval    = nanmean(RTleftpostStimCond(:,(1:end-1),sacor,cV));
        SE      = nanstd(RTleftpostStimCond(:,(1:end-1),sacor,cV))./sqrt(Ns);
       % jbfill(binspostStimT(1:end-1)+25,-(Mval+SE),-(Mval-SE),cmap2(sacor,:),cmap2(sacor,:),1,.5);,hold on
        h   = plot(binspostStimT(1:end-1)+binsize/2,-Mval,'.:','LineWidth',.5,'Color',cmap1(sacor,:));
          
       % hline(0)
       % vline(200:200:800,':k')
        vline([0],'-k')
        box on
        view([90 90])
        axis([-180 1000 -plim plim])
        
        if ismember(px,[1,5])
            set(gca,'XTick',0:100:1000,'XTickLabel',{'0','','200','','400','','600','','800','','1000'})
        else
           set(gca,'XTick',0:100:1000,'XTickLabel',{})
        end
        
   end
   if ismember(px,[5:8])
    set(gca,'YTick',-0.25:.125:.25,'YTickLabel',{'0.25','','0','','0.25'})
    title(SUMMARY.condLabels{px}(1:2))
    ylabel('Relative Frequency')
    else
         set(gca,'YTick',-0.25:.125:.25,'YTickLabel',{})
    end
    for sacor = 1:3
        Mval    = nanmean(RTrightpostStimCond(:,(1:end-1),sacor,cV)-RTleftpostStimCond(:,(1:end-1),sacor,cV));
        SE      = nanstd(RTrightpostStimCond(:,(1:end-1),sacor,cV)-RTleftpostStimCond(:,(1:end-1),sacor,cV))./sqrt(Ns);
        jbfill(binspostStimT(1:end-1)+25,(Mval+SE),(Mval-SE),cmap2(sacor,:),cmap2(sacor,:),1,.5,.75);,hold on
       h   = plot(binspostStimT(1:end-1)+25,Mval,'LineWidth',.5,'Color',cmap1(sacor,:));
      
    end
    if px==1 || px==5
        xlabel('Reaction Time (ms)','FontSize',10)
    end
   px=px+1; 
end
figsize     = [4*4.6 4*4.6*fh.Position(4)/fh.Position(3)];
 doimage(gcf,fullfile(patheye,'figures'),'pdf',['latToMoveperSidePerCond_' namegr],'600','painters',figsize,1)
% axis([0 1000 0 1])
% ylabel('p(right)','FontSize',12)
% xlabel('time (ms)','FontSize',12)




%%
% proportion of movements

% order of movement (from trial start) after stimulation
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
 doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['OrderMovAfterStim_' namegr],[],1)

%%
% movement to the right by order after stim
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

%  doimage(gcf,fullfile(patheye,'figures'),...
 %              'tiff',['probRight_' namegr],[],1)

%%
% movement to the right by time after stimulation

%info
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
doimage(gcf,fullfile(patheye,'figures'),...
               'tiff',['probRightAfterStimInf_' namegr],[],1)
% not info
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
doimage(gcf,fullfile(patheye,'figures'),...
               'tiff',['probRightAfterStimUni_' namegr],[],1)
           
           
           
           %%
%%
% logistic model ALL
clc
numbi = 13;
bfmc = (numbi-1)*2;
probh = @(x) 1./(1+exp(-x));
        
for unin = [0,1]
     for nb = 1:numbi
             indx = find(YY(:,1)==nb & XY(:,2)==unin);
             if length(indx)>10
                 %auxXY = [XY(indx,1:2),XY(indx,1).*XY(indx,2)];
                 %auxY = YY(indx,2);
                 %auxXY(auxXY==0) = -1;
                 if unin==1
                    auxY = (YY(indx,2)==0 & XY(indx,1)==1) | (YY(indx,2)==1 & XY(indx,1)==0);
                 else
                    auxY = (YY(indx,2)==0 & XY(indx,1)==0) | (YY(indx,2)==1 & XY(indx,1)==1); 
                 end
                 sum(auxY)/length(indx);
                 %auxXY = [XY(indx,3)];
                 %nb
                 %auxres = fitglm(auxXY,auxY,'Distribution','Binomial')
                 auxXY = [XY(indx,[3,4])];
                 T = table(auxY,auxXY(:,1),auxXY(:,2),'VariableNames',{'p_ext','cross','subject'});
                 auxres = fitglme(T,'p_ext ~ cross + (1|subject)','Distribution','Binomial')
             else

             end
             coeffs     = auxres.Coefficients.Estimate;
             coeffsSE   = auxres.Coefficients.SE;
             df         = auxres.Coefficients.DF;
             critT      = [tinv(.05/2/bfmc,df),tinv(1-.05/2/bfmc,df)];
             coeffsCI   = [coeffs+coeffsSE.*critT(:,1) coeffs coeffs+coeffsSE.*critT(:,2)]; 
             p1(nb,:,unin+1)         = probh(coeffsCI(1,:));
             p2(nb,:,unin+1)         = probh(sum(coeffsCI));
             OR(nb,:,unin+1)         = exp(coeffsCI(2,:));
             fprintf('bin %d, p_external uninfo %1.3f\n',nb,1/(1+exp(-([1 0]*auxres.Coefficients.Estimate))))
             fprintf('bin %d, p_external info %1.3f\n',nb,1/(1+exp(-([1 1]*auxres.Coefficients.Estimate))))
             fprintf('bin %d, OR external info vs uninfo %1.2f\n\n',nb,exp(auxres.Coefficients.Estimate(2)))
     end
     
end
%%
cmap2 = cbrewer('qual','Pastel1',4);
labs = {'Uncross','Cross'};
ORfix = 60;
for unin = [0,1]
fh = figure; hold on
 plot(binspostStimT(1:numbi)+binsize/2,OR(:,2,unin+1)/ORfix,'.--','MarkerSize',10,'Color',[0 0 0],'LineWidth',1)
    
     jbfill(binspostStimT(1:numbi)+binsize/2,p1(:,1,unin+1)',p1(:,3,unin+1)',cmap2(1,:),cmap2(1,:),1,.5),hold on
     plot(binspostStimT(1:numbi)+binsize/2,p1(:,2,unin+1),'.-','MarkerSize',10,'Color',cmap1(1,:),'LineWidth',1)
     sigind = sum(p1(:,:,unin+1)>.5,2)==3 | sum(p1(:,:,unin+1)<.5,2)==3;
     if any(sigind)
         plot(binspostStimT(sigind)+binsize/2,p1(sigind,2,unin+1),'*k','MarkerSize',4,'Color',cmap1(1,:),'LineWidth',1)
     end
     jbfill(binspostStimT(1:numbi)+binsize/2,p2(:,1,unin+1)',p2(:,3,unin+1)',cmap2(2,:),cmap2(2,:),1,.5),hold on
     plot(binspostStimT(1:numbi)+binsize/2,p2(:,2,unin+1),'.-','MarkerSize',10,'Color',cmap1(2,:),'LineWidth',1)
     sigind = sum(p2(:,:,unin+1)>.5,2)==3 | sum(p2(:,:,unin+1)<.5,2)==3;
     
     if any(sigind)
         plot(binspostStimT(sigind)+binsize/2,p2(sigind,2,unin+1),'*k','MarkerSize',4,'Color',cmap1(2,:),'LineWidth',1)
     end
     set(gca,'XTick',0:100:600,'XTickLabel',{'0','','200','','400','','600'},'YTick',[0:.25:1],'FontSize',6)
     xlabel('Time (ms)','FontSize',10)
     if unin==0
        ylabel('p(cue side)','FontSize',8)
     else
         ylabel('p(cue side & external)','FontSize',8)
     end
     box off
     hline(.5,':k') ,hold on
     line([600 600],[0 .5],'Color',[0 0 0],'LineWidth',.5)
     line([595 600],[.25 .25],'Color',[0 0 0],'LineWidth',.5)
     line([595 600],[.5 .5],'Color',[0 0 0],'LineWidth',.5)
     hline(1/60,':r')
     text(610,.25,num2str(ORfix/4),'FontSize',6)
     text(610,.5,num2str(ORfix/2),'FontSize',6)
   %  text(690,.19,'OR','Rotation',90,'FontSize',10)
     axis([0 600 0 1])
     figsize     = [4.6 4.6*fh.Position(4)/fh.Position(3)];
     
 doimage(fh,fullfile(patheye,'figures'),'pdf',[labs{unin+1} '_pext'],'opengl',figsize,1)

end
