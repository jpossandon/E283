%%
% PROBABILITY TO MOVE TO TARGET FIRST
binsize     = 2.5;
bind        = [0:binsize:15,+inf];
result.firstsubjN       = zeros(length(bind)-1,length(subjects));
result.firstNN          = zeros(length(bind)-1,length(subjects));
result.firstsubjNall    = zeros(length(bind)-1,length(subjects));
result.firstNNall       = zeros(length(bind)-1,length(subjects));
for ss = 1:length(subjects)
    indxS       = find(data.subject==subjects(ss));
    indxTfix    = find(data.type==1 & data.subject==subjects(ss) & data.start<0 & data.end>90);
%     indxSfix    = find(data.type==2 & data.subject==subjects(ss) & data.start<0 & data.end>0);
    if ~isempty(indxTfix)
        auxevidence      = (data.DToTarg(indxTfix)./posVec.pixxdeg);
        [N,BIN]         = histc(auxevidence,bind);
        indxTfix(BIN==0)= [];
        BIN(BIN==0)     = [];
        auxonT          = data.onTarg([indxTfix+2]);  
        result.firstsubjN(:,ss) = result.firstsubjN(:,ss)+accumarray(BIN',auxonT',[length(bind)-1,1]); 
        result.firstNN(:,ss)    = result.firstNN(:,ss)+accumarray(BIN',1,[length(bind)-1,1]); 
    end
    indxTfix(find(data.onTarg(indxTfix+2))) = [];       % remove the ones that actually end on target
    for symb = 1:48
        data.auxD   = sqrt((data.posx-posVec.scr(1,symb)).^2+(data.posy-posVec.scr(2,symb)).^2)./posVec.pixxdeg;
        if ~isempty(indxTfix)
           auxevidence      = data.auxD(indxTfix);
%            undist(symb) = mean(auxevidence)
            [N,BIN]         = histc(auxevidence,bind);
            indxTfix(BIN==0)= [];
            BIN(BIN==0)     = [];
            auxonSymb       = data.elfix(indxTfix+2)==symb;  % BEFORE IT WAS WITHOUT THE +2, RECHECK AGAIN
            result.firstsubjNall(:,ss) = result.firstsubjNall(:,ss)+accumarray(BIN',auxonSymb',[length(bind)-1,1]); 
            result.firstNNall(:,ss)    = result.firstNNall(:,ss)+accumarray(BIN',1,[length(bind)-1,1]); 
        end
    end
end
probto1 = result.firstsubjN./result.firstNN;
probtoall1 = result.firstsubjNall./result.firstNNall;
cmap1 = cbrewer('qual','Set2',3);
cmap2 = cbrewer('qual','Pastel2',3);
figure,hold on
plot(bind(1:end-1)+binsize/2,probto1,'-','Color',cmap2(1,:),'LineWidth',.5)
plot(bind(1:end-1)+binsize/2,probtoall1,'-','Color',cmap2(2,:),'LineWidth',.5)
plot(bind(1:end-1)+binsize/2,probto1-probtoall1,'-','Color',cmap2(3,:),'LineWidth',.5)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto1,2),nanstd(probto1,0,2)./sqrt(sum(~isnan(probto1),2)),'.-','Color',cmap1(1,:),'LineWidth',2)
errorbar(bind(1:end-1)+binsize/2,nanmean(probtoall1,2),nanstd(probtoall1,0,2)./sqrt(sum(~isnan(probtoall1),2)),'.-','Color',cmap1(2,:),'LineWidth',2)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto1-probtoall1,2),nanstd(probto1-probtoall1,0,2)./sqrt(sum(~isnan(probto1-probtoall1),2)),'.-','Color',cmap1(3,:),'LineWidth',2)
set(gca,'Xtick',bind)
xlabel('Distance to Target (degree)'),ylabel('p first mov to target')
xlim([0 bind(end-1)+binsize])
 doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['Ptotargfirstmov'],'600','painters',[],0)

% DO CONTROL
%%
% PROBABILITY TO MOVE TO TARGET IN N-MOVS BY DISTANCE
 binsize     = 1.5;
 bind        = [0:binsize:15,+inf];
nmovs       = 10;
result.subjN       = zeros(length(bind)-1,nmovs,length(subjects));
result.NN          = zeros(length(bind)-1,nmovs,length(subjects));
result.subjN2d     = zeros(length(bind)-1,length(bind)-1,length(subjects));
result.NN2d        = zeros(length(bind)-1,length(bind)-1,length(subjects));
for ss = 1:length(subjects)
    indxS   = find(data.subject==subjects(ss));
    for tt = unique(data.trial(indxS))
        indxTfix    = find(data.trial==tt & data.type==1 & data.subject==subjects(ss) &data.start>0);
        if ~isempty(indxTfix)
           auxevidence      = (data.DToTarg(indxTfix)./posVec.pixxdeg);
            [N,BIN]         = histc(auxevidence,bind);
            indxTfix(BIN==0)= [];
            BIN(BIN==0)     = [];
            auxonT          = data.onTarg([indxTfix]);
            if length(indxTfix)>3                                           %-1 fix          -2 fix
                result.subjN2d(:,:,ss) = result.subjN2d(:,:,ss)+accumarray([BIN(2:end-1)' BIN(1:end-2)'] ,auxonT(3:end)',[length(bind)-1,length(bind)-1]); 
                result.NN2d(:,:,ss)    = result.NN2d(:,:,ss)+accumarray([BIN(2:end-1)' BIN(1:end-2)'],1,[length(bind)-1,length(bind)-1]);
            end
            for nmov = 1:nmovs
                if length(BIN)>=nmov
                    result.subjN(:,nmov,ss) = result.subjN(:,nmov,ss)+accumarray(BIN(1:end-nmov)',auxonT(nmov+1:end)',[length(bind)-1,1]); 
                    result.NN(:,nmov,ss) = result.NN(:,nmov,ss)+accumarray(BIN(1:end-nmov)',1,[length(bind)-1,1]);
                end
            end

        end
    end
end

lineColors     = cmocean('thermal',nmovs); 
probto         = result.subjN./result.NN;
figure,
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
for rr = 1:size(probto,2)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto(:,rr,:),3),nanstd(probto(:,rr,:),0,3)./sqrt(sum(~isnan(probto(:,rr,:)),3)),'.-','LineWidth',1),hold  on
end
% plot(bind(1:end-1)+binsize/2,nanmean(probto,3),'.-'),xlabel('Distance to target (degree)'),ylabel('p #mov to target')
xlim([0 bind(end-1)+binsize])
colormap(lineColors), hc=colorbar;
hc.Limits = [0.5 nmovs+.5];
caxis([0.5 nmovs+.5])
hc.Ticks = [1 ceil(nmovs/2) nmovs];
hc.Label.String = '# mov';
xlabel('Distance to target (degree)'),ylabel('p(fix target in #mov)')
doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['PtotargperDistandOrdnotNorm'],'600','painters',[],0)
            
figure,
lineColors     = cmocean('thermal',length(bind)-1); 
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
for rr = 1:size(probto,1)
errorbar(1:nmovs,nanmean(probto(rr,:,:),3),nanstd(probto(rr,:,:),0,3)./sqrt(sum(~isnan(probto(rr,:,:)),3)),'.-','LineWidth',1),hold  on
end
xlabel('order'),ylabel('p(fix target)')
xlim([0.5 nmovs+.5])
colormap(lineColors), hc=colorbar;
hc.Limits = [0 bind(end-1)+binsize];
caxis([0 bind(end-1)+binsize])
hc.Ticks = [0 binsize*4 binsize*(length(bind)-3)]+binsize/2;
hc.Label.String = 'dist. to target (degree)';
doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['PtotargperOrd'],'600','painters',[],0)

           
figure,imagesc(nanmean(probto,3)')
axis xy 
xlabel('Distance to target (degree)'),ylabel('mov order pre target')
caxis([0 1]), colorbar     
 doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['PtotargperDistOrdim'],'600','painters',[],0)
           
probto2d = result.subjN2d./result.NN2d;
probto2dmar = sum(result.subjN2d,2)./sum(result.NN2d,2);


% plot(bind(1:end-1)+binsize/2,nanmean(probto2d,3),'.-'),xlabel('-1 Distance to Target (degree)'),ylabel('p next mov to target')
% xlim([0 bind(end-1)+binsize])
% colormap(lineColors), hc=colorbar;
% hc.Limits = [0 bind(end-1)+binsize];
% caxis([0 bind(end-1)+binsize])
% hc.Ticks = [0 bind(round(length(bind)/2)) bind(end-1)+binsize];
% hc.Label.String = '-2 dist. to target (degree)';
% 
%   doimage(gcf,fullfile(patheye,'figures'),...
%                  'tiff',['PnextMovtoTargbyDistsfirstandall'],'600','painters',[],1)
             
lineColors     = cmocean('thermal',size(probto2d,2)); 
figure,hold  on
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
for rr = 1:size(probto2d,1)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto2d(:,rr,:),3),nanstd(probto2d(:,rr,:),0,3)./sqrt(sum(~isnan(probto2d(:,rr,:)),3)),'.-','LineWidth',1),hold  on
end
errorbar(bind(1:end-1)+binsize/2,nanmean(probto2dmar,3),nanstd(probto2dmar,0,3)./sqrt(sum(~isnan(probto2dmar),3)),'k--','LineWidth',2),hold  on
xlabel('-1 Distance to Target (degree)'),ylabel('p next mov to target')
% plot(bind(1:end-1)+binsize/2,nanmean(probto2d,3),'.-'),
xlim([0 bind(end-1)+binsize])
colormap(lineColors), hc=colorbar;
hc.Limits = [0 bind(end-1)+binsize];
caxis([0 bind(end-1)+binsize])
hc.Ticks = [0 bind(round(length(bind)/2)) bind(end-1)+binsize];
hc.Label.String = '-2 dist. to target (degree)';

   doimage(gcf,fullfile(patheye,'figures'),...
                  'tiff',['PnextMovtoTargbyDists'],'600','painters',[],0)

%%
% PROBABILITY TO MOVE TO ANY ELEMENT IN N-MOVS BY DISTANCE
% binsize         = 1.5;
% bind            = [0:binsize:30,+inf];
% nmovs           = 11;
result.subjNall        = zeros(length(bind)-1,nmovs,length(subjects));
result.NNall           = zeros(length(bind)-1,nmovs,length(subjects));
result.subjN2dall      = zeros(length(bind)-1,length(bind)-1,length(subjects));
result.NN2dall         = zeros(length(bind)-1,length(bind)-1,length(subjects));

for ss = 1:length(subjects)
    ss
    indxS           = find(data.subject==subjects(ss));
    for symb = 1:48
        data.auxD   = sqrt((data.posx-posVec.scr(1,symb)).^2+(data.posy-posVec.scr(2,symb)).^2)./posVec.pixxdeg;
        indxSymb    = find(data.type==1 & data.subject==subjects(ss) & data.tpos ~= symb &data.start>0);
        utrials     = unique(data.trial(indxSymb));
        for tt = utrials
            indxTfix    = find(data.trial==tt & data.type==1 & data.subject==subjects(ss));
            if ~isempty(indxTfix)
                auxevidence = data.auxD(indxTfix);
                [N,BIN]     = histc(auxevidence,bind);
                indxTfix(BIN==0) = [];
                BIN(BIN==0) = [];
                auxonSymb   = data.elfix(indxTfix)==symb;
                if length(auxonSymb)>2
                    result.subjN2dall(:,:,ss) = result.subjN2dall(:,:,ss)+accumarray([BIN(2:end-1)' BIN(1:end-2)'] ,auxonSymb(3:end)',[length(bind)-1,length(bind)-1]); 
                    result.NN2dall(:,:,ss) = result.NN2dall(:,:,ss)+accumarray([BIN(2:end-1)' BIN(1:end-2)'],1,[length(bind)-1,length(bind)-1]);
                end
                for nmov = 1:nmovs
                    if length(BIN)>=nmov
                        result.subjNall(:,nmov,ss) = result.subjNall(:,nmov,ss)+accumarray(BIN(1:end-nmov)',auxonSymb(nmov+1:end)',[length(bind)-1,1]); 
                        result.NNall(:,nmov,ss) = result.NNall(:,nmov,ss)+accumarray(BIN(1:end-nmov)',1,[length(bind)-1,1]);
                    end
                end
            end
        end
    end 
  
end

lineColors     = cmocean('thermal',nmovs); 
probtoAll = result.subjNall./result.NNall;
save(fullfile(patheye,'eyedata','pbyorderdist.mat'),'result') 
figure,
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
for rr = 1:size(probto,2)
errorbar(bind(1:end-1)+binsize/2,nanmean(probtoAll(:,rr,:),3),nanstd(probtoAll(:,rr,:),0,3)./sqrt(sum(~isnan(probtoAll(:,rr,:)),3)),'.-','LineWidth',1),hold  on
end
xlim([0 bind(end-1)+binsize])
colormap(lineColors), hc=colorbar;
hc.Limits = [0.5 nmovs+.5];
caxis([0.5 nmovs+.5])
hc.Ticks = [1 ceil(nmovs/2) nmovs];
hc.Label.String = '# mov';
xlabel('Distance to element (degree)'),ylabel('p(fix element in #mov)')

  doimage(gcf,fullfile(patheye,'figures'),...
                 'tiff',['PtoElperDistOrd'],'600','painters',[],1)
            
figure,
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
plot(bind(1:end-1)+binsize/2,nanmean(probto-probtoAll,3),'.-'),xlabel('Distance to element (degree)'),ylabel('p #mov to element')
xlim([0 30])
colormap(lineColors), hc=colorbar;
hc.Limits = [0.5 nmovs+.5];
caxis([0.5 nmovs+.5])
hc.Ticks = [1 ceil(nmovs/2) nmovs];
hc.Label.String = '# fix';
%  doimage(gcf,fullfile(patheye,'figures'),...
%                 'tiff',['PtoTgtminElperDistOrd'],'600','painters',[],1)           
figure,
lineColors     = cmocean('thermal',length(bind-1)); 
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
plot(nanmean(probtoAll,3)','.-'),xlabel('order'),ylabel('p #mov to target')
xlim([0 nmovs+1])
colormap(lineColors), hc=colorbar;
hc.Limits = [binsize/2 length(bind)+.5];
caxis([0.5 bind(end-1)+1.5])
hc.Ticks = [1 bind(end-1)/2 bind(end-1) bind(end)];
hc.Label.String = 'dist. to target (degree)'

%  doimage(gcf,fullfile(patheye,'figures'),...
%                 'tiff',['PtoelperOrd'],'600','painters',[],1)


figure,
lineColors     = cmocean('thermal',length(bind-1)); 
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
plot(nanmean(probto-probtoAll,3)','.-'),xlabel('order'),ylabel('p #mov to target')
xlim([0 nmovs+1])
colormap(lineColors), hc=colorbar;
hc.Limits = [binsize/2 length(bind)+.5];
caxis([0.5 bind(end-1)+1.5])
hc.Ticks = [1 bind(end-1)/2 bind(end-1) bind(end)];
hc.Label.String = 'dist. to target (degree)'

%  doimage(gcf,fullfile(patheye,'figures'),...
%                 'tiff',['PtoTgtminElperOrd'],'600','painters',[],1)
            
probto2dall = result.subjN2dall./result.NN2dall;
probto2dallmar = sum(result.subjN2dall,2)./sum(result.NN2dall,2);

figure,hold  on
plot(bind(1:end-1)+binsize/2,probto1-probtoall1,'-','Color',cmap2(3,:),'LineWidth',.5)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto1-probtoall1,2),nanstd(probto1-probtoall1,0,2)./sqrt(sum(~isnan(probto1-probtoall1),2)),'.-','Color',cmap1(3,:),'LineWidth',2)
plot(bind(1:end-1)+binsize/2,squeeze(probto2dmar-probto2dallmar),':','Color',cmap2(3,:),'LineWidth',.5)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto2dmar-probto2dallmar,3),nanstd(probto2dmar,0,3)./sqrt(sum(~isnan(probto2dmar),3)),'.:','Color',cmap1(3,:),'LineWidth',2),hold  on


lineColors     = cmocean('thermal',size(probto2dall,2)); 
figure
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
plot(bind(1:end-1)+binsize/2,nanmean(probto2dall,3),'.-'),xlabel('-1 Distance to element (degree)'),ylabel('p next mov to element')
xlim([0 bind(end-1)+binsize])
colormap(lineColors), hc=colorbar;
hc.Limits = [0 bind(end-1)+binsize];
caxis([0 bind(end-1)+binsize])
hc.Ticks = [0 bind(round(length(bind)/2)) bind(end-1)+binsize];
hc.Label.String = '-2 dist. to element (degree)';
 

lineColors     = cmocean('thermal',size(probto2dall,2)); 
figure
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
for rr = 1:size(probto2d,1)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto2d(:,rr,:)-probto2dall(:,rr,:),3),nanstd(probto2d(:,rr,:)-probto2dall(:,rr,:),0,3)./sqrt(sum(~isnan(probto2d(:,rr,:)-probto2dall(:,rr,:)),3)),'.-','LineWidth',1),hold  on
end
% errorbar(bind(1:end-1)+binsize/2,nanmean(probto1-probtoall1,2),nanstd(probto1-probtoall1,0,2)./sqrt(sum(~isnan(probto1-probtoall1),2)),'k.--','LineWidth',1.5)
errorbar(bind(1:end-1)+binsize/2,nanmean(probto2dmar-probto2dallmar,3),nanstd(probto2dmar-probto2dallmar,0,3)./sqrt(sum(~isnan(probto2dmar-probto2dallmar),3)),'k.--','LineWidth',1.5)

% plot(bind(1:end-1)+binsize/2,nanmean(probto2d,3)-nanmean(probto2dall,3),'.-'),
xlabel('-1 Distance to Target (degree)'),ylabel('p next mov (to target - to element)')
xlim([0 bind(end-1)+binsize])
colormap(lineColors), hc=colorbar;
hc.Limits = [0 bind(end-1)+binsize];
caxis([0 bind(end-1)+binsize])
hc.Ticks = [0 bind(round(length(bind)/2)) bind(end-1)+binsize];
hc.Label.String = '-2 dist. to target (degree)';
  doimage(gcf,fullfile(patheye,'figures'),...
                 'tiff',['PnextMovtoTargbyDistsnrom'],'600','painters',[],1)


%%
% FIT pool
datafit = [bind(1:end-1)'+binsize/2 round(100*[mean(squeeze(probto2dmar)-squeeze(probto2dallmar),2) ones(size(probto2dmar,1),1)])];
datafit(:,2)=datafit(:,3)-datafit(:,2);
options.sigmoidName = 'norm';   % choose a cumulative Gauss as the sigmoid  
options.expType     = 'YesNo';
resultfit = psignifit(datafit,options);
 figure
 plotPsych(resultfit);
 fcn = @(b,x) normcdf(x, b(1), b(2));                    % Objective Function
NRCF = @(b) norm(datafit(:,2)/100 - fcn(b,datafit(:,1)));                     % Norm Residual Cost Function
B = fminsearch(NRCF, [0; 10])

    x         = linspace(min(resultfit.data(:,1)),max(resultfit.data(:,1)),100);
 fitValues = (1-resultfit.Fit(3)-resultfit.Fit(4))*arrayfun(@(x) resultfit.options.sigmoidHandle(x,resultfit.Fit(1),resultfit.Fit(2)),x)+resultfit.Fit(4);
figure,
plot(x,1-fitValues), hold on
%  plot(resultfit.conf_Intervals(1,:,1),repmat(resultfit.Fit(4)+resultfit.options.threshPC*(1-resultfit.Fit(3)-resultfit.Fit(4)),1,2))
 errorbar(bind(1:end-1)+binsize/2,nanmean(probto2dmar-probto2dallmar,3),nanstd(probto2dmar-probto2dallmar,0,3)./sqrt(sum(~isnan(probto2dmar-probto2dallmar),3)),'k.--','LineWidth',1.5)

 %%
 for ss = 1:size(result.NN2d,3)
    datafit = [bind(1:end-1)'+binsize/2 round(100*[squeeze(probto2dmar(:,:,ss))-squeeze(probto2dallmar(:,:,ss)) ones(size(probto2dmar,1),1)])];
    datafit(:,2)=datafit(:,3)-datafit(:,2);
    options.sigmoidName = 'norm';   % choose a cumulative Gauss as the sigmoid  
    options.expType     = 'YesNo';
    resultfitsuj(ss) = psignifit(datafit,options);
ss
 end
figure,plot(allfit(1,:),SUMMARY.rT,'.k')
box off

xlabel('fit 50% threshold')
ylabel('RT (s)')
resultfit = rmfield(resultfit,{'Posterior','weight'})
resultfitsuj = rmfield(resultfitsuj,{'Posterior','weight'})
save(fullfile(patheye,'eyedata','pbyorderfit.mat'),'resultfit','resultfitsuj') 

%%
auxdata                 = struct_select(data,{'type','value'},{'==1','>0'},2);
fh = figure;
% fh.Position = [30 312 1411 383];
cmap1 = cbrewer('qual','Set1',9);
for nback = 1:5;
    bla=[];
    N = [];
    thrsdur= prctile(auxdata.dur,[99.9]);
    xl = 0:.15:3;%[0 prctile(auxdata.evidence(nback,auxdata.evidence(nback,:)>0),10:10:100)]
    for suj = 1:length(subjects)
        for sp =1:length(xl)-1,
            auxindx     = find(auxdata.subject == subjects(suj) & auxdata.dur<thrsdur & auxdata.evidence(nback,:)>xl(sp) & auxdata.evidence(nback,:)<xl(sp+1));
            bla(sp,suj)  =nanmedian(auxdata.dur(auxindx));
            N(sp,suj)    =length(auxindx);
            a=a+1;
        end
    end
%     subplot(1,3,1),hold on
    h(nback) = plot(xl(1:end-1)+diff(xl)/2,nanmean(bla,2),'.-','Color',cmap1(nback,:)); hold on
    errorbar(xl(1:end-1)+diff(xl)/2,nanmean(bla,2),nanstd(bla,1,2)./sqrt(length(subjects)),'Color',cmap1(nback,:))
    xlabel('cumulative evidence'),ylabel('fix dur')
    ylim([100 250])
   % subplot(1,3,2),hold on
   % plot(nanmedian(bla,1),'.-','Color',cmap1(nback,:))
   % xlabel('preT order'),ylabel('fix dur')
   % ylim([100 250])
end
legend(h,{'same','1-back','2-','3-','4-'})
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTANCE TO TARGET DISTRIBUTION PER ORDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize = 1;
bind    = [0:binsize:30,+inf];
auxN    = [];
for ss = 1:length(subjects)
    indxS           = find(data.subject==subjects(ss));
    for oF = 1:15
        indxTfix        = find(data.type==1 & data.subject==subjects(ss) & data.orderPreT==oF);
        onT             = data.onTarg(indxTfix+2*oF);
        auxD            = data.DToTarg(indxTfix)./posVec.pixxdeg;
        auxN(:,oF,ss)   = histc(auxD(logical(onT)),bind)./sum(onT);
    end
end

figure,
lineColors     = cmocean('thermal',oF);
set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
plot(bind+binsize/2,mean(auxN,3),'.-')
set(gca,'XTick',3:3:30)
xlim([0 30])
xlabel('Distance to target (degree)'),ylabel('Relative Frequency')
colormap(lineColors), hc=colorbar;
hc.Limits = [0.5 oF+.5];
caxis([0.5 oF+.5])
hc.Ticks = [1 ceil(oF/2) oF];
hc.Label.String = '# fix'
 doimage(gcf,fullfile(patheye,'figures'),...
                'tiff',['PtotargperDistOrdnorm'],'600','painters',[],1)
           
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORDER DISTRIBUTION PER DISTANCE TO TARGET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize = 1;
bind    = [0:binsize:10,+inf];
auxN    = [];
for ss = 1:length(subjects)
    indxS           = find(data.subject==subjects(ss));
    for oD = 1:length(bind)-1
        indxTfix        = find(data.type==1 & data.subject==subjects(ss) & data.DToTarg./posVec.pixxdeg>bind(oD)  & data.DToTarg./posVec.pixxdeg<bind(oD+1) & data.tCorrect==1);
        auxOrder        = data.orderPreT(indxTfix)';
        auxOrder(auxOrder>20) = 20;
        auxN(:,oD,ss)   = accumarray(auxOrder+1,1,[21 1])./length(indxTfix);
    end
end

 figure,
 lineColors     = cmocean('thermal',oD);
 set(gca, 'ColorOrder',lineColors, 'NextPlot', 'replacechildren');
 plot(0:size(auxN,1)-1,nanmean(auxN,3),'.-')
 xlabel('Order'),ylabel('Relative Frequency')
 colormap(lineColors), hc=colorbar;
 hc.Limits = [-0.5 oD+1-0.5];
 caxis([-0.5 oD+1-0.5])
 hc.Label.String = 'dist to targ (degree)';
  doimage(gcf,fullfile(patheye,'figures'),...
                 'tiff',['distOrderperDist'],'600','painters',[],1)