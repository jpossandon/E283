%% 
% pooled fixation duration according to 2 positions previous target 
% fixation (within a distance of two elements)
xbins = -3*posVec.dsX+posVec.dsX/2:posVec.dsX:+3*posVec.dsX-posVec.dsX/2;
ybins = -3*posVec.dsY+posVec.dsY/2:posVec.dsY:+3*posVec.dsY-posVec.dsY/2;
nxb    = length(xbins)-1;
nyb    = length(ybins)-1;
auxdat =data; 
fix1 = auxdat.type==1 & auxdat.orderPreT==1& auxdat.event_order>1;
fixdur1 = auxdat.dur(fix1); 
[~,xb1]=histc(auxdat.xToTarg(fix1),xbins);
[~,yb1]=histc(auxdat.yToTarg(fix1),ybins);
fix2 = auxdat.type==1 & auxdat.orderPreT==2;
[~,xb2]=histc(auxdat.xToTarg(fix2),xbins);
[~,yb2]=histc(auxdat.yToTarg(fix2),ybins);
elimfix = xb1==0 | yb1==0 | xb2==0 | yb2==0;
xb1(elimfix) = [];
yb1(elimfix) = [];
xb2(elimfix) = [];
yb2(elimfix) = [];
fixdur1(elimfix) = [];

xx = nxb*(xb1-1)+xb2;
yy = nyb*(yb1-1)+yb2;
N =accumarray([xx',yy'],ones(length(xx),1),[nxb*nxb,nyb*nyb]);
M =accumarray([xx',yy'],fixdur1,[nxb*nxb,nyb*nyb],@mean);

fh = figure;
fh.Position = [4 408 1200 390];
subplot(1,2,1)
imagesc(.5:1:24.5,.5:1:24.5,N)
vline(5:5:20,'r')
hline(5:5:20,'r')
colorbar
title('N')
subplot(1,2,2)
imagesc(.5:1:24.5,.5:1:24.5,M)
vline(5:5:20,'r')
hline(5:5:20,'r')
colorbar
title('M')
caxis([100 300])
cmap = cmocean('thermal');
colormap(cmap)
doimage(gcf,fullfile(patheye,'figures'),...
              'tiff',['T1durperT2pos_' namegr],[],1)

%%
% DISTANCE TO TARGET VS TIME TO TARGET PER ORDER
bins = [0:25:1250,Inf];
cmap = cbrewer('qual','Set1',18);
figure,hold on
for opT = 1:10
    auxfix = auxdat.type==1 & auxdat.orderPreT==opT& auxdat.event_order>1 & auxdat.DToTarg<1408 & auxdat.tToTarg>-5000;
     plot(auxdat.DToTarg(auxfix)/45,-auxdat.tToTarg(auxfix),'.','Color',cmap(opT,:),'MarkerSize',6)
%     bihist = accumarray(ceil([-auxdat.tToTarg(auxfix)'/25,auxdat.DToTarg(auxfix)'/(posVec.pixxdeg/2)]),1,[200,62]);
% [n,c] = hist3(ceil([-auxdat.tToTarg(auxfix)'/25,auxdat.DToTarg(auxfix)'/(posVec.pixxdeg/2)]),{0:25:5000,0:.5:31})
% contour(c{1},c{2},n)
end
% axis([0 1500 0 10000])
xlabel('Distance to target')
ylabel('Time to target')
%%
% order pre targer per distance to target normalized per distancepooled
degbins     = 0:binsiz:30;
resDtotarg  = [];
auxdat     = struct_select(data,{'type'},{'==1'},2);
auxDtotarg = auxdat.DToTarg/posVec.pixxdeg;
[n,bins] = histc(auxDtotarg,degbins);
for dbi = 1:length(degbins)
    auxdist = auxdat.orderPreT(find(bins==dbi));
    auxdist(isnan(auxdist)) = [];
    resDtotarg(:,dbi) = accumarray(auxdist'+1,1,[63,1])./length(auxdist);
end
figure,imagesc(degbins,0:62,resDtotarg)
axis xy
 axis([0 30 0 30])
xlabel('Distance to target (degree)')
ylabel('Order pre Target')
% end
%%
% fixation duration with the respect to the relationship between the two
% fixations previous target fixation in respect to their distance to the target
binsiz = 1;
edgesrho = 0:binsiz:16;
for ss = 1:length(subjects)

    auxdat    = struct_select(data,{'subject'},{['==' num2str(subjects(ss))]},2);
    fix1 = find(auxdat.type==1 & auxdat.orderPreT==1 & auxdat.event_order>3 );
    fix2 = find(auxdat.type==1 & auxdat.orderPreT==2 & auxdat.event_order>2);
    fix3 = find(auxdat.type==1 & auxdat.orderPreT==3 & auxdat.event_order>1);
    fix4 = find(auxdat.type==1 & auxdat.orderPreT==4);
    fix1dur = auxdat.dur(fix1);
    fix2dur = auxdat.dur(fix2);
    fix3dur = auxdat.dur(fix3);
    fix4dur = auxdat.dur(fix4);

    fix1fix2dur = auxdat.end(fix1)-auxdat.start(fix2);
    fix1fix3dur = auxdat.end(fix1)-auxdat.start(fix3);
    [~,rho1] = cart2pol(auxdat.xToTarg(fix1),auxdat.yToTarg(fix1));
    [~,rho2] = cart2pol(auxdat.xToTarg(fix2),auxdat.yToTarg(fix2));
    [~,rho3] = cart2pol(auxdat.xToTarg(fix3),auxdat.yToTarg(fix3));
    [~,rho4] = cart2pol(auxdat.xToTarg(fix4),auxdat.yToTarg(fix4));

    [N1,BIN1]   = histc(rho1./posVec.pixxdeg,edgesrho);
    [N2,BIN2]   = histc(rho2./posVec.pixxdeg,edgesrho);
    [N3,BIN3]   = histc(rho3./posVec.pixxdeg,edgesrho);
    [N4,BIN4]   = histc(rho4./posVec.pixxdeg,edgesrho);
    indxuse1    =  find(BIN1>0);
    indxuse2    = find(BIN2>0);
    indxuse3    = find(BIN3>0);
    indxuse4    = find(BIN4>0);
    indxuse12   =  find(BIN1>0 & BIN2>0);
    indxuse23   =  find(BIN2>0 & BIN3>0);
    indxuse123  =  find(BIN1>0 & BIN2>0 & BIN3>0);
    dur1dist1(:,ss)  = accumarray(BIN1(indxuse1)',fix1dur(indxuse1),[length(N1)-1 1],@mean);
    dur1dist2(:,ss)  = accumarray(BIN2(indxuse2)',fix1dur(indxuse2),[length(N1)-1 1],@mean);
    dur1dist3(:,ss)  = accumarray(BIN3(indxuse3)',fix1dur(indxuse3),[length(N1)-1 1],@mean);
    dur1dist4(:,ss)  = accumarray(BIN4(indxuse4)',fix1dur(indxuse4),[length(N1)-1 1],@mean);
    dur12dist2(:,ss)  = accumarray(BIN2(indxuse2)',fix1fix2dur(indxuse2),[length(N1)-1 1],@mean);
    dur1dist12(:,:,ss)  = accumarray([BIN1(indxuse12)',BIN2(indxuse12)'],fix1dur(indxuse12),[length(N1)-1 length(N2)-1],@mean);
    dur2dist12(:,:,ss)  = accumarray([BIN1(indxuse12)',BIN2(indxuse12)'],fix2dur(indxuse12),[length(N1)-1 length(N2)-1],@mean);
    
    dur12dist12(:,:,ss)  = accumarray([BIN1(indxuse12)',BIN2(indxuse12)'],fix1fix2dur(indxuse12),[length(N1)-1 length(N2)-1],@mean);
    Ndist12(:,:,ss)  = accumarray([BIN1(indxuse12)',BIN2(indxuse12)'],1,[length(N1)-1 length(N2)-1]);
    Ndist23(:,:,ss)  = accumarray([BIN2(indxuse23)',BIN3(indxuse23)'],1,[length(N2)-1 length(N3)-1]);
    
    
    dur2dist2(:,ss)  = accumarray(BIN2(indxuse2)',fix2dur(indxuse2),[length(N1)-1 1],@mean);
    dur2dist3(:,ss)  = accumarray(BIN3(indxuse3)',fix2dur(indxuse3),[length(N1)-1 1],@mean);
    dur3dist3(:,ss)  = accumarray(BIN3(indxuse3)',fix1dur(indxuse3),[length(N1)-1 1],@mean);
end
dur1dist1(dur1dist1==0) = NaN;
dur1dist2(dur1dist2==0) = NaN;
dur1dist3(dur1dist3==0) = NaN;
dur2dist2(dur2dist2==0) = NaN;
dur3dist3(dur3dist3==0) = NaN;
dur1dist4(dur1dist4==0) = NaN;
dur1dist12( dur1dist12==0)= NaN;
dur2dist12( dur2dist12==0)= NaN;
dur12dist12( dur12dist12==0)= NaN;

% FIXATION DURATION PER DISTANCE TO TARGET AT DIFFERENT PRETARGET ORDER
  figure, clear h, hold on
h(1) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur1dist1,2),nanstd(dur1dist1,0,2)./sqrt(sum(~isnan(dur1dist1),2)));
h(2) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur1dist2,2),nanstd(dur1dist2,0,2)./sqrt(sum(~isnan(dur1dist2),2)));
h(3) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur1dist3,2),nanstd(dur1dist3,0,2)./sqrt(sum(~isnan(dur1dist3),2)));
h(4) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur1dist4,2),nanstd(dur1dist4,0,2)./sqrt(sum(~isnan(dur1dist4),2)));
h(5) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur12dist2,2),nanstd(dur12dist2,0,2)./sqrt(sum(~isnan(dur12dist2),2)));
h(6) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur2dist2,2),nanstd(dur2dist2,0,2)./sqrt(sum(~isnan(dur2dist2),2)));
h(7) = errorbar(edgesrho(1:end-1)+binsiz/2,nanmean(dur2dist3,2),nanstd(dur2dist3,0,2)./sqrt(sum(~isnan(dur2dist3),2)));
xlabel('Distance to target')
ylabel('Fixation duration')
ll =legend(h,{'dist -1 dur -1','dist -2 dur -1','dist -3 dur -1','dist -4 dur -1','dist -2 dur -1+2','dist -2 dur -2','dist -3 dur -2'},'Location','East')
ll.Position = [0.6 0.5 0.25 0.2];
% doimage(gcf,fullfile(patheye,'figures'),...
%               'tiff',['fixdurperTdist_' namegr],'600','painters',[],1)

% 2D COLORPLOT OF FIXATION DURATION AGAINST DISTANCE TO TARGET -1 AND -2
figure,imagesc(edgesrho(1:end-1)+binsiz/2,edgesrho(1:end-1)+binsiz/2,nanmean(dur1dist12,3))
axis xy

ylabel('-1 Distance to target (degree)')
xlabel('-2 Distance to target (degree)')
caxis([0 350])
colormap jet,colorbar
title('Duration 1')
doimage(gcf,fullfile(patheye,'figures'),...
           'tiff',['Dur1perDtoTatg12'],'600','painters',[],1)

% 2D COLORPLOT OF FIXATION DURATION -2 AGAINST DISTANCE TO TARGET -1 AND -2
figure,imagesc(edgesrho(1:end-1)+binsiz/2,edgesrho(1:end-1)+binsiz/2,nanmean(dur2dist12,3))
axis xy
ylabel('-1 Distance to target (degree)')
xlabel('-2 Distance to target (degree)')
caxis([0 350])
colormap jet,colorbar
title('Duration 2')
doimage(gcf,fullfile(patheye,'figures'),...
       'tiff',['Dur2perDtoTatg12'],'600','painters',[],1)
   
% 2D COLORPLOT OF FIXATION DURATION -1 + -2 AGAINST DISTANCE TO TARGET -1 AND -2
figure,imagesc(edgesrho(1:end-1)+binsiz/2,edgesrho(1:end-1)+binsiz/2,nanmean(dur12dist12,3))
axis xy
ylabel('-1 Distance to target (degree)')
xlabel('-2 Distance to target (degree)')
caxis([0 600])
colormap jet,colorbar
title('Duration 1+2')
doimage(gcf,fullfile(patheye,'figures'),...
               'tiff',['Dur12perDtoTatg12'],'600','painters',[],1)

% 2D COLORPLOT OF AVERAGE N PER DISTANCE TO TARGET -1 AND -2

figure,imagesc(edgesrho(1:end-1)+binsiz/2,edgesrho(1:end-1)+binsiz/2,nanmean(Ndist12,3))
axis xy
ylabel('-1 Distance to target (degree)')
xlabel('-2 Distance to target (degree)')
colormap jet,colorbar
title('N')
doimage(gcf,fullfile(patheye,'figures'),...
   'tiff',['NperDtoTatg12'],'600','painters',[],1)

%%
% all fixation duration in color scale agains distance to target -1 and -2
auxdat    = struct_select(data,{'type'},{'==1'},2);
fix1      = find(auxdat.type==1 & auxdat.orderPreT==1 & auxdat.event_order>2 );
fix2      = find(auxdat.type==1 & auxdat.orderPreT==2 & auxdat.event_order>1 );

fix1dur = auxdat.dur(fix1);
[~,rho1] = cart2pol(auxdat.xToTarg(fix1),auxdat.yToTarg(fix1));
[~,rho2] = cart2pol(auxdat.xToTarg(fix2),auxdat.yToTarg(fix2));
    
cmap    = flipud(cmocean('thermal'));
fh      = figure;,hold on
fh.Position = [2 164 1440 440];
clim = 400;
for ff = 1:length(fix1dur)
    auxdur = fix1dur(ff);
    if auxdur>clim,auxdur=clim;,end
    subplot(1,2,1),hold on
    plot(rho1(ff)./posVec.pixxdeg,rho2(ff)./posVec.pixxdeg,'.','Color',cmap(round(auxdur/clim*256),:),'MarkerSize',4)
    subplot(1,2,2),hold on
     plot(rho1(ff)./posVec.pixxdeg,rho1(ff)./rho2(ff),'.','Color',cmap(round(auxdur/clim*256),:),'MarkerSize',4)
end
% axis([-pi-.1 pi+.1 0 3])
sp1 = subplot(1,2,1);
sp1.Position = [.05 .13 .42 .8];
hc = colorbar('colormap',cmap,'XTick',0:1/(clim/100):1,'XTickLabel',[0:100:clim]);
axis([0 15 0 30])
 xlabel('Distance to Target fix(-1)')
ylabel('Distance to Target fix(-2)')
sp2 = subplot(1,2,2);
sp2.Position = [.55 .13 .42 .8];
hc = colorbar('colormap',cmap,'XTick',0:1/(clim/100):1,'XTickLabel',[0:100:clim]);
axis([0 15 0 5])
 xlabel('Distance to Target fix(-1)')
ylabel('Ratio Distance to Target fix(-1)/fix(-2)')
 doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['T1durperT1T2dist_' namegr],'600','painters',[],1)
%%
% fixation duration t1 and t2 according to position in the grid with
% respcet to the target
labels = {'T2onT_T1onT','T2nextT_T1onT','T2far_T1onT',...
        'T2onT_T1nextT','T2nextT_T1nextT','T2far_T1nextT',...
        'T2onT_T1far','T2nextT_T1far','T2far_T1far'};
codeonT     = [1 0 0 1 0 0 1 0 0;...
               1 1 1 0 0 0 0 0 0];
codenextT   = [0 1 0 0 1 0 0 1 0;...
               0 0 0 1 1 1 0 0 0];
for ss = 1:length(subjects)

    auxdat    = struct_select(data,{'subject'},{['==' num2str(subjects(ss))]},2);
    auxdat.nextToTarg = any([auxdat.nextToTargH;auxdat.nextToTargV;auxdat.nextToTargD]);
    
    for cc = 1:length(labels)
    auxindx         = find(auxdat.type(1:end-2)==1 & auxdat.orderPreT(1:end-2)==2 & auxdat.onTarg(1:end-2) == codeonT(1,cc) & auxdat.nextToTarg(1:end-2) == codenextT(1,cc)  &...
                        auxdat.type(3:end)==1 & auxdat.orderPreT(3:end)==1 & auxdat.onTarg(3:end) == codeonT(2,cc) & auxdat.nextToTarg(3:end) == codenextT(2,cc));
    nn(ss,cc)       = length(auxindx);
    mmT2(ss,cc)      = mean(auxdat.dur(auxindx));
    mmT1(ss,cc)      = mean(auxdat.dur(auxindx+2));
    end
   
end

%
M1      = nanmean(mmT1);
SE1     = nanstd(mmT1)./sqrt(length(subjects));
M2      = nanmean(mmT2);
SE2     = nanstd(mmT2)./sqrt(length(subjects));
xpos    =[1:3,5:7,9:11];

fh=figure;
hold on
hh(1) = plot(xpos,M1,'sk','MarkerFaceColor',[1 0 0]);
errorbar(xpos,M1,SE1,'sk')
hh(2) = plot(xpos,M2,'sk','MarkerFaceColor',[0 1 0]);
errorbar(xpos,M2,SE2,'sk')
text([2,6,10],[260,260,260],{'T1ont','T1next','T1far'},'HorizontalAlignment','center')
text(xpos,250*ones(1,9),arrayfun(@num2str,round(mean(nn)), 'UniformOutput', false),'HorizontalAlignment','center')
ylim([100 300])
set(gca,'XTick',xpos,'XTickLabel',{'T2onT','T2nextT','T2far'})
ylabel('Fixation duration (ms)')    
legend(hh,{'T1','T2'})
 doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['T1T2durperT1T2dist_' namegr],[],1)