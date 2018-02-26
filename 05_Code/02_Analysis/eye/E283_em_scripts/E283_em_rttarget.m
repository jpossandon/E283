%% 
% pooled fixation duration according to 2 posiutions previous target 
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
% fixation duration wit the respect to the relationship between the two
% fixations previous target fixation in respect to their rleatives angles
% and distance to the target
fix1 = find(auxdat.type==1 & auxdat.orderPreT==1 & auxdat.event_order>1 );
 fix2 = find(auxdat.type==1 & auxdat.orderPreT==2);
%fix2 = fix1-2;
fix1dur = auxdat.dur(fix1);
fix2dur = auxdat.dur(fix2);
[theta1,rho1] = cart2pol(auxdat.xToTarg(fix1),auxdat.yToTarg(fix1));
[theta2,rho2] = cart2pol(auxdat.xToTarg(fix2),auxdat.yToTarg(fix2));
angledif = mod(theta1-theta2+pi,2*pi)-pi;
ampdif = rho1./rho2;
cmap    = flipud(cmocean('thermal'));
fh = figure,hold on
fh.Position = [2 164 1440 440];
clim = 400;
for ff = 1:length(fix1dur)
    auxdur = fix1dur(ff);
    if auxdur>clim,auxdur=clim;,end
    subplot(1,2,1),hold on
%     plot(angledif(ff),ampdif(ff),'.','Color',cmap(round(auxdur/clim*256),:),'MarkerSize',4)
    plot(rho1(ff),rho2(ff),'.','Color',cmap(round(auxdur/clim*256),:),'MarkerSize',4)
    subplot(1,2,2),hold on
     plot(rho1(ff),rho1(ff)./rho2(ff),'.','Color',cmap(round(auxdur/clim*256),:),'MarkerSize',4)
end
% axis([-pi-.1 pi+.1 0 3])
sp1 = subplot(1,2,1);
sp1.Position = [.05 .13 .42 .8];
hc = colorbar('colormap',cmap,'XTick',0:1/(clim/100):1,'XTickLabel',[0:100:clim]);
axis([0 1000 0 1000])
 xlabel('Distance to Target fix(-1)')
ylabel('Distance to Target fix(-2)')
sp2 = subplot(1,2,2);
sp2.Position = [.55 .13 .42 .8];
hc = colorbar('colormap',cmap,'XTick',0:1/(clim/100):1,'XTickLabel',[0:100:clim]);
axis([0 1000 0 5])
 xlabel('Distance to Target fix(-1)')
ylabel('Ratio Distance to Target fix(-1)/fix(-2)')
 doimage(gcf,fullfile(patheye,'figures'),...
             'tiff',['T1durperT1T2dist_' namegr],[],1)
%%

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

%%
M1      = mean(mmT1);
SE1     = std(mmT1)./sqrt(length(subjects));
M2      = mean(mmT2);
SE2     = std(mmT2)./sqrt(length(subjects));
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