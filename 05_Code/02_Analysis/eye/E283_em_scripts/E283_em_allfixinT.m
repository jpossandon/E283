%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple check that eyedata correspond with target positions and border
% between elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, hold on
set(gcf,'Position',[0 0 1280 700])
axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
plot(data.posx(data.type==1 & data.start>0),data.posy(data.type==1 & data.start>0),'.k','MarkerSize',3)
plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',16)
axis ij
vline(posVec.hLims)
hline(posVec.vLims,'r:')
tightfig
% doimage(gcf,fullfile(patheye,'figures'),...
%             'pdf',['allFixs_' namegr],[],1)
        
        %%
% this is for the different types of cue
cuevalues = [1 2 5 6;9 10 13 14];
infLabel  = {'Info','unInfo'};
for ii = 1:2


    for sp = [1,3,2,4]
            fh = figure;, hold on
set(fh,'Position',[20 20 posVec.stimRes/2])
        %subplot(2,2,sp), hold on
        axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
        indxtoplot = find(data.type==1 & data.value==cuevalues(ii,sp) & data.start>0);
        length(indxtoplot)
         indx50     = indxtoplot(randsample(length(indxtoplot),round(length(indxtoplot)/2)));
            length(indx50)
%         indx50     = indxtoplot;
        plot(data.posx(indx50),...
            data.posy(indx50),'.k','MarkerSize',1,'LineWidth',.1)
        for el = 1:size(posVec.scr,2)
            rectangle('Position',[posVec.scr(1,el)-posVec.pixxdeg*.75/2 posVec.scr(2,el)-posVec.pixxdeg*.75/2 posVec.pixxdeg*.75 posVec.pixxdeg*.75],'Curvature',[1 1],'EdgeColor',[0 0.447 0.741],'LineWidth',.4)
        end
%         plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',6,'LineWidth',.5)
        axis ij
       % vline(posVec.hLims)
       % hline(posVec.vLims,'r:')
%         title(SUMMARY.condLabels{sp+(ii-1)*4})
        axis([posVec.hLims(1)-85 posVec.hLims(end)+85 posVec.vLims(1)-70 posVec.vLims(end)+70])
        set(gca,'XTick',posVec.hLims,'XTickLabels',[-15,-11.3,-7.5,-3.7,0,3.7,7.5,11.3,15],...
            'YTick',posVec.vLims,'YTickLabels',round((posVec.vLims-posVec.stimRes(2)/2+1)./posVec.pixxdeg*10)/10,...
            'FontSize',6)
        figsize     = [4.6 4.6*fh.Position(4)/fh.Position(3)];
        axis off
        tightfig
         doimage(gcf,fullfile(patheye,'figures'),...
                     'pdf',['allFixs_' namegr '_' SUMMARY.condLabels{sp+(ii-1)*4}],'900','opengl',figsize,1)

    end
   % tightfig
   % doimage(gcf,fullfile(patheye,'figures'),...
    %                'pdf',['allFixs_' namegr '_' infLabel{ii}],[],1)
end

%%
%%
% 4 subplots for left/right uninfo/info
cuevalues = [1 6;
    2 5;
    9 13;
    10 14];
infLabel  = {'Info','unInfo'};

fh  = figure; hold on
set(fh,'Position',[20 20 posVec.stimRes])
auxdata = data;
auxdata.indx = zeros(1,length(data.start));
auxdata.indx(randsample(length(data.start),length(data.start)/2)) = 1;
auxdata = struct_select(auxdata,{'indx'},{'==1'},2);
spcords = [.06          .53;  
           .53 .53;
           .06          .06;
           .53 .06];
for sp = 1:4
    subplot('position',[spcords(sp,1),spcords(sp,2),.41,.41]), hold on
    axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
     plot(auxdata.posx(auxdata.type==1 & (auxdata.value==cuevalues(sp,1) | auxdata.value==cuevalues(sp,2)) & auxdata.start>0),...
        auxdata.posy(auxdata.type==1 & (auxdata.value==cuevalues(sp,1) | auxdata.value==cuevalues(sp,2)) & auxdata.start>0),'.k','MarkerSize',1,'LineWidth',.1)
    for el = 1:size(posVec.scr,2)
        rectangle('Position',[posVec.scr(1,el)-posVec.pixxdeg*.75/2 posVec.scr(2,el)-posVec.pixxdeg*.75/2 posVec.pixxdeg*.75 posVec.pixxdeg*.75],'Curvature',[1 1],'EdgeColor',[0 0.447 0.741],'LineWidth',.4)
    end
    %    plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',6,'LineWidth',.5)
    axis ij
    line([posVec.hLims' posVec.hLims']',[ones(9,1)*[0 posVec.stimRes(2)]]','LineStyle',':','LineWidth',.2,'Color',[0.85 0.325 0.098]);
    
    line([ones(7,1)*[0 posVec.stimRes(1)]]',[posVec.vLims' posVec.vLims']','LineStyle',':','LineWidth',.3,'Color',[0.85 0.325 0.098]);
    
%     hv = vline(posVec.hLims);
%     hh = hline(posVec.vLims,'r:');
    
    axis([posVec.hLims(1)-85 posVec.hLims(end)+85 posVec.vLims(1)-70 posVec.vLims(end)+70])
%     set(gca,'XTick',posVec.hLims,'XTickLabels',[-15,-11.3,-7.5,-3.7,0,3.7,7.5,11.3,15],...
%         'YTick',posVec.vLims,'YTickLabels',round((posVec.vLims-posVec.stimRes(2)/2+1)./posVec.pixxdeg*10)/10,...
%         'FontSize',6)
   axis off 
    
end
figsize     = [2*4.6 2*4.6*fh.Position(4)/fh.Position(3)];
    
%      doimage(gcf,fullfile(patheye,'figures'),...
%          'pdf',['allFixs_perextside'],figsize,1)
% tightfig
% doimage(gcf,fullfile(patheye,'figures'),...
%                'pdf',['allFixs_' namegr '_' infLabel{ii}],[],1)
