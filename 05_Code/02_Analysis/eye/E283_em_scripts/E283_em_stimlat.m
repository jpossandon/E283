s=1;
nsacs = 5;
load(fullfile(patheye,'eyedata','allsampledata.mat')) 
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
  figure
        set(gcf,'Position',[30 53 1251 652])
        N = zeros(length(allgazes),8,3);
for ss = 1:length(allgazes)  % loop through subjects
    for gOr = 1:2            % loop thoough saccade order after stim
      
        
        for ev = 1:length(allgazes(ss).gaze)  % loop thorough saccades
            if allgazes(ss).order(ev)==gOr
                wsP = find(cV==allgazes(ss).cond(ev));  % conditiom
                subplot(2,4,subpltmap(wsP)),hold on
                gazetoplot      = allgazes(ss).gazesac{ev};
                if ~isempty(gazetoplot) & allgazes(ss).latposStim(ev)<600
                    gazetoplot      = (gazetoplot-repmat(gazetoplot(:,1),1,size(gazetoplot,2)))/posVec.pixxdeg;
                    gazetoplot(2,:) = gazetoplot(2,:)+allgazes(ss).latposStim(ev);
                    if any(gazetoplot(1,:)<0) & any(gazetoplot(1,:)>0)  % reversal from baseline stast
                        if gazetoplot(1,end)<0
                        plot(gazetoplot(1,:),gazetoplot(2,:),'-','Color',[1 0 0],'MarkerSize',1,'LineWidth',.5)
                        N(ss,wsP,1) = N(ss,wsP,1)+1;
                        end
                        if gazetoplot(1,end)>0
                            plot(gazetoplot(1,:),gazetoplot(2,:),'-','Color',[0 0 1],'MarkerSize',1,'LineWidth',.5)
                            N(ss,wsP,2) = N(ss,wsP,2)+1;
                        end
                        
                        
                    else
                        N(ss,wsP,3) = N(ss,wsP,3)+1; 
%                         plot(gazetoplot(1,:),gazetoplot(2,:),'-','Color',[.7 .7 .7],'MarkerSize',1,'LineWidth',.1)
                    end
                end
            end
        end
    
        for wsP=1:8
            subplot(2,4,subpltmap(wsP))
            axis([-10 10 0 600])
            vline(0)
            hline(0)
            title(SUMMARY.condLabels(wsP))
        end
%          doimage(gcf,fullfile(patheye,'figures','sacpersubj'),...
%                    'tiff',sprintf('s%d_postStim%d_%s',subjects(ss),gOr,namegr),[],1)
    end
end