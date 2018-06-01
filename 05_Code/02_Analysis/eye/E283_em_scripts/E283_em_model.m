%%
% modelling

%%
% pooled transition matrices and target missesper subject
tmax_bypos = zeros(6,8,48,8,length(subjects));
tmax_cent  = nan(11,15,48,8,length(subjects));
conditions = [1 2 5 6 9 10 13 14];
for ss = 1:length(subjects)
    for cc = 1:length(conditions)
        auxdata    = struct_select(data,{'subject','value'},{['==' num2str(subjects(ss))],['==' num2str(conditions(cc))]},2);

        for el = 1:48
            % fixations after start
            auxindx                 = find(auxdata.type==1 & auxdata.elfix==el & auxdata.orderPreT>0 & auxdata.start>0);
            auxindx(find(auxdata.type(auxindx+2)~=1))=[];
            auxnext                 = auxdata.elfix(auxindx+2);
            auxnext(isnan(auxnext)) = [];
            tmax_bypos(:,:,el,cc,ss)      = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
            [Iel,Jel]               = ind2sub([6 8],el);
            tmax_cent([6:6+6-1]-Iel+1,[8:8+8-1]-Jel+1,el,cc,ss) = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
        end
        auxindx                 = find(auxdata.type==1 & auxdata.orderPreT>0 & auxdata.start<0 & auxdata.end>0);
        auxindx(find(auxdata.type(auxindx+2)~=1))=[];    
        auxnext                 = auxdata.elfix(auxindx+2);
        auxnext(isnan(auxnext)) = [];
        tmax_first(:,:,cc,ss)      = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
    end
       % target miss 
    auxdata    = struct_select(data,{'subject'},{['==' num2str(subjects(ss))]},2);
    orderOnTarg = auxdata.orderPreT(logical(auxdata.onTarg));
    orderOnTarg(isnan(orderOnTarg)) = [];
    orderOnTarg(orderOnTarg>5) =5;
    ordOnTarg(ss,:) = accumarray(orderOnTarg'+1,1,[6 1])./length(orderOnTarg);
   
   % fixation and saccade duration
   fixdurs{ss}          =  auxdata.dur(find(auxdata.type==1 & auxdata.orderPreT>0 & auxdata.start>0 & ~auxdata.revisit & ~auxdata.refix));
   fixdursrefix{ss}     =  auxdata.dur(find(auxdata.type==1 & auxdata.orderPreT>0 & auxdata.start>0 & auxdata.refix));
   fixdursrevisit{ss}   =  auxdata.dur(find(auxdata.type==1 & auxdata.orderPreT>0 & auxdata.start>0 & auxdata.revisit));
   sacdurs{ss}          =  auxdata.dur(find(auxdata.type==2 & auxdata.start>0));
   latfmov{ss}          = auxdata.end(find(auxdata.type==1 & auxdata.orderPreT>0 & auxdata.start<0 & auxdata.end>0));
   
   % nrefix, nrevisit, lags
   nrefix(ss)           = sum(auxdata.refix)./length(unique(auxdata.trial));
   nrevisits(ss)        = sum(auxdata.revisit>0)./length(unique(auxdata.trial));
   nnewfixs(ss)            = sum(~auxdata.refix & auxdata.revisit==0)./length(unique(auxdata.trial));
   nlags{ss}            = auxdata.revisit(auxdata.revisit>0);
end
%%
reps    = 1;
trans   = 0;
prev    = 1;
tlimit  = 1;
trefix  = 1;
for rep= 1:reps
    for ss = 1:length(subjects)
        fprintf('\nRepetition %d/%d Subject %d',rep,reps,ss)
        auxdata    = struct_select(data,{'subject'},{['==' num2str(subjects(ss))]},2);
        lagsTP = [];
        for tr=1:592
            fprintf('.');
            tgtIndx = unique(auxdata.tpos(auxdata.trial==tr));
            lags    = [];
            fI      = 1;
            trDur   = 0;
            findx   = [];
            refix   = 0;
            revisit    = 0;
            nextFixOnt = 0;
            while fI>0
                if fI ==1                                                       % first movement
                    trDur = trDur + double(randsample(latfmov{ss},1));          % latency first fixation
                    if trans
                        findx(fI) = randsample(1:48,1,true,reshape(tmax_first(:,:,ss),[1, 48]));                                 % first fixation
                    else
                        findx(fI) = randsample(1:48,1);
                    end       
                else
                    nhF = 1;
                    if nextFixOnt
                        findx(fI) = tgtIndx;
                        nhF = 0;
                    end
                    if nhF
                        if fI==2
                            if trans
                                findx(fI) = randsample(1:48,1,true,reshape(tmax_bypos(:,:,findx(fI-1),ss),[1, 48]));                                 % first fixation
                            else
                                findx(fI) = randsample(1:48,1);
                            end  
                        else
                            if trans
                                findx(fI) = randsample(1:48,1,true,reshape(tmax_bypos(:,:,findx(fI-1),ss),[1, 48]));                                 % first fixation
                            else
                                findx(fI) = randsample(1:48,1);
                            end
                        end
                    end
                    if findx(end)==findx(end-1)
                        refix = refix+1;
                    elseif ismember(findx(end),findx(1:end-2))
                        revisit = revisit+1;
                        lags = [lags, length(findx)-find(findx(1:end-1)==findx(end),1,'last')];
                    end
                end
                if findx(fI) == tgtIndx
                    % if random number is above the miss probability the
                    % target is found, refix and 
                    randT = rand(1);
                    if randT<ordOnTarg(ss,1)
                        if fI>1
                            distfI = sqrt(sum((posVec.scr(:,findx(fI))-...
                                    posVec.scr(:,findx(fI-1))).^2))./posVec.pixxdeg;       % distance between elements
%                             trDur  = trDur + round(10^(sacfit(1) + sacfit(2)*log10(distfI)));
                        end
                        sruns(ss,rep).hit(tr)  = 1;
                        sruns(ss,rep).rT(tr)   = double(trDur)/1000;
                        sruns(ss,rep).nFix(tr) = fI-1;
                        c1 = 0; %why?
                        if fI>2
                            if findx(end)==findx(end-1)
                                refix = refix-1;
                                c1    = c1+1;
                            end
                            if ismember(findx(end),findx(1:end-2))
                                revisit = revisit-1;
                                c1   = c1+1;
                            end
                        end
                        sruns(ss,rep).nrevisit(tr)  = revisit;
                        sruns(ss,rep).nrefix(tr) = refix;
                        sruns(ss,rep).nnewFix(tr)  = fI-revisit-refix+c1-1;    % doublecheked this is correct
                        lagsTP = [lagsTP,lags]; 
                        fI = 0;             % found!!
                    elseif randT>ordOnTarg(ss,1) & randT<sum(ordOnTarg(ss,1:2))
                        nextFixOnt = 1;
                    else
                        nextFixOnt = 0;
                    end
                end
                if fI>1 % here we add fixation duration (according to fixation type) and saccade duration time (according to subject main sequence)
                    distfI = sqrt(sum((posVec.scr(:,findx(fI))-...
                                    posVec.scr(:,findx(fI-1))).^2))./posVec.pixxdeg;       % distance between elements
%                         trDur  = trDur + round(10^(sacfit(1) + sacfit(2)*log10(distfI)));                     % expected saccade duration
                    if ismember(findx(end),findx(1:end-1))
                        if findx(end)-findx(end-1) == 0             % refixation
                                 trDur = trDur + double(randsample(fixdursrefix{ss},1));
                        else                                     % revisit
                                 trDur = trDur + double(randsample(fixdursrevisit{ss},1));
                        end
                    else               % new fixation
                             trDur = trDur + double(randsample(fixdurs{ss},1));
                    end
                end
                if fI >0
                    if tlimit ==1
                        dyndur = 12;
                    else
                        dyndur = 30;
                    end
                    if trDur/1000>dyndur
                        sruns(ss,rep).hit(tr)      = 0;
                        sruns(ss,rep).rT(tr)       = NaN;
                        sruns(ss,rep).nFix(tr)     = fI;
                        sruns(ss,rep).nrevisit(tr)  = revisit;
                        sruns(ss,rep).nrefix(tr) = refix;
                        sruns(ss,rep).nnewFix(tr)  = fI;
%                         lagsmiss = [lagsmiss,lags];
                        fI = 0;
                    else
                        fI = fI+1;
%                         [i j] = ind2sub([6 6],findx(fI-1));
                    end
                end
            end
        end
         sruns(ss,rep).lagsTP     = lagsTP;
           
    end
end
fprintf('\n');
%%
 for ss = 1:length(subjects)
     nrefix_sim(ss)   = sum(sruns(ss).nrefix)./592;
     nrevisit_sim(ss) = sum(sruns(ss).nrevisit)./592;
     nnewfixs_sim(ss) = sum(sruns(ss).nnewFix)./592;
 end