%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pooled transition matrices per pair location cell grid/ target per
% condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
tmax_bypos = zeros(6,8,48,8);
tmax_bypos2mov = zeros(6,8,48,8);
tmax_cent  = nan(11,15,48,8);
conditions = [1 2 5 6 9 10 13 14];

for cc = 1:length(conditions)
    auxdata    = struct_select(data,{'value'},{['==' num2str(conditions(cc))]},2);

    for el = 1:48
        % fixations after start
        auxindx                 = find(auxdata.type==1 & auxdata.elfix==el & auxdata.orderPreT>0 & auxdata.start>0);
        auxindx(find(auxdata.type(auxindx+2)~=1))=[];
        auxnext                 = auxdata.elfix(auxindx+2);
        auxindx                 = find(auxdata.type==1 & auxdata.elfix==el & auxdata.orderPreT>1 & auxdata.start>0);
        auxnext2                = auxdata.elfix(auxindx+4);
        auxnext(isnan(auxnext)) = [];
        auxnext2(isnan(auxnext2)) = [];
        tmax_bypos(:,:,el,cc)      = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
        tmax_bypos2mov(:,:,el,cc)  = reshape(accumarray(auxnext2',1,[48,1]),[6 8])./length(auxnext2);
        [Iel,Jel]               = ind2sub([6 8],el);
        tmax_cent([6:6+6-1]-Iel+1,[8:8+8-1]-Jel+1,el,cc) = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
    end
    auxindx                 = find(auxdata.type==1 & auxdata.orderPreT>0 & auxdata.start<0 & auxdata.end>0);
    auxindx(find(auxdata.type(auxindx+2)~=1))=[];    
    auxnext                 = auxdata.elfix(auxindx+2);
    auxnext(isnan(auxnext)) = [];
    tmax_first(:,:,cc)      = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
end
       % target miss 
    auxdata                         = data;
    orderOnTarg                     = auxdata.orderPreT(logical(auxdata.onTarg));
    orderOnTarg(isnan(orderOnTarg)) = [];
    orderOnTarg(orderOnTarg>5)      = 5;
    ordOnTarg                       = accumarray(orderOnTarg'+1,1,[6 1])./length(orderOnTarg);
   resultTM.tmax_cent = tmax_cent;
    %    % nrefix, nrevisit, lags
%    nrefix(ss)           = sum(auxdata.refix)./length(unique(auxdata.trial));
%    nrevisits(ss)        = sum(auxdata.revisit>0)./length(unique(auxdata.trial));
%    nnewfixs(ss)            = sum(~auxdata.refix & auxdata.revisit==0)./length(unique(auxdata.trial));
%    nlags{ss}            = auxdata.revisit(auxdata.revisit>0);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% probability to fixate target 1-2-3 back by pos
% pooled data analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auxdata                 = struct_select(data,{'type','value'},{'==1','>0'},2);
auxdata.next1elfix      = [auxdata.elfix(2:end) NaN];
auxdata.next1elfix(auxdata.orderPreT==0) = NaN;
auxdata.next2elfix      = [auxdata.elfix(3:end) NaN];
auxdata.next2elfix(auxdata.orderPreT==0 | auxdata.orderPreT==1) = NaN;
auxdata.next3elfix      = [auxdata.elfix(4:end) NaN];
auxdata.next3elfix(auxdata.orderPreT==0 | auxdata.orderPreT==1 | auxdata.orderPreT==2) = NaN;

% one back target and not in target
clear excess ptotarg expptotarg
nbs = {'1bnT','2b1bnT','3b2b1bnT'};
for nb = 1:length(nbs)
    tmax_cent_totgt         = zeros(11,15);
    tmax_cent_nottotgt      = zeros(11,15);
    tmax_exp                = zeros(11,15);
    tmax_exp_fix1b          = zeros(11,15);
    tmax_exp_fix2b          = zeros(11,15);
    tmax_exp_2mov           = zeros(11,15);
    % selecting data according to type of n-back analysis
    if strcmp(nbs{nb},'1bnT')
        auxindx             = find(auxdata.orderPreT==1);
    elseif strcmp(nbs{nb},'2b1bnT')
        auxindx             = find(auxdata.orderPreT==2);
        %auxindx(find(auxdata.onTarg(auxindx+1))) = [];                      % do not take in account 2-back that have 1-back in target
    elseif strcmp(nbs{nb},'3b2b1bnT')
        auxindx             = find(auxdata.orderPreT==3);
       % auxindx(find(auxdata.onTarg(auxindx+1) | auxdata.onTarg(auxindx+2))) = [];  
    end
    auxel           = auxdata.elfix(auxindx);
    auxtgt          = auxdata.tpos(auxindx);
    auxcond         = auxdata.value(auxindx);
    
    mtx             = [auxel',auxtgt',auxcond'];                        % which element is two back, which elements is the trial target, wha tcondition
    
    if strcmp(nbs{nb},'1bnT')
        auxall          = [auxel',auxtgt',auxcond'];
    elseif strcmp(nbs{nb},'2b1bnT')
        auxel1b         = auxdata.elfix(auxindx+1);
        auxall          = [auxel',auxtgt',auxcond',auxel1b'];
    elseif strcmp(nbs{nb},'3b2b1bnT')
        auxel2b         = auxdata.elfix(auxindx+1);
        auxel1b         = auxdata.elfix(auxindx+2);
        auxall          = [auxel',auxtgt',auxcond',auxel1b',auxel2b'];
    end
    % remove NaNs
    auxall(any(isnan(auxall),2),:) = [];
    [~,IA]  = unique(auxall(:,1:3),'rows');
    auxall  = auxall(IA,:);
    length(auxall)
    % loop trough all empirical n-back position to target position in every given condition
    for aa = 1:size(auxall,1)
   
        % find starting and target position in transition matrix
        [Iel,Jel]    = ind2sub([6 8],auxall(aa,1));
        [Itgt,Jtgt]  = ind2sub([6 8],auxall(aa,2));
%         if Iel==Itgt & Jel==Jtgt
%             aa
%         end
        % # of movement from a given element position to a given target (in a given
        % condition)
        if strcmp(nbs{nb},'1bnT')
            totgt = sum(auxdata.orderPreT==1  &  auxdata.elfix==auxall(aa,1) & auxdata.next1elfix==auxall(aa,2) & auxdata.value==auxall(aa,3));
            tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + totgt;
        elseif strcmp(nbs{nb},'2b1bnT')
              % this is not considering the specific path
             totgt = sum(auxdata.orderPreT==2  & auxdata.elfix==auxall(aa,1) & auxdata.next2elfix==auxall(aa,2) & auxdata.value==auxall(aa,3));
              tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + totgt;
        elseif strcmp(nbs{nb},'3b2b1bnT')
             totgt = sum(auxdata.orderPreT==3  & auxdata.elfix==auxall(aa,1) & auxdata.next3elfix==auxall(aa,2) & auxdata.value==auxall(aa,3));
              tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + totgt;
        end
                
        % # of events from a given element that are not directed to the target (in a
        % given condition)
        if strcmp(nbs{nb},'1bnT')
            nottotgt = sum(auxdata.elfix==auxall(aa,1) & ~isnan(auxdata.next1elfix) & auxdata.next1elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
            tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgt;
        elseif strcmp(nbs{nb},'2b1bnT')
            % this is not considering the specific path
             nottotgt = sum(auxdata.elfix==auxall(aa,1) & ~isnan(auxdata.next2elfix) & auxdata.next2elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
             tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgt;
             % this is considering the same 1b
%              nottotgtpath = sum(auxdata.elfix==auxall(aa,1) & auxdata.next1elfix==auxall(aa,4) & auxdata.next2elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
%              tmax_cent_nottotgt_path(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt_path(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgtpath;
        elseif strcmp(nbs{nb},'3b2b1bnT')
               nottotgt = sum(auxdata.elfix==auxall(aa,1) & ~isnan(auxdata.next3elfix) & auxdata.next3elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
             tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgt;
         end
        
        % expected transitions of the same type given the transition matrix
        % (that takes in accoun starting position, vector of movement and condition)
        if strcmp(nbs{nb},'1bnT')
            tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
           (nottotgt+totgt).*tmax_bypos(Itgt,Jtgt,auxall(aa,1),find(conditions==auxall(aa,3)),:);
        elseif strcmp(nbs{nb},'2b1bnT')
            % when more than 1-back, we need to sum the probabilities 
            % according to all possibles path from starting element to
            % target. This is done in two ways: (1) using a one-movement 
            % transition matrix, multiplying the probabilities of a given path,
            % and summing the probabilities of different paths.
            % (2) using a combination of the one-movement transition path
            % and the empirical probabilities obserbed for n-back movements
           
            p2b1btgt    = 0;           % this takes all possible paths, always taking the  
            p2b1bfixtgt = 0;
            for e = 1:48   %loop trough possible destination of 2b to 1b saccade
                [I1b,J1b]       = ind2sub([6 8],e);
                p2bto1b         = tmax_bypos(I1b,J1b,auxall(aa,1),find(conditions==auxall(aa,3)),:);
                if p2bto1b>0
                   p2b1btgt     = p2b1btgt+p2bto1b.*tmax_bypos(Itgt,Jtgt,e,find(conditions==auxall(aa,3)),:); 
                   p2b1bfixtgt  = nansum([p2b1bfixtgt,p2bto1b.*resultTM.ptotarg{1}(6+I1b-Itgt,8+J1b-Jtgt)]); 
                end
            end
                tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
            (nottotgt+totgt).*p2b1btgt;
                tmax_exp_fix1b(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp_fix1b(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
            (nottotgt+totgt).*p2b1bfixtgt;
                tmax_exp_2mov(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp_2mov(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
           (nottotgt+totgt).*tmax_bypos2mov(Itgt,Jtgt,auxall(aa,1),find(conditions==auxall(aa,3)),:);
        
        elseif strcmp(nbs{nb},'3b2b1bnT')
            
            p3b2b1btgt    = 0;           % this takes all possible paths, always taking the  
            p3b2bfixtgt   = 0;
            p3b2b1bfixtgt = 0;
            for e2 = 1:48   %loop trough possible 3b to 2b transition
                [I2b,J2b]       = ind2sub([6 8],e2);
                 p3bto2b         = tmax_bypos(I2b,J2b,auxall(aa,1),find(conditions==auxall(aa,3)),:);
                if p3bto2b>0 % If 3b to 2b exist
                    p3b2bfixtgt  = nansum([p3b2bfixtgt,p3bto2b.*resultTM.ptotarg{2}(6+I2b-Itgt,8+J2b-Jtgt)]); 
                    for e1 = 1:48  %loop trough possible 2b to 1b location
                        [I1b,J1b]       = ind2sub([6 8],e1);
                        p2bto1b         = tmax_bypos(I1b,J1b,e2,find(conditions==auxall(aa,3)),:);
                        if p2bto1b>0 % If 2b to 1b exist
                            p3b2b1btgt     = p3b2b1btgt+p3bto2b.*p2bto1b.*tmax_bypos(Itgt,Jtgt,e1,find(conditions==auxall(aa,3)),:); 
                            p3b2b1bfixtgt  = nansum([p3b2b1bfixtgt,p3bto2b*p2bto1b.*resultTM.ptotarg{1}(6+I1b-Itgt,8+J1b-Jtgt)]);
                        end
                    end
                end
            end
                tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
            (nottotgt+totgt).*p3b2b1btgt;
                tmax_exp_fix2b(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp_fix2b(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
             (nottotgt+totgt).*p3b2bfixtgt;
            tmax_exp_fix1b(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp_fix1b(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
             (nottotgt+totgt).*p3b2b1bfixtgt;
        end
    end
     resultTM.num_excess{nb}         = tmax_cent_totgt-tmax_exp;
     resultTM.num_excess_fix1b{nb}   = tmax_cent_totgt-tmax_exp_fix1b;
     resultTM.ptotarg{nb}            = tmax_cent_totgt./(tmax_cent_totgt+tmax_cent_nottotgt); %p to move to the target from a given location and n back
     resultTM.expptotarg{nb}         = tmax_exp./(tmax_cent_totgt+tmax_cent_nottotgt);
     resultTM.abs_pdiff{nb}          = resultTM.ptotarg{nb}-resultTM.expptotarg{nb}; 
     resultTM.OR{nb}                 = (resultTM.ptotarg{nb}./(1-resultTM.ptotarg{nb}))./(resultTM.expptotarg{nb}./(1-resultTM.expptotarg{nb}));
     resultTM.expptotarg_fix1b{nb}   = tmax_exp_fix1b./(tmax_cent_totgt+tmax_cent_nottotgt);
     resultTM.abs_pdiff1b{nb}        = resultTM.ptotarg{nb}-resultTM.expptotarg_fix1b{nb}; 
     resultTM.OR1b{nb}               = (resultTM.ptotarg{nb}./(1-resultTM.ptotarg{nb}))./(resultTM.expptotarg_fix1b{nb}./(1-resultTM.expptotarg_fix1b{nb}));
     resultTM.expptotarg_fix2b{nb}   = tmax_exp_fix2b./(tmax_cent_totgt+tmax_cent_nottotgt);
     resultTM.abs_pdiff2b{nb}        = resultTM.ptotarg{nb}-resultTM.expptotarg_fix2b{nb}; 
     resultTM.OR2b{nb}               = (resultTM.ptotarg{nb}./(1-resultTM.ptotarg{nb}))./(resultTM.expptotarg_fix2b{nb}./(1-resultTM.expptotarg_fix2b{nb}));
     
   %      expptotarg_2mov{nb}    = tmax_exp_2mov./(tmax_cent_totgt+tmax_cent_nottotgt);
     resultTM.observado{nb}          = tmax_cent_totgt;
     resultTM.noobservado{nb}        = tmax_cent_nottotgt;
     resultTM.esperado{nb}           = tmax_exp;
     resultTM.esperado_fix1b{nb}     = tmax_exp_fix1b;
     resultTM.esperado_fix2b{nb}     = tmax_exp_fix2b;
     resultTM.esperado_2mov{nb}      = tmax_exp_2mov;
     resultTM.pout{nb}               = reshape(myBinomTest(resultTM.observado{nb}(:),resultTM.observado{nb}(:)+resultTM.noobservado{nb}(:),resultTM.expptotarg{nb}(:),'two'),[11 15]);
     resultTM.pout1b{nb}             = reshape(myBinomTest(resultTM.observado{nb}(:),resultTM.observado{nb}(:)+resultTM.noobservado{nb}(:),resultTM.expptotarg_fix1b{nb}(:),'two'),[11 15]);
     resultTM.pout2b{nb}             = reshape(myBinomTest(resultTM.observado{nb}(:),resultTM.observado{nb}(:)+resultTM.noobservado{nb}(:),resultTM.expptotarg_fix2b{nb}(:),'two'),[11 15]);

end
save(fullfile(patheye,'eyedata','transitionP.mat'),'resultTM') 
%%}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evidence and cumulative evidence and p to fix and fix duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(patheye,'eyedata','transitionP.mat'),'resultTM')
data.evidence           = zeros(5,length(data.trial));
proxyev                 =  resultTM.abs_pdiff{1} ;
blip = resultTM.observado{1}+resultTM.noobservado{1};

proxyev(isnan(proxyev)) = 0;        
proxyev(blip<62) = 0;        
YY = [];XY=[];
for ss = 1:length(subjects)
    indxS           = find(data.subject==subjects(ss));
    for tt = unique(data.trial(indxS))
        indxTfix    = find(data.trial==tt & data.type==1 & data.subject==subjects(ss) & data.orderPreT>0);
        [Itgt,Jtgt] = ind2sub([6 8],unique(data.tpos(indxTfix)));
        [Iel,Jel]   = ind2sub([6 8],data.elfix(indxTfix));
        tmindx      = sub2ind([11,15],6+(Iel-Itgt),8+(Jel-Jtgt));   
        tmindxnan   = find(isnan(tmindx));
        tmindx(tmindxnan)       = 1;
        auxevidence             = proxyev(tmindx);
        auxevidence(tmindxnan)  = 0;
        for tor = 2:5
            auxevidence(tor,tor:end) = auxevidence(tor-1,tor-1:end-1);%+auxevidence(1,tor:end);
        end
        data.evidence(:,indxTfix) = auxevidence; 
%         if ~isempty(indxTfix)
%             if length(indxTfix)>1
%                 YY = [YY;data.onTarg([indxTfix(2:end) indxTfix(end)+2])'];
%             else
%                 YY = [YY;data.onTarg([indxTfix(end)+2])'];
%             end
%             if size(auxevidence,2)>4
%                 XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN;auxevidence(1,1:end-2)'],[NaN;NaN;NaN;auxevidence(1,1:end-3)'],[NaN;NaN;NaN;NaN;auxevidence(1,1:end-4)'],[ss.*ones(size(auxevidence,2),1)]]];
%             elseif size(auxevidence,2)==4
%                 XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN;auxevidence(1,1:end-2)'],[NaN;NaN;NaN;auxevidence(1,1:end-3)'],[NaN;NaN;NaN;NaN],[ss;ss;ss;ss]]];
%             elseif size(auxevidence,2)==3
%                 XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN;auxevidence(1,1:end-2)'],[NaN;NaN;NaN],[NaN;NaN;NaN],[ss;ss;ss]]];
%             elseif size(auxevidence,2)==2
%                 XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN],[NaN;NaN],[NaN;NaN],[ss;ss]]];
%             elseif size(auxevidence,2)==1
%                 XY = [XY;[[auxevidence(1,:)'],[NaN],[NaN],[NaN],[NaN],[ss]]];
%             end
%         end
    end
end
save(fullfile(patheye,'eyedata','alleyedataFULLevidence.mat'),'data')    % info from subjects EDF files, augments with augmentinfALLdata.m 

% T = table(YY,XY(:,1),XY(:,2),XY(:,3),XY(:,4),XY(:,5),XY(:,6),'VariableNames',{'p_totarg','pre1','pre2','pre3','pre4','pre5','subject'});
%                  auxres = fitglme(T,'p_totarg ~ pre1 + pre2 + pre3 + pre4 + pre5 + (1|subject)','Distribution','Binomial')
%1/(1+exp(-(auxres.Coefficients.Estimate(1) + auxres.Coefficients.Estimate(2))))
%exp(auxres.Coefficients.Estimate(1) + auxres.Coefficients.Estimate(5))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance evidence and cumulative evidence and p to fix and fix duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(patheye,'eyedata','transitionP.mat'),'resultTM')
YY = [];XY=[];
for ss = 1:length(subjects)
    indxS           = find(data.subject==subjects(ss));
    for tt = unique(data.trial(indxS))
        indxTfix    = find(data.trial==tt & data.type==1 & data.subject==subjects(ss) & data.orderPreT>0);
        auxevidence             = (data.DToTarg(indxTfix)./posVec.pixxdeg);
%         auxevidence(auxevidence>10)          = 10;
%         for tor = 2:5
%             auxevidence(tor,tor:end) = auxevidence(tor-1,tor-1:end-1);%+auxevidence(1,tor:end);
%         end
%         data.evidence(:,indxTfix) = auxevidence; 
        if ~isempty(indxTfix)
            if length(indxTfix)>1
                YY = [YY;data.onTarg([indxTfix(2:end) indxTfix(end)+2])'];
            else
                YY = [YY;data.onTarg([indxTfix(end)+2])'];
            end
            if size(auxevidence,2)>4
                XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN;auxevidence(1,1:end-2)'],[NaN;NaN;NaN;auxevidence(1,1:end-3)'],[NaN;NaN;NaN;NaN;auxevidence(1,1:end-4)'],[ss.*ones(size(auxevidence,2),1)]]];
            elseif size(auxevidence,2)==4
                XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN;auxevidence(1,1:end-2)'],[NaN;NaN;NaN;auxevidence(1,1:end-3)'],[NaN;NaN;NaN;NaN],[ss;ss;ss;ss]]];
            elseif size(auxevidence,2)==3
                XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN;auxevidence(1,1:end-2)'],[NaN;NaN;NaN],[NaN;NaN;NaN],[ss;ss;ss]]];
            elseif size(auxevidence,2)==2
                XY = [XY;[[auxevidence(1,:)'],[NaN;auxevidence(1,1:end-1)'],[NaN;NaN],[NaN;NaN],[NaN;NaN],[ss;ss]]];
            elseif size(auxevidence,2)==1
                XY = [XY;[[auxevidence(1,end)'],[NaN],[NaN],[NaN],[NaN],[ss]]];
            end
        end
    end
end
% save(fullfile(patheye,'eyedata','alleyedataFULLevidence.mat'),'data')    % info from subjects EDF files, augments with augmentinfALLdata.m 
 T = table(YY,XY(:,1),XY(:,2),XY(:,3),XY(:,4),XY(:,5),XY(:,6),'VariableNames',{'p_totarg','pre1','pre2','pre3','pre4','pre5','subject'});
                  auxres = fitglme(T,'p_totarg ~ pre1 + pre2 + pre3 + pre4 + pre5 + (1|subject)','Distribution','Binomial')
% 1/(1+exp(-(auxres.Coefficients.Estimate(1) + auxres.Coefficients.Estimate(2))))
%exp(auxres.Coefficients.Estimate(1) + auxres.Coefficients.Estimate(5))
figure,,hold on
xx = .1:.1:20;
for preT = 1:5
   h(preT) = plot(xx,1./(1+exp(-(auxres.Coefficients.Estimate(1) + xx.*auxres.Coefficients.Estimate(preT+1)))),'.-') 
end

T = table(YY,XY(:,1),XY(:,6),'VariableNames',{'p_totarg','pre1','subject'});
                  auxres1 = fitglme(T,'p_totarg ~ pre1 + (1|subject)','Distribution','Binomial')
T = table(YY,XY(:,1)+XY(:,2),XY(:,6),'VariableNames',{'p_totarg','pre1maspre2','subject'});
                  auxres12 = fitglme(T,'p_totarg ~ pre1maspre2 + (1|subject)','Distribution','Binomial')
T = table(YY,XY(:,1)+XY(:,2)+XY(:,3),XY(:,6),'VariableNames',{'p_totarg','pre1maspre2maspre3','subject'});
auxres123 = fitglme(T,'p_totarg ~ pre1maspre2maspre3 + (1|subject)','Distribution','Binomial')
clear h
figure,hold on
h(1) = plot(xx,1./(1+exp(-(auxres1.Coefficients.Estimate(1) + xx.*auxres1.Coefficients.Estimate(2)))),'.-') 
h(2) = plot(xx,1./(1+exp(-(auxres12.Coefficients.Estimate(1) + xx.*auxres12.Coefficients.Estimate(2)))),'.-') 
h(3) = plot(xx,1./(1+exp(-(auxres123.Coefficients.Estimate(1) + xx.*auxres123.Coefficients.Estimate(2)))),'.-') 


% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FIGURES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setAbsoluteFigureSize

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSITION MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fignames = {'transitionMatrix','transitionMatrixTOtgt',...
            'ptotag','expptotarg','ptotag2','expptotarg1b',...
            'ptotag3','expptotarg2b','ptorarg-exptotarg',...
            'ptorarg2-exptotarg1b','ptorarg3-exptotarg2b'};
matrices = {nansum(nanmean(resultTM.tmax_cent,4),3)/48,resultTM.observado{1}./sum(resultTM.observado{1}(:)),...
            resultTM.ptotarg{1},resultTM.expptotarg{1},resultTM.ptotarg{2},resultTM.expptotarg_fix1b{2},...
            resultTM.ptotarg{3},resultTM.expptotarg{3},resultTM.abs_pdiff{1},...
            resultTM.abs_pdiff1b{2},resultTM.abs_pdiff2b{3}};
matricesp = {[],[],[],[],[],[],[],[],resultTM.pout{1}, resultTM.pout1b{2}, resultTM.pout2b{3}};        
clims    = [-3 1;-3 1;
            0 60;0 60;0 60;0 60;
            0 60;0 60;0 30;
            0 30;0 30];

for ff=1:length(matrices)
    thisTM          = matrices{ff};
    cmap            = brighten(flipud(cbrewer('seq','Blues',128)),-.4);
    cmap(110:128,:) = [];
%     avgTM           = nansum(nanmean(tmax_cent,4),3)/48;
    clim            = clims(ff,:);
    rectcolor       = [.8 0 0];
    fh              = figure;, hold on
    fg.Position(4)  = 560;
    fh.Name         = fignames{ff};
    if strcmp(fignames{ff},'transitionMatrix') || strcmp(fignames{ff},'transitionMatrixTOtgt')
        imagesc(-7:7,-5:5,log10(thisTM*100)), hold on  
    else
        imagesc(-7:7,-5:5,thisTM*100), hold on  
    end
    
    axis ij
    line([[-6.5:6.5]' [-6.5:6.5]']',[ones(14,1)*[-5.5 5.5]]','LineWidth',.1,'Color','k');
    line([ones(12,1)*[-7.5 7.5]]',[[-5.5:5.5]' [-5.5:5.5]']','LineWidth',.1,'Color','k');
    xlabel('Rel. horiz. cell position')
    ylabel('Rel. vert. cell position')
    if strcmp(fignames{ff},'transitionMatrix') || strcmp(fignames{ff},'transitionMatrixTOtgt')
        rectangle('Position',[-1.5,-1.5,3,3],'EdgeColor',rectcolor)
        text(1.5,-1.5 ,sprintf('%2.1f%%',sum(sum(thisTM(5:6,7:9)))*100),'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','top','Color',rectcolor)
        rectangle('Position',[-2.5,-2.5,5,5],'EdgeColor',rectcolor)
        text(2.5,-2.5 ,sprintf('%2.1f%%',sum(sum(thisTM(4:8,6:10)))*100),'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','top','Color',rectcolor)
        rectangle('Position',[-3.5,-3.5,7,7],'EdgeColor',rectcolor)
        text(3.5,-3.5 ,sprintf('%2.1f%%',sum(sum(thisTM(3:9,5:11)))*100),'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','top','Color',rectcolor)
    end
    if any(strfind(fignames{ff},'-'))
        auxp = matricesp{ff}<0.05/11/15/3;
        [I,J] = ind2sub([11 15],find(auxp)');
        for rr = 1:length(I)
            rectangle('Position',[(J(rr)-8)-.5,(I(rr)-6)-.5,1,1],'EdgeColor',rectcolor)
        end
    end
    set(gca,'Position',[.14 .17 .68 .8],'XTick',[-7,0,7],'YTick',[-5,0,5],'FontSize',10)
    colormap(cmap)
    set(gca)
    axis([-7.5 7.5 -5.5 5.5])
    caxis(clim)

    hc              = colorbar;
    hc.Position     = [.85 .17+(.8)/4 .04 (.8)/2];
    hc.Limits       = clim;
    if strcmp(fignames{ff},'transitionMatrix') || strcmp(fignames{ff},'transitionMatrixTOtgt')
        hc.Ticks        = clim(1):1:clim(end);
        hc.TickLabels   = 10.^hc.Ticks;
    else
        hc.Ticks        = clim(1):10:clim(end);
    end
    %    hc.TickLabels = {'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}'};%[10.^[-6:-1]]
    hc.FontSize     = 6;
    figsize         = [4.6 4.6*fh.Position(4)/fh.Position(3)];
    doimage(gcf,fullfile(patheye,'figures'),...
               'pdf',fignames{ff},figsize  ,1)
end