%%
% pooled transition matrices and target missesper subject
tmax_bypos = zeros(6,8,48,8);
tmax_cent  = nan(11,15,48,8);
conditions = [1 2 5 6 9 10 13 14];

for cc = 1:length(conditions)
    auxdata    = struct_select(data,{'value'},{['==' num2str(conditions(cc))]},2);

    for el = 1:48
        % fixations after start
        auxindx                 = find(auxdata.type==1 & auxdata.elfix==el & auxdata.orderPreT>0 & auxdata.start>0);
        auxindx(find(auxdata.type(auxindx+2)~=1))=[];
        auxnext                 = auxdata.elfix(auxindx+2);
        auxnext(isnan(auxnext)) = [];
        tmax_bypos(:,:,el,cc)      = reshape(accumarray(auxnext',1,[48,1]),[6 8])./length(auxnext);
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
   
  
%    % nrefix, nrevisit, lags
%    nrefix(ss)           = sum(auxdata.refix)./length(unique(auxdata.trial));
%    nrevisits(ss)        = sum(auxdata.revisit>0)./length(unique(auxdata.trial));
%    nnewfixs(ss)            = sum(~auxdata.refix & auxdata.revisit==0)./length(unique(auxdata.trial));
%    nlags{ss}            = auxdata.revisit(auxdata.revisit>0);

%%
% probability to fixate target 1-2-3 back by pos
% pooled data analysis
tmax_cent_totgt         = zeros(11,15);
tmax_cent_nottotgt      = zeros(11,15);
tmax_exp                = zeros(11,15);
tmax_exp_fix1b          = zeros(11,15);


auxdata                 = struct_select(data,{'type','value'},{'==1','>0'},2);
auxdata.next1elfix      = [auxdata.elfix(2:end) NaN];
auxdata.next1elfix(auxdata.orderPreT==0) = NaN;
auxdata.next2elfix      = [auxdata.elfix(3:end) NaN];
auxdata.next2elfix(auxdata.orderPreT==0 | auxdata.orderPreT==1) = NaN;

% one back target and not in target
clear excess ptotarg expptotarg
nbs = {'1bnT','2b1bnT'};
for nb = 1:length(nbs)
    
    % selecting data according to type of n-back analysis
    if strcmp(nbs{nb},'1bnT')
        auxindx             = find(auxdata.orderPreT==1 & auxdata.onTarg==0);
    elseif strcmp(nbs{nb},'2b1bnT')
        auxindx             = find(auxdata.orderPreT==2);
        auxindx(find(auxdata.onTarg(auxindx+1))) = [];                      % do not take in account 2-back that have 1-back in target
    end
    auxel           = auxdata.elfix(auxindx);
    auxel           = auxdata.elfix(auxindx);
    auxtgt          = auxdata.tpos(auxindx);
    auxcond         = auxdata.value(auxindx);
    
    mtx             = [auxel',auxtgt',auxcond'];
    
    if strcmp(nbs{nb},'1bnT')
        auxall          = [auxel',auxtgt',auxcond'];
    elseif strcmp(nbs{nb},'2b1bnT')
        auxel1b         = auxdata.elfix(auxindx+1);
        auxall          = [auxel',auxtgt',auxcond',auxel1b'];
    end
    % remove NaNs
    auxall(any(isnan(auxall),2),:) = [];
    [~,IA]  = unique(auxall(:,1:3),'rows');
    auxall  = auxall(IA,:);
    
    % loop trough all empirical n-back position to target position in every given condition
    for aa = 1:size(auxall,1)
   
        % find starting and target position in transition matrix
        [Iel,Jel]    = ind2sub([6 8],auxall(aa,1));
        [Itgt,Jtgt]  = ind2sub([6 8],auxall(aa,2));
        % # of movement from a given element position to a given target (in a given
        % condition)
        if strcmp(nbs{nb},'1bnT')
            totgt = sum(auxdata.orderPreT==1  & auxdata.onTarg==0 & auxdata.elfix==auxall(aa,1) & auxdata.next1elfix==auxall(aa,2) & auxdata.value==auxall(aa,3));
            tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + totgt;
        elseif strcmp(nbs{nb},'2b1bnT')
              % this is not considering the specific path
             totgt = sum(auxdata.orderPreT==2  & auxdata.elfix==auxall(aa,1) & auxdata.next2elfix==auxall(aa,2) & auxdata.value==auxall(aa,3));
              tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_totgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + totgt;
        end
                
        % # of events from a given element that are not directed to the target (in a
        % given condition)
        if strcmp(nbs{nb},'1bnT')
            nottotgt = sum(auxdata.elfix==auxall(aa,1) & auxdata.next1elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
            tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgt;
        elseif strcmp(nbs{nb},'2b1bnT')
            % this is not considering the specific path
             nottotgt = sum(auxdata.elfix==auxall(aa,1) & auxdata.next2elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
             tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgt;
             % this is considering the same 1b
%              nottotgtpath = sum(auxdata.elfix==auxall(aa,1) & auxdata.next1elfix==auxall(aa,4) & auxdata.next2elfix~=auxall(aa,2) & auxdata.tpos==auxall(aa,2) & auxdata.value==auxall(aa,3));
%              tmax_cent_nottotgt_path(6+(Iel-Itgt),8+(Jel-Jtgt)) =  tmax_cent_nottotgt_path(6+(Iel-Itgt),8+(Jel-Jtgt)) + nottotgtpath;
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
                   p2b1bfixtgt  = nansum([p2b1bfixtgt,p2bto1b.*ptotarg{1}(6+I1b-Itgt,8+J1b-Jtgt)]); 
                end
            end
                tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
            (nottotgt+totgt).*p2b1btgt;
                tmax_exp_fix1b(6+(Iel-Itgt),8+(Jel-Jtgt)) = tmax_exp_fix1b(6+(Iel-Itgt),8+(Jel-Jtgt)) +...
            (nottotgt+totgt).*p2b1bfixtgt;
        end
    end
     morethanX{nb}          = [tmax_cent_totgt+tmax_cent_nottotgt]>500;
     excess{nb}             = tmax_cent_totgt./tmax_exp;
     excess_fix1b{nb}       = tmax_cent_totgt./tmax_exp_fix1b;
     ptotarg{nb}            = tmax_cent_totgt./(tmax_cent_totgt+tmax_cent_nottotgt);
     expptotarg{nb}         = tmax_exp./(tmax_cent_totgt+tmax_cent_nottotgt);
     expptotarg_fix1b{nb}   = tmax_exp_fix1b./(tmax_cent_totgt+tmax_cent_nottotgt);
end
% two back target, one back in target
auxindx             = find(auxdata.orderPreT==2);


% end