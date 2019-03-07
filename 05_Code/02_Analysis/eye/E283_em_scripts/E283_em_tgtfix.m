%%
result = struct('nmovs',zeros(70,31),'nnt',zeros(31,1),'nntdur',zeros(31,3),'nntDtoTarg',zeros(31,3),...
    'ntt',zeros(31,1),'nttdur',zeros(31,3),'nttDtoTarg',zeros(31,3),...
    'tntt',zeros(31,1),'tnttdur',zeros(31,4),'tnttDtoTarg',zeros(31,4),...
    'ttt',zeros(31,1),'tttdur',zeros(31,3),'tttDtoTarg',zeros(31,3),...
    'tnt',zeros(31,1),'tntdur',zeros(31,3),'tntDtoTarg',zeros(31,3),...
    'nntnt',zeros(31,1),'nntntdur',zeros(31,5),'nntntDtoTarg',zeros(31,5),...
    'tt',zeros(31,1),'ttdur',zeros(31,2),'ttDtoTarg',zeros(31,2),...
    'nt',zeros(31,1),'ntdur',zeros(31,2),'ntDtoTarg',zeros(31,2),...
    't',zeros(31,1),'Mntn',zeros(31,1),'Mntndur',zeros(31,3),'MntnDtoTarg',zeros(31,3),...
    'trlmissed',zeros(31,1),'tmissed',zeros(31,1))

for ss = 1:length(subjects)
    ss
    indxS   = find(data.subject==subjects(ss));
    for tt = unique(data.trial(indxS))
        indxfix = find(data.subject==subjects(ss) & data.trial==tt & data.type==1 & data.end>0);
        auxOnT  = data.onTarg(indxfix);
        if any(auxOnT)
            if auxOnT(end)==1
                result.nmovs(length(auxOnT),ss) = result.nmovs(length(auxOnT),ss)+1;
                if length(auxOnT)>=3 
                    if sum(auxOnT(end-2:end) == [0 0 1])==3
                        result.nnt(ss) = result.nnt(ss)+1;
                        result.nntdur(ss,:)  = result.nntdur(ss,:)+data.dur(indxfix(end-2:end));
                        result.nntDtoTarg(ss,:)  = result.nntDtoTarg(ss,:)+data.DToTarg(indxfix(end-2:end));
                    elseif sum(auxOnT(end-2:end) == [0 1 1])==3
                        result.ntt(ss) = result.ntt(ss)+1;
                        result.nttdur(ss,:)  = result.nttdur(ss,:)+data.dur(indxfix(end-2:end));
                        result.nttDtoTarg(ss,:)  = result.nttDtoTarg(ss,:)+data.DToTarg(indxfix(end-2:end));
                        if length(auxOnT)>=4
                            if sum(auxOnT(end-3:end) == [1 0 1 1])==4
                            result.tntt(ss) = result.tntt(ss)+1;
                            result.tnttdur(ss,:)  = result.tnttdur(ss,:)+data.dur(indxfix(end-3:end));
                            result.tnttDtoTarg(ss,:)  = result.tnttDtoTarg(ss,:)+data.DToTarg(indxfix(end-3:end));
                            end
                        end
                    elseif sum(auxOnT(end-2:end) == [1 1 1])==3
                        result.ttt(ss) = result.ttt(ss)+1;
                        result.tttdur(ss,:)  = result.tttdur(ss,:)+data.dur(indxfix(end-2:end));
                        result.tttDtoTarg(ss,:)  = result.tttDtoTarg(ss,:)+data.DToTarg(indxfix(end-2:end));
                    
                    elseif sum(auxOnT(end-2:end) == [1 0 1])==3
                        result.tnt(ss) = result.tnt(ss)+1;
                        result.tntdur(ss,:)  = result.tntdur(ss,:)+data.dur(indxfix(end-2:end));
                        result.tntDtoTarg(ss,:)  = result.tntDtoTarg(ss,:)+data.DToTarg(indxfix(end-2:end));
                        if length(auxOnT)>=5
                            if sum(auxOnT(end-4:end) == [0 0 1 0 1])==5
                                result.nntnt(ss) = result.nntnt(ss)+1;
                                result.nntntdur(ss,:)  = result.nntntdur(ss,:)+data.dur(indxfix(end-4:end));
                                result.nntntDtoTarg(ss,:)  = result.nntntDtoTarg(ss,:)+data.DToTarg(indxfix(end-4:end));
                            end
                        end
                    end
                elseif length(auxOnT)==2 
                    if sum(auxOnT(end-1:end) == [0 1])==2
                        result.nt(ss) = result.nt(ss)+1;
                    elseif sum(auxOnT(end-1:end) == [1 1])==2
                        result.tt(ss) = result.tt(ss)+1;
                    end
                elseif length(auxOnT)==1  
                    result.t(ss) = result.t(ss)+1;
                end
            else
                result.trlmissed(ss) = result.trlmissed(ss)+1;
            end
            if length(auxOnT)>3 
                if any(auxOnT(1:end-3))
                    result.tmissed(ss) = result.tmissed(ss)+sum(auxOnT(1:end-3));
                    auxmOnT = find(auxOnT(1:end-3));
                    for aot = 1:length(auxmOnT)
                        if auxmOnT(aot)>1 & auxmOnT(aot)<length(auxOnT)-2
                           result.Mntn(ss) = result.Mntn(ss)+1;
                            result.Mntndur(ss,:)  = result.Mntndur(ss,:)+data.dur(indxfix(auxmOnT(aot)-1:auxmOnT(aot)+1));
                                result.MntnDtoTarg(ss,:)  = result.MntnDtoTarg(ss,:)+data.DToTarg(indxfix(auxmOnT(aot)-1:auxmOnT(aot)+1));
                           
                        end
                    end
                end
            end
        else
            result.trlmissed(ss) = result.trlmissed(ss)+1;
        end
    end
end
result.nntdur = result.nntdur./repmat(result.nnt,1,3);
result.nttdur = result.nttdur./repmat(result.ntt,1,3);
result.tnttdur = result.tnttdur./repmat(result.tntt,1,4);
result.tntdur = result.tntdur./repmat(result.tnt,1,3);
result.tttdur = result.tttdur./repmat(result.ttt,1,3);
result.nntntdur = result.nntntdur./repmat(result.nntnt,1,5);
result.Mntndur = result.Mntndur./repmat(result.Mntn,1,3);
result.nntDtoTarg = result.nntDtoTarg./repmat(result.nnt,1,3);
result.nttDtoTarg = result.nttDtoTarg./repmat(result.ntt,1,3);
result.tntDtoTarg = result.tntDtoTarg./repmat(result.tnt,1,3);
result.tttDtoTarg = result.tttDtoTarg./repmat(result.ttt,1,3);
result.nntntDtoTarg = result.nntntDtoTarg./repmat(result.nntnt,1,5);
result.MntnDtoTarg = result.MntnDtoTarg./repmat(result.Mntn,1,3);