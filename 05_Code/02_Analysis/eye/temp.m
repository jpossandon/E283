%%
load('/Users/jossando/trabajo/E283/07_Analysis/03_Eye/eyedata/alleyedataFULLevidence')
edgest = 0:.5:30;
indfix = data.type==1;
indsac = data.type==2;

figure, hold on
for opt = 0:10
    indxopt = data.orderPreT==opt;
    [n,c]=hist(data.DToTarg(indfix & indxopt)/45,edgest);
    h(opt+1) = plot(c,n,'.-');
end
xlabel('DtoTarg')
doimage(gcf,fullfile(patheye,'figures'),...
               'tiff',['DtoTatg_perorder'],'600','painters',[],1)

edgest = 0:50:800;
figure, hold on
for opt = 0:10
    indxopt = data.orderPreT==opt;
    [n,c]=hist(data.dur(indfix & indxopt),edgest);
    h(opt+1) = plot(c,n,'.-')
end
xlabel('fixdur')
