function Events_Formation4_2(path,P)
% Rainfall events identification
load([path 'BF_Swift'])
LMG=0;
r0=0;
p5=LMG;
p6=r0;

PRE=cell(1,length(A));
for i=1:length(A)
  PRE{i}=CPM_RE(P(:,i),p5,p6);
end

clear LMG Pks i p1 p2 p3 p4 p5 p6 r0
% save([path 'RFE_TarP'],'A','BF','FE','IBFIm','LSP','P','PRE','Q')
S=matfile([path 'RFE_Swift'],'Writable',true);
S.P=P;
S.PRE=PRE;
% save([path 'RFE_Swift'],'P','PRE')
end