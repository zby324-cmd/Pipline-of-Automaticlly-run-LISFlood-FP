function Baseflow_Separation3(path)
% load([path 'RA_TarP'])
dataFile = [path 'RA_Swift.mat'];

%% check file exist - Lzq25
if exist(dataFile, 'file')
    load(dataFile);
else
    error('File %s not exist，check the pathway.', dataFile);
end

% Initialization for baseflow separation
RNC=5e-5*ones(1,length(A));
LSP=0.827*24*A.^.2;

p1=LRW;
p2=RNC;
p3=r2MA;
p4=BFIU; % Initial BFI
p5=LSP;
p6=exp(-fp_var(1,:));
p7=uer;
p8=ler;

BF=cell(1,length(A));
for i=1:length(A)
  BF{i}=FRCK_RCK(Q(:,i),p1(i),p2(i),p3(i),p4(i),p5(i),p6(i),p7(i),p8(i));
end

BFIm=zeros(1,length(A));
for i=1:length(A);
  BFIm(i)=nansum(BF{i}.f)/nansum(Q(:,i));
end

% Iterate to find parameters for baseflow separation
it=6; % times of iteration
IBFIm=[BFIU;BFIm;zeros(it,length(A))];

for y=1:it
  p4=BFIm;

  BF=cell(1,length(A));
  for i=1:length(A)
    BF{i}=FRCK_RCK(Q(:,i),p1(i),p2(i),p3(i),p4(i),p5(i),p6(i),p7(i),p8(i));
  end

  for i=1:length(A);
    BFIm(i)=nansum(BF{i}.f)/nansum(Q(:,i));
  end

  IBFIm(y+2,:)=BFIm;
end
dBFI=diff(IBFIm,1);

clear BFIU BFIm LRW RNC dBFI fp i it ler p1 p2 p3 p4 p5 p6 p7 p8 r2MA uer y
save([path 'BF_Swift'],'A','BF','IBFIm','LSP','Q','-v7.3');
% S=matfile([path 'BF_Swift'],'Writable',true);
% S.BF=BF;
% S.IBFIm=IBFIm;
% S.LSP=LSP;
end

