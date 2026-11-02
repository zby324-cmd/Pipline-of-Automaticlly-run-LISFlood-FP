function Recession_Analysis2(Q,P,A,dirWork)
%% input
% Q stream flow (m^3/s)
% P precipitaton (mm/h)
% A area (km^2)
P(P<0)=NaN;
for i=1:length(A)
  Q(:,i)=3600000*Q(:,i)/(1000000*A(i));
end
% Baseflow index under UKIH method
a=1.11;
Lc=5; % Length of searching step (h)
Qd=zeros(size(Q,1)/24,length(A));
for i=1:length(A)
  for d=1:size(Q,1)/24
    bk=Q(1+24*(d-1):24*d,i);

    k=~isnan(bk);
    Qd(d,i)=sum(bk(k))/length(k);
  end
end
Qd(Qd==0)=NaN;

Qdi=zeros(size(Q,1)/24,length(A)); % Interpolated daily streamflow
for i=1:length(A)
  Qe=[(1:size(Qd,1))' Qd(:,i)];
  Qe(isnan(Qe(:,2)),:)=[];
  Qdi(:,i)=interp1(Qe(:,1),Qe(:,2),1:size(Qd,1));
end

Qbs=zeros(size(Qd,1),size(Qd,2));
BFIU=zeros(1,length(A));
for i=1:length(A)
  Qbs(:,i)=UKIH(Qdi(:,i),Lc,a);
  BFIU(i)=nansum(Qbs(:,i))/nansum(Qdi(:,i));
end

% Empirical CDF of flow
ugp=.06;
bgp=.02;

ulr=zeros(1,length(A));
llr=zeros(1,length(A));
uer=zeros(1,length(A));
ler=zeros(1,length(A));

CDF=zeros(size(Q,1),size(Q,2));
for i=1:length(A)
  CDF(:,2*i-1)=sort(Q(:,i));
  for n=1:size(Q,1)
    CDF(n,2*i)=n/size(Q,1);
  end
  [~,c]=min(abs(nanmean(Q(:,i))-CDF(:,2*i-1)));
  ulr(i)=CDF(c,2*i)+ugp;
  llr(i)=CDF(c,2*i)-bgp;

  [~,c]=min(abs(ulr(i)-CDF(:,2*i)));
  uer(i)=CDF(c,2*i-1);
  [~,c]=min(abs(llr(i)-CDF(:,2*i)));
  ler(i)=CDF(c,2*i-1);
end

% Determine the recession coefficient
LRW=9;
r2MA=.88;
p1=LRW;
p2=r2MA;

K=cell(1,length(A));
fp=zeros(2,length(A));
for i=1:length(A);
  k=NaN*ones(length(Q(:,i)),2);
  for t=p1(i)+1:length(Q(:,i))-p1(i)
    x=(t-p1(i):t+p1(i))';
    X=[ones(length(x),1) x];
    y=Q(x,i);
    if length(find(~isnan(y)))>round(1.5*p1(i)+1) % Avoid the NaN
%       [b,~,~,~,r]=regress(y,X);
      [b,~,r]=regressNoLic(y,X); %% revised
      k(t,1)=-b(2); % Q_t slope for every bin
      k(t,2)=Q(t,i);
      k(t,3)=r(1); % coefficient of determination
    end
  end
  k(k(:,1)<=1e-9 | k(:,2)==0 | k(:,3)<p2(i),:)=NaN; % 1e-7 is estimated from the plots
  k(isnan(k(:,1)) | isnan(k(:,2)),:)=[];
  K{i}=k;
  [f,~,g]=regressNoLic(k(:,1),[ones(length(k(:,2)),1) k(:,2)]);%revised
%   [f,~,~,~,g]=regress(k(:,1),[ones(length(k(:,2)),1) k(:,2)]);
  fp(:,i)=[f(2);g(1)];
end

% rename fp to avoid conflict - Lzq25
fp_var = fp;
clear fp

clear CDF K Lc Qbs Qd Qdi Qe X a b bgp bk c d f g i k llr n p1 p2 r t ugp
clear ulr x y
save([dirWork 'RA_Swift'],'A','BFIU','LRW','Q','fp_var','ler','r2MA','uer', '-v7.3')
