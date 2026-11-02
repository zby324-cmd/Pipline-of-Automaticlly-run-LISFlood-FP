% Yiwen Mei
% University of Connecticut
% 11/17/2012

% Description:
% The code a) locate the turning points of the hydrograph; b) construct
% envelope flow time series for turning points and c) separate baseflow
% from streamflow by two methods; 

% Inputs:
% flow: Streamflow time series
%  p1 : Half-regression period (half length of the period to calculate the
%       recession constant. 6h)
%  p2 : Constant slope threshold (1e-7min-2 is suggested in Blume et al.
%       2007. 1e-5h-2 is using)
%  p3 : Recession restriction (r2 statistics of the Q-t corellation. r2
%       lower than 0.9 is not considered as recession period)
%  p4 : BFImax (following Eckhardt, 2005)
%  p5 : Length of blocking period (N=0.827*A^0.2 is suggested by Chow, 1988)
%  p6 : Recession constant K (K=exp(-t/tc))
%  p7 : Upper limit of envelope baseflow rate
%  p8 : Lower limit of envelope baseflow rate

% Outputs:
% Baseflow.RiPs: Location of points in rising limb
% Baseflow.RePs: Location of points in recession limb
% Baseflow.Qelp: Flow envelope hydrograph
% Baseflow.belp: Baseflow envelope hydrograph
%   Baseflow.k : Information of the recession coefficients
%  Baseflow.sl : Baseflow time series from Method 1
%   Baseflow.f : Baseflow time series from Method 2

function Baseflow=FRCK_RCK(flow,p1,p2,p3,p4,p5,p6,p7,p8)

k=NaN*ones(length(flow),3);
for t=p1+1:length(flow)-p1
  x=(t-p1:t+p1)';
  X=[ones(size(x)) x];
  y=log(flow(x));
  if length(find(~isnan(y)))>round(1.5*p1+1)
%     [b,~,~,~,r]=regress(y,X);
    [b,~,r]=regressNoLic(y,X); %% revised
    k(t,1)=-b(2); % Slope for every time step
    k(t,2)=r(1); % r2 stats for every time step
  end
end

for t=2:length(flow)
  k(t,3)=k(t,1)-k(t-1,1); % Change of slope in h-2
end

mb=NaN*ones(length(flow)-round(p5)+1,1); % Maxima of each block
for t=1:length(flow)-round(p5)+1
  if ~isempty(find(~isnan(flow(t:t+round(p5)-1)),1))
    mb(t)=find(flow(t:t+round(p5)-1)>=max(flow(t:t+round(p5)-1)),1)+t-1;
  end
end

for t=2:length(flow)-round(p5)
  if ~isnan(mb(t))
    if mb(t)~=mb(t-1) && mb(t)~=mb(t+1) % Find the plateau
      mb(t)=NaN;
    end
  end
end
mb(isnan(mb))=[];
mb=unique(mb);

mb1=mb; % Compared to the two outliers
for bk=2:length(mb)-1
  if flow(mb(bk))>=flow(mb(bk-1)) && flow(mb(bk))>=flow(mb(bk+1))
    mb1(bk)=mb(bk);
  else
    mb1(bk)=NaN;
  end
end
mb1(isnan(mb1))=[];
mb=cat(1,0,mb1);

Qelp=flow; % Envelope flow rate series (form by the maxima)
for bk=2:length(mb)-1
  Qelp(round((mb(bk)+mb(bk-1))/2):round((mb(bk)+mb(bk+1))/2))=flow(mb(bk));
end
Qelp(1:round((mb(2)+mb(1))/2))=flow(mb(2));
Qelp(round((mb(length(mb))+mb(length(mb)-1))/2):length(flow))=flow(mb(length(mb)-1));
Qelp(Qelp<mean(flow))=mean(flow);

belp=p4*Qelp;
belp(belp<p8)=p8;
belp(belp>p7)=p7;

Qck=flow;
Qck(abs(k(:,3))>p2 | (k(:,2)<p3 & flow>belp) | flow>belp | k(:,1)<0)=NaN; % Select the points
Qck([1:p1 length(flow)-p1+1:length(flow)])=NaN;

RePs=find(Qck>=0);

rp=find(flow<=belp & k(:,2)<p3);
idx=0;
for t=2:length(rp)-1
  if flow(rp(t))<=flow(rp(t-1)) && flow(rp(t))<flow(rp(t+1))
    idx=idx+1;
    RiPs(idx,1)=rp(t);
  end
end

TPs=union([1;RiPs;RePs;length(flow)],[1;RiPs;RePs;length(flow)]);
TPs=[TPs flow(TPs)]; % Turning points

CK_sl=zeros(length(flow),1); % Baseflow time series (Method 1)
for tao=1:size(TPs,1)-1;
  for t=TPs(tao,1):TPs(tao+1,1)
    CK_sl(t)=((TPs(tao+1,2)-TPs(tao,2))/(TPs(tao+1,1)-TPs(tao,1)))*t+(TPs(tao+1,1)*TPs(tao,2)-TPs(tao,1)*TPs(tao+1,2))/(TPs(tao+1,1)-TPs(tao,1));
    if CK_sl(t)>flow(t) || CK_sl(t)<0;
      CK_sl(t)=flow(t);
    end
  end
end

CK_f=flow;
for t=2:length(flow)
  CK_f(t)=((1-p4)*p6*CK_sl(t-1)+(1-p6)*p4*flow(t))/(1-p6*p4);
  if CK_f(t)<0
    CK_f(t)=0;
  elseif CK_f(t)>flow(t)
    CK_f(t)=flow(t);
  end
end

Baseflow.RiPs=RiPs;
Baseflow.RePs=RePs;
Baseflow.Qelp=Qelp;
Baseflow.belp=belp;
Baseflow.k=k;
Baseflow.sl=CK_sl;
Baseflow.f=CK_f;
