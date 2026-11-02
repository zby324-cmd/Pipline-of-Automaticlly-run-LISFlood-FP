% Yiwen Mei
% University of Connecticut
% 5/28/2013

% Description:
% The code calculate the baseflow time series based on UKIH method (Gutard
% 1992). This method is for daily flow series.

% Inputs:
% flow: flow time series;
%  p1 : Block periods
%  p2 : minima ratio

% Outputs:
% Qb: Baseflow time series

function Qb=UKIH(flow,p1,p2)

bm=floor(length(flow)/p1);

TP=ones(bm,1);
for b=1:bm
  [~,d]=min(flow(p1*b-(p1-1):p1*b));
  TP(b)=p1*(b-1)+d;
end

TP1=TP;
for t=2:length(TP)-1
  if flow(TP(t+1))/flow(TP(t))<=p2 && flow(TP(t-1))/flow(TP(t))<=p2
    TP1(t)=NaN;
  end
end
TP(isnan(TP1))=[];
TP=(union([1 length(flow)],TP))';

Qbs=zeros(length(flow),1);
for b=1:length(TP)-1
  for t=TP(b):TP(b+1)
    Qbs(t)=((flow(TP(b+1))-flow(TP(b)))/(TP(b+1)-TP(b)))*t+(TP(b+1)*flow(TP(b))-TP(b)*flow(TP(b+1)))/(TP(b+1)-TP(b));
    
    if Qbs(t)>flow(t)
      Qbs(t)=flow(t);
    end
  end
end

Qb=Qbs;
