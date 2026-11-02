% Yiwen Mei
% University of Connecticut
% 1/22/2012

% Description:
% The code locate the beginning and end of all the rainfall event

% Inputs:
% rain: Rainfall time series
%  p1: Maximum pause of rainfall
%  p2: Low-rain rate threshold (typically 0)

% Outputs:
%  RE: Beginnings, ends, duration, cumulative volume and centroids of
%      raifnall event

function RE=CPM_RE(rain,p1,p2)

rm=nanmean(rain);

P=rain;
fP=flipud(P);
for t=1:length(rain)-p1
  if nansum(P(t:t+p1))<=p2*(p1+1)*rm
    P(t)=NaN;
  end

  if nansum(fP(t:t+p1))<=p2*(p1+1)*rm
    fP(t)=NaN;
  end
end
fP=flipud(fP);
P(isnan(fP))=NaN;

Pa=find(P>=0);
Pb=ones(length(Pa),1);
for x=2:length(Pa)
  if Pa(x)>Pa(x-1)+1
    Pb(x)=Pb(x-1)+1;
  else
    Pb(x)=Pb(x-1);
  end
end

Pc=zeros(max(Pb),5);
for x=1:max(Pb)
  Pc(x,1)=Pa(find(Pb==x,1,'first'));
  Pc(x,2)=Pa(find(Pb==x,1,'last'));
  Pc(x,3)=length(Pc(x,1):Pc(x,2));
  Pc(x,4)=nansum(rain(Pc(x,1):Pc(x,2)));
  Pc(x,5)=max(rain(Pc(x,1):Pc(x,2)));
end

tcr=zeros(size(Pc,1),1);
for x=1:size(Pc,1)
  t=(Pc(x,1):Pc(x,2))';
  tcr(x)=nansum(t.*rain(t))/nansum(rain(t));
end

RE=[Pc tcr];
