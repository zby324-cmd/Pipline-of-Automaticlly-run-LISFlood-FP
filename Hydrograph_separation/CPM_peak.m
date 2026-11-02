% Yiwen Mei
% University of Connecticut
% 11/17/2012

% Description:
% The code find the location of peaks

% Inputs:
% flow: Streamflow time series
%  p1 : BFImax (following Eckhardt, 2005)
%  p2 : Length of blocking period / general length of events
%  p3 : Basin response ratio (Norbiato et al., 2011)
%  p4 : Ratio for peak selection 
%  pc : Location of points in recession limb
%  pc1: Baseflow time series from Method 2 of CK_sl_f
%  pc2: Envelope flow hydrograph

% Outputs:
% Qma: Location of the peaks

function Pks=CPM_peak(flow,p1,p2,p3,p4,pc,pc1,pc2)

Qc=flow; % Initial and final Characteristic flow points
for t=2:length(flow)-1
  if flow(t)==flow(t-1) && flow(t)==flow(t+1)
    Qc(t)=NaN;
  end
end
Qc1(:,1:2)=[find(~isnan(Qc)) flow(~isnan(Qc))];

skm=0; % Maximum number of skiped time step
for x=1:size(Qc1,1)-1
  if Qc1(x,2)==Qc1(x+1,2)
    skm=skm+1;
  end
end

sk=0;
for x=1:size(Qc1,1)-skm-1
  if Qc1(x+sk,2)==Qc1(x+sk+1,2)
    Qc2(x,1)=floor((Qc1(x+sk,1)+Qc1(x+1+sk,1))/2); % Characteristic flow points
    Qc2(x,2)=Qc1(x+sk,2);
    sk=sk+1;
  else
    Qc2(x,:)=Qc1(x+sk,:);
  end
end

Qma=Qc2(:,1); % Maxima flow Characteristic points
for x=2:size(Qc2,1)-1
  if Qc2(x,2)>Qc2(x-1,2) && Qc2(x,2)>Qc2(x+1,2)
    Qma(x)=Qc2(x,1);
  else
    Qma(x)=NaN;
  end
end
Qma(isnan(Qma))=[];
Qma(length(Qma))=[];

tp_l=max(pc);
Qma(flow(Qma)./pc1(Qma)<=p3 | flow(Qma)./(p1*nanmean(flow))<p3 | Qma>tp_l | flow(Qma)./pc2(Qma)<p4,:)=NaN;
Qma(isnan(Qma))=[];

Qma(:,2)=1; % Exclude the adjadcent Maxima
for x=2:length(Qma)
  if Qma(x,1)-Qma(x-1,1)>p2
    Qma(x,2)=Qma(x-1,2)+1;
  else
    Qma(x,2)=Qma(x-1,2);
  end
end

Qma1=Qma;
for j=1:max(Qma(:,2))
  X=find(Qma(:,2)==j);
  if length(X)>2
    for x=2:length(X)-1
      if flow(Qma(X(x),1))<flow(Qma(X(x)-1,1)) && flow(Qma(X(x),1))<flow(Qma(X(x)+1,1))
        Qma1(X(x),:)=NaN;
      end
    end
  end
end
Qma=Qma1(:,1);
Qma(isnan(Qma))=[];

Qma(:,2)=1;
for x=2:length(Qma)
  if Qma(x,1)-Qma(x-1,1)>p2
    Qma(x,2)=Qma(x-1,2)+1;
  else
    Qma(x,2)=Qma(x-1,2);
  end
end

Qma1=Qma;
sk=0;
for j=1:max(Qma(:,2))
  X=Qma(Qma(:,2)==j,1);
  x=find(flow(X)<max(flow(X)));
  Qma1(x+sk,:)=NaN;
  sk=sk+length(X);
  x=find(Qma1(:,2)==j);
  if length(x)>1;
    x=x(x~=round(mean(x)));
    Qma1(x,:)=NaN;
  end
end
Qma=Qma1(:,1);
Qma(isnan(Qma))=[];

Pks=Qma;
