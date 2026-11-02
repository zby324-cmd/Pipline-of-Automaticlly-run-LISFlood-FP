% Yiwen Mei
% University of Connecticut
% 11/17/2012

% Description:
% The code identify potential hydrologic events

% Inputs:
% peak: Location of peak flow
% RiPs: Points in rising limb
% RePs: Points in recession limb

% Outputs:
% SPEs: Locations of single peak hydrologic events
% MPEs: Locations of multiple peak hydrologic events
% MPEP: Locations of peaks of multiple peak hydrologic events

function FE=CPM_FE(peak,RiPs,RePs)

PE=zeros(length(peak)-2,3);
for x=1:length(peak)-1
  PE(x,1)=RiPs(find(peak(x)-RiPs>0,1,'last'));
  PE(x,2)=peak(x);
  PE(x,3)=RePs(find(peak(x)-RePs<0,1,'first'));
end

J1=zeros(size(PE,1),1);
J=zeros(size(PE,1)-1,1);
for x=1:size(PE,1)
  X=setdiff(1:size(PE,1),x);

  for y=1:length(X)
    if ~isempty(max(intersect(PE(x,1):PE(x,3),PE(X(y),1):PE(X(y),3))))
      J(y)=max(intersect(PE(x,1):PE(x,3),PE(X(y),1):PE(X(y),3)));
    else
      J(y)=0;
    end
  end
  if ~isempty(min(J(J~=0)))
    J1(x)=min(J(J~=0));
  else
    J1(x)=0;
  end
end
J=J1(J1~=0);

J=intersect(J,J);
MPEs=zeros(length(J),1);
MPEP=cell(length(J),1);
for y=1:length(J)
  MPEs(y,1)=min(PE(J1==J(y),1));
  MPEs(y,2)=max(PE(J1==J(y),3));
  MPEP{y}=PE(J1==J(y),2);
end

SPEs=PE(J1==0,:);

FE.SPEs=SPEs;
FE.MPEs=MPEs;
FE.MPEP=MPEP;
