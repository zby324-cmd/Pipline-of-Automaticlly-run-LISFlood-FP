% Yiwen Mei
% University of Connecticut
% 7/25/2013

% Description:
% The code associates the rainfall events with runoff events and calculate
% the correspoinding event properties

% Inputs:
%   flow  : Streamflow time series
%   rain  : Rainfall time series
% baseflow: Baseflow time series
%    sh   : Starting points of flow events
%    ph   : Peaks of flow events
%    eh   : End points of flow events
%    sr   : Starting points of rainfall events
%    cr   : Centroid of rainfall events
%    er   : End points of rainfall events
%    p1   : Length of searching window (h)
%    p2   : Mean event time lag (h)

% Outputs:
% Mat_RE.D : Durations of combined-rainfall event, flow event and event
%            time lag
% Mat_RE.V : Volume of event rainfall, flow and baseflow
% Mat_RE.tE: Timing of the matched rainfall and flow events (beginnings of
%            combined-rainfall & flow, ends of combined-rainfall & flow,
%            centroids of combined-rainfall & flow and maximum flow peak)
% Mat_RE.R : Event-based runoff coefficient and baseflow index

function Mat_RE=CPM(flow,rain,bFRCK,sh,ph,eh,sr,cr,er,p1,p2)

eh=[1;eh];
rain(isnan(rain))=0;

ch=nan(size(sh,1),1);
ccr=nan(size(sh,1),1);

pbE=cell(size(sh,1),1);
pcE=cell(size(sh,1),1);
precE=cell(size(sh,1),1);

GAEs=NaN*zeros(size(sh,1),2);
GBEs=NaN*zeros(size(sh,1),2);
for n=1:size(sh,1)
  th=sh(n):eh(n+1);
  ch(n)=nansum(flow(th)'.*th)/nansum(flow(th)); % Centroid of flow event

% Start of hydrograph after rainfall within a range
% Range is the mean time lag for wet period or the gap between two flow
% event, LSP is used for quality control
  pbE{n}=find(sh(n)-sr>=0 & sh(n)-sr<min([p1,sh(n)-eh(n),p2]));
  sr_e=findGAREvt(sr,er,sh(n),eh(n),p1,flow,rain);
  if ~isempty(sr_e)
      pbE{n}=sr_e;
  end
% Centroid of hyetograph before hydrograph but after the last pbE{n} 
  if ~isempty(pbE{n})
    pcE{n}=find(cr>cr(pbE{n}(length(pbE{n}))) & ch(n)>=cr);
  else % Centroid of hyetograph before hydrograph but after start of hydrograph
    pcE{n}=find(cr>sh(n) & ch(n)>=cr);
  end
  precE{n}=union(pcE{n},pbE{n});

  if ~isempty(precE{n})
    tcr=sr(precE{n}(1)):er(precE{n}(length(precE{n})));
    ccr(n)=nansum(rain(tcr)'.*tcr)/nansum(rain(tcr));
  end

  if ~isempty(pbE{n})
    GAEs(n,:)=[sr(min(pbE{n})) er(max(pbE{n}))];
  end
  if ~isempty(pcE{n})
    GBEs(n,:)=[sr(min(pcE{n})) er(max(pcE{n}))];
  end
end

lag=ch-ccr;

postE=cell(size(sh,1),1);
aRE=cell(size(sh,1),1);

Dh=nan(size(sh,1),1);
Dcr=nan(size(sh,1),1);

Vcr=nan(size(sh,1),1);
Vt=nan(size(sh,1),1);
Vb=nan(size(sh,1),1);

tE=nan(size(sh,1),7);

GCEs=nan(size(sh,1),2);
for n=1:size(sh,1)
  postE{n}=find(ch(n)<cr & eh(n+1)>er & eh(n+1)-cr>lag(n));

  if ~isempty(postE{n})
    GCEs(n,:)=[sr(min(postE{n})) er(max(postE{n}))];
  end

  aRE{n}=union(precE{n},postE{n});

  if ~isempty(aRE{n})
    tcr=sr(aRE{n}(1)):er(aRE{n}(length(aRE{n}))); % Period of combined-
    th=sh(n):eh(n+1); % Period of flow event      % rainfall event

    Dh(n)=length(~isnan(th)); % Duration of flow event
    Dcr(n)=length(tcr); % Duration of rainfall event

    Vcr(n)=nansum(rain(tcr)); % Combined-event rainfall volume
    Vt(n)=nansum(flow(th)); % Event flow volume
    Vb(n)=nansum(bFRCK(th)); % Event baseflow volume by FRCK

    ccr(n)=nansum(rain(tcr)'.*tcr)/Vcr(n); % Centroid of combined-rainfall event
    tE(n,:)=[min(tcr) min(th) max(tcr) max(th) ccr(n) ch(n) ph(n)]; % Event timing

  else
    Dh(n)=NaN;
    Dcr(n)=NaN;
    Vcr(n)=NaN;
    Vt(n)=NaN;
    Vb(n)=NaN;
    ccr(n)=NaN;
    %% revised to keep events of wrong precip
    tE(n,:)=[NaN sh(n) NaN eh(n+1) NaN ch(n) ph(n)];
  end
end

lag=ch-ccr; % Event time lag
BFI=Vb./Vt; % Event-based baseflow index by FRCK
RC=(Vt-Vb)./Vcr; % Event-based runoff coefficient by FRCK

Mat_RE.D=[Dcr Dh lag];
Mat_RE.V=[Vcr Vt Vb];
Mat_RE.tE=tE;
Mat_RE.R=[RC BFI];
Mat_RE.GEs=[GAEs GBEs GCEs];
