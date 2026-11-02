function [sh,eh]=Events_Formation4_1(path,perctile)
%% extract flood events without precpitation
load([path 'BF_Swift'])

% Peak flow identification
p1=IBFIm(size(IBFIm,1),:);
p2=LSP;
p3=1./IBFIm(size(IBFIm,1),:);
p4=zeros(1,length(A));
% p4=[0 IBFIm(size(IBFIm,1),2) IBFIm(size(IBFIm,1),3) IBFIm(size(IBFIm,1),4) IBFIm(size(IBFIm,1),5) IBFIm(size(IBFIm,1),6) 0]; % Dont need to report

Pks=cell(1,length(A));
FE=cell(1,length(A));

for i=1:length(A)
  Pks{i}=CPM_peak(Q(:,i),p1(i),p2(i),p3(i),p4(i),BF{i}.RePs,BF{i}.f,BF{i}.Qelp);
  minPeak = prctile(Q(:,i),perctile);
  Pks{i}(Q(Pks{i})<minPeak)=[];% revised 1, to remove small floods
% Formation of flow hydrologic events
  Pks{i}(Pks{i}<min(BF{i}.RiPs))=[];% revised 2, remove peaks has no rising part due to the limited time span.
  FE{i}=CPM_FE(Pks{i},BF{i}.RiPs,BF{i}.RePs);
end
% S=matfile([path 'RFE_Swift'],'Writable',true);
% S.FE=FE;
sh=cell(1,length(A));
eh=cell(1,length(A));
for i=1:length(A)
  %% moved from Events_Matching5
  % rules to remove events with high start or end flow
  if ~isempty(FE{i}.MPEs)
      mph=zeros(size(FE{i}.MPEs,1),1);% peak locations of multi-peak events
      for ie=1:length(mph)
          [~,mph(ie)]=max(Q(FE{i}.MPEs(ie,1):FE{i}.MPEs(ie,2)));
          mph(ie)=mph(ie)+FE{i}.MPEs(ie,1)-1;
      end
      indexRemove=Q(FE{i}.MPEs(:,1))./Q(mph)>IBFIm(size(IBFIm,1),i) | Q(FE{i}.MPEs(:,2))./Q(mph)>IBFIm(size(IBFIm,1),i);
      FE{i}.MPEs(indexRemove,:)=[];
      FE{i}.MPEP(indexRemove,:)=[];
      MPEs1=FE{i}.MPEs(:,1);
      MPEs2=FE{i}.MPEs(:,2);
  else
      MPEs1=[];
      MPEs2=[];
  end
  if ~isempty(FE{i}.SPEs)
      sph=zeros(size(FE{i}.SPEs,1),1); %peak locations of single-peak events
      for ie=1:length(sph)
          [~,sph(ie)]=max(Q(FE{i}.SPEs(ie,1):FE{i}.SPEs(ie,3)));
          sph(ie)=sph(ie)+FE{i}.SPEs(ie,1)-1;
      end
      indexRemove=Q(FE{i}.SPEs(:,1))./Q(sph)>IBFIm(size(IBFIm,1),i) | Q(FE{i}.SPEs(:,3))./Q(sph)>IBFIm(size(IBFIm,1),i);
      FE{i}.SPEs(indexRemove,:)=[];
      SPEs1=FE{i}.SPEs(:,1);
      SPEs3=FE{i}.SPEs(:,3);
  else
      SPEs1=[];
      SPEs3=[];
  end
  ish=sort([SPEs1;MPEs1]);
  ieh=sort([SPEs3;MPEs2]);
  sh{i}=ish;
  eh{i}=ieh;
end

save([path 'RFE_Swift'],'A','BF','FE','IBFIm','LSP','Q','-v7.3');

end
