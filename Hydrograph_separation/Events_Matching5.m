function Events_Matching5(path)
% load([path 'RFE_TarP'])
load([path 'RFE_Swift'])

% Rainfall-runoff event and RC
mlag=LSP;
p1=LSP;
p2=mlag;

evts=cell(1,length(A));
lag=cell(1,length(A));
for i=1:length(A)
      if isempty(FE{i}.MPEs(:,1))
        MPEs1=[];
        MPEs2=[];
      else
          MPEs1=FE{i}.MPEs(:,1);
          MPEs2=FE{i}.MPEs(:,2);
      end
      if isempty(FE{i}.SPEs(:,1))
          SPEs1=[];
          SPEs3=[];
      else
          SPEs1=FE{i}.SPEs(:,1);
          SPEs3=FE{i}.SPEs(:,3);
      end
      sh=sort([SPEs1;MPEs1]);
      eh=sort([SPEs3;MPEs2]);   
    
%   sh=sort([FE{i}.SPEs(:,1);FE{i}.MPEs(:,1)]);
%   eh=sort([FE{i}.SPEs(:,3);FE{i}.MPEs(:,2)]);

  mph=zeros(size(FE{i}.MPEP,1),1);
  for m=1:size(FE{i}.MPEP,1)
    [~,b]=max(Q(FE{i}.MPEP{m},i));
    mph(m)=FE{i}.MPEP{m}(b);
  end
  ph=sort(cat(1,FE{i}.SPEs(:,2),mph));

  sr=PRE{i}(:,1);
  er=PRE{i}(:,2);
  cr=PRE{i}(:,6);

  evts{i}=CPM(Q(:,i),P(:,i),BF{i}.f,sh,ph,eh,sr,cr,er,p1(i),p2(i));

  lag{i}=evts{i}.D(:,3);
  lag{i}((evts{i}.D(:,3)<0 | isnan(evts{i}.D(:,3)) | evts{i}.D(:,3)==max(evts{i}.D(:,3)) | evts{i}.R(:,1)>1))=[];
  mlag(i)=mean(lag{i});
end

it=5;
Imlag=[p2;mlag;zeros(it,length(A))];
for y=1:it
  p2=mlag;

  for i=1:length(A)
      if isempty(FE{i}.MPEs(:,1))
        MPEs1=[];
        MPEs2=[];
      else
          MPEs1=FE{i}.MPEs(:,1);
          MPEs2=FE{i}.MPEs(:,2);
      end
      if isempty(FE{i}.SPEs(:,1))
          SPEs1=[];
          SPEs3=[];
      else
          SPEs1=FE{i}.SPEs(:,1);
          SPEs3=FE{i}.SPEs(:,3);
      end
      sh=sort([SPEs1;MPEs1]);
      eh=sort([SPEs3;MPEs2]);
%     sh=sort([FE{i}.SPEs(:,1);FE{i}.MPEs(:,1)]);
%     eh=sort([FE{i}.SPEs(:,3);FE{i}.MPEs(:,2)]);

    mph=zeros(size(FE{i}.MPEP,1),1);
    for m=1:size(FE{i}.MPEP,1)
      [~,b]=max(Q(FE{i}.MPEP{m},i));
      mph(m)=FE{i}.MPEP{m}(b);
    end
    ph=sort(cat(1,FE{i}.SPEs(:,2),mph));

    sr=PRE{i}(:,1);
    er=PRE{i}(:,2);
    cr=PRE{i}(:,6);

    evts{i}=CPM(Q(:,i),P(:,i),BF{i}.f,sh,ph,eh,sr,cr,er,p1(i),p2(i));

    lag{i}=evts{i}.D(:,3);
    lag{i}((evts{i}.D(:,3)<0 | isnan(evts{i}.D(:,3)) | evts{i}.D(:,3)==max(evts{i}.D(:,3)) | evts{i}.R(:,1)>1))=[];
    mlag(i)=mean(lag{i});
  end

  Imlag(y+2,:)=mlag;
end

EP=cell(1,length(A)); % (1) tb_r, (2) tb_f, (3) te_r, (4) te_f, (5) tc_r, (6) tc_f, (7) tp, (8) Dr,

for i=1:length(A)     % (9) Df, (10) tlag, (11) Vr, (12) Vf, (13) Vb, (14) RC, (15) BFI, (16) Np,
                      % (17) tb_rA, (18) te_rA, (19) tb_rB, (20) te_rB, (21) tb_rC, (22) te_rC
  
  %% category of the flood event. 
  % 1: event of erroneous timing and volume
  % 2: event of correct volume but wrong timing
  % 3: event of both correct volume and timing
  eType=ones(size(evts{i}.tE,1),1);
  EP{i}=[evts{i}.tE evts{i}.D evts{i}.V evts{i}.R ones(size(evts{i}.D,1),1) evts{i}.GEs,eType];
%   EP{i}(evts{i}.D(:,3)<0 | evts{i}.R(:,1)>.9 | evts{i}.V(:,2)./evts{i}.D(:,2)<nanmean(Q(:,i)),:)=NaN;
%   EP{i}(isnan(EP{i}(:,1)),:)=[];
  EP{i}(evts{i}.R(:,1)<.9 & evts{i}.D(:,3)<0,end)=2;
  EP{i}(evts{i}.R(:,1)<.9 & evts{i}.D(:,3)>=0,end)=3;
%   for n=1:size(EP{i},1)
%     if isnan(EP{i}(n,2)) || isnan(EP{i}(n,7)) || isnan(EP{i}(n,4))
%         EP{i}(n,:)=NaN;
%     else
%         if Q(EP{i}(n,2),i)/Q(EP{i}(n,7),i)>IBFIm(size(IBFIm,1),i) || Q(EP{i}(n,4),i)/Q(EP{i}(n,7),i)>IBFIm(size(IBFIm,1),i)
%             EP{i}(n,:)=NaN;
%         end
%     end
%   end
%   EP{i}(isnan(EP{i}(:,1)),:)=[];
end

dmlag=diff(Imlag,1);

MPEP=cell(1,length(A));
for i=1:length(A)
  l=0;

  for o=1:size(FE{i}.MPEs,1)
    me=find(EP{i}(:,2)==FE{i}.MPEs(o,1));

    EP{i}(me,16)=length(FE{i}.MPEP{o});

    if ~isempty(me)
      l=l+1;

      MPEP{i}{l}=FE{i}.MPEP{o};
    end
  end
end

% The baseflow time series
FRCK=zeros(size(Q,1),length(A));
for i=1:length(A)
  FRCK(:,i)=BF{i}.f;
end

clear BF FE IBFIm Imlag LSP PRE b cr dmlag eh er evts i it k l lag m me
clear mlag mph n o p1 p2 ph sh sr tE y
% save([path 'Results_TarP'],'A','EP','FRCK','MPEP','P','Q')
save([path 'Results_Swift'],'A','EP','FRCK','MPEP','P','Q','-v7.3');
