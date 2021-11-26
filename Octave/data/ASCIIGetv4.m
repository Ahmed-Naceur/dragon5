% ASCSIIGetv4
% Get the record data in an ASCII file created by DRAGON
% AUTHOR   :      R. CHAMBON
% date     :      nov 2009
%
function data = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,RecToGet)
nbdata=0;
siz=size(RecordPos);
MaxLvl=length(siz);
if MaxLvl==2
    if siz(2)==1
        MaxLvl=1;
    end
end
for i=1:siz(ilvl)-1
  CurRecInd(ilvl)=i;
  strpos='';
  for i=1:ilvl
    strpos=[strpos int2str(CurRecInd(i)+1) ','];
  end
  for i=ilvl+1:MaxLvl
    strpos=[strpos int2str(CurRecInd(i)) ','];
  end
  strpos=[strpos '1' ];
  pos=eval(['RecordPos(' strpos ')']);
  lgn=eval(['RecordNumb(' strpos ')']);
  if lgn~=0
    if(RecordName(lgn,:)==RecToGet)
      fseek(pfin,pos,'bof');
      Line=fgetl(pfin);
      typdata=str2num(Line(19:26));
      nbdata=str2num(Line(27:34));
      if(strcmp(RecToGet(1:4),'elt#')==0)
        Line=fgetl(pfin);
      end
      break;
    end
  else
    break;
  end
end
if nbdata==0
  error(['ASCIIGetv4: record ''' RecToGet ''' not found'])
end

if nbdata>0
  if typdata==3
    tmp='';
  else
    tmp=[ ];
  end
  if ~feof(pfin)
    if typdata==1
      nline=ceil(nbdata/8);
    elseif typdata==2
      nline=ceil(nbdata/5);
    elseif typdata==3
      nline=ceil(nbdata/20);
      nskip=ceil(nbdata/8);
      for i=1:nskip
        fgetl(pfin); % skip '         4' lines
      end
    elseif typdata==4
      nline=ceil(nbdata/4);
    end
    for i=1:nline
      Line=fgetl(pfin);
      if typdata==3
        tmp=[tmp Line];
      else
        tmp=horzcat(tmp,str2num(Line));
      end
    end
  end
  if typdata==3
    data=tmp;
  else
    data=tmp';
  end
  CurRecInd(ilvl)=1;
end
