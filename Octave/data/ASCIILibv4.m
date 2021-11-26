% ASCSIILib
% Get the table of content of active directory
% AUTHOR   :      R. CHAMBON
% date     :      nov 2009
%
function ASCIILibv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd)
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
    fseek(pfin,pos,'bof');
    Line=fgetl(pfin);
    typdata=str2num(Line(19:26));
    nbdata=str2num(Line(27:34));
    if(typdata==3)
      nbdata=nbdata*4;
    end
    disp(['-> '  RecordName(lgn,:) '  len= ' num2str(nbdata)]);
  else
    break;
  end
end
