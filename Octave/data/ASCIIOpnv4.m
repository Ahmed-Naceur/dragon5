% ASCIIOpnv4
% Open file, read it, extract the structure of the records
% AUTHOR   :      R. CHAMBON
% date     :      nov 2009

function [pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd] = ASCIIOpnv4(FileName);
pfin=fopen( FileName );
indrec=1;
CurRecInd=0;
MaxLvl=1;
if pfin>0
  ilvl=1;
  Line=fgetl(pfin);
  while (~feof(pfin))
    lnotbypass=1;
    curlvl=str2num(Line(3:10));
    currecnamelen=str2num(Line(11:18));
    curtype=str2num(Line(19:26));
    curnum=str2num(Line(27:34));
    pos=ftell(pfin)-73;
    %73 to point at the begining of the record line '-> ... <-
    if (currecnamelen~=0)
      Line=fgetl(pfin);
      currecname=Line(1:12);
    elseif curlvl>0 && curtype ~=99
      currecname=['elt#' Line(73:80)];
      pos=pos-8;
      % -8 for list element records are longer:
      % number of elt added at the end of line (12345678)
    end
    if (curlvl==ilvl)
      CurRecInd(ilvl)=CurRecInd(ilvl)+1;
    elseif (curlvl==ilvl+1)
      MaxLvl=max(curlvl,MaxLvl);
      ilvl=curlvl;
      CurRecInd(ilvl)=1;
    elseif curlvl>0 & curtype ~=99
      for i=curlvl+1:ilvl
        CurRecInd(i)=1;
      end
      CurRecInd(curlvl)=CurRecInd(curlvl)+1;
      ilvl=curlvl;
    else
      lnotbypass=0;
    end
    if lnotbypass==1
      strpos='';
      for i=1:ilvl
        strpos=[strpos int2str(CurRecInd(i)+1) ','];
      end
      for i=ilvl+1:MaxLvl
        strpos=[strpos int2str(CurRecInd(i)) ','];
      end
      strpos=[strpos '1' ];
      eval(['RecordPos(' strpos ')=pos;']);
      eval(['RecordNumb(' strpos ')=indrec;']);
      RecordName(indrec,:)=currecname;
      indrec=indrec+1;
    end
    % change position of cursor to next record line
    if (curtype==0)|(curtype==10)
      nskip=0;
    elseif (curtype==1)
      nl=floor((curnum-1)/8)+1; %8 int / line
      nskip=nl+curnum*10; % 10 bytes / int + 1 byte / line
    elseif (curtype==2)
      nl=floor((curnum-1)/5)+1; %5 real / line
      nskip=nl+curnum*16; % 16 bytes / real + 1 byte / line
    elseif (curtype==3)
      nl=floor((curnum-1)/8)+1; %10 char / '         4'
      nl=nl+floor((curnum-1)/20)+1; %20 char / line
      nskip=nl+curnum*14; % 10+4 bytes / char + 1 byte / line
    elseif (curtype==4)
      nl=floor((curnum-1)/4)+1; %4 double / line
      nskip=nl+curnum*20; % 20 bytes / double + 1 byte / line
    end
    fseek(pfin,nskip,'cof');
    Line=fgetl(pfin);
  end
  CurRecInd=ones(1,MaxLvl);
  ilvl=1;
else
  ilvl=0
end
% ASCIILibv4(pfin,ilvl,RecordPos,RecordNumb,RecordName)
