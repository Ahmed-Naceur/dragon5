% ASCIISixv4
% Change the current index to the one of the specified record
%  used for ASCII file created by DRAGON read with ASCIIOpn.m
% AUTHOR   :      R. CHAMBON
% date     :      nov 2009
%
% note:
%   move=0 to reset pfin index
%   move=1 to reset pfin to the line corresponding to 'RecToMoveTo' at THIS
%           level. If the record does not exist nbdata = 0
%   move=2 to reset pfin to FIRST record with level= ilvl-1
%
function [pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,RecToMoveTo,move)
nbdata=0;
len=length(CurRecInd);
% Return to first line
if (move==0 || (move==2 && ilvl==1)|| (move==2 && ilvl==2))
  fseek(pfin, 0, 'bof');
  ilvl=1;
  CurRecInd=ones(1,len);
  % Return to FIRST record with level= ilvl-1
  %   for ilvl=1 or 2 it is the bof,
  %   for ilvl>2 the CurRecInd should not correspond to the bof
elseif move==2
  CurRecInd(ilvl-1)=1;
  ilvl=ilvl-1;
elseif move==1
  tmpRecInd=CurRecInd;
  siz=size(RecordPos);
  MaxLvl=length(siz);
  for i=1:siz(ilvl)-1
    tmpRecInd(ilvl)=i;
    strpos='';
    for i=1:ilvl
      strpos=[strpos int2str(tmpRecInd(i)+1) ','];
    end
    for i=ilvl+1:MaxLvl
      strpos=[strpos int2str(tmpRecInd(i)) ','];
    end
    strpos=[strpos '1' ];
    lgn=eval(['RecordNumb(' strpos ')']);
    if lgn~=0
      if(RecordName(lgn,:)==RecToMoveTo)
        CurRecInd=tmpRecInd;
        ilvl=ilvl+1;
        nbdata=-1;
        break;
      end
    else
      break;
    end
  end
  if nbdata==0
    disp(['WARNING: record ''' RecToMoveTo ''' not found'])
    %     fseek(pfin,position,'bof');
  end
else
  disp('WRONG move number');
end
