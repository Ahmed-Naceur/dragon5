% Validation of ASCIIOpnv4, ASCIILenv4, ASCIIGetv4 and ASCIISixv4 functions
% subroutine: ASCIIDmpv4
%-----------------
clear

FileName = 'TemplateVerifv4_data/fmap' ;

disp(['open file : ' FileName ])
[pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd] = ASCIIOpnv4(FileName);

disp('Directory Content')
ASCIILibv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd)

disp('Complete data Structure Content')
ASCIIDmpv4(pfin,ilvl,RecordPos,RecordNumb,RecordName)

disp('get real data level 1')
VOLUME = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'VOLUME      ');
disp('2.85292793E+04=?=VOLUME(1)')
VOLUME(1)
disp('1=?=leng')
leng=ASCIILenv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'VOLUME      ')
% 
disp('get integer data level 1')
STV = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'STATE-VECTOR');
disp('520=?=STV(3)')
STV(3)
disp('40=?=leng')
leng=ASCIILenv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'STATE-VECTOR')
% 
disp('get character data level 1')
XNAME = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'XNAME       ');
disp('-   1   2   3=?=XNAME(1:13)')
XNAME(1:13)
disp('28*4=112=?=leng')
leng=ASCIILenv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'XNAME       ')
% 
disp('move to a directory')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'GEOMAP      ',1);
CurRecInd
disp('get real data level 2')
MESHX = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'MESHX       ');
disp('6.00000000E+01=?=MESHX(2)')
MESHX(2)
disp('get integer data level 2')
STV2 = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'STATE-VECTOR');
disp('9408=?=STV2(6)')
STV2(6)
disp('move to the root')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,' ',0);
CurRecInd
% 
disp('move to a directory which is not existing lvl=1')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'NOT-HERE    ',1);
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,' ',0);
% 
disp('move to a directory which is not existing lvl=2')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'GEOMAP      ',1);
CurRecInd
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'NOT-HERE    ',1);
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,' ',0);
CurRecInd

disp('move up from a directory level 2')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'GEOMAP      ',1);
CurRecInd
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,' ',2);
CurRecInd

disp('move to a directory level 1')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'FUEL        ',1);
CurRecInd
disp('move to a directory level 2')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'elt#00000002',1);
CurRecInd
disp('move up from a directory level 3')
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,' ',2);
CurRecInd


fclose(pfin);



