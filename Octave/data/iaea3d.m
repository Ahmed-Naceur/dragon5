% 
% Plot the interpolated flux computed by the Trivac/data/iaea3d.x2m
%   (get file 'AIFLUiaea_33' from results and copy it in 'iaea3d_data' folder)
%

clear

FileName = 'iaea3d_data/AIFLUiaea_33' ;
tit=['iaea3d-DUAL-33'] ;

savefig=1;
ext='eps';

disp(['open file : ' FileName ])
[pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd] = ASCIIOpnv4(FileName);
ista= ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'STATE-VECTOR');
mxi = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'MXI         ');
myi = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'MYI         ');
mzi = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'MZI         ');
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'FLUX        ',1);
flup(:,1)   = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'elt#00000001');
flup(:,2)   = ASCIIGetv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'elt#00000002');
[pfin,ilvl,nbdata,CurRecInd] = ASCIISixv4(pfin,ilvl,RecordPos,RecordNumb,RecordName,CurRecInd,'            ',0);
fclose(pfin);

ng=ista(1); nx=ista(2); ny=ista(3); nz=ista(4);

flui=zeros(nx,ny,nz,ng);
for ig=1:ng
id=0;
for iz=1:nz
for iy=1:ny
for ix=1:nx
id=id+1;
flui(ix,iy,iz,ig)=flup(id,ig);
end
end
end
end

figure
for k=1:20
subplot(4,5,k)
surf(mxi,myi,flui(:,:,k,2))
view(46,53)
axis([0 180 0 180 -1e-6 5e-6])
caxis([-1e-6, 5e-6])
title(num2str(k))
end

if(savefig==1)
  saveas(gcf,tit,ext);
end
