% ASCIIDmpv4
% Print the data-structure architecture
% AUTHOR   :      R. CHAMBON
% date     :      nov 2009

function ASCIIDmpv4(pfin,ilvl,RecordPos,RecordNumb,RecordName)
siz=size(RecordPos);
nel=prod(siz);
MaxLvl=length(siz);
CurRecInd=ones(1,MaxLvl);
CurRecInd(MaxLvl)=0;
spacemax=MaxLvl*3+2;

disp('''RecordPos'' index   |-> ASCII file')
for i=1:nel
    tmplvl=MaxLvl;
    CurRecInd(MaxLvl)=CurRecInd(MaxLvl)+1;
    if CurRecInd(MaxLvl)>siz(MaxLvl)
        while (CurRecInd(tmplvl)>siz(tmplvl))
            CurRecInd(tmplvl)=1;
            tmplvl=tmplvl-1;
            if tmplvl==0
                break;
            end
            CurRecInd(tmplvl)=CurRecInd(tmplvl)+1;
        end
    end
    strpos='';
    for il=1:MaxLvl
        strpos=[strpos int2str(CurRecInd(il)) ','];
    end
    strpos=[strpos '1' ];
    eval(['indrec=RecordNumb(' strpos ');']);
    if indrec>0
        leng=length(strpos);
        for j=leng-1:spacemax
            strpos(j)=' ';
        end
        spacelvl='|';
        for il=2:MaxLvl
            if CurRecInd(il)>1
                spacelvl=[spacelvl int2str(il)];
            end
        end
        disp([strpos spacelvl '-> '  RecordName(indrec,:)]);
    elseif indrec==0
        bypass=CurRecInd(tmplvl)-siz(tmplvl)+1;
        byplvl=1;
        if tmplvl+1<=MaxLvl
            for il=tmplvl+1:MaxLvl
                byplvl=byplvl*siz(il);
            end
        end
        bypasstot=bypass*byplvl-1;
        i=i+bypasstot;
        CurRecInd(tmplvl)=1;
        tmplvl=tmplvl-1;
        if tmplvl==0
            break;
        end
        CurRecInd(tmplvl)=CurRecInd(tmplvl)+1;
        if CurRecInd(tmplvl)>siz(tmplvl)
            while (CurRecInd(tmplvl)>siz(tmplvl))
                CurRecInd(tmplvl)=1;
                tmplvl=tmplvl-1;
                if tmplvl==0
                    break;
                end
                CurRecInd(tmplvl)=CurRecInd(tmplvl)+1;
            end
        end
        CurRecInd(MaxLvl)=0;
    end
    if tmplvl==0
        break;
    end
end
