#!/bin/env python
""" generation of module dependances for a xxx.f90 file """

import string,os
from subprocess import Popen, PIPE

def noComment(mot):
    """ remove comment characters in a word """
    ind= string.find(mot,'!')
    if ind != -1:
        mot= mot[:ind]
    return mot

class md9:
    """ a file, its module and its dependances """
    def __init__(self,fic):
        self.name= fic
        self.used= []
        #find all module dependances
        listeLignes= open(self.name,'r').readlines()
        self.module= "dummy"
        for ligne in listeLignes:
            ligne= string.lower(ligne)
            if string.count(ligne,'use'):
                ligne= string.replace(ligne,',',' ')
                listeMots= string.split(ligne)
                if listeMots[0] == 'use':
                    mod= noComment(listeMots[1])
                    if mod not in self.used and mod != 'intrinsic':
                        self.used.append(mod)
            elif string.count(ligne,'module'):
                listeMots= string.split(ligne)
                if listeMots[0] == 'module' and \
                   noComment(listeMots[1]) != 'procedure':
                    self.module= listeMots[1]
    def utilise(self,other):
        return (other.module in self.used)

def listdir_shell(path):
    p = Popen(('ls', path), shell=False, stdout=PIPE, close_fds=True)
    return [path.rstrip('\n') for path in p.stdout.readlines()]

listeFichier= filter(lambda x: string.find(x,'.f90')!=-1,listdir_shell(os.getcwd()))
listeClasse= map(md9,listeFichier)
#sort the file list
listeRes= []
while listeClasse:
    lgav= len(listeClasse)
    for fic in listeClasse:
        if len(filter(lambda x,y=fic:y.utilise(x),listeClasse)) == 0:
            listeRes.append(fic)
            listeClasse.remove(fic)
    lgap= len(listeClasse)
    if lgav == lgap and lgap != 0:
        raise RuntimeError("make_depend: cross-references found")
#result output
resstr=''
for fic in listeRes:
    resstr=resstr+fic.name+' '
print resstr
