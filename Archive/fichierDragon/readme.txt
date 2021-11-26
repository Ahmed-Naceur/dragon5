Voici les fichiers de ma version officieuse de DRAGON5.

Je liste ici les modifications qui sont apportées à chaque fichier et leur qualité sur 
l'échelle (--:médiocre,-,+,++:Implémentation quasi définitive).

---------------------------------------------------------------------------------------------
Dans le dossier src/ de DRAGON/ :
---------------------------------------------------------------------------------------------
 
FMAC01.f :    (--) J'ai procédé à une petite permutation des lettres associés à proton et 
                   positron pour permettre la lecture des fichiers avec électron. Temporaire
                   seulement
              (-)  J'ai appliqué une correction pour la lecture des frontières énergétiques
                   et des stopping power associés.
-- Naceur, Ahmed: Taken

ASMDRV.f :    (++) Transfert de E_{cutoff} dans une macrolib
              (++) Transfert du terme T=alpha/2 dans une macrolib seulement si au moins une
                   valeur non-nulle
- Naceur, Ahmed: Taken 

SNFLUX.f :    (++) Récupération des valeurs de E_{cutoff} et de T=alpha/2
              (+) Calcul de l'énergie sous le cutoff (++ pour le 1D, j'ai pas encore testé en 2D/3D)
              (++) Ajout de DCUTOFF comme argument de SNFLUX.f
              (+)  Transfert de alpha à la subroutine SNFP1P.f
- Naceur, Ahmed: Taken

SNFP1P.f :    (+)  Traitement direct du terme alpha
              (--) Source monodirectionnelle implémenté directement dans cette sous-routine,
                   à décommenter pour utiliser. Problème oscillatoire avec celle-ci à bas
                   nombre de groupe. Normalisation non fonctionnelle, donc manuelle seulement.
- Naceur, Ahmed: Taken

SNFT1P.f :    (--) Source monodirectionnelle implémenté directement dans cette sous-routine,
                   à décommenter pour utiliser.Normalisation non fonctionnelle, donc manuelle
                   seulement.
- Naceur, Ahmed: Taken

SNF.f :       (++) Ajout de DCUTOFF comme argument
- Naceur, Ahmed: Taken

DOORFVR.f :   (++) Ajout de DCUTOFF comme argument
- Naceur, Ahmed: Taken

FLU2DR.f :    (++) Transfert de DCUTOFF dans une macrolib
- Naceur, Ahmed: Taken

HEADRV.f :    (++) Ajout du cutoff au calcul de l'énergie de déposition
- Naceur, Ahmed: Taken

---------------------------------------------------------------------------------------------
Dans le dossier src/ de Trivac/ :
---------------------------------------------------------------------------------------------

OUTDRV.f :    (++) Transfert de DCUTOFF d'une macrolib à une autre
= Naceur, Ahmed: Do you want to say "OUTAUX.f" 
                 since OUTDRV.f is not provided in the current folder.


- Naceur, Ahmed: Taken OUTAUX.f




------
NOTES, FEEDBACK: Naceur, Ahmed
------
SNFP1P.f warning: unused ZCODE
