
/*****************************************/
/*             CLE-2000 API              */
/*   AUTHOR OF FORTRAN VERSION: R. Roy   */
/*     AUTHOR: A. Hebert ; 31/07/10      */
/*****************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cle2000.h"
#define maxent 1000  /* maximum number of module arguments */

int_32 cle2000_c(int_32 ilevel,
                 int_32 (*dummod)(char *, int_32, char (*)[13], int_32 *, int_32 *, lcm **, char (*)[73]),
                 char *filenm, int_32 iprint, lifo *my_param)
{
   char *nomsub = "cle2000_c";
   int_32 ret_val = 0;
   static char *cdclkw[] = {"PROCEDURE", "MODULE", "LINKED_LIST", "XSM_FILE",
                            "SEQ_BINARY", "SEQ_ASCII", "DIR_ACCESS", "HDF5_FILE",
                            "PARAMETER"};
   int ldatav, lobnew;
   kdi_file *iKDI, *icobj;
   FILE *icinp = NULL, *icout = NULL, *icfile;
   lifo *my_iptrun;
   char hwrite[73] = " ";
   int_32 nusec2, jrecin;
   lcm *kparam = NULL;
   int_32 iretcd, iloop1, iparam, jparam, lparam;
   int_32 jdispe, nentry, nmodul, ilogin, ityp, nitma;
   float_32 flott;
   double_64 dflot;
   char hentry[maxent][13];
   int_32 ientry[maxent], jentry[maxent];
   lcm *kentry[maxent];
   char cmodul[13], hparam[73], text[73], cproce[73] = " ";
   char filinp[73], filobj[73], filout[73], filkdi[73];
   
/* COMPILE MAIN INPUT INTO OBJECT FILE */
   if (strcmp(filenm, " ") == 0) {
      icinp = stdin;
      icfile = NULL;
      icout = stdout;
      strcpy(filobj, "_DUMMY");
   } else {
      sprintf(filinp, "%s.c2m",filenm);
      sprintf(filobj, "%s.o2m",filenm);
      sprintf(filout, "%s.l2m",filenm);
      icfile = fopen(filobj, "r");
   }
   if (icfile == NULL) {
/*    OPEN SOURCE FILE '.c2m' */
      if (icinp != stdin) {
         icinp = fopen(filinp, "r");
         if (icinp == NULL) goto L9003;
      }
      if (icout != stdout) {
         icout = fopen(filout, "w+");
         if (icout == NULL) goto L9005;
      }

/*    CREATE OBJECT FILE SUFFIX IS '.o2m' */
      icobj = kdiop_c(filobj, 0);
      if (icobj == NULL) goto L9001;

/*    COMPILE NEW FUNCTION */
      iretcd = clepil(icinp, icout, icobj, clecst);
      if (iretcd != 0) {
         printf("%s: COMPILING _MAIN.c2m FILE (ERROR CODE) IRC=%d\n", nomsub,(int)iretcd);
         goto L666;
      }

/*    ADD OBJECTS/MODULES TO OBJECT FILE */
      iretcd = objpil(icobj, icout, 0);
      if (iretcd != 0) {
         printf("%s: BAD OBJECTS _MAIN.c2m FILE (ERROR CODE) IRC=%d\n", nomsub,(int)iretcd);
         goto L666;
      }

/*    CLOSE & KEEP SOURCE & OUTPUT FILES */
      if (icout != stdout) {
         iretcd = fclose(icout);
         if (iretcd != 0) goto L9006;
      }
      if (icinp != stdin) {
         iretcd = fclose(icinp);
         if (iretcd != 0) goto L9004;
      }
   } else {
      iretcd = fclose(icfile);
      if (iretcd != 0) goto L9002;
      icobj = kdiop_c(filobj, 1);
      if (icobj->fd == NULL) {
         printf("%s: DID YOU FORGET TO COMPILE *%s*?\n", nomsub, filobj);
         goto L9001;
      }
   }

/* NOW, MAKE A COPY OF OBJECT FILE */
   if (strcmp(filenm," ") == 0) {
     sprintf(filkdi,"_main%.3d", (int)ilevel);
   } else {
     sprintf(filkdi,"_%s%.3d", filenm,(int)ilevel);
   }
   iKDI = kdiop_c(filkdi, 0);
   if (iKDI == NULL) goto L9007;
   iretcd = clecop(icobj, iKDI);
   if (iretcd != 0) {
      printf("%s: COPYING PREVIOUS FILE  *%s* IRC=%d\n", nomsub, cproce, (int)iretcd);
      goto L666;
   }

/* CLOSE AND KEEP ORIGINAL OBJECT FILE */
   iretcd = kdicl_c(icobj, 1);
   if (iretcd != 0) goto L9002;
   if (strcmp(filobj, "_DUMMY") == 0) {
      iretcd = remove(filobj);
      if (iretcd != 0) goto L9006;
   }

   if (iprint > 0) printf("%s: STARTING EXECUTION ON _MAIN.o2m FILE\n", nomsub);
   redopn_c(iKDI, stdout, hwrite, 0);
   cleopn(&my_iptrun);
   nusec2 = 0;

L10:
/* GET SENTENCE */
   jdispe = 0, nentry = 0, nmodul = 0;
   redget_c(&ityp, &nitma, &flott, text, &dflot);
/* TREAT FIRST WORD */
   if (ityp == 3) {
      ilogin = 0;
      for (iloop1 = 0; iloop1 < 9; ++iloop1) {
         if (strcmp(cdclkw[iloop1], text) == 0) ilogin = iloop1+1;
      }
/*    OUTSIDE THE DATA SECTION ( *HERE*  ::  ... ) */
L30:
      if (strcmp(text, ":=") == 0) {
         jdispe = 2;
      } else if (strcmp(text, "::") == 0) {
/*       FOR PROCEDURE/MODULE WITHOUT DATA, BRANCH NOW */
         ldatav = 1;
         goto L40;
      } else if (strcmp(text, ";") == 0) {
/*       FOR PROCEDURE/MODULE WITH    DATA, BRANCH NOW */
         ldatav = 0;
         goto L40;
      } else {
         lifo_node *my_node;
         my_node = clenode(&my_iptrun, text);
         if (my_node == NULL) {
            my_node = (lifo_node *) malloc(sizeof(lifo_node));
            strcpy(my_node->name, text);
            strcpy(my_node->OSname, " ");
            clepush(&my_iptrun, my_node);
         }
         if (ilogin != 0) {
            if (strcmp(text, cdclkw[ilogin-1]) == 0) {
/*             DECLARATION ITSELF IS A MODULE */
               iparam = 2;
            } else {
/*             TYPE IS SET TO VALUE < 0 (UNDEFINED) */
               iparam = -ilogin;
            }
            strcpy(hparam, text);
            my_node->type = iparam;
            jparam = -1; my_node->access = jparam;
            kparam = NULL;  my_node->value.mylcm = kparam;
            lparam = 0;  my_node->lparam = lparam;
            if (iparam < 0) strcpy(my_node->OSname, hparam);
         } else {
            if(my_node == NULL) {
               printf("%s: NODE DOES NOT EXIST\n", nomsub);
               goto L666;
            }
            iparam = my_node->type;
            kparam = my_node->value.mylcm;
            lparam = my_node->lparam;
            strcpy(hparam, my_node->OSname);
         }
         if (nmodul == 0 && (abs(iparam) == 1 || abs(iparam) == 2)) {
/*          ONCE MODULE/PROCEDURE FOUND, RESET *JDISPE=2, (READ-ONLY MODE) */
            strcpy(cmodul, text);
            jdispe = 2;
            if (abs(iparam) == 2) {
               nmodul = 1;
            } else {
               nmodul = -1;
               strcpy(cproce, hparam);
               if (iparam == -1) {
                  printf("%s: FILE *%s* DOES NOT EXIST\n", nomsub, hparam);
                  goto L666;
               }
            }
         } else {
            lobnew = 1;
            if (nentry != 0) {
               for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
                  if (strcmp(hentry[iloop1], text) == 0) {
                     if (jentry[iloop1] == 0) {
/*                      OBJECT GOES TO (MODIFICATION) MODE */
                        jentry[iloop1] = 1;
                     } else {
                        printf("%s: INCONSISTENT CALL (text=%s)\n", nomsub,text);
                        goto L666;
                     }
                     lobnew = 0;
                  }
               }
            }
            if (lobnew) {
               ++(nentry);
               if (nentry > maxent) {
                  printf("%s: TOO MANY OBJECTS\n", nomsub);
                  goto L666;
               }
               strcpy(hentry[nentry-1], text);
               jentry[nentry-1] = jdispe;
            }
         }
      }
      redget_c(&ityp, &nitma, &flott, text, &dflot);
      goto L30;
L40:
      if (nmodul == 0) {
         if (nentry == 0) {
            strcpy(cmodul, "IOX:");
         } else {
            strcpy(cmodul, "EQU:");
         }
         nmodul = 1;
      }
      if (nmodul == 1) {
/*       FOR MODULES */
/*       IF NOT LDATAV DISCONNECT READER */
         if (!ldatav) redcls_c(&iKDI, &icout, hwrite, &jrecin);
         if (ilogin == 0) {
            if (strcmp(cmodul, "IOX:") == 0) {
               int_32 minput = 0;
               if (nentry != 0) {
                  printf("%s: MODULE *IOX:* WITH INVALID PARAMETERS\n", nomsub);
                  goto L666;
               }
               drviox(my_param, minput, &nusec2);
            } else {
               char hparam_c[maxent][73];
               for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
                  lifo_node *my_node;
                  my_node = clenode(&my_iptrun, hentry[iloop1]);
                  if (my_node == NULL) {
                     printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
                     goto L666;
                  }
                  iparam = my_node->type;
                  jparam = my_node->access;
                  if (iparam == 3) kparam = my_node->value.mylcm;
                  if (iparam == 7) lparam = my_node->lparam;
#if defined(HDF5_LIB)
                  if (iparam == 8) kparam = (lcm*)my_node->value.myhdf5;
#endif
                  strcpy(hparam, my_node->OSname);
                  strcpy(hparam_c[iloop1], hparam);
/*                CONSISTENCY TESTS */
                  if (jentry[iloop1] == 0) {
                     if (jparam == 1 || jparam == 2 || iparam >= 0) {
                        printf("%s: %s *%s* ALREADY EXISTS\n", nomsub, cdclkw[abs(iparam)-1], hentry[iloop1]);
                        ret_val = -2;
                        goto L666;
                     }
                     iparam = -iparam;
                     my_node->type = iparam;
                  } else if (jentry[iloop1] == 1) {
                     if (jparam == 2) {
                        printf("%s: %s *%s* IS PROTECTED\n", nomsub, cdclkw[abs(iparam)-1], hentry[iloop1]);
                        ret_val = -2;
                        goto L666;
                     } else if (iparam <= 0) {
/*                      ALLOW ACCESS TO ANY FILE AT MAIN LEVEL */
                        if (ilevel == 1 && (iparam == -4 || iparam == -5 || iparam == -6 || iparam == -7 || iparam == -8)) {
                           iparam = -iparam;
                           my_node->type = iparam;
                        } else {
                           printf("%s: %s *%s* IS NOT DEFINED(1)\n", nomsub, cdclkw[-iparam-1], hentry[iloop1]);
                           ret_val = -2;
                           goto L666;
                        }
                     }
                  } else if (jentry[iloop1] == 2) {
                     if (iparam <= 0) {
/*                      ALLOW ACCESS TO ANY FILE AT MAIN LEVEL */
                        if (ilevel == 1 && (iparam == -4 || iparam == -5 || iparam == -6 || iparam == -7 || iparam == -8)) {
                           iparam = -iparam;
                           my_node->type = iparam;
                        } else {
                           printf("%s: %s *%s* IS NOT DEFINED(2)\n", nomsub, cdclkw[-iparam-1], hentry[iloop1]);
                           ret_val = -2;
                           goto L666;
                        }
                     }
                  } else {
                     printf("%s: INVALID JENTRY=%d\n", nomsub,(int)jentry[iloop1]);
                     goto L666;
                  }
                  jdispe = jentry[iloop1];
                  if (iprint > 1) {
                     if (jdispe == 0) {
                        printf("%s: OPEN %s *%s* IN CREATION MODE\n", nomsub, cdclkw[iparam-1], hparam);
                     } else if (jdispe == 1) {
                        printf("%s: OPEN %s *%s* IN MODIFICATION MODE\n", nomsub, cdclkw[iparam-1], hparam);
                     } else if (jdispe == 2) {
                        printf("%s: OPEN %s *%s* IN READ/ONLY MODE\n", nomsub, cdclkw[iparam-1], hparam);
                     }
                  }
                  if (iparam == 3 || iparam == 4) {
                     if (jdispe > 0) kentry[iloop1] = kparam;
                     lcmop_c(&kentry[iloop1], hparam, jdispe, iparam-2, 0);
                     kparam = kentry[iloop1];
                  } else if (iparam == 5 || iparam == 6 || iparam == 7) {
                     kentry[iloop1] = NULL;
#if defined(HDF5_LIB)
                  } else if (iparam == 8) {
                     if (jdispe > 0) kentry[iloop1] = kparam;
                     hid_t myhdf5 = 0;
                     if (jdispe == 0) {
                       myhdf5 = H5Fcreate(hparam, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
                       if (iprint > 1) {
                         printf("%s: create HDF5 file at address=%lld\n",nomsub, (long long int)myhdf5);
                       }
                     } else if (jdispe == 1) {
                       myhdf5 = H5Fopen(hparam, H5F_ACC_RDWR, H5P_DEFAULT);
                       if (iprint > 1) {
                         printf("%s: open HDF5 file in read-write mode at address=%lld\n",nomsub, (long long int)myhdf5);
                       }
                     } else if (jdispe == 2) {
                       myhdf5 = H5Fopen(hparam, H5F_ACC_RDONLY, H5P_DEFAULT);
                       if (iprint > 1) {
                         printf("%s: open HDF5 file in read-only mode at address=%lld\n",nomsub, (long long int)myhdf5);
                       }
                     }
                     if (myhdf5 < 0) {
                       printf("%s: H5Fopen failure on HDF5 file '%s'.\n",nomsub,hparam);
                       goto L666;
                     }
                     kentry[iloop1] = (lcm*)myhdf5;
                     kparam = kentry[iloop1];
#endif
                  } else {
                     printf("%s: USE %s *%s* IS IMPOSSIBLE. INVALID IPARAM (%d)\n", nomsub,
                            cdclkw[iparam-1], hparam, (int)iparam);
                     goto L666;
                  }
                  ientry[iloop1] = iparam - 2;
                  if (jdispe == 0) {
                     my_node->value.mylcm = kparam;
                  } else if (strcmp(cmodul, "DELETE:") == 0) {
                     if (abs(jparam) != 1) {
                        printf("%s: KIL %s *%s* IS IMPOSSIBLE. INVALID DELETE\n", nomsub,
                               cdclkw[iparam-1], hparam);
                        goto L666;
                     }
                     kparam = 0;
                     iparam = -iparam;
                     my_node->type = iparam;
                     if (iparam == 3) my_node->value.mylcm = kparam;
                  }
               }

/*             CALLING MODULES */
               jdispe = 1;
               if (strcmp(cmodul, "END:") == 0) {
                  if (ldatav) {
                     printf("%s: *END:* HAS NO DATA\n", nomsub);
                     goto L666;
                  } else if (nentry != 0) {
                     printf("%s: *END:* HAS NO OBJECT\n", nomsub);
                     goto L666;
                  }
                  iretcd = 0;
               } else if (strcmp(cmodul, "DELETE:") == 0) {
/*                STANDARD DELETE MODULE (SEE BELOW). */
                  jdispe = 2;
                  iretcd = 0;
               } else if (dummod != NULL) {
/*                CALLING ANOTHER STANDARD UTILITY MODULE in ANSI-C. */
                  fflush(stdout);
                  iretcd = (*dummod)(cmodul, nentry, hentry, ientry, jentry, kentry, hparam_c);
               } else if (dummod == NULL) {
                  printf("%s: MODULE *%s* NOT FOUND; DUMMOD NOT SET\n", nomsub, cmodul);
                  goto L666;
               }
               if (iretcd != 0) {
                  printf("%s: MODULE *%s* NOT FOUND\n", nomsub, cmodul);
                  goto L666;
               }

/*             CLOSE EVERYTHING */
               for (iloop1 = 0; iloop1 < nentry; ++iloop1) {
                  iparam = ientry[iloop1];
                  if (iprint > 1) {
                     lifo_node *my_node;
                     my_node = clenode(&my_iptrun, hentry[iloop1]);
                     if (my_node == NULL) {
                        printf("%s: UNABLE TO FIND NODE FOR %s\n", nomsub, hentry[iloop1]);
                        goto L666;
                     }
                     strcpy(hparam, my_node->OSname);
                     if (jdispe == 1) {
                        printf("%s: CLS %s *%s*\n", nomsub, cdclkw[iparam-1], hparam);
                     } else {
                        printf("%s: KIL %s *%s* (AS IF IT NEVER EXISTED)\n", nomsub,
                               cdclkw[iparam-1], hparam);
                     }
                  }
                  if (iparam == 1 || iparam == 2) {
                     lcmcl_c(&kentry[iloop1], jdispe);
                  } else if (iparam == 3 || iparam == 4 || iparam == 5) {
                     if (jdispe == 2) {
                        iretcd = remove(hentry[iloop1]);
                        if (iretcd != 0) {
                           printf("%s: REMOVE FAILURE. file=*%s*\n", nomsub, hentry[iloop1]);
                           goto L666;
                        }
                        kparam = 0;
                     }
#if defined(HDF5_LIB)
                  } else if (iparam == 6) {
                     iretcd = H5Fclose((hid_t)kentry[iloop1]);
                     if (iretcd != 0) {
                        printf("%s: HDF5 CLOSE FAILURE. file=*%s* iretcd=%d\n", nomsub, hentry[iloop1],
                        iretcd);
                        goto L666;
                     }
                     if (jdispe == 2) {
                        iretcd = remove(hentry[iloop1]);
                        if (iretcd != 0) {
                           printf("%s: REMOVE FAILURE. file=*%s*\n", nomsub, hentry[iloop1]);
                           goto L666;
                        }
                        kparam = 0;
                     }
#endif
                  } else {
                     printf("%s: UNABLE TO CLOSE *%s*\n", nomsub, hentry[iloop1]);
                     goto L666;
                  }
               }
            }
         } else {
/*          CHARGE ENTRIES FOR DECLARATION MODULES */
            if (strcmp(cmodul, "PARAMETER") == 0) {
/*             *PARAMETER   * DECLARATION MODULE */
               int_32 minput = 0;
               iretcd = kdrprm(&my_iptrun, &my_param, minput, nentry, jentry, hentry);
               if (iretcd != 0) {
                  printf("%s: PROBLEM ACCEPTING PARAMETERS IER=%d\n", nomsub, (int)iretcd);
                  goto L666;
               }
               nusec2 = nentry;
            } else if (strcmp(cmodul, "PROCEDURE") == 0) {
/*             *PROCEDURE   * DECLARATION MODULE */
               iretcd = kdrdpr(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "MODULE") == 0) {
/*             *MODULE      * DECLARATION MODULE */
               iretcd = kdrdmd(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "LINKED_LIST") == 0) {
/*             *LINKED_LIST * DECLARATION MODULE */
               iretcd = kdrdll(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "XSM_FILE") == 0) {
/*             *XSM_FILE    * DECLARATION MODULE */
               iretcd = kdrdxf(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "SEQ_BINARY") == 0) {
/*             *SEQ_BINARY  * DECLARATION MODULE */
               iretcd = kdrdsb(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "SEQ_ASCII") == 0) {
/*             *SEQ_ASCII   * DECLARATION MODULE */
               iretcd = kdrdsa(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "DIR_ACCESS") == 0) {
/*             *DIR_ACCESS  * DECLARATION MODULE */
               iretcd = kdrdda(&my_iptrun, nentry, hentry);
            } else if (strcmp(cmodul, "HDF5_FILE") == 0) {
/*             *HDF5_FILE  * DECLARATION MODULE */
               iretcd = kdrdh5(&my_iptrun, nentry, hentry);
            } else {
/*             OTHERWISE, DECLARATION MODULE IS NOT AVAILABLE */
               printf("%s: DECLARATION MODULE *%s* NOT AVAILABLE IN THIS CODE\n", nomsub, cmodul);
               goto L666;
            }
            if (iretcd != 0) {
               printf("%s: PROBLEM WITH MODULE *%s*\n", nomsub, cmodul);
               goto L666;
            }
         }
         if (iprint > 0) printf("%s: END OF MODULE *%s*\n", nomsub, cmodul);
         if (!ldatav && strcmp(cmodul, "END:") != 0) {
/*          RECONNECT READER IF DISCONNECTED OUTSIDE END: */
            redopn_c(iKDI, stdout, hwrite, jrecin);
         }
      } else {
/*       FOR PROCEDURES */
         int_32 minput;
         lifo *my_param_daughter;

         minput = -1;
         iretcd = kdrprm(&my_iptrun, &my_param_daughter, minput, nentry, jentry, hentry);
         if (iretcd != 0) {
            printf("%s: PROBLEM PASSING PARAMETERS\n", nomsub);
            goto L666;
         } else if (my_param_daughter == NULL) {
            printf("%s: MISSING call_daughter SUB-STRUCTURE\n", nomsub);
            goto L666;
         }

/*       HERE, ONE STEP UP */
         ++(ilevel);
         if (ldatav) {
/*          FOR PROCEDURES, READ DATA SECTION ( ...   :: *HERE* ) */
            drviox(my_param_daughter, minput, &nusec2);
         }
         if (iprint > 0) printf("%s: BEG PROCEDURE *%s* (LEVEL: STEP UP)\n", nomsub, cmodul);

/*       CLOSE THE READER AT CURRENT LEVEL */
         redcls_c(&iKDI, &icout, hwrite, &jrecin);
         
/*       RECURSIVE CALL TO cle2000_c.c */
         iretcd = cle2000_c(ilevel, dummod, cmodul, iprint, my_param_daughter);
         if (iretcd != 0) return iretcd;

/*       REOPEN THE READER AT CURRENT LEVEL */
         redopn_c(iKDI, stdout, hwrite, jrecin);
         
/*       RETURN BACK TO PREVIOUS LEVEL AT SELECTED RECORD */
         minput = 1;
         iretcd = kdrprm(&my_iptrun, &my_param_daughter, minput, nentry, jentry, hentry);
         if (iretcd != 0) {
            printf("%s: PROBLEM RETURNING PARAMETERS\n", nomsub);
            goto L666;
         }

/*       RECOVERING OUTPUT IN DATA FIELD */
         drviox(my_param_daughter, minput, &nusec2);

/*       CLEANING NON-DUMMY OBJECTS */
         iretcd = kdrcln(my_param_daughter, iprint);
         if (iretcd != 0) {
            printf("%s: PROBLEM CLEANING NON-DUMMY OBJECTS\n", nomsub);
            goto L666;
         }

/*       HERE, ONE STEP DOWN */
         --(ilevel);
         if (iprint > 0) {
            if (ityp == 9) {
               printf("%s: END PROCEDURE  NO MORE DATA (LEVEL: STEP DOWN). LEV=%d\n", nomsub, (int)ilevel);
            } else {
               printf("%s: END PROCEDURE  RETURN       (LEVEL: STEP DOWN). LEV=%d\n", nomsub, (int)ilevel);
            }
         }
      }
   } else if (ityp == 9 || ityp == 10) {
      goto L100;
   } else {
      printf("%s: INVALID TYPE\n", nomsub);
      goto L666;
   }
   goto L10;

L100:
/* RECOVER NEW OBJECTS IN LIFO STACK */
   if (my_param != NULL) {
      for (iloop1 = 0; iloop1 < my_param->nup; ++iloop1) {
         char dparam[13];
         lifo_node *my_node, *my_node_daughter;
         my_node = clepos(&my_param, iloop1);
         if ((my_node->type <= -10) || (my_node->type > 0)) continue;
         strcpy(dparam, my_node->name_daughter);
         my_node_daughter = clenode(&my_iptrun, dparam);
         if (my_node_daughter == NULL) {
            printf("%s: UNABLE TO FIND NODE FOR %s<-->%s at position %d\n", nomsub, dparam, my_node->name, (int)iloop1);
            goto L666;
      }
         if (my_node_daughter->type != -my_node->type) {
            printf("%s: INCONSISTENT TYPE IN NODES %s<-->%s at position %d\n", nomsub, dparam, my_node->name, (int)iloop1);
            goto L666;
      }
         if (my_node->type == -3) my_node->value.mylcm = my_node_daughter->value.mylcm;
         my_node->type = -my_node->type;
      }
   }

/* DESTROY THE LIFO STACK */
   iretcd = kdrcln(my_iptrun, iprint);
   if (iretcd != 0) {
      printf("%s: PROBLEM CLEANING NON-DUMMY OBJECTS\n", nomsub);
      goto L666;
   }

/* CLOSE AND DESTROY MAIN OBJECT FILE */
   iretcd = kdicl_c(iKDI, 2);
   if (iretcd != 0) goto L9008;
   if ((iprint > 0) && (ret_val == 0)) {
      printf("%s: SUCCESSFUL EXECUTION AT LEVEL %d\n", nomsub, (int)ilevel);
   }
   return ret_val;
L666:
   return 666;
L9001:
   printf("%s: ERROR WHEN OPENING OBJECT FILE *%s*\n", nomsub, filenm);
   goto L666;
L9002:
   printf("%s: ERROR WHEN CLOSING OBJECT FILE *%s*\n", nomsub, filenm);
   goto L666;
L9003:
   printf("%s: ERROR WHEN OPENING SOURCE FILE *.c2m\n", nomsub);
   goto L666;
L9004:
   printf("%s: ERROR WHEN CLOSING SOURCE FILE *.c2m\n", nomsub);
   goto L666;
L9005:
   printf("%s: ERROR WHEN OPENING OUTPUT FILE *.l2m\n", nomsub);
   goto L666;
L9006:
   printf("%s: ERROR WHEN CLOSING OUTPUT FILE *.l2m; IRC=%d\n", nomsub, (int)iretcd);
   goto L666;
L9007:
   if (strcmp(filenm," ") == 0) {
     printf("%s: ERROR WHEN OPENING OBJECT FILE _main%.3d\n", nomsub, (int)ilevel);
   } else {
     printf("%s: ERROR WHEN OPENING OBJECT FILE _%s%.3d\n", nomsub, filenm, (int)ilevel);
   }
   goto L666;
L9008:
   if (strcmp(filenm," ") == 0) {
     printf("%s: ERROR WHEN CLOSING OBJECT FILE _main%.3d\n", nomsub, (int)ilevel);
   } else {
     printf("%s: ERROR WHEN CLOSING OBJECT FILE _%s%.3d\n", nomsub, filenm, (int)ilevel);
   }
   goto L666;
}

int_32 kdrcln(lifo *my_iptrun, int_32 iprint)
{
   char *nomsub = "kdrcln";
   int_32 ret_val = 0;
   char *cdclkw[] = {"PROCEDURE", "MODULE", "LINKED_LIST", "XSM_FILE",
                     "SEQ_BINARY", "SEQ_ASCII", "DIR_ACCESS", "PARAMETER"};

   while (my_iptrun->nup > 0) {
      lifo_node *my_node;
      int_32 iparam,jparam;

      my_node = clepop(&my_iptrun);
      if (my_node == NULL) {
         printf("%s: POP FAILURE IN LIFO STACK.\n",nomsub);
         ret_val = -1;
         goto L20;
      }
      if (iprint > 0) printf("%s: CLEANING FOR NODE %d (*%s*).\n", nomsub, (int)my_iptrun->nup, my_node->name);

      iparam = my_node->type;
      jparam = my_node->access;
      if (iparam <= 10) {
         if (jparam == -1 && (abs(iparam) >= 3 && abs(iparam) < 8)) {
            if (iparam > 0) {
               int_32 jdispe = 1;
/*             DESTROY OBJECT */
               if (iparam == 3) {
                  char recnam[73];
                  strcpy(recnam, my_node->OSname);
                  lcmop_c(&my_node->value.mylcm, recnam, jdispe, iparam-2, 0);
                  lcmcl_c(&my_node->value.mylcm, 2);
               } else if (iparam == 4 || iparam == 5 || iparam == 6 || iparam == 7) {
                  ret_val = (int_32)remove(my_node->OSname);
                  if (ret_val != 0) {
                     printf("%s: CANNOT DESTROY %s FILE %s (irc=%d).\n",nomsub,cdclkw[iparam-1],my_node->name,(int)ret_val);
                     ret_val = -2;
                     goto L20;
                  }
               }
               if (iprint > 0) printf("%s: DEL %s %s (WILL NEVER EXIST ANYMORE).\n",nomsub,cdclkw[iparam-1],my_node->name);
               my_node->type = -iparam;
            } else {
               if (iprint > 1) printf("%s: DEL %s %s (WAS NOT DEFINED ANYWAY).\n",nomsub,cdclkw[-iparam-1],my_node->name);
            }
         }
      }
      free(my_node);
   }
   ret_val = clecls(&my_iptrun);
   if (ret_val != 0) {
      printf("%s: LIFO STACK NOT EMPTY (irc=%d).\n",nomsub, (int)ret_val);
      ret_val = -3;
   }
L20:
   return ret_val;
}
