
/*****************************************/
/*      C++ CLE-2000 OBJECT WRAPPER      */
/*    AUTHOR: A. HEBERT ; 2012/10/07     */
/*****************************************/

/*
Copyright (C) 2012 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include "Cle2000.hxx"

using namespace std;
using namespace ganlib;
static ostringstream hsmg;

ganlib::Cle2000::Cle2000(string sname) {
  this->edit = 0;
  this->procName = sname;
  this->stack.reset();
}
ganlib::Cle2000::Cle2000(string sname, int_32 edit) {
  this->edit = edit;
  this->procName = sname;
  this->stack.reset();
}
ganlib::Cle2000::Cle2000(string sname, LifoPtr jstack) {
  this->edit = 0;
  this->procName = sname;
  this->stack = jstack;
}
ganlib::Cle2000::Cle2000(string sname, int_32 edit, LifoPtr jstack) {
  this->edit = edit;
  this->procName = sname;
  this->stack = jstack;
}

ganlib::Cle2000::~Cle2000() throw(Cle2000Exception) {
  cout << "Cle2000 destructor called." << endl;
}

void ganlib::Cle2000::setLifo(LifoPtr myLifo) {
  this->stack = myLifo;
}

void ganlib::Cle2000::exec() throw(Cle2000Exception) {
  int_32 ier, ilevel = 1;

// close the LCM objects
  for (int_32 ipos=0; ipos<this->stack->getMax(); ++ipos) {
    int_32 myTypeNode = this->stack->typeNode(ipos);
    if ((myTypeNode == 3) || (myTypeNode == 4)) {
      if (this->stack->accessNode(ipos) > 0) {
        ClcmPtr myClcm; this->stack->node(ipos, myClcm);
        try {
          lcm *myLcmStructure = myClcm->extract();
          lifo *myLifo = this->stack->extract();
          lifo_node *myNode = clepos(&myLifo, ipos);
          strcpy(myNode->OSname, myLcmStructure->hname);
          lcmcl_c(&myLcmStructure, 1);
        } catch(...) {
          throw Cle2000Exception("Exception catched by lcmcl_c");
        }
      }
    }
  }

// call the parametrized procedure
  try {
    ier = cle2000_c(ilevel, &donmod, (char *)this->procName.c_str(), this->edit,
                   (lifo*)this->stack->extract());
                   
  } catch(...) {
    throw Cle2000Exception("Exception catched by cle2000_c");
  }
  if (ier != 0) {
    hsmg << "Cle2000: cle2000 failure (" << this->procName << ".c2m). ier=" << ier;
    throw Cle2000Exception(hsmg.str());
  }

// reopen the LCM objects
  for (int ipos=0; ipos<this->stack->getMax(); ++ipos) {
    int_32 myTypeNode = this->stack->typeNode(ipos);
    if ((myTypeNode == 3) || (myTypeNode == 4)) {
      int_32 access = this->stack->accessNode(ipos);
      if (access == 0) access=1;
      ClcmPtr myClcm; this->stack->node(ipos, myClcm);
      try {
        lcm *myLcmStructure = myClcm->extract();
        int_32 myTypeNode = this->stack->typeNode(ipos);
        char *myOSname = (char *)this->stack->OSName(ipos).c_str();
        lcmop_c(&myLcmStructure, myOSname, access, myTypeNode-2, 0);
      } catch(...) {
        throw Cle2000Exception("Exception catched by lcmop_c");
      }
    }
  }
}
