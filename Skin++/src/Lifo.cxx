
/*****************************************/
/*    C++ LIFO STACK OBJECT WRAPPER      */
/*    AUTHOR: A. HEBERT ; 2012/10/07     */
/*****************************************/

/*
Copyright (C) 2012 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include "Lifo.hxx"

using namespace std;
using namespace ganlib;
static ostringstream hsmg;

ganlib::Lifo::Lifo() throw(LifoException) {
  cout << "New Lifo object constructed.'" << endl;
  try {
    cleopn(&(this->addr));
  } catch(...) {
    throw LifoException("Exception catched by cleopn");
  }
  this->global_list = new vector<ClcmPtr>;
}

ganlib::Lifo::~Lifo() throw(LifoException) {
  cout << "Lifo destructor called (" << this->getMax() << " nodes remaining)." << endl;
  int_32 nitma = clecls(&(this->addr));
  if(nitma != 0) throw LifoException("Lifo destruction failure: Lifo not empty");
  delete this->global_list;
}

//----------------------------------------------------------------------
//--                               pop                                --
//----------------------------------------------------------------------

void ganlib::Lifo::pop() throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(0)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(1)");
  }
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(int_32 &myInteger) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(1)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(1)");
  }
  if (myNode->type != 11) {
    hsmg << "pop() expecting an integer value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myInteger = myNode->value.ival;
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(float_32 &myFloat) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(2)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(2)");
  }
  if (myNode->type != 12) {
    hsmg << "pop() expecting a real value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myFloat = myNode->value.fval;
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(string &myString) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(3)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(3)");
  }
  if (myNode->type != 13) {
    hsmg << "pop() expecting a string value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myString = string(myNode->value.hval);
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(double &myDouble) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(4)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(4)");
  }
  if (myNode->type != 14) {
    hsmg << "pop() expecting a double precision value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myDouble = myNode->value.dval;
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(bool &myBool) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(5)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(5)");
  }
  if (myNode->type != 15) {
    hsmg << "pop() expecting a boolean value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myBool = (myNode->value.ival == 1);
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(ClcmPtr &myClcm) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(6)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(6)");
  }
  if ((myNode->type != 3) && (myNode->type != 4)) {
    hsmg << "pop() expecting a LCM or XSM object (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  lcm *myLcm = myNode->value.mylcm;
  myClcm = ClcmPtr(new Clcm(myLcm, myNode->type, myNode->access, myNode->OSname));
  clepop(&(this->addr));
  free(myNode);
}

void ganlib::Lifo::pop(string myFile, string stype) throw(LifoException){
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(7)");
  try {
    myNode = clepos(&(this->addr), this->addr->nup-1);
  } catch(...) {
    throw LifoException("Exception catched by clepos(7)");
  }
  if (myNode->type == 5) {
    stype = "BINARY";
  } else if (myNode->type == 6) {
    stype = "ASCII";
  } else if (myNode->type == 7) {
    stype = "DA";
  } else {
    hsmg << "pop(): unknown file type";
    throw LifoException(hsmg.str());
  }
  myFile = string(myNode->value.hval);
  clepop(&(this->addr));
  free(myNode);
}

//----------------------------------------------------------------------
//--                              push                                --
//----------------------------------------------------------------------

void ganlib::Lifo::push(string sname, const int_32 myInteger) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str()); myNode->type = 11; myNode->value.ival = myInteger;
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(1)");
  }
}

void ganlib::Lifo::push(string sname, const float_32 myFloat) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str()); myNode->type = 12; myNode->value.fval = myFloat;
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(2)");
  }
}

void ganlib::Lifo::push(string sname, const string myString) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str()); myNode->type = 13;
  strcpy(myNode->value.hval, (char *)myString.c_str());
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(3)");
  }
}

void ganlib::Lifo::push(string sname, const double myDouble) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str()); myNode->type = 14; myNode->value.dval = myDouble;
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(4)");
  }
}

void ganlib::Lifo::push(string sname, bool myBool) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  int_32 ibool = 0; if (myBool) ibool = 1;
  strcpy(myNode->name, (char *)sname.c_str()); myNode->type = 15; myNode->value.ival = ibool;
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(5)");
  }
}

void ganlib::Lifo::push(string sname, ClcmPtr myClcm) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str());
  myNode->type = myClcm->getType()+2;
  myNode->access = myClcm->getAccess();
  lcm *myLcm = myClcm->extract();
  myNode->value.mylcm = myLcm;
  strcpy(myNode->OSname, (char *)myClcm->getName().c_str());
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(6)");
  }
  (this->global_list)->push_back(myClcm);
}

void ganlib::Lifo::push(string sname, string myFile, string stype) throw(LifoException){
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str());
  if (stype == "BINARY") {
    myNode->type = 5;
  } else if (stype == "ASCII") {
    myNode->type = 6;
  } else if (stype == "DA") {
    myNode->type = 7;
  } else {
    hsmg << "push() unknown file type (" << stype << ")";
    throw LifoException(hsmg.str());
  }
  strcpy(myNode->value.hval, (char *)myFile.c_str());
  strcpy(myNode->OSname, (char *)myFile.c_str());
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(7)");
  }
}

void ganlib::Lifo::push(string sname, string myFile, string stype, string OSname) throw(LifoException){
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node));
  strcpy(myNode->name, (char *)sname.c_str());
  if (stype == "BINARY") {
    myNode->type = 5;
  } else if (stype == "ASCII") {
    myNode->type = 6;
  } else if (stype == "DA") {
    myNode->type = 7;
  } else {
    hsmg << "push(): unknown file type (" << stype << ")";
    throw LifoException(hsmg.str());
  }
  strcpy(myNode->value.hval, (char *)myFile.c_str());
  strcpy(myNode->OSname, (char *)OSname.c_str());
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(8)");
  }
}

//----------------------------------------------------------------------
//--                            pushEmpty                             --
//----------------------------------------------------------------------

void ganlib::Lifo::pushEmpty(string sname, string nodeType) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node)); 
  strcpy(myNode->name, (char *)sname.c_str());
  if (nodeType == "I") {
    myNode->type = -11 ;
  } else if (nodeType == "F") {
    myNode->type = -12 ;
  } else if (nodeType == "D") {
    myNode->type = -13 ;
  } else if (nodeType == "S") {
    myNode->type = -14 ;
  } else if (nodeType == "B") {
    myNode->type = -15 ;
  } else if (nodeType == "LCM") {
    myNode->type = -3 ;
    myNode->access = 0 ;
    strcpy(myNode->OSname, (char *)sname.c_str());
  } else if (nodeType == "XSM") {
    myNode->type = -4 ;
    myNode->access = 0 ;
    strcpy(myNode->OSname, (char *)sname.c_str());
  } else if (nodeType == "BINARY") {
    myNode->type = -5 ;
    myNode->access = 0 ;
  } else if (nodeType == "ASCII") {
    myNode->type = -6 ;
    myNode->access = 0 ;
  } else if (nodeType == "DA") {
    myNode->type = -7 ;
    myNode->access = 0 ;
  } else {
    hsmg << nodeType << "is an invalid type in pushEmpty";
    throw LifoException(hsmg.str());
  }
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(9)");
  }
}

void ganlib::Lifo::pushEmpty(string sname, string nodeType, string OSname) throw(LifoException) {
  lifo_node * myNode = (lifo_node *) malloc(sizeof(lifo_node)); 
  strcpy(myNode->name, (char *)sname.c_str());
  if (nodeType == "XSM") {
    myNode->type = -4 ;
    myNode->access = 0 ;
  } else if (nodeType == "BINARY") {
    myNode->type = -5 ;
    myNode->access = 0 ;
  } else if (nodeType == "ASCII") {
    myNode->type = -6 ;
    myNode->access = 0 ;
  } else if (nodeType == "DA") {
    myNode->type = -7 ;
    myNode->access = 0 ;
  } else {
    hsmg << nodeType << "is an invalid type in pushEmpty";
    throw LifoException(hsmg.str());
  }
  strcpy(myNode->OSname, (char *)OSname.c_str());
  try {
    clepush(&(this->addr), myNode);
  } catch(...) {
    throw LifoException("Exception catched by clepush(9)");
  }
}

//----------------------------------------------------------------------
//--                        node by position                          --
//----------------------------------------------------------------------

void ganlib::Lifo::node(int_32 ipos, int_32 &myInteger) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(11)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(11)");
  }
  if (myNode->type != 11) {
    hsmg << "node() expecting an integer value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myInteger = myNode->value.ival;
}

void ganlib::Lifo::node(int_32 ipos, float_32 &myFloat) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(12)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(12)");
  }
  if (myNode->type != 12) {
    hsmg << "node() expecting a real value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myFloat = myNode->value.fval;
}

void ganlib::Lifo::node(int_32 ipos, string &myString) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(13)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(13)");
  }
  if (myNode->type != 13) {
    hsmg << "node() expecting a string value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myString = string(myNode->value.hval);
}

void ganlib::Lifo::node(int_32 ipos, double &myDouble) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(14)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(14)");
  }
  if (myNode->type != 14) {
    hsmg << "node() expecting a double precision value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myDouble = myNode->value.dval;
}

void ganlib::Lifo::node(int_32 ipos, bool &myBool) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(15)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(15)");
  }
  if (myNode->type != 11) {
    hsmg << "node() expecting a boolean value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myBool = (myNode->value.ival == 1);
}

void ganlib::Lifo::node(int_32 ipos, ClcmPtr &myClcm) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(16)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(16)");
  }
  if ((myNode->type != 3) && (myNode->type != 4)) {
    hsmg << "node() expecting a LCM or XSM object (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  lcm *myLcm = myNode->value.mylcm;
  myClcm = ClcmPtr(new Clcm(myLcm, myNode->type, myNode->access, myNode->OSname));
}

void ganlib::Lifo::node(int_32 ipos, string myFile, string stype) throw(LifoException){
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(17)");
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(17)");
  }
  if (myNode->type == 5) {
    stype = "BINARY";
  } else if (myNode->type == 6) {
    stype = "ASCII";
  } else if (myNode->type == 7) {
    stype = "DA";
  } else {
    hsmg << "node(): unknown file type";
    throw LifoException(hsmg.str());
  }
  myFile = string(myNode->value.hval);
}

//----------------------------------------------------------------------
//--                          node by name                            --
//----------------------------------------------------------------------

void ganlib::Lifo::node(string sname, int_32 &myInteger) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(21)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(21)");
  }
  if (myNode->type != 11) {
    hsmg << "node() expecting an integer value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myInteger = myNode->value.ival;
}

void ganlib::Lifo::node(string sname, float_32 &myFloat) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(22)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(22)");
  }
  if (myNode->type != 12) {
    hsmg << "node() expecting a real value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myFloat = myNode->value.fval;
}

void ganlib::Lifo::node(string sname, string &myString) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(23)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(23)");
  }
  if (myNode->type != 13) {
    hsmg << "node() expecting a string value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myString = string(myNode->value.hval);
}

void ganlib::Lifo::node(string sname, double &myDouble) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(24)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(24)");
  }
  if (myNode->type != 14) {
    hsmg << "node() expecting a double precision value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myDouble = myNode->value.dval;
}

void ganlib::Lifo::node(string sname, bool &myBool) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(25)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(25)");
  }
  if (myNode->type != 11) {
    hsmg << "node() expecting a boolean value (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  myBool = (myNode->value.ival == 1);
}

void ganlib::Lifo::node(string sname, ClcmPtr &myClcm) throw(LifoException) {
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to recover from stack(26)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(26)");
  }
  if ((myNode->type != 3) && (myNode->type != 4)) {
    hsmg << "node() expecting a LCM or XSM object (" << myNode->type << " found)";
    throw LifoException(hsmg.str());
  }
  lcm *myLcm = myNode->value.mylcm;
  myClcm = ClcmPtr(new Clcm(myLcm, myNode->type, myNode->access, myNode->OSname));
}

void ganlib::Lifo::node(string sname, string myFile, string stype) throw(LifoException){
  lifo_node *myNode;
  if (this->addr->nup == 0)  throw LifoException("no nodes to pop from stack(27)");
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(27)");
  }
  if (myNode->type == 5) {
    stype = "BINARY";
  } else if (myNode->type == 6) {
    stype = "ASCII";
  } else if (myNode->type == 7) {
    stype = "DA";
  } else {
    hsmg << "node(): unknown file type";
    throw LifoException(hsmg.str());
  }
  myFile = string(myNode->value.hval);
}

//----------------------------------------------------------------------

int_32 ganlib::Lifo::getMax() {
  return this->addr->nup;
}

int_32 ganlib::Lifo::typeNode(int_32 ipos) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(51)");
  }
  return myNode->type;
}

int_32 ganlib::Lifo::accessNode(int_32 ipos) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(52)");
  }
  return myNode->access;
}

string ganlib::Lifo::OSName(int_32 ipos) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(53)");
  }
  return myNode->OSname;
}

string ganlib::Lifo::Name(int_32 ipos) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clepos(&(this->addr), ipos);
  } catch(...) {
    throw LifoException("Exception catched by clepos(54)");
  }
  return myNode->name;
}

int_32 ganlib::Lifo::typeNode(string sname) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(55)");
  }
  return myNode->type;
}

int_32 ganlib::Lifo::accessNode(string sname) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(56)");
  }
  return myNode->access;
}

string ganlib::Lifo::OSName(string sname) throw(LifoException) {
  lifo_node *myNode;
  try {
    myNode = clenode(&(this->addr), sname.c_str());
  } catch(...) {
    throw LifoException("Exception catched by clenode(57)");
  }
  return myNode->OSname;
}

void ganlib::Lifo::lib() throw(LifoException) {
  try {
    clelib(&(this->addr));
  } catch(...) {
    throw LifoException("Exception catched by clelib");
  }
}

lifo *ganlib::Lifo::extract() {
  return this->addr;
}
