
/*****************************************/
/*        C++ LCM OBJECT WRAPPER         */
/*    AUTHOR: A. HEBERT ; 2010/12/31     */
/*****************************************/

/*
Copyright (C) 2010 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include "Clcm.hxx"

using namespace std;
using namespace ganlib;
static ostringstream hsmg;

/* abort condition setting */
void xabort_c(char *msg){
  cout << "Clcm xabort called: " << string(msg) << endl;
  throw LCMexception(msg);
}

ganlib::Clcm::Clcm(const string stype, const string myName, const string myPath) throw(LCMexception) {
  cout << "New Clcm object called '" << myName << "'" << endl;
  this->root = this;
  this->global_list = new vector<AnyPtr>;
  this->isroot = true;
  this->version = 4;
  bool l_import = false;
  bool l_exists = false;
  int_32 imode=0;
  if (stype == "LCM") {
    this->objtype = 1;
  } else if (stype == "XSM") {
    DIR *dp;
    struct dirent *dirp;
    this->objtype = 2;
    vector<string> names = vector<string>();
    if((dp  = opendir(myPath.c_str())) == NULL) {
      hsmg << "no XSM file found (missing directory name: '" << myPath << ")'" << endl;
      throw LCMexception(hsmg.str());
    }
    while ((dirp = readdir(dp)) != NULL) {
      names.push_back(string(dirp->d_name));
    }
    closedir(dp);
    for (size_t i=0; i<names.size(); i++) {
      if (names[i] == myName) l_exists = true;
    }
  } else if (stype == "BINARY") {
    this->objtype = 3;
  } else if (stype == "ASCII") {
    this->objtype = 4;
  } else if (stype == "DA") {
    this->objtype = 5;
  } else if (stype == "LCM_IMP_BINARY") {
    imode = 1;
    this->objtype = 1;
    l_import = true;
  } else if (stype == "XSM_IMP_BINARY") {
    imode = 1;
    this->objtype = 2;
    l_import = true;
  } else if (stype == "LCM_IMP_ASCII") {
    imode = 2;
    this->objtype = 1;
    l_import = true;
  } else if (stype == "XSM_IMP_ASCII") {
    imode = 2;
    this->objtype = 2;
    l_import = true;
  } else if (stype == "LCM_IMP_ASCII_V3") {
    imode = 2;
    this->objtype = 1;
    this->version = 3;
    l_import = true;
  } else {
    throw LCMexception("invalid type: '" + stype + "'");
  }
  this->access = 0;
  this->name = myName;
  this->path = myPath;
  this->addr = NULL;
  this->lrda = 128;
  try {
    if (l_exists) {
      this->addr = NULL;
      cout << "XSM file " << path << myName << " already exists\n";
      this->open("READ-ONLY");
      this->access = 2;
    } else {
      this->open("NEW");
      this->access = 1;
    }
  } catch(...) {
    throw LCMexception("Exception catched by new(1)");
  }
  if (this->objtype <= 2) {
    this->nammy = "/";
    this->empty = true;
    this->ilong = -1;
  }
  if (l_import) {
    this->empty = false;
    this->untoutched = false;
    try {
      DIR *dp;
      struct dirent *dirp;
      vector<string> names = vector<string>();
      if((dp  = opendir(myPath.c_str())) == NULL) {
        hsmg << "no ASCII file found (missing directory name: '" << myPath << ")'" << endl;
        throw LCMexception(hsmg.str());
      }
      while ((dirp = readdir(dp)) != NULL) {
        names.push_back(string(dirp->d_name));
      }
      closedir(dp);
      for (size_t i=0; i<names.size(); i++) {
        string::size_type loc = names[i].find("_");
        if( loc != string::npos ) {
          string fileName = names[i].substr(loc, string::npos);
          if (fileName == "_"+this->name) {
            FILE *file;
            file = fopen((myPath+names[i]).c_str(), "r");
            if (this->version == 4) {
              lcmexp_c(&(this->addr), 0, file, imode, 2);
            }
            else {
              lcmexpv3_c(&(this->addr), 0, file, imode, 2);
            }
            fclose(file);
            return;
          }
        }
      }
      // Modified R. Chambon 09/2011
      for (size_t i=0; i<names.size(); i++) {
        if (myName.compare(names[i]) == 0) {
          FILE *file;
          file = fopen((myPath+names[i]).c_str(), "r");
          if (this->version == 4) {
            lcmexp_c(&(this->addr), 0, file, imode, 2);
          }
          else {
            lcmexpv3_c(&(this->addr), 0, file, imode, 2);
          }
          fclose(file);
          return;
        }
      }
      //
      hsmg << "File '" << "_" << this->name << "' not found in '" << myPath <<"'";
      throw LCMexception(hsmg.str());
    } catch(...) {
      throw LCMexception("Exception catched by new(2)");
    }
  }
}

ganlib::Clcm::Clcm(lcm *daughter, Clcm *root, const int_32 type, const string name, const int_32 ilong, const string key) throw(LCMexception) {
  this->global_list = root->global_list;
  this->root = root;
  this->addr = daughter;
  this->access = 1;
  this->name = name;
  this->path = "./";
  this->nammy = key;
  this->empty = true;
  this->ilong = ilong;
  this->objtype = type;
  this->untoutched = false;
  this->isroot = false;
  this->version = root->version;
}

ganlib::Clcm::Clcm(lcm *myLcm, const int_32 type, const int_32 access, const string OSname) {
  this->root = this;
  this->global_list = new vector<AnyPtr>;
  this->isroot = false;
  this->version = 4;
  this->addr = myLcm;
  this->access = access;
  this->ilong = -1;
  this->objtype = type-2;
  this->empty = false;
  this->name = OSname;
  this->nammy = "/";
  this->path = "";
  this->untoutched = true;
}

ganlib::Clcm::Clcm(const string stype, ClcmPtr myClcm) throw(LCMexception) {
  if (myClcm->objtype > 2) throw LCMexception("copy constructor not available for file objects");
  if (stype == "LCM") {
    this->objtype = 1;
  } else if (stype == "XSM") {
    this->objtype = 2;
  }
  bool lopen = true;
  if (myClcm->access == 0) {
    lopen = false;
    myClcm->open("READ-ONLY");
  }
  this->version = root->version;
  this->name = myClcm->name + ".bis";
  cout << "New Clcm object (clone) called '" << this->name << "'" << endl;
  this->open("NEW");
  this->access = 1;
  lcm *myLcm = myClcm->addr;
  lcmequ_c(&myLcm, &(this->addr));
  this->root = this;
  this->global_list = new vector<AnyPtr>;
  this->isroot = true;
  this->path = myClcm->path;
  this->lrda = myClcm->lrda;
  if (!lopen) {
    myClcm->close("KEEP");
    this->close("KEEP");
  }
}

ganlib::Clcm::~Clcm() {
  if (this->isroot) {
    cout << "Clcm root destructor called. LCM object=" << this->name << endl;
    if ((this->addr != NULL) & (this->access == 1)) {
      cout << "Destroy object " << this->name << endl;
      this->close("DESTROY");
    } else if ((this->addr != NULL) & (this->access == 2)) {
      cout << "Keep object " << this->name << endl;
      this->close("KEEP");
    }
    if (this->global_list != NULL) delete this->global_list;
  }
}

ostream & ganlib::operator<<(ostream &s, ClcmPtr myClcm) throw(LCMexception) {
  if (myClcm->objtype > 2) throw LCMexception("operator<< not available for file objects");
  FILE *file = fopen("_dummy.txt", "w");
  if (file == NULL) {
    throw LCMexception("fopen failure in operator<<. file=_dummy.txt");
  }
  lcmexp_c(&(myClcm->addr), 0, file, 2, 1);
  fclose(file);
  ifstream ifs ( "_dummy.txt" , ifstream::in );
  while (ifs.good()) s << (char) ifs.get();
  ifs.close();
  if( remove("_dummy.txt") != 0 ) {
    throw LCMexception("Error deleting file after <<: _dummy.txt");
  }
  return s;
}

void ganlib::Clcm::expor(const string stype) {
  this->expor(stype, "_"+this->name);
}

void ganlib::Clcm::expor(const string stype, const string new_name) {
  if (this->objtype > 2) throw LCMexception("expor not available for file objects");
  int_32 imode;
  if (stype == "BINARY") {
    imode = 1;
  } else if (stype == "ASCII") {
    imode = 2;
  } else {
    throw LCMexception("invalid type: '" + stype + "'");
  }
  size_t found=new_name.find("_");
  if (found == string::npos) {
    throw LCMexception("Export file name ("+new_name+") must contains _");
  }
  FILE *file = fopen(new_name.c_str(), "w");
  if (file == NULL) {
    throw LCMexception("fopen failure in expor(). file="+new_name);
  }
  lcmexp_c(&(this->addr), 0, file, imode, 1);
  fclose(file);
}

void ganlib::Clcm::except() {
  throw LCMexception("Programmer-triggered exception");
}

string ganlib::Clcm::getName() throw() {
  return this->name;
}

string ganlib::Clcm::getPath() throw() {
  return this->path;
}

int_32 ganlib::Clcm::getType() throw() {
  return this->objtype;
}

int_32 ganlib::Clcm::getLength() throw(LCMexception) {
  char namlcm[73], nammy[13];
  int_32 empty, ilong, lcml, access;
  if (this->objtype > 2) throw LCMexception("getLength not available for file objects");
  lcminf_c(&(this->addr), namlcm, nammy, &empty, &ilong, &lcml, &access);
  return ilong;
}

string ganlib::Clcm::getDirectory() throw(LCMexception) {
  char namlcm[73], nammy[13];
  int_32 empty, ilong, lcml, access;
  if (this->objtype > 2) throw LCMexception("getDirectory not available for file objects");
  lcminf_c(&(this->addr), namlcm, nammy, &empty, &ilong, &lcml, &access);
  return nammy;
}

int_32 ganlib::Clcm::getAccess() throw(LCMexception) {
  char namlcm[73], nammy[13];
  int_32 empty, ilong, lcml;
  if (this->objtype <= 2) {
    lcminf_c(&(this->addr), namlcm, nammy, &empty, &ilong, &lcml, &(this->access));
  }
  return this->access;
}

int_32 ganlib::Clcm::getLrda() throw() {
  return this->lrda;
}

int_32 ganlib::Clcm::getVersion() throw() {
  return this->version;
}

bool ganlib::Clcm::isEmpty() throw(LCMexception) {
  char namlcm[73], nammy[13];
  int_32 empty, ilong, lcml, access;
  if (this->objtype > 2) throw LCMexception("isEmpty not available for file objects");
  lcminf_c(&(this->addr), namlcm, nammy, &empty, &ilong, &lcml, &access);
  return empty;
}

bool ganlib::Clcm::isNew() throw() {
  return this->untoutched;
}

void ganlib::Clcm::open(const string saccess) throw(LCMexception) {
  if (this->access != 0) {
    throw LCMexception("the object " + this->name + " is already open");
  }
  if (saccess == "NEW") {
    this->untoutched = true;
  } else if (saccess == "READ/WRITE") {
    if ((this->addr == NULL)&(this->objtype == 1)) {
      throw LCMexception("the object " + this->name + " is destroyed(1)");
    }
    this->access = 1;
    this->untoutched = false;
  } else if (saccess == "READ-ONLY") {
    if ((this->addr == NULL)&(this->objtype == 1)) {
      throw LCMexception("the object " + this->name + " is destroyed(2)");
    }
    this->access = 2;
    this->untoutched = false;
  } else {
    throw LCMexception("invalid access: '" + saccess + "'");
  }
  if (this->objtype == 1) {
    try {
      lcmop_c(&(this->addr), (char *)(this->name).c_str(), this->access, this->objtype, 99);
    } catch(...) {
      throw LCMexception("Exception catched by open (1)");
    }
  } else if (this->objtype == 2) {
    try {
      lcmop_c(&(this->addr), (char *)(this->path+this->name).c_str(), this->access, this->objtype, 0);
    } catch(...) {
      throw LCMexception("Exception catched by open (1)");
    }
  } else if (this->objtype <= 5) {
    this->addr = NULL;
    bool exs = false;
    DIR *dp;
    struct dirent *dirp;
    vector<string> names = vector<string>();
    if((dp  = opendir(this->path.c_str())) == NULL) {
      hsmg << "no file found (missing directory name: '" << this->path << ")'" << endl;
      throw LCMexception(hsmg.str());
    }
    while ((dirp = readdir(dp)) != NULL) {
      names.push_back(string(dirp->d_name));
    }
    closedir(dp);
    for (size_t i=0; i<names.size(); i++) {
      if (names[i] == this->name) exs = true;
    }
    if (this->access == 0) {
      if (exs) {
        hsmg << "New file " << this->path << this->name << " exists";;
        throw LCMexception(hsmg.str());
      }
    } else if (this->access == 1) {
      if (!exs) {
        hsmg << "Old file " << this->path << this->name << " does not exist(1)";;
        throw LCMexception(hsmg.str());
      }
    } else if (this->access == 2) {
      if (!exs) {
        hsmg << "Old file " << this->path << this->name << " does not exist(2)";;
        throw LCMexception(hsmg.str());
      }
    }
  } else {
    hsmg << "invalid type in open (" << this->objtype << ")";
    throw LCMexception(hsmg.str());
  }
}

void ganlib::Clcm::close(const string saccess) throw(LCMexception) {
  int iact;
  if (saccess == "KEEP") {
    iact = 1;
  } else if (saccess == "DESTROY") {
    iact = 2;
  } else if (saccess == "ERASE") {
    iact = 3;
  } else {
    throw LCMexception("invalid condition: '" + saccess + "'");
  }
  if (this->access != 0) {
    if (this->objtype <= 2) {
      try {
        lcmcl_c(&(this->addr), iact);
      } catch(...) {
        throw LCMexception("Exception catched by close(1)");
      }
    } else if (this->objtype <= 5) {
      if (iact >= 2) {
        try {
          this->aFile.close();
          if( remove((this->path+this->name).c_str()) != 0 ) {
            throw LCMexception("Error deleting file: " + this->path + this->name);
          }
          cout << "File deleted:" << this->path  << this->name << endl;
        } catch(...) {
          throw LCMexception("Exception catched by close(2)");
        }
      }
    }
  } else {
    throw LCMexception("the object " + this->name + " is already closed");
  }
  if (iact == 2) {
    this->addr = NULL;
  }
  this->access = 0;
}

int_32 ganlib::Clcm::length(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2) throw LCMexception("length not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (itylcm == 3) mylength=4*mylength;
  return mylength;
}

int_32 ganlib::Clcm::length(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2) throw LCMexception("length not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (itylcm == 3) mylength=4*mylength;
  return mylength;
}

int_32 ganlib::Clcm::type(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2) throw LCMexception("type not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (itylcm == 99) {
    hsmg << "length: no block named '" << key << "'";
    throw LCMexception(hsmg.str());
  }
  return itylcm;
}

int_32 ganlib::Clcm::type(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2) throw LCMexception("type not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (itylcm == 99) {
    hsmg << "length: no block at position " << iset << ".";
    throw LCMexception(hsmg.str());
  }
  return itylcm;
}

void ganlib::Clcm::lib() throw(LCMexception) {
  if (this->objtype > 2) throw LCMexception("lib not available for file objects");
  lcmlib_c(&(this->addr));
}

void ganlib::Clcm::val() throw(LCMexception) {
  if (this->objtype > 2) throw LCMexception("val not available for file objects");
  lcmval_c(&(this->addr), " ");
}

//----------------------------------------------------------------------
//--                               set                                --
//----------------------------------------------------------------------

ClcmPtr ganlib::Clcm::setDictionary(const std::string key) throw(LCMexception) {
  if (this->objtype > 2) throw LCMexception("setDictionary not available for file objects");
  lcm *daughter = lcmdid_c(&(this->addr), (char *)key.c_str());
  return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, -1, key));
}

ClcmPtr ganlib::Clcm::setDictionary(const int_32 iset) throw(LCMexception) {
  if (this->objtype > 2) throw LCMexception("setDictionary not available for file objects");
  lcm *daughter = lcmdil_c(&(this->addr), iset);
  return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, -1, " "));
}

ClcmPtr ganlib::Clcm::setList(const std::string key, const int_32 ilong) throw(LCMexception) {
  if (this->objtype > 2) throw LCMexception("setList not available for file objects");
  lcm *daughter = lcmlid_c(&(this->addr), (char *)key.c_str(), ilong);
  return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, 0, key));
}

ClcmPtr ganlib::Clcm::setList(const int_32 iset, const int_32 ilong) throw(LCMexception) {
  if (this->objtype > 2) throw LCMexception("setList not available for file objects");
  lcm *daughter = lcmlil_c(&(this->addr), iset, ilong);
  return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, 0, " "));
}

//----------------------------------------------------------------------
//--                               get                                --
//----------------------------------------------------------------------

ClcmPtr ganlib::Clcm::getClcm(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2) throw LCMexception("getClcm not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getClcm(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  }
  if (itylcm == 0) {
    lcm *daughter = lcmgid_c(&(this->addr), (char *)key.c_str());
    return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, -1, key));
  } else if (itylcm == 10) {
    lcm *daughter = lcmgid_c(&(this->addr), (char *)key.c_str());
    return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, 0, key));
  } else {
    hsmg << "getClcm(" << key << ") expecting a Clcm object(1)";
    throw LCMexception(hsmg.str());
  }
}

ClcmPtr ganlib::Clcm::getClcm(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2) throw LCMexception("getClcm not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getClcm(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  }
  if (itylcm == 0) {
    lcm *daughter = lcmgil_c(&(this->addr), iset);
    return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, -1, " "));
  } else if (itylcm == 10) {
    lcm *daughter = lcmgil_c(&(this->addr), iset);
    return ClcmPtr(new Clcm(daughter, this->root, this->objtype, this->name, 0, " "));
  } else {
    hsmg << "getClcm(" << iset << ") expecting a Clcm object(2)";
    throw LCMexception(hsmg.str());
  }
}

/// @cond DEV
// Use this destructor to avoid deleting LCM array pointed by the shared_array object
template<typename T>
struct my_deleter
{
   void operator()(T* p) { }
};
/// @endcond

IntPtrConst ganlib::Clcm::getInt(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getInt(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 1) {
    hsmg << "getInt(" << key << ") expecting an integer array(1)";
    throw LCMexception(hsmg.str());
  }
  if (this->objtype == 1) {
    // Use a reference without copy (pinning).
    int_32 *array;
    lcmgpd_c(&(this->addr), (char *)key.c_str(), &array);
    IntPtrConst iarray(array, my_deleter<const int_32>());
    return iarray;
  } else if (this->objtype == 2) {
    // Use a reference with copy.
    IntPtrConst iarray(new int_32[mylength]);
    lcmget_c(&(this->addr), (char *)key.c_str(), (int_32 *)iarray.get());
    return iarray;
  } else {
    throw LCMexception("getInt not available for file objects");
  }
}

IntPtrConst ganlib::Clcm::getInt(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getInt(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 1) {
    hsmg << "getInt(" << iset << ") expecting an integer array(2)";
    throw LCMexception(hsmg.str());
  }
  if (this->objtype == 1) {
    // Use a reference without copy (pinning).
    int_32 *array;
    lcmgpl_c(&(this->addr), iset, &array);
    IntPtrConst iarray(array, my_deleter<const int_32>());
    return iarray;
  } else if (this->objtype == 2) {
    // Use a reference with copy.
    IntPtrConst iarray(new int_32[mylength]);
    lcmgdl_c(&(this->addr), iset, (int_32 *)iarray.get());
    return iarray;
  } else {
    throw LCMexception("getInt not available for file objects");
  }
}

FloatPtrConst ganlib::Clcm::getFloat(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getFloat(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 2) {
    hsmg << "getFloat(" << key << ") expecting a single precision real array(1)";
    throw LCMexception(hsmg.str());
  }
  if (this->objtype == 1) {
    // Use a reference without copy (pinning).
    float_32 *array;
    lcmgpd_c(&(this->addr), (char *)key.c_str(), (int_32 **)&array);
    FloatPtrConst iarray(array, my_deleter<const float_32>());
    return iarray;
  } else if (this->objtype == 2) {
    // Use a reference with copy.
    FloatPtrConst iarray(new float_32[mylength]);
    lcmget_c(&(this->addr), (char *)key.c_str(), (int_32 *)iarray.get());
    return iarray;
  } else {
    throw LCMexception("getFloat not available for file objects");
  }
}

FloatPtrConst ganlib::Clcm::getFloat(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getFloat(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 2) {
    hsmg << "getFloat(" << iset << ") expecting a single precision real array(2)";
    throw LCMexception(hsmg.str());
  }
  if (this->objtype == 1) {
    // Use a reference without copy (pinning).
    float_32 *array;
    lcmgpl_c(&(this->addr), iset, (int_32 **)&array);
    FloatPtrConst iarray(array, my_deleter<const float_32>());
    return iarray;
  } else if (this->objtype == 2) {
    // Use a reference with copy.
    FloatPtrConst iarray(new float_32[mylength]);
    lcmgdl_c(&(this->addr), iset, (int_32 *)iarray.get());
    return iarray;
  } else {
    throw LCMexception("getFloat not available for file objects");
  }
}

StringPtrConst ganlib::Clcm::getString(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getString not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getString(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 3) {
    hsmg << "getString(" << key << ") expecting a character array(1)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpd_c(&(this->addr), (char *)key.c_str(), &iarray);
  return StringPtrConst(new string((char*)iarray, 4*mylength));
}

StringPtrConst ganlib::Clcm::getString(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getString not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getString(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 3) {
    hsmg << "getString(" << iset << ") expecting a character array(2)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpl_c(&(this->addr), iset, &iarray);
  return StringPtrConst(new string((char*)iarray, 4*mylength));
}

StringVecPtr ganlib::Clcm::getVecString(const std::string key, const int_32 size) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getVecString not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getVecString(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 3) {
    hsmg << "getVecString(" << key << ") expecting a character array(1)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpd_c(&(this->addr), (char *)key.c_str(), &iarray);
  StringVecPtr ivecPtr = StringVecPtr(new vector<string>);
  for (int i=0; i<size; i++) {
    int_32 nbchar = ( mylength*4 + size - 1) / size;
    int_32 istart = nbchar*i;
    if (i == size-1) nbchar = mylength*4 - nbchar*i;
    if (nbchar != 0) {
      string myString = string((char *)iarray, istart, nbchar);
      ivecPtr->push_back(myString.substr(0, myString.find_last_not_of(" ")+1));
    }
  }
  return ivecPtr;
}

StringVecPtr ganlib::Clcm::getVecString(const int_32 iset, const int_32 size) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getVecString not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getVecString(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 3) {
    hsmg << "getVecString(" << iset << ") expecting a character array(2)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpl_c(&(this->addr), iset, &iarray);
  StringVecPtr ivecPtr = StringVecPtr(new vector<string>);
  for (int i=0; i<size; i++) {
    int_32 nbchar = ( mylength*4 + size - 1) / size;
    int_32 istart = nbchar*i;
    if (i == size-1) nbchar = mylength*4 - nbchar*i;
    if (nbchar != 0) {
      string myString = string((char *)iarray, istart, nbchar);
      ivecPtr->push_back(myString.substr(0, myString.find_last_not_of(" ") + 1));
    }
  }
  return ivecPtr;
}

DoublePtrConst ganlib::Clcm::getDouble(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getDouble(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 4) {
    hsmg << "getDouble(" << key << ") expecting a double precision real array(1)";
    throw LCMexception(hsmg.str());
  }
  if (this->objtype == 1) {
    // Use a reference without copy (pinning).
    double_64 *array;
    lcmgpd_c(&(this->addr), (char *)key.c_str(), (int_32 **)&array);
    DoublePtrConst iarray(array, my_deleter<const double_64>());
    return iarray;
  } else if (this->objtype == 2) {
    // Use a reference with copy.
    DoublePtrConst iarray(new double_64[mylength]);
    lcmget_c(&(this->addr), (char *)key.c_str(), (int_32 *)iarray.get());
    return iarray;
  } else {
    throw LCMexception("getDouble not available for file objects");
  }
}

DoublePtrConst ganlib::Clcm::getDouble(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getDouble(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 4) {
    hsmg << "getDouble(" << iset << ") expecting a double precision real array(2)";
    throw LCMexception(hsmg.str());
  }
  if (this->objtype == 1) {
    // Use a reference without copy (pinning).
    double_64 *array;
    lcmgpl_c(&(this->addr), iset, (int_32 **)&array);
    DoublePtrConst iarray(array, my_deleter<const double_64>());
    return iarray;
  } else if (this->objtype == 2) {
    // Use a reference with copy.
    DoublePtrConst iarray(new double_64[mylength]);
    lcmgdl_c(&(this->addr), iset, (int_32 *)iarray.get());
    return iarray;
  } else {
    throw LCMexception("getDouble not available for file objects");
  }
}

BoolPtrConst ganlib::Clcm::getBool(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getBool not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getBool(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 5) {
    hsmg << "getBool(" << key << ") expecting a boolean array(1)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpd_c(&(this->addr), (char *)key.c_str(), &iarray);
  BoolPtrConst ibool(new bool[mylength]);
  bool *myBool = (bool *)ibool.get();
  for (int i=0; i<mylength; i++) myBool[i] = (iarray[i] == 1);
  return ibool;
}

BoolPtrConst ganlib::Clcm::getBool(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getBool not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getBool(" << iset << ") entry is missing(2)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 5) {
    hsmg << "getBool(" << iset << ") expecting a boolean array(2)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpl_c(&(this->addr), iset, &iarray);
  BoolPtrConst ibool(new bool[mylength]);
  bool *myBool = (bool *)ibool.get();
  for (int i=0; i<mylength; i++) myBool[i] = (iarray[i] == 1);
  return ibool;
}

ComplexPtrConst ganlib::Clcm::getComplex(const std::string key) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getComplex not available for file objects");
  lcmlen_c(&(this->addr), (char *)key.c_str() ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getComplex(" << key << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 6) {
    hsmg << "getComplex(" << key << ") expecting a Complex array(1)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpd_c(&(this->addr), (char *)key.c_str(), &iarray);
  complex<float_32> *icomplex = new complex<float_32>[mylength];
  for (int i=0; i<mylength; i++) {
    icomplex[i] = complex<float_32>((float_32)iarray[2*i],(float_32)iarray[1+2*i]);
  }
  return ComplexPtrConst(icomplex);
}

ComplexPtrConst ganlib::Clcm::getComplex(const int_32 iset) throw(LCMexception) {
  int_32 mylength, itylcm;
  if (this->objtype > 2)  throw LCMexception("getComplex not available for file objects");
  lcmlel_c(&(this->addr), iset ,&mylength, &itylcm);
  if (mylength == 0) {
    hsmg << "getComplex(" << iset << ") entry is missing(1)";
    throw LCMexception(hsmg.str());
  } else if (itylcm != 6) {
    hsmg << "getComplex(" << iset << ") expecting a Complex array(2)";
    throw LCMexception(hsmg.str());
  }
  // Use a reference with copy.
  int_32 *iarray;
  lcmgpl_c(&(this->addr), iset, &iarray);
  complex<float_32> *icomplex = new complex<float_32>[mylength];
  for (int i=0; i<mylength; i++) {
    icomplex[i] = complex<float_32>((float_32)iarray[2*i],(float_32)iarray[1+2*i]);
  }
  return ComplexPtrConst(icomplex);
}

//----------------------------------------------------------------------
//--                               keys                               --
//----------------------------------------------------------------------

StringVecPtr ganlib::Clcm::keys() throw(LCMexception) {
  if (this->objtype > 2)  throw LCMexception("keys not available for file objects");
  if (this->getLength() != -1) throw LCMexception("keys only available for dictionaries");
  StringVecPtr ivecPtr = StringVecPtr(new vector<string>);
  char first[13], key[13];
  strcpy(first, " ");
  lcmnxt_c(&(this->addr), first);
  strcpy(key, first);
  while (true) {
    lcmnxt_c(&(this->addr), key);
    ivecPtr->push_back(string(key));
    if (strcmp(key, first) == 0) break;
  }
  return ivecPtr;
}

//----------------------------------------------------------------------
//--                               put                                --
//----------------------------------------------------------------------

void ganlib::Clcm::put(const std::string key, IntPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  int_32 *myArrayPointer = (int_32 *)myArray.get();
  if ((this->objtype == 1) && (this->global_list != NULL)) {
    // Use a reference without copy (pinning) for memory-resident lcm objects.
    lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 1, myArrayPointer);
    refpush(&(this->addr), myArrayPointer);
    (this->global_list)->push_back(myArray);
  } else if (this->objtype == 1) {
    // Use a reference with copy for memory-resident lcm objects (the lcm object is
    // coming from outside C++).
    int_32 *iarray = setara_c(myLength);
    for (int i=0; i<myLength; i++) iarray[i] = myArrayPointer[i];
    lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 1, iarray);
  } else if (this->objtype == 2) {
    // Use a reference with copy for xsm file objects.
    xsm *ipxsm = (xsm *)(this->addr);
    xsmput_c(&ipxsm, (char *)key.c_str(), myLength, 1, myArrayPointer);
  }
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, IntPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  int_32 *myArrayPointer = (int_32 *)myArray.get();
  if ((this->objtype == 1) && (this->global_list != NULL)) {
    // Use a reference without copy (pinning) for memory-resident lcm objects.
    lcmppl_c(&(this->addr), iset, myLength, 1, myArrayPointer);
    refpush(&(this->addr), myArrayPointer);
    (this->global_list)->push_back(myArray);
  } else if (this->objtype == 1) {
    // Use a reference with copy for memory-resident lcm objects (the lcm object is
    // coming from outside C++).
    int_32 *iarray = setara_c(myLength);
    for (int i=0; i<myLength; i++) iarray[i] = myArrayPointer[i];
    lcmppl_c(&(this->addr), iset, myLength, 1, iarray);
  } else if (this->objtype == 2) {
    // Use a reference with copy for xsm file objects.
    xsm *ipxsm = (xsm *)(this->addr) + iset;
    xsmput_c(&ipxsm, " ", myLength, 1, myArrayPointer);
  }
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, FloatPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  float_32 *myArrayPointer = (float_32 *)myArray.get();
  if ((this->objtype == 1) && (this->global_list != NULL)) {
    // Use a reference without copy (pinning) for memory-resident lcm objects.
    lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 2, (int_32 *)myArrayPointer);
    refpush(&(this->addr), (int_32 *)myArrayPointer);
    (this->global_list)->push_back(myArray);
  } else if (this->objtype == 1) {
    // Use a reference with copy for memory-resident lcm objects (the lcm object is
    // coming from outside C++).
    int_32 *iarray = setara_c(myLength);
    for (int i=0; i<myLength; i++) iarray[i] = (int_32)myArrayPointer[i];
    lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 2, iarray);
  } else if (this->objtype == 2) {
    // Use a reference with copy for xsm file objects.
    xsm *ipxsm = (xsm *)(this->addr);
    xsmput_c(&ipxsm, (char *)key.c_str(), myLength, 2, (int_32 *)myArrayPointer);
  }
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, FloatPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  float_32 *myArrayPointer = (float_32 *)myArray.get();
  if ((this->objtype == 1) && (this->global_list != NULL)) {
    // Use a reference without copy (pinning) for memory-resident lcm objects.
    lcmppl_c(&(this->addr), iset, myLength, 2, (int_32 *)myArrayPointer);
    refpush(&(this->addr), (int_32 *)myArrayPointer);
    (this->global_list)->push_back(myArray);
  } else if (this->objtype == 1) {
    // Use a reference with copy for memory-resident lcm objects (the lcm object is
    // coming from outside C++).
    int_32 *iarray = setara_c(myLength);
    for (int i=0; i<myLength; i++) iarray[i] = (int_32)myArrayPointer[i];
    lcmppl_c(&(this->addr), iset, myLength, 2, iarray);
  } else if (this->objtype == 2) {
    // Use a reference with copy for xsm file objects.
    xsm *ipxsm = (xsm *)(this->addr) + iset;
    xsmput_c(&ipxsm, " ", myLength, 2, (int_32 *)myArrayPointer);
  }
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, StringPtr myArray) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if ((myArray.get())->length() <= 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  int_32 ilen = ((myArray.get())->length()+3)/4;
  int_32 *iarray = (int_32 *)setara_c(ilen);
  for (int i=0; i<ilen; i++) strncpy((char *)&iarray[i], "    ", 4);
  string myString = *myArray.get();
  strncpy((char *)iarray, myString.c_str(), myString.length());
  lcmppd_c(&(this->addr), (char *)key.c_str(), ilen, 3, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, StringPtr myArray) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if ((myArray.get())->length() <= 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  int_32 ilen = ((myArray.get())->length()+3)/4;
  int_32 *iarray = (int_32 *)setara_c(ilen);
  for (int i=0; i<ilen; i++) strncpy((char *)&iarray[i], "    ", 4);
  string myString = *myArray.get();
  strncpy((char *)iarray, myString.c_str(), myString.length());
  lcmppl_c(&(this->addr), iset, ilen, 3, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, StringVecPtr myArray) throw(LCMexception) {
  const vector<string> *ivec = (vector<string>*)myArray.get();
  const int_32 mySize = ivec->size();
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if (mySize <= 0) {
    hsmg << "put(" << key << ") has zero size(1)";
    throw LCMexception(hsmg.str());
  }
  int_32 maxLength = 0;
  for (int i=0; i<mySize; i++) {
    maxLength = max(maxLength, (int_32)(*ivec)[i].length());
  }
  if (maxLength == 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  maxLength = (maxLength+3)/4;
  int_32 *iarray = (int_32 *)setara_c(maxLength*mySize);
  for (int i=0; i<maxLength*mySize; i++) strncpy((char *)&iarray[i], "    ", 4);
  for (int i=0; i<mySize; i++) {
    string myString = (*ivec)[i];
    strncpy((char *)&iarray[maxLength*i], myString.c_str(), myString.length());
  }
  lcmppd_c(&(this->addr), (char *)key.c_str(), maxLength*mySize, 3, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, StringVecPtr myArray) throw(LCMexception) {
  const vector<string> *ivec = (vector<string>*)myArray.get();
  const int_32 mySize = ivec->size();
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if (mySize <= 0) {
    hsmg << "put(" << iset << ") has zero size(2)";
    throw LCMexception(hsmg.str());
  }
  int_32 maxLength = 0;
  for (int i=0; i<mySize; i++) {
    maxLength = max(maxLength, (int_32)(*ivec)[i].length());
  }
  if (maxLength == 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  maxLength = (maxLength+3)/4;
  int_32 *iarray = (int_32 *)setara_c(maxLength*mySize);
  for (int i=0; i<maxLength*mySize; i++) strncpy((char *)(&iarray[i]), "    ", 4);
  for (int i=0; i<mySize; i++) {
    string myString = (*ivec)[i];
    strncpy((char *)(&iarray[maxLength*i]), myString.c_str(), myString.length());
  }
  lcmppl_c(&(this->addr), iset, maxLength*mySize, 3, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, DoublePtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  double_64 *myArrayPointer = (double_64 *)myArray.get();
  if ((this->objtype == 1) && (this->global_list != NULL)) {
    // Use a reference without copy (pinning) for memory-resident lcm objects.
    lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 4, (int_32 *)myArrayPointer);
    refpush(&(this->addr), (int_32 *)myArrayPointer);
    (this->global_list)->push_back(myArray);
  } else if (this->objtype == 1) {
    // Use a reference with copy for memory-resident lcm objects (the lcm object is
    // coming from outside C++).
    int_32 *iarray = setara_c(2*myLength);
    for (int i=0; i<2*myLength; i++) iarray[i] = ((int_32 *)myArrayPointer)[i];
    lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 4, iarray);
  } else if (this->objtype == 2) {
    // Use a reference with copy for xsm file objects.
    xsm *ipxsm = (xsm *)(this->addr);
    xsmput_c(&ipxsm, (char *)key.c_str(), myLength, 4, (int_32 *)myArrayPointer);
  }
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, DoublePtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  double_64 *myArrayPointer = (double_64 *)myArray.get();
  if ((this->objtype == 1) && (this->global_list != NULL)) {
    // Use a reference without copy (pinning) for memory-resident lcm objects.
    lcmppl_c(&(this->addr), iset, myLength, 4, (int_32 *)myArrayPointer);
    refpush(&(this->addr), (int_32 *)myArrayPointer);
    (this->global_list)->push_back(myArray);
  } else if (this->objtype == 1) {
    // Use a reference with copy for memory-resident lcm objects (the lcm object is
    // coming from outside C++).
    int_32 *iarray = setara_c(2*myLength);
    for (int i=0; i<2*myLength; i++) iarray[i] = ((int_32 *)myArrayPointer)[i];
    lcmppl_c(&(this->addr), iset, myLength, 4, iarray);
  } else if (this->objtype == 2) {
    // Use a reference with copy for xsm file objects.
    xsm *ipxsm = (xsm *)(this->addr) + iset;
    xsmput_c(&ipxsm, " ", myLength, 4, (int_32 *)myArrayPointer);
  }
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, BoolPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  bool *myArrayPointer = (bool *)myArray.get();
  // Use a reference with copy.
  int_32 *iarray = setara_c(myLength);
  for (int i=0; i<myLength; i++) {
    myArrayPointer[i] ? iarray[2*i] = 1 : iarray[2*i] = 0;
  }
  lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 5, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, BoolPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  bool *myArrayPointer = (bool *)myArray.get();
  // Use a reference with copy.
  int_32 *iarray = setara_c(myLength);
  for (int i=0; i<myLength; i++) {
    myArrayPointer[i] ? iarray[2*i] = 1 : iarray[2*i] = 0;
  }
  lcmppl_c(&(this->addr), iset, myLength, 5, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, ComplexPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << key << ") has zero length(1)";
    throw LCMexception(hsmg.str());
  }
  complex<float_32> *myArrayPointer = (complex<float_32> *)myArray.get();
  // Use a reference with copy.
  int_32 *iarray = setara_c(2*myLength);
  for (int i=0; i<myLength; i++) {
    iarray[2*i] = (int_32)myArrayPointer[i].real();
    iarray[1+2*i] = (int_32)myArrayPointer[i].imag();
  }
  lcmppd_c(&(this->addr), (char *)key.c_str(), myLength, 6, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const int_32 iset, ComplexPtr myArray, const int_32 myLength) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  } else if (myLength <= 0) {
    hsmg << "put(" << iset << ") has zero length(2)";
    throw LCMexception(hsmg.str());
  }
  complex<float_32> *myArrayPointer = (complex<float_32> *)myArray.get();
  // Use a reference with copy.
  int_32 *iarray = setara_c(2*myLength);
  for (int i=0; i<myLength; i++) {
    iarray[2*i] = (int_32)myArrayPointer[i].real();
    iarray[1+2*i] = (int_32)myArrayPointer[i].imag();
  }
  lcmppl_c(&(this->addr), iset, myLength, 6, iarray);
  myArray.reset();
  return;
}

void ganlib::Clcm::put(const std::string key, ClcmPtr myClcm) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << key << ") forbidden if the object is not open in modification mode(1)";
    throw LCMexception(hsmg.str());
  }
  lcm *daughter = lcmgid_c(&(this->addr), (char *)key.c_str());
  lcm *myLcm = myClcm->addr;
  lcmequ_c(&myLcm, &daughter);
  return;
}

void ganlib::Clcm::put(const int_32 iset, ClcmPtr myClcm) throw(LCMexception) {
  if (this->getAccess() != 1) {
    hsmg << "put(" << iset << ") forbidden if the object is not open in modification mode(2)";
    throw LCMexception(hsmg.str());
  }
  lcm *daughter = lcmgil_c(&(this->addr), iset);
  lcm *myLcm = myClcm->addr;
  lcmequ_c(&myLcm, &daughter);
  return;
}

lcm *ganlib::Clcm::extract() {
  return this->addr;
}
