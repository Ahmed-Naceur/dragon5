
/*****************************************/
/*       EXCEPTION CLASS FOR LCM         */
/*    AUTHOR: A. HEBERT ; 2010/12/31     */
/*****************************************/

/*
Copyright (C) 2010 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include "LCMexception.hxx"
using namespace std;
using namespace ganlib;

ganlib::LCMexception::LCMexception(const string& msg) : runtime_error(msg) {
  cout << "A LCMexception was thrown with message: " << msg << endl;
}
