
/*****************************************/
/*      EXCEPTION CLASS FOR Lifo         */
/*    AUTHOR: A. HEBERT ; 2012/10/07     */
/*****************************************/

/*
Copyright (C) 2012 Ecole Polytechnique de Montreal

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
*/

#include "LifoException.hxx"
using namespace std;
using namespace ganlib;

ganlib::LifoException::LifoException(const string& msg) : runtime_error(msg) {
  cout << "A LifoException was thrown with message: " << msg << endl;
}
