/**
 * Exception class for LCM uses in C++
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2010)
 */
#ifndef LCMexception_HXX
#define LCMexception_HXX

#include <string>
#include <iostream>
#include <stdexcept>

namespace ganlib {

class LCMexception : public std::runtime_error {
public:
  LCMexception(const std::string& msg = "");
}; // class LCMexception */

} // namespace ganlib
#endif
