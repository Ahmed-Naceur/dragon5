/**
 * Exception class for Lifo stack uses in C++
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2012)
 */
#ifndef LifoException_HXX
#define LifoException_HXX

#include <string>
#include <iostream>
#include <stdexcept>

namespace ganlib {

class LifoException : public std::runtime_error {
public:
  LifoException(const std::string& msg = "");
}; // class LifoException */

} // namespace ganlib
#endif
