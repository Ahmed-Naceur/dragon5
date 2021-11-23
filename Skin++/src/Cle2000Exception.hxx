/**
 * Exception class for Cle2000 wrapper uses in C++
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2012)
 */
#ifndef Cle2000Exception_HXX
#define Cle2000Exception_HXX

#include <string>
#include <iostream>
#include <stdexcept>

namespace ganlib {

class Cle2000Exception : public std::runtime_error {
public:
  Cle2000Exception(const std::string& msg = "");
}; // class Cle2000Exception */

} // namespace ganlib
#endif
