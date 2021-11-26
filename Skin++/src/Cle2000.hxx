/**
 * Call a parametrized Cle2000 procedure. A Lifo object is used to manage
 * input/output parameters for this CLE-2000 procedure.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2012)
 */
#ifndef Cle2000_HXX
#define Cle2000_HXX

#include <cstddef> // for __GLIBCXX__
 
#ifdef __GLIBCXX__ // for memory management with shared_ptr
#  include <tr1/memory>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <memory>
#endif

#include "Cle2000Exception.hxx"
#include "Lifo.hxx"
extern "C" {
#include <dirent.h>
int_32 donmod(char *cmodul, int_32 nentry, char (*hentry)[13], int_32 *ientry,
       int_32 *jentry, lcm **kentry, char (*hparam)[73]);
}
#define Cle2000Ptr boost::shared_ptr<Cle2000>

namespace ganlib {

/**
 * Call a parametrized CLE-2000 procedure. A CLE-2000 procedure is implemented
 * as a 80-column ascii file. The name of that file must have a <tt>".c2m"</tt>
 * suffix. A Lifo object is used to manage input/output parameters for this
 * CLE-2000 procedure.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2012)
 */
class Cle2000 {
public:
  /** use this constructor to create a new CLE-2000 object
   * @param sname string containing the name of the ascii file containing the CLE-2000
   * procedure (without the <tt>".c2m"</tt> suffix)
   */
  Cle2000(std::string sname);

  /** use this constructor to create a new Cle2000 object
   * @param sname string containing the name of the ascii file containing the CLE-2000
   * procedure (without the <tt>".c2m"</tt> suffix)
   * @param edit user-defined edition index. Increasing value of <tt>edit</tt> will
   * cause increasing amount of listing information. Set <tt>edit</tt> to zero to avoid
   * listing information.
   */
  Cle2000(std::string sname, int_32 edit);

  /** use this constructor to create a new Cle2000 object with an embedded Lifo stack
   * @param sname string containing the name of the ascii file containing the CLE-2000
   * procedure (without the <tt>".c2m"</tt> suffix)
   * @param jstack Lifo stack to include in Cle2000 object
   */
  Cle2000(std::string sname, LifoPtr jstack);

  /** use this constructor to create a new Cle2000 object with an embedded Lifo stack
   * @param sname string containing the name of the ascii file containing the CLE-2000
   * procedure (without the <tt>".c2m"</tt> suffix)
   * @param edit user-defined edition index. Increasing value of <tt>edit</tt> will
   * cause increasing amount of listing information. Set <tt>edit</tt> to zero to avoid
   * listing information.
   * @param jstack Lifo stack to include in Cle2000 object
   */
  Cle2000(std::string sname, int_32 edit, LifoPtr jstack);

  /** Close and destroy a Cle2000 object.
   */
  ~Cle2000() throw(Cle2000Exception);

  /** attach a lifo stack to the Cle2000 object
   * @param myLifo Lifo stack containing input/output parameters
   */
  void setLifo(LifoPtr myLifo);

  /** call the native CLE-2000 procedure
   */
  void exec() throw(Cle2000Exception);

private:
  std::string procName;
  int_32 edit;
  LifoPtr stack;
}; // class Cle2000

} // namespace ganlib
#endif
