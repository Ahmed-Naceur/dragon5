/**
 * This class is a C++ wrapper for a while loop calling a
 * Cle-2000 procedure named THMComponent.c2m at each iteration.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2013)
 */
#ifndef __THM_HXX__
#define __THM_HXX__

#include "Cle2000.hxx"
#include "Communication.hxx"

class THM {
public:
  /** use this constructor to create a new THM object
   */
  THM();

  /** use this method to assign a Calcium component to the THM object
   * @param component Calcium component reference
   */
  void initialize(void* component);

  /** use this method to execute the THM object
   */
  void run();

private:
  Communication communicator_;
}; // class THM

#endif
