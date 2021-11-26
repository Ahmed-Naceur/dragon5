/**
 * This class is a C++ wrapper for a while loop calling a
 * Cle-2000 procedure named PowComponent.c2m at each iteration.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2013)
 */
#ifndef __POW_HXX__
#define __POW_HXX__

#include "Cle2000.hxx"
#include "Communication.hxx"

class POW {
public:
  /** use this constructor to create a new POW object
   */
  POW();

  /** use this method to assign a Calcium component to the POW object
   * @param power thermal reactor power in MW
   * @param component Calcium component reference
   */
  void initialize(double power, void* component);

  /** use this method to execute the POW object
   */
  void run();

private:
  double power_;
  Communication communicator_;
}; // class POW

#endif
