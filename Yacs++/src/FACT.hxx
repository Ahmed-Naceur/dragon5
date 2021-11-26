/**
 * This class is a C++ wrapper for calling the FACT.c2m CLE-2000
 procedure from YACS.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2013)
 */
#ifndef __FACT_HXX__
#define __FACT_HXX__

class FACT {
public:
  /** use this constructor to create a new FACT object
   */
  FACT();

  /** use this method to execute the FACT object
   * @param a input integer
   * @param b output integer containing factorial of a
   */
  void run(long a, long& b);
}; // class FACT

#endif
