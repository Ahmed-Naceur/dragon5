/**
 * This class is an implementation of the C++ bindings for Calcium.
 * Calcium capabilities with iteration synchronization are available for a program written in C++
 * by using methods belonging to the Communication class.
 *
 * The Communication class uses predefined declarations for some datatypes:
 * <table border="0">
 * <tr> <td><tt>#define IntPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<int_32></tt></td> </tr>
 * <tr> <td><tt>#define FloatPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<float_32></tt></td> </tr>
 * <tr> <td><tt>#define DoublePtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<double_64></tt></td> </tr>
 * <tr> <td><tt>#define BoolPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<bool></tt></td> </tr>
 * </table>
 * <P>
 * Moreover, send methods are constructing shared arrays that <i>cannot</i> be modified in the calling C++
 * code. To prevent this, they are declared <tt>const</tt> using the following predefined declarations:
 * <table border="0">
 * <tr> <td><tt>#define IntPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const int_32></tt></td> </tr>
 * <tr> <td><tt>#define FloatPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const float_32></tt></td> </tr>
 * <tr> <td><tt>#define DoublePtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const double_64></tt></td> </tr>
 * <tr> <td><tt>#define BoolPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const bool></tt></td> </tr>
 * </table>
 * <P>
 *
 * @author Hadrien Leroyer, EDF and Alain Hebert, Ecole Polytechnique
 * de Montreal (2013)
 */
#if ! defined( __Communication_hxx__ )
#define __Communication_hxx__

#include <iostream>
#include "CalciumException.hxx"
extern "C" {
#include <unistd.h>
#include "calcium.h"
}
#include "Cle2000.hxx"

class Communication {

private:

    void* compo_;
     
public:
  /** use this constructor to create a new empty Communication object
   */
    Communication();

  /** use this constructor to create a new Communication object with an embedded component
   * @param compo embedded component
   */
    Communication(const Communication& compo);

virtual ~Communication();

  /** initialize the Communication object and set the embedded component
   * @param compo embedded Calcium component
   */
    int initialize(void* compo);

  /** close the Communication object
   */
    int terminate();

  /** send a single integer value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val integer value
   */
    int send(const int iteration, const std::string portName, const int& val );

  /** send an integer array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the integer array
   * @param tab integer array
   */
    int send(const int iteration, const std::string portName, const int size, IntPtrConst& tab );

  /** send a single real value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val real value
   */
    int send(const int iteration, const std::string portName, const float& val );

  /** send an real array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the real array
   * @param tab real array
   */
    int send(const int iteration, const std::string portName, const int size, FloatPtrConst& tab );

  /** send a single double precision real value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val double precision real value
   */
    int send(const int iteration, const std::string portName, const double& val );

  /** send an double precision real array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the double precision real array
   * @param tab double precision real array
   */
    int send(const int iteration, const std::string portName, const int size, DoublePtrConst& tab );

  /** send a single boolean value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val boolean value
   */
    int send(const int iteration, const std::string portName, const bool& val );

  /** send an boolean array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the boolean array
   * @param tab boolean array
   */
    int send(const int iteration, const std::string portName, const int size, BoolPtrConst& tab );

  /** receive a single integer value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val integer value
   */
    int recv(int& iteration, const std::string portName, int& val );

  /** receive an integer array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the integer array
   * @param tab integer array
   */
    int recv(int& iteration, const std::string portName, const int size, IntPtr& tab );

  /** receive a single real value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val real value
   */
    int recv(int& iteration, const std::string portName, float& val );

  /** receive an real array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the real array
   * @param tab real array
   */
    int recv(int& iteration, const std::string portName, const int size, FloatPtr& tab );

  /** receive a single double precision real value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val double precision real value
   */
    int recv(int& iteration, const std::string portName, double& val );

  /** receive an double precision real array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the double precision real array
   * @param tab double precision real array
   */
    int recv(int& iteration, const std::string portName, const int size, DoublePtr& tab );

  /** receive a single boolean value
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param val boolean value
   */
    int recv(int& iteration, const std::string portName, bool& val );

  /** receive a boolean array
   * @param iteration iteration index
   * @param portName name of the datastream
   * @param size size of the oolean array
   * @param tab boolean array
   */
    int recv(int& iteration, const std::string portName, const int size, BoolPtr& tab );

};

#endif //#if ! defined( __Communication_hxx__ )
