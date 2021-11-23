/**
 * This class is an implementation of the C++ bindings for LCM.
 * LCM capabilities are available for a program written in C++
 * by using methods belonging to the Clcm class.
 * <P> A Clcm object in C++ can encapsulate a native LCM object
 * which can be memory-resident or persistent (using the XSM C API)
 * or a file. A LCM or XSM native object can contains dictionaries
 * (aka associative tables or hash tables) and lists.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2010)
 */
#ifndef Clcm_HXX
#define Clcm_HXX

#include <cstddef> // for __GLIBCXX__
 
#ifdef __GLIBCXX__ // for memory management with shared_ptr
#  include <tr1/memory>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <memory>
#endif

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <complex>
#include<boost/shared_ptr.hpp>
#include<boost/shared_array.hpp>
#include<boost/variant.hpp>
#include "LCMexception.hxx"
extern "C" {
#include <dirent.h>
#include "lcm.h"
}

#define IntPtr boost::shared_array<int_32>
#define FloatPtr boost::shared_array<float_32>
#define StringPtr boost::shared_ptr<std::string>
#define DoublePtr boost::shared_array<double_64>
#define BoolPtr boost::shared_array<bool>
#define ComplexPtr boost::shared_array< std::complex<float_32> >
#define StringVecPtr boost::shared_ptr< std::vector<std::string> >
#define ClcmPtr boost::shared_ptr<Clcm>

// support for pinning. 
#define AnyPtr boost::variant<IntPtr, FloatPtr, StringPtr, DoublePtr, BoolPtr, ComplexPtr>
#define IntPtrConst boost::shared_array<const int_32>
#define FloatPtrConst boost::shared_array<const float_32>
#define StringPtrConst boost::shared_ptr<const std::string>
#define DoublePtrConst boost::shared_array<const double_64>
#define BoolPtrConst boost::shared_array<const bool>
#define ComplexPtrConst boost::shared_array<const std::complex<float_32> >

namespace ganlib {

/**
 * This class is an implementation of the C++/Boost bindings for LCM.
 * LCM capabilities are available for a program written in C++
 * by using methods belonging to the Clcm class. These methods
 * encapsulate the LCM C API calls used as "extern"C" functions.
 * <P> A Clcm object in C++ can encapsulate a native LCM object
 * which can be memory-resident or persistent (using the XSM C API)
 * or a file. A LCM or XSM native object can contains dictionaries
 * (aka associative tables or hash tables) and lists.
 * <P>
 * The Clcm class uses predefined declarations for some datatypes:
 * <table border="0">
 * <tr> <td><tt>#define IntPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<int_32></tt></td> </tr>
 * <tr> <td><tt>#define FloatPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<float_32></tt></td> </tr>
 * <tr> <td><tt>#define StringPtr</tt>:</td> &nbsp; <td><tt>boost::shared_ptr<std::string></tt></td> </tr>
 * <tr> <td><tt>#define DoublePtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<double_64></tt></td> </tr>
 * <tr> <td><tt>#define BoolPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array<bool></tt></td> </tr>
 * <tr> <td><tt>#define ComplexPtr</tt>:</td> &nbsp; <td><tt>boost::shared_array< complex<float_32> ></tt></td> </tr>
 * <tr> <td><tt>#define StringVecPtr</tt>:</td> &nbsp; <td><tt>boost::shared_ptr< std::vector<std::string> ></tt></td> </tr>
 * <tr> <td><tt>#define ClcmPtr</tt>:</td> &nbsp; <td><tt>boost::shared_ptr<Clcm></tt></td> </tr>
 * </table>
 * <P>
 * Moreover, get methods are constructing shared arrays that <i>cannot</i> be modified in the calling C++
 * code. To prevent this, they are declared <tt>const</tt> using the following predefined declarations:
 * <table border="0">
 * <tr> <td><tt>#define IntPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const int_32></tt></td> </tr>
 * <tr> <td><tt>#define FloatPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const float_32></tt></td> </tr>
 * <tr> <td><tt>#define StringPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_ptr<const std::string></tt></td> </tr>
 * <tr> <td><tt>#define DoublePtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const double_64></tt></td> </tr>
 * <tr> <td><tt>#define BoolPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const bool></tt></td> </tr>
 * <tr> <td><tt>#define ComplexPtrConst</tt>:</td> &nbsp; <td><tt>boost::shared_array<const complex<float_32> ></tt></td> </tr>
 * </table>
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2010)
 */

class Clcm {
public:
  /**
   * Use this constructor to create a Clcm object of a specified type.
   * The new LCM object is open in modification mode and is in a state
   * such that the <tt>isNew()</tt> method will return <tt>true</tt>.
   * If a XSM file <tt>name</tt> exists, the Clcm object is pointing to
   * this file open in <tt>"READ-ONLY"</tt> mode and is in a state
   * such that the <tt>isNew()</tt> method will return <tt>false</tt>.
   * If the type is <tt>"LCM_IMP"</tt> or <tt>"XSM_IMP"</tt>, the new
   * LCM object is open in modification mode and the <tt>isNew()</tt>
   * method will return <tt>false</tt>.
   * <p>
   * A Clcm object is created using a construct similar to the following one:
   * <pre>
   * ClcmPtr myNewClcm = ClcmPtr(new Clcm("XSM_IMP_ASCII", "Multicompo", "./"));
   * </pre>
   * where a XSM file is imported from ASCII file <tt>./Multicompo</tt>.
   *
   * @param stype type of the new Clcm object. This variable is chosen among
   * the following values:
   * <ul><dl>
   * <dt> <tt>"LCM"</tt> <dd> memory-resident LCM object
   * <dt> <tt>"XSM"</tt> <dd> XSM object (or persistent LCM object)
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * <dt> <tt>"DA"</tt> <dd> binary direct-access file (with 128-word records)
   * <dt> <tt>"LCM_IMP_BINARY"</tt> <dd> memory-resident LCM object constructed by
   * importing from a binary sequential file named <tt>"_"</tt>+name
   * <dt> <tt>"XSM_IMP_BINARY"</tt> <dd> XSM object constructed by importing
   * from a binary sequential file named <tt>"_"</tt>+name
   * <dt> <tt>"LCM_IMP_ASCII"</tt> <dd> memory-resident LCM object constructed by
   * importing from an ASCII sequential file named <tt>"_"</tt>+name
   * <dt> <tt>"XSM_IMP_ASCII"</tt> <dd> XSM object constructed by importing
   * from an ASCII sequential file named <tt>"_"</tt>+name
   * </dl></ul>
   * @param myName creation name of the new Clcm object
   * @param myPath path to access the associated file. This can be a sequential import
   *        file (with <tt>"_"</tt> prefix), an <tt>"XSM"</tt>, <tt>"BINARY"</tt>,
   *        <tt>"ASCII"</tt>  or <tt>"DA"</tt> file. Set to <tt>"./"</tt> to work
   *        in the current directory.
   */
  Clcm(const std::string stype, const std::string myName, const std::string myPath) throw(LCMexception);

  /// @cond DEV
  Clcm(lcm *, Clcm *, const int_32, const std::string, const int_32, const std::string) throw(LCMexception);
  /// @endcond

  /**
   * Use this constructor to create a clone of an existing Clcm object. This is a deep-copy
   * operation. A clone is created using the following construct:
   * <pre>
   * ClcmPtr myClcmClone = ClcmPtr(new Clcm("LCM", myClcm));
   * </pre>
   * @param stype type of the new Clcm object. This variable is chosen among
   * the following values:
   * <ul><dl>
   * <dt> <tt>"LCM"</tt> <dd> memory-resident LCM object
   * <dt> <tt>"XSM"</tt> <dd> XSM object (or persistent LCM object)
   * </dl></ul>
   * @param myClcm existing ClcmPtr object accessed in read-only mode. This object must me of
   * LCM or XSM type.
   */
  Clcm(const std::string stype, ClcmPtr myClcm) throw(LCMexception);

  /**
   * Use this constructor to encapsulate an open LCM or XSM object into a Clcm object. The
   * existing LCM or XSM object is not garbage collected (it may belong to another Clcm object).
   * @param mylcm existing LCM or XSM object.
   * @param type type of object (=1: LCM; =2:XSM).
   * @param access access mode of object (=1: ; =2:).
   * @param OSname operating system name of object.
   */
  Clcm(lcm *mylcm, const int_32 type, const int_32 access, const std::string OSname);

  /** Close and destroy a Clcm object if in modification mode; close and
   * keep a Clcm object if in read-only mode.
   */
  ~Clcm();

  /** Export a Clcm object into an ostream (ascii stream). This method is not available for
   * file-type Clcm objects. A Clcm object can be dumped to the standard output using:
   * <pre>
   * cout << myClcm;
   * </pre>
    * @param s initial std::ostream objet.
    * @param myClcm ClcmPtr object to export.
    * @return final std::ostream objet including the ClcmPtr object.
   */
  friend std::ostream & operator<<(std::ostream &s, ClcmPtr myClcm) throw(LCMexception);

  /** Serialize and save the object content on a sequential file. This method is not
   * available for file-type Clcm objects. The name of the sequential file is
   * the catenation of <tt>"_"</tt> with the name of the Clcm object.
   * @param stype type of the export file. This variable is chosen among
   * the following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * </dl></ul>
   */
  void expor(const std::string stype);

 /** Serialize and save the object content on a sequential file. This method is not
   * available for file-type Clcm objects.
   * @param stype type of the export file. This variable is chosen among
   * the following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * </dl></ul>
   * @param new_name name of the sequential file. This name must begin by
   * character <tt>"_"</tt>.
   */
  void expor(const std::string stype, const std::string new_name);

  /** Cause an exception (used to debug XABORT in C++).
   */
  void except();

  /** Return the name of the Clcm object.
   */
  std::string getName() throw();

  /** Return the path to access the associated file.
   */
  std::string getPath() throw();

  /** Return the type of the Clcm object.
   * @return =1: memory-resident LCM object; =2: persistent LCM object;
   * =3: binary sequential file; =4: ASCII sequential file; =5 binary
   * direct-access file.
   */
  int_32 getType() throw();

  /** Return the length of a list-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * @return =-1: dictionary; >=1: length of the list-type Clcm object.
   */
  int_32 getLength() throw(LCMexception);

  /** Return the name of the accessible directory of a dictionary-type
   * Clcm object. This method is not available for file-type Clcm objects.
   * @return =<tt>"/"</tt> for a dictionary on root or name of the accessible
   * directory.
   */
  std::string getDirectory() throw(LCMexception);

  /** Return the access type of a Clcm object.
   * @return =0: the object is closed; =1: the object is open in modification
   * mode; =2: the object is open in read-only mode.
   */
  int_32 getAccess() throw(LCMexception);

  /** Return the number of words in a record of a direct access-type Clcm object.
   */
  int_32 getLrda() throw();
	
  /** Return the original version of the imported object.
   */
  int_32 getVersion() throw();
	
  /** Return true if the a dictionary-type Clcm object is empty. This method
   * is not available for file-type Clcm objects.
   */
  bool isEmpty() throw(LCMexception);

  /** Return true if the dictionary-type Clcm object is new. This method
   * is not available for file-type Clcm objects.
   */
  bool isNew() throw();

  /** Open a new Clcm object or reopen an existing Clcm object already closed by
   * the close(<tt>"KEEP"</tt>) method. In this case, the Clcm object is reopen
   * in a state such as the <tt>isNew()</tt> method will return false.
   * @param saccess type of open. This variable is chosen among the following
   * values:
   * <ul><dl>
   * <dt> <tt>"NEW"</tt> <dd> open a new Clcm object.
   * <dt> <tt>"READ/WRITE"</tt> <dd> open in modification mode.
   * <dt> <tt>"READ-ONLY"</tt> <dd> open in read-only mode.
   * </dl></ul>
   */
  void open(const std::string saccess) throw(LCMexception);

  /** Close a Clcm object.
   * @param saccess =<tt>"KEEP"</tt>: close without destruction of the object
   * content; =<tt>"DESTROY"</tt>: close with destruction of the object content.
   */
  void close(const std::string saccess) throw(LCMexception);

  /**
   * Return the length of a block of information in a dictionary-type
   * Clcm object. This method is not available for file-type Clcm objects.
   * @param key key identification of the block in the dictionary
   * @return the number of components for an elementary block; -1 for
   * a daughter dictionary; 0: the block does't exist; length of the list
   * for a daughter list; length in characters for a string array.
   */
  int_32 length(const std::string key) throw(LCMexception);

  /**
   * Return the length of a block of information in a list-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param iset index of the block in the list. The first list element is
   * stored at index 0.
   * @return the number of components for an elementary block; -1 for
   * a daughter dictionary; 0: the block does't exist;  length of the list
   * for a daughter list; length in characters for a string array.
   */
  int_32 length(const int_32 iset) throw(LCMexception);

  /**
   * Return the type of a block of information in a dictionary-type
   * Clcm object. This method is not available for file-type Clcm objects.
   * @param key key identification of the block in the dictionary
   * @return =0: dictionary; =1: integer (int_32); =2: real number (float);
   * =3: string; =4: real number (double); =5: boolean; =6: Complex object;
   * =10: list.
   */
  int_32 type(const std::string key) throw(LCMexception);

  /**
   * Return the type of a block of information in a list-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param iset index of the block in the list. The first list element is
   * stored at index 0.
   * @return =0: dictionary; =1: integer (int_32); =2: real number (float);
   * =3: string; =4: real number (double); =5: boolean; =6: Complex object;
   * =10: list.
   */
  int_32 type(const int_32 iset) throw(LCMexception);

  /** Print the table of contents of a dictionary- or list-type Clcm object.
   * This method is not available for file-type Clcm objects.
   */
  void lib() throw(LCMexception);

  /** Validate a dictionary- or list-type Clcm object. Detect possible memory
   * corruption.
   * This method is not available for file-type Clcm objects.
   */
  void val() throw(LCMexception);

  /** Set a daughter dictionary-type Clcm object from a dictionary-type Clcm
   * object. This method is not available for file-type Clcm objects.
   * @param key key identification of the daughter Clcm object in the dictionary
   * @return daughter ClcmPtr object of dictionary or list type.
   */
  ClcmPtr setDictionary(const std::string key) throw(LCMexception);

  /** Set a daughter dictionary-type Clcm object from a list-type Clcm
   * object. This method is not available for file-type Clcm objects.
   * @param iset index of the daughter Clcm object in the list. The first list
   * element is stored at index 0.
   * @return daughter ClcmPtr object of dictionary or list type.
   */
  ClcmPtr setDictionary(const int_32 iset) throw(LCMexception);

  /** Set a daughter list-type Clcm object from a dictionary-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param key key identification of the daughter Clcm object in the dictionary
   * @param ilong initial length of the heterogeneous list
   * @return daughter ClcmPtr object of dictionary or list type.
   */
  ClcmPtr setList(const std::string key, const int_32 ilong) throw(LCMexception);

  /** Set a daughter list-type Clcm object from a list-type Clcm object
   * This method is not available for file-type Clcm objects.
   * @param iset index of the daughter Clcm object in the list. The first list
   * element is stored at index 0.
   * @param ilong initial length of the heterogeneous list
   * @return daughter ClcmPtr object of dictionary or list type.
   */
  ClcmPtr setList(const int_32 iset, const int_32 ilong) throw(LCMexception);

  /** Recover a daughter Clcm object from an existing dictionary-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param key key identification of the daughter Clcm object in the dictionary
   * @return daughter ClcmPtr object of dictionary or list type.
   */
  ClcmPtr getClcm(const std::string key) throw(LCMexception);

  /** Recover a daughter Clcm object from an existing list-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param iset index of the daughter Clcm object in the list. The first list
   * element is stored at index 0.
   * @return daughter ClcmPtr object of dictionary or list type.
   */
  ClcmPtr getClcm(const int_32 iset) throw(LCMexception);

  /** Recover an integer array from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getInt</tt>. This method is not available for file-type Clcm objects.
   * <p>
   * <b>Example</b>: An integer array named <tt>"myArray"</tt> is read from
   * Clcm object named <tt>multicompo</tt> using
   * <pre>
   *  IntPtrConst ia = multicompo->getInt("myArray");
   *  for(int i = 0; i < multicompo->length("myArray"); ++i) 
   *     cout << "ia(" << i << ")=" << ia[i] << endl;
   * </pre>
   * @param key key identification of the integer array in the dictionary
   * @return array of integer values stored as <tt>IntPtr</tt> object.
   */
  IntPtrConst getInt(const std::string key) throw(LCMexception);

  /** Recover an integer array from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getInt</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the integer array in the list. The first list element
   * is stored at index 0.
   * @return array of integer values stored as <tt>IntPtr</tt> object.
   */
  IntPtrConst getInt(const int_32 iset) throw(LCMexception);

  /** Recover a single precision real array from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getFloat</tt>. This method is not available for file-type Clcm objects.
   * @param key key identification of the real array in the dictionary
   * @return array of single precision real values stored as <tt>FloatPtr</tt> object.
   */
  FloatPtrConst getFloat(const std::string key) throw(LCMexception);

  /** Recover a single precision real array from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getFloat</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the real array in the list. The first list element
   * is stored at index 0.
   * @return array of single precision real values stored as <tt>FloatPtr</tt> object.
   */
  FloatPtrConst getFloat(const int_32 iset) throw(LCMexception);

  /** Recover a string pointer from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the string returned
   * by <tt>getString</tt>. This method is not available for file-type Clcm objects.
   * @param key key identification of the character array in the dictionary
   * @return character information stored as <tt>StringPtr</tt> object.
   */
  StringPtrConst getString(const std::string key) throw(LCMexception);

  /** Recover a string pointer from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the string returned
   * by <tt>getString</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the integer array in the list. The first list element
   * is stored at index 0.
   * @return character information stored as <tt>StringPtr</tt> object.
   */
  StringPtrConst getString(const int_32 iset) throw(LCMexception);

  /** Recover a vector-of-string pointer from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the vector-of-string returned
   * by <tt>getVecString</tt>. This method is not available for file-type Clcm objects.
   * @param key key identification of the character array in the dictionary
   * @param size number of components in the vector-of-string
   * @return vector-of-string containing the character information stored as
   * <tt>StringVecPtr</tt> object.
   */
  StringVecPtr getVecString(const std::string key, const int_32 size) throw(LCMexception);

  /** Recover a vector-of-string pointer from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the vector-of-string returned
   * by <tt>getVecString</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the integer array in the list. The first list element
   * is stored at index 0.
   * @param size number of components in the vector-of-string
   * @return vector-of-string containing the character information stored as
   * <tt>StringVecPtr</tt> object.
   */
  StringVecPtr getVecString(const int_32 iset, const int_32 size) throw(LCMexception);

  /** Recover a double precision real array from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getDouble</tt>. This method is not available for file-type Clcm objects.
   * @param key key identification of the double precision array in the dictionary
   * @return array of double precision real values stored as <tt>DoublePtr</tt> object.
   */
  DoublePtrConst getDouble(const std::string key) throw(LCMexception);

  /** Recover a double precision real array from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getDouble</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the double precision array in the list. The first list element
   * is stored at index 0.
   * @return array of double precision real values stored as <tt>DoublePtr</tt> object.
   */
  DoublePtrConst getDouble(const int_32 iset) throw(LCMexception);

  /** Recover a boolean array from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getBool</tt>. This method is not available for file-type Clcm objects.
   * @param key key identification of the boolean array in the dictionary
   * @return array of boolean values stored as <tt>BoolPtr</tt> object.
   */
  BoolPtrConst getBool(const std::string key) throw(LCMexception);

  /** Recover a boolean array from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getBool</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the boolean array in the list. The first list element
   * is stored at index 0.
   * @return array of boolean values stored as <tt>BoolPtr</tt> object.
   */
  BoolPtrConst getBool(const int_32 iset) throw(LCMexception);

  /** Recover a complex array from a dictionary-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getComplex</tt>. This method is not available for file-type Clcm objects.
   * @param key key identification of the complex array in the dictionary
   * @return array of complex values stored as <tt>ComplexPtr</tt> object.
   */
  ComplexPtrConst getComplex(const std::string key) throw(LCMexception);

  /** Recover a complex array from a list-type Clcm object.
   * <b>General rule:</b> Never try to modify, deallocate or free the object returned
   * by <tt>getComplex</tt>. This method is not available for file-type Clcm objects.
   * @param iset index of the complex array in the list. The first list element
   * is stored at index 0.
   * @return array of complex values stored as <tt>ComplexPtr</tt> object.
   */
  ComplexPtrConst getComplex(const int_32 iset) throw(LCMexception);

  /** Recover a vector-of-string pointer containing the keys of a dictionary-type Clcm
   * object. <b>General rule:</b> Never try to modify, deallocate or free the vector-of-string
   * returned by <tt>keys</tt>. This method is not available for file-type or list-type
   * Clcm objects.
   * @return <tt>StringVecPtr</tt> of a vector-of-string containing the dictionary keys
   */
  StringVecPtr keys() throw(LCMexception);

  /** Store an integer array in a dictionary-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * <p>
   * <b>Example</b>: An integer array named <tt>"myIntArray"</tt> is inserted in
   * Clcm object named <tt>multicompo</tt> using
   * <pre>
   *  static int_32 myStaticArray[] = {33, 22, 11, 89};
   *  IntPtr myIntPtr = IntPtr(new int_32[4]);
   *  for(int i = 0; i < 4; ++i) myIntPtr[i] = myStaticArray[i];
   *  multicompo->put("myIntArray", myIntPtr, 4);
   * </pre>
   * @param key key identification of the block in the dictionary
   * @param myArray integer array stored as <tt>IntPtr</tt> object.
   * @param myLength number of components in integer array.
   */
  void put(const std::string key, IntPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store an integer array in a list-type  Clcm object. This method is not available
   * for file-type Clcm objects.
   * @param iset index of the integer array in the list. The first list element
   * is stored at index 0.
   * @param myArray integer array stored as <tt>IntPtr</tt> object.
   * @param myLength number of components in integer array.
   */
  void put(const int_32 iset, IntPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a single precision real array in a dictionary-type Clcm object. This method
   * is not available for file-type Clcm objects.
   * @param key key identification of the block in the dictionary
   * @param myArray single precision real array stored as <tt>FloatPtr</tt> object.
   * @param myLength number of components in real array.
   */
  void put(const std::string key, FloatPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a single precision real array in a list-type Clcm object. This method is
   * not available for file-type Clcm objects.
   * @param iset index of the real array in the list. The first list element
   * is stored at index 0.
   * @param myArray single precision real array stored as <tt>FloatPtr</tt> object.
   * @param myLength number of components in real array.
   */
  void put(const int_32 iset, FloatPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a string pointer in a dictionary-type Clcm object. This method is not available
   * for file-type Clcm objects.
   * @param key key identification of the character array in the dictionary
   * @param myArray character information stored as <tt>StringPtr</tt> object.
   */
  void put(const std::string key, StringPtr myArray) throw(LCMexception);

  /** Store a string pointer in a list-type Clcm object. This method is not available for
   * file-type Clcm objects.
   * @param iset index of the string in the list. The first list element is stored at index 0.
   * @param myArray character information stored as <tt>StringPtr</tt> object.
   */
  void put(const int_32 iset, StringPtr myArray) throw(LCMexception);

  /** Store a vector-of-string pointer in a dictionary-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param key key identification of the character array in the dictionary
   * @param myArray vector-of-string containing the character information stored as
   * <tt>StringVecPtr</tt> object.
   */
  void put(const std::string key, StringVecPtr myArray) throw(LCMexception);

  /** Store a vector-of-string pointer in a list-type Clcm object.
   * This method is not available for file-type Clcm objects.
   * @param iset index of the character array in the list. The first list element
   * is stored at index 0.
   * @param myArray vector-of-string containing the character information stored as
   * <tt>StringVecPtr</tt> object.
   */
  void put(const int_32 iset, StringVecPtr myArray) throw(LCMexception);

  /** Store a double precision real array in a dictionary-type Clcm object. This method is
   * not available for file-type Clcm objects.
   * @param key key identification of the double precision array in the dictionary
   * @param myArray double precision array stored as <tt>DoublePtr</tt> object.
   * @param myLength number of components in double precision array.
   */
  void put(const std::string key, DoublePtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a double precision real array in a list-type Clcm object. This method is
   * not available for file-type Clcm objects.
   * @param iset index of the double precision array in the list. The first list element
   * is stored at index 0.
   * @param myArray double precision array stored as <tt>DoublePtr</tt> object.
   * @param myLength number of components in double precision array.
   */
  void put(const int_32 iset, DoublePtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a boolean array in a dictionary-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * @param key key identification of the block in the dictionary
   * @param myArray boolean array stored as <tt>BoolPtr</tt> object.
   * @param myLength number of components in boolean array.
   */
  void put(const std::string key, BoolPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a boolean array in a list-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * @param iset index of the boolean array in the list. The first list element
   * is stored at index 0.
   * @param myArray boolean array stored as <tt>BoolPtr</tt> object.
   * @param myLength number of components in boolean array.
   */
  void put(const int_32 iset, BoolPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a complex array in a dictionary-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * <p>
   * <b>Example</b>: A complex array named <tt>"myComplArray"</tt> is inserted in
   * Clcm object named <tt>multicompo</tt> using
   * <pre>
   *  ComplexPtr myComplexPtr = ComplexPtr(new complex<float_32>[2]);
   *  myComplexPtr[0] = complex<float_32>(3.3, 8.9);
   *  myComplexPtr[1] = complex<float_32>(-3.3, 7.9);
   *  multicompo->put("myComplArray", myComplexPtr, 2);
   * </pre>
   * @param key key identification of the block in the dictionary
   * @param myArray complex array stored as <tt>ComplexPtr</tt> object.
   * @param myLength number of components in complex array.
   */
  void put(const std::string key, ComplexPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a complex array in a list-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * @param iset index of the complex array in the list. The first list element
   * is stored at index 0.
   * @param myArray complex array stored as <tt>ComplexPtr</tt> object.
   * @param myLength number of components in complex array.
   */
  void put(const int_32 iset, ComplexPtr myArray, const int_32 myLength) throw(LCMexception);

  /** Store a daughter Clcm object in a dictionary-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * @param key key identification of the block in the dictionary
   * @param myClcm daughter Clcm object.
   */
  void put(const std::string key, ClcmPtr myClcm) throw(LCMexception);

  /** Store a daughter Clcm object in a list-type Clcm object. This method is not
   * available for file-type Clcm objects.
   * @param iset index of the complex array in the list. The first list element
   * is stored at index 0.
   * @param myClcm daughter Clcm object.
   */
  void put(const int_32 iset, ClcmPtr myClcm) throw(LCMexception);

  /** Extract the LCM structure.
   * @return ANSI C pointer of the embedded LCM structure.
   */
  lcm *extract();

private:
  std::vector<AnyPtr> *global_list; // container for the global references
  lcm *addr;                  // address of the LCM object
  Clcm *root;                 // address of Clcm root object
  std::ofstream aFile;        // embedded file if type >= 3
  int_32 lrda;                // used with direct access files only
  std::string name;           // object name
  std::string path;           // path to access the Clcm file to import
  std::string nammy;          // directory name inside LCM object
  bool empty;                 // empty flag
  int_32 ilong;               // -1 (dictionary) or list length
  int_32 objtype;             // type of the object =1: LCM; =2: XSM; =3: binary Fortran file
                              // =4: ASCII Fortran file; =5: direct access Fortran file
  int_32 access;              // =0: closed object; =1: object in read/write mode
                              // =2: object in read-only mode
  int_32 version;             // original ganlib version of object
  bool untoutched;            // brand new object flag
  bool isroot;                // true if addr is pointing on the root object
}; // class Clcm
  std::ostream & operator<<(std::ostream &, ClcmPtr) throw(LCMexception);

} // namespace ganlib
#endif
