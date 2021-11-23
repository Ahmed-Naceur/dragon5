/**
 * This class is an implementation of the C++ bindings for a Lifo
 * stack, as required by the Cle-2000 class.
 * Lifo capabilities are available for a program written in C++
 * by using methods belonging to the Lifo class.
 * <P> A Lifo object in C++ can encapsulate a native Lifo stack
 * using the GANLIB5 API in ANSI C
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2012)
 */
#ifndef Lifo_HXX
#define Lifo_HXX

#include <cstddef> // for __GLIBCXX__
 
#ifdef __GLIBCXX__ // for memory management with shared_ptr
#  include <tr1/memory>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <memory>
#endif

#include "LifoException.hxx"
#include "Clcm.hxx"
extern "C" {
#include <dirent.h>
#include "cle2000.h"
}
#define LifoPtr boost::shared_ptr<Lifo>

namespace ganlib {

/**
 * This class is an implementation of the C++/Boost bindings for a last-in-first-out
 * (lifo) stack used with CLE-2000. Lifo management capabilities for a program
 * written in C++ are available by using methods belonging to the Lifo class.
 * These methods encapsulate the lifo API calls used as "extern"C" functions.
 * <P> A Lifo object in C++ can encapsulate a native lifo stack used to
 * manage CLE-2000 parameters.
 * <P> A lifo stack can contain defined or undefined (empty) nodes; used to represent
 * known or unknown parameters, respectively.
 * <P> <I>Note:</I> There is a constraint_32 in CLE-2000. LCM (or XSM) objects and files
 * must be pushed before single-value nodes in the stack.
 * <P>
 *
 * @author Alain Hebert, Ecole Polytechnique de Montreal (2012)
 */
class Lifo {
public:
  /**
   * Use this constructor to create an empty Lifo object.
   */
  Lifo() throw(LifoException);

  /** Close and destroy a Lifo object.
   */
  ~Lifo() throw(LifoException);

  /** Pop the value on top of the Lifo stack without recovering it. If the node is empty, an
   * exception is thrown
   */
  void pop() throw(LifoException);

  /** Pop the integer value on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myInteger integer value node to pop from the lifo stack.
   */
  void pop(int_32 &myInteger) throw(LifoException);

  /** Pop the real value on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myFloat real value node to pop from the lifo stack.
   */
  void pop(float_32 &myFloat) throw(LifoException);

  /** Pop the character string on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myString character string node (limited to 72 characters) to pop from the lifo stack.
   */
  void pop(std::string &myString) throw(LifoException);

  /** Pop the double precision value on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myDouble double precision value node to pop from the lifo stack.
   */
  void pop(double &myDouble) throw(LifoException);

  /** Pop the boolean value on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myBool boolean value node to pop from the lifo stack.
   */
  void pop(bool &myBool) throw(LifoException);

  /** Pop the ClcmPtr object on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myClcm LCM object node to pop from the lifo stack.
   */
  void pop(ClcmPtr &myClcm) throw(LifoException);

  /** Pop the file object on top of the Lifo stack. If the node is empty or if the argument type
   * is wrong, an exception is thrown
   * @param myFile operating system (OS) name (limited to 72 characters) associated to the file object
   * node to pop from the lifo stack.
   * @param stype type of file. This variable is selected by <tt>pop</tt> among the
   * following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * <dt> <tt>"DA"</tt> <dd> binary direct-access file (with 128-word records)</ul>
   */
  void pop(std::string myFile, std::string stype) throw(LifoException);

  /** Return the integer value node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param myInteger integer value node.
   */
  void node(int_32 ipos, int_32 &myInteger) throw(LifoException);

  /** Return the real value node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param myFloat real value node.
   */
  void node(int_32 ipos, float_32 &myFloat) throw(LifoException);

  /** Return the character string node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param myString character string node (limited to 72 characters).
   */
  void node(int_32 ipos, std::string &myString) throw(LifoException);

  /** Return the double precision value node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param myDouble double precision value node.
   */
  void node(int_32 ipos, double &myDouble) throw(LifoException);

  /** Return the boolean value node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param myBool integer value node.
   */
  void node(int_32 ipos, bool &myBool) throw(LifoException);

  /** Return the ClcmPtr object node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param myClcm ClcmPtr object node.
   */
  void node(int_32 ipos, ClcmPtr &myClcm) throw(LifoException);

  /** Return the file object node with a given position. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @param node operating system (OS) name (limited to 72 characters) associated to the file object
   * node to pop from the lifo stack.
   * @param stype type of LCM object or file. This variable is selected by <tt>pop</tt> among the
   * following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * <dt> <tt>"DA"</tt> <dd> binary direct-access file (with 128-word records)</ul>
   */
  void node(int_32 ipos, std::string node, std::string stype) throw(LifoException);

  /** Return the integer value node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param myInteger integer value node.
   */
  void node(std::string sname, int_32 &myInteger) throw(LifoException);

  /** Return the real value node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param myFloat real value node.
   */
  void node(std::string sname, float_32 &myFloat) throw(LifoException);

  /** Return the character string node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param myString character string node (limited to 72 characters).
   */
  void node(std::string sname, std::string &myString) throw(LifoException);

  /** Return the double precision value node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param myDouble double precision value node.
   */
  void node(std::string sname, double &myDouble) throw(LifoException);

  /** Return the boolean value node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param myBool integer value node.
   */
  void node(std::string sname, bool &myBool) throw(LifoException);

  /** Return the ClcmPtr object node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param myClcm ClcmPtr object node.
   */
  void node(std::string sname, ClcmPtr &myClcm) throw(LifoException);

  /** Return the file object node with a given name. If the node is empty or if the argument type
   * is wrong, an exception is thrown. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @param node operating system (OS) name (limited to 72 characters) associated to the file object
   * node to pop from the lifo stack.
   * @param stype type of LCM object or file. This variable is selected by <tt>pop</tt> among the
   * following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * <dt> <tt>"DA"</tt> <dd> binary direct-access file (with 128-word records)</ul>
   */
  void node(std::string sname, std::string node, std::string stype) throw(LifoException);

  /** Return the OSname of the node associated with a given name. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @return OSName of the file
   */
  std::string OSName(std::string sname) throw(LifoException);

  /** Return the OSname of the node at a given position. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @return OSName of the file
   */
  std::string OSName(int_32 ipos) throw(LifoException);

  /** Return the name of the node at a given position. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @return name of the node
   */
  std::string Name(int_32 ipos) throw(LifoException);

  /** Return the type of the node associated with a given name. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @return node type (3= LCM object; 4= XSM file; 5= seq binary; 6= seq ascii; 7= DA binary;
   11= integer value; 12= real value; 13= character std::string; 14= double precision value;
   15= logical value). A negative value indicates an empty node.
   */
  int_32 typeNode(std::string sname) throw(LifoException);

  /** Return the type of the node at a given position. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @return node type (3= LCM object; 4= XSM file; 5= seq binary; 6= seq ascii; 7= DA binary;
   11= integer value; 12= real value; 13= character std::string; 14= double precision value;
   15= logical value). A negative value indicates an empty node.
   */
  int_32 typeNode(int_32 ipos) throw(LifoException);

  /** Return the access of the node associated with a given name. The lifo stack is not modified
   * @param sname node name (limited to 12 characters)
   * @return node access (0= creation; 1:modification; 2=read-only).
   */
  int_32 accessNode(std::string sname) throw(LifoException);

  /** Return the access of the node at a given position. The lifo stack is not modified
   * @param ipos node position in lifo stack
   * @return node access (0= creation; 1:modification; 2=read-only).
   */
  int_32 accessNode(int_32 ipos) throw(LifoException);

  /** Push a new integer on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myInteger integer value to push into the lifo stack
   */
  void push(std::string sname, const int_32 myInteger) throw(LifoException);

  /** Push a new real on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myFloat real value to push into the lifo stack
   */
  void push(std::string sname, const float myFloat) throw(LifoException);

  /** Push a new node object on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myString string object to push into the lifo stack.
   */
  void push(std::string sname, const std::string myString) throw(LifoException);

  /** Push a new double precision value on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myDouble double precision value to push into the lifo stack
   */
  void push(std::string sname, const double_64 myDouble) throw(LifoException);

  /** Push a new boolean value on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myBool boolean value to push into the lifo stack
   */
  void push(std::string sname, const bool myBool) throw(LifoException);

  /** Push a new ClcmPtr node object on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myClcm ClcmPtr LCM/XSM object to push into the lifo stack.
   * result
   */
  void push(std::string sname, const ClcmPtr myClcm) throw(LifoException);

  /** Push a new node object with a type on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myFile file object to push into the lifo stack.
   * @param stype type of LCM object or file. This variable is chosen among
   * the following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * <dt> <tt>"DA"</tt> <dd> binary direct-access file (with 128-word records)
   */
  void push(std::string sname, std::string myFile, std::string stype) throw(LifoException);

  /** Push a new node object with a type and an <tt>"OSname"</tt> on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param myFile file object to push into the lifo stack.
   * @param stype type of LCM object or file. This variable is chosen among
   * the following values:
   * <ul><dl>
   * <dt> <tt>"BINARY"</tt> <dd> binary sequential file
   * <dt> <tt>"ASCII"</tt> <dd> ACSII sequential file
   * <dt> <tt>"DA"</tt> <dd> binary direct-access file (with 128-word records)</ul>
   * @param OSname operating system (OS) name associated to the file (limited to 72 characters).
   */
  void push(std::string sname, std::string myFile, std::string stype,
            std::string OSname) throw(LifoException);

  /** Push an empty node on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param nodeType node type (='I': integer value; 'F': real value; 'D': double precision value;
   * 'B': boolean value;'S': character string value; ='LCM': memory-resident LCM object; ='XSM': XSM file;
   * ='BINARY': binary file; ='ASCII': ascii file; ='DA': direct-access file (with 128-word records))
   */
  void pushEmpty(std::string sname, std::string nodeType) throw(LifoException);

  /** Push an empty node on top of the Lifo stack
   * @param sname node name (limited to 12 characters)
   * @param nodeType node type (='XSM': XSM file; ='BINARY': binary file; ='ASCII': ascii file;
   * ='DA': direct-access file (with 128-word records))
   * @param OSname operating system (OS) name associated to the file (limited to 72 characters).
   */
  void pushEmpty(std::string sname, std::string nodeType, std::string OSname) throw(LifoException);

  /** Gives the number of nodes in the lifo stack
   * @return number of nodes in lifo stack
   */
  int_32 getMax();

  /** Print_32 the table-of-content of a lifo stack
   */
  void lib() throw(LifoException);

  /** Extract the lifo structure.
   * @return ANSI C pointer of the embedded lifo structure.
   */
  lifo *extract();

private:
  std::vector<ClcmPtr> *global_list; // container for the global references
  lifo *addr;                  // address of the Lifo object
}; // class Lifo

} // namespace ganlib
#endif
