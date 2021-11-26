#include "FACT.hxx"
#include "Cle2000.hxx"
using namespace boost;  using namespace std;  using namespace ganlib;

FACT::FACT() {
  cout << "New FACT object constructed.'" << endl;
}

void FACT::run(long a, long& b) {
   // construct the lifo stack
   LifoPtr ipLifo = LifoPtr(new Lifo());
   ipLifo->push("input_val",int_32(a));
   ipLifo->pushEmpty("output_val","I");

   // call the procedure with in-out CLE-2000 variables
   Cle2000Ptr ipCle2000 = Cle2000Ptr(new Cle2000("fact", 0, ipLifo));
   ipCle2000->exec();

   // recover and erase the lifo stack
   int_32 output_val; ipLifo->pop(output_val);
   b = (long)output_val;

   // empty the lifo stack
   while (ipLifo->getMax() > 0) ipLifo->pop();
   
   // end of case -----------------------------------------------------------
   cout << "FACT::run: successful end of execution" << endl;
}
