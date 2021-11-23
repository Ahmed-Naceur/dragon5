#include "Cle2000.hxx"
using namespace boost;  using namespace std;  using namespace ganlib;

int main(int argc, char **argv) {
   // construct the lifo stack
   LifoPtr ipLifo = LifoPtr(new Lifo());
   ipLifo->push("input_val",5);
   ipLifo->pushEmpty("output_val","I");

   // call the procedure with in-out CLE-2000 variables
   Cle2000Ptr ipCle2000 = Cle2000Ptr(new Cle2000("fact", 0, ipLifo));
   ipCle2000->exec();

   // recover and erase the lifo stack
   int_32 output_val; ipLifo->pop(output_val);
   cout << endl << "factorial(5)=" << output_val << endl << endl;
   if (output_val != 120) throw Cle2000Exception("Test failure 1");

   // empty the lifo stack
   while (ipLifo->getMax() > 0) ipLifo->pop();
   
   // end of case -----------------------------------------------------------

   // construct a new lcm object
   ClcmPtr iplistgen = ClcmPtr(new Clcm("LCM_IMP_ASCII", "Geo3D.txt", "./"));
   iplistgen->lib();
   StringPtrConst myString = iplistgen->getString("SIGNATURE");
   cout << "signature=" << *myString << "<---" << endl;
   iplistgen->close("KEEP");

   // construct a new lifo stack
   LifoPtr ipLifo2 = LifoPtr(new Lifo());
   ipLifo2->push("Geo3D",iplistgen);
   ipLifo2->pushEmpty("nb_mix","I");
   ipLifo2->lib();

   // call the procedure with in-out CLE-2000 variables
   Cle2000Ptr ipCle2000_bis = Cle2000Ptr(new Cle2000("grep", 0, ipLifo2));
   ipCle2000_bis->exec();

   // recover and erase the lifo stack
   int_32 nb_mix; ipLifo2->pop(nb_mix);
   cout << endl << "number of mixtures = " << nb_mix << endl << endl;
   if (nb_mix != 17) throw Cle2000Exception("Test failure 2");

   // empty the lifo stack
   while (ipLifo2->getMax() > 0) ipLifo2->pop();

   // end of case -----------------------------------------------------------

   // construct the lifo stack
   LifoPtr ipLifo3 = LifoPtr(new Lifo());
   ipLifo3->pushEmpty("saphyb", "XSM", "./saphyb_UOX_AICN");
   ipLifo3->push("draglib", string("draglibendfb7r1SHEM295"));
   ipLifo3->push("type", string("FRA"));
   ipLifo3->push("combustible", string("UOX"));
   ipLifo3->push("barres", string("AICN"));
   ipLifo3->push("homoge", string("HOMOGENE"));
   ipLifo3->push("u235", float_32(3.7));
   ipLifo3->push("burnup", float_32(85000.));
   ipLifo3->push("nrgoup", int_32(2));
   ipLifo3->lib();

   // call the procedure with in-out CLE-2000 variables
   Cle2000Ptr ipCle2000_ter = Cle2000Ptr(new Cle2000("pwr2010", 0, ipLifo3));
   ipCle2000_ter->exec();

   // recover the saphyb
   ClcmPtr xsm_object; ipLifo3->node("saphyb", xsm_object);
   xsm_object->open("READ-ONLY");
   xsm_object->lib();
   xsm_object->close("KEEP");

   // empty the lifo stack
   while (ipLifo3->getMax() > 0) ipLifo3->pop();

   // end of case -----------------------------------------------------------
   cout << "Skin++: normal end of execution" << endl;
   return 0;
}
