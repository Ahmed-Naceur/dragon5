#include "THM.hxx"
using namespace boost;  using namespace std;  using namespace ganlib;

THM::THM() {
  cout << "New THM object constructed.'" << endl;
}

void THM::initialize(void* component)
{
communicator_.initialize(component);
}

void THM::run()
{
   // construct the Lifo stack for IniThmCompo
   LifoPtr ipLifo1 = LifoPtr(new Lifo());
   ipLifo1->pushEmpty("Fmap", "LCM");
   ipLifo1->pushEmpty("Matex", "LCM");

   // call IniComponent Cle-2000 procedure
   Cle2000Ptr IniComponent = Cle2000Ptr(new Cle2000("IniThmCompo", 0, ipLifo1));
   IniComponent->exec();
   cout << "IniThmCompo execution completed" << endl;

   // recover the output LCM objects
   ClcmPtr Fmap; ipLifo1->node("Fmap", Fmap);
   ClcmPtr Matex; ipLifo1->node("Matex", Matex);
   IntPtrConst stateVector = Fmap->getInt("STATE-VECTOR");
   long mylength = stateVector[0] * stateVector[1];
   long npar = stateVector[7];

   // empty the Lifo stack
   while (ipLifo1->getMax() > 0) ipLifo1->pop();
   
   //iteration loop
   ClcmPtr Flux, Thm;
   float_32 densB = 2000.;
   int iter = 0;
   int continueLoop = 1;
   LifoPtr ipLifo3 = LifoPtr(new Lifo());
   Cle2000Ptr ThmComponent = Cle2000Ptr(new Cle2000("ThmComponent", 0, ipLifo3));

   while (continueLoop == 1) {
     ++ iter;
     cout << "THM: ITERATION NUMBER:" << iter << " continueLoop=" << continueLoop << endl;

     FloatPtr powerTab(new float_32[mylength]);
     communicator_.recv(iter, "powerTab", mylength, powerTab );
     Fmap->put("BUND-PW", powerTab, mylength);
     FloatPtr irradiationTab(new float_32[mylength]);
     communicator_.recv(iter, "irradiationTab", mylength, irradiationTab );
     Fmap->put("BURN-INST", irradiationTab, mylength);
     Fmap->val();
	
     // construct the Lifo stack for ThmComponent
     ipLifo3->push("Fmap", Fmap);
     if (iter == 1) {
       ipLifo3->pushEmpty("Thm", "LCM");
     } else {
       ipLifo3->push("Thm", Thm);
     }
     ipLifo3->push("iter", iter);
     ipLifo3->push("densB", densB);
     ipLifo3->pushEmpty("CONV", "B");

     // call ThmComponent Cle-2000 procedure
     ThmComponent->exec();
     cout << "ThmComponent execution completed." << endl;
	
     // recover and send thermo-hydraulics information
     ipLifo3->node("Thm", Thm);
     ClcmPtr Jpmap = Fmap->getClcm("PARAM");
     for (int ipar=0; ipar<npar; ipar++) {
       ClcmPtr Kpmap = Jpmap->getClcm(ipar);
       StringPtrConst pname = Kpmap->getString("P-NAME");
       if (*pname == "C-BORE      ") continue;
       IntPtrConst ptype = Kpmap->getInt("P-TYPE");
       if (ptype[0] != 2) throw Cle2000Exception("THM failure");
       FloatPtrConst myArray = Kpmap->getFloat("P-VALUE");
       if (*pname == "T-FUEL      ") {
         communicator_.send(iter, "fuel_temperature", mylength, myArray);
       } else if (*pname == "D-COOL      ") {
         communicator_.send(iter, "water_density", mylength, myArray);
       } else if (*pname == "T-COOL      ") {
         communicator_.send(iter, "water_temperature", mylength, myArray);
       }
     }

     // recover convergence flag on top of Lifo object
     bool CONV;
     ipLifo3->pop(CONV);
     if (CONV) continueLoop = 0;

     // empty the Lifo stack
     while (ipLifo3->getMax() > 0) ipLifo3->pop();
	
     communicator_.send(iter, "continueLoop", continueLoop);
     cout << "THM: Value of continueLoop : "<< continueLoop << " at iteration " << iter << endl;
   }
   cout << "THM: close the Calcium communicator" << endl ;
   communicator_.terminate();
}
