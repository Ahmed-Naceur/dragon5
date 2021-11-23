#include "POW.hxx"
using namespace boost;  using namespace std;  using namespace ganlib;

POW::POW() {
  cout << "New POW object constructed.'" << endl;
}

void POW::initialize(double power, void* component)
{
power_ = power;
communicator_.initialize(component);
}

void POW::run()
{
   // construct the Lifo stack for IniPowCompo
   cout << "POW::run" << endl;
   LifoPtr ipLifo1 = LifoPtr(new Lifo());
   ipLifo1->pushEmpty("Fmap", "LCM");
   ipLifo1->pushEmpty("Matex", "LCM");
   ipLifo1->pushEmpty("Cpo", "LCM");
   ipLifo1->pushEmpty("Track", "LCM");

   // call IniPowCompo Cle-2000 procedure
   Cle2000Ptr IniPowCompo = Cle2000Ptr(new Cle2000("IniPowCompo", 0, ipLifo1));
   IniPowCompo->exec();
   cout << "IniPowCompo execution completed" << endl;

   // recover the output LCM objects
   ClcmPtr Fmap; ipLifo1->node("Fmap", Fmap);
   ClcmPtr Matex; ipLifo1->node("Matex", Matex);
   ClcmPtr Cpo; ipLifo1->node("Cpo", Cpo);
   ClcmPtr Track; ipLifo1->node("Track", Track);
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
   float Keff_conv = 1.0;
   LifoPtr ipLifo2 = LifoPtr(new Lifo());
   Cle2000Ptr PowComponent = Cle2000Ptr(new Cle2000("PowComponent", 0, ipLifo2));
   
   while (continueLoop == 1) {
     ++ iter;
     cout << "POW: ITERATION NUMBER:" << iter << " continueLoop=" << continueLoop << endl;

     // construct the Lifo stack for PowComponent
     ipLifo2->push("Fmap", Fmap);
     ipLifo2->push("Matex", Matex);
     if (iter == 1) {
       ipLifo2->pushEmpty("Flux", "LCM");
     } else {
       ipLifo2->push("Flux", Flux);
     }
     ipLifo2->push("Cpo", Cpo);
     ipLifo2->push("Track", Track);
     ipLifo2->push("iter", iter);
     ipLifo2->push("powi", float_32(power_));
     ipLifo2->push("densB", densB);

     // call PowComponent Cle-2000 procedure
     cout << "call PowComponent->exec" << endl;
     PowComponent->exec();
     cout << "PowComponent execution completed" << endl;
     ipLifo2->node("Flux", Flux);
     FloatPtrConst Keff = Flux->getFloat("K-EFFECTIVE");
     cout << "POW: iter=" << iter << " ------------- Keffective=" << Keff[0] << endl;
     Keff_conv = Keff[0];

     // send reactor physics information
     FloatPtrConst powerTab = Fmap->getFloat("BUND-PW");
     communicator_.send(iter, "powerTab", mylength, powerTab);
     FloatPtrConst irradiationTab = Fmap->getFloat("BURN-INST");
     communicator_.send(iter, "irradiationTab", mylength, irradiationTab);

     // receive thermo-hydraulics information
     ClcmPtr Jpmap = Fmap->getClcm("PARAM");
     IntPtr myIntPtr(new int_32[1]); myIntPtr[0] = 2;
     for (int ipar=0; ipar<npar; ipar++) {
       ClcmPtr Kpmap = Jpmap->getClcm(ipar);
       StringPtrConst pname = Kpmap->getString("P-NAME");
       if (*pname == "T-FUEL      ") {
         FloatPtr myArray(new float_32[mylength]);
         communicator_.recv(iter, "fuel_temperature", mylength, myArray);
         Kpmap->put("P-VALUE", myArray, mylength);
         Kpmap->put("P-TYPE", myIntPtr, 1);
       } else if (*pname == "D-COOL      ") {
         FloatPtr myArray(new float_32[mylength]);
         communicator_.recv(iter, "water_density", mylength, myArray);
         Kpmap->put("P-VALUE", myArray, mylength);
         Kpmap->put("P-TYPE", myIntPtr, 1);
       } else if (*pname == "T-COOL      ") {
         FloatPtr myArray(new float_32[mylength]);
         communicator_.recv(iter, "water_temperature", mylength, myArray);
         Kpmap->put("P-VALUE", myArray, mylength);
         Kpmap->put("P-TYPE", myIntPtr, 1);
       }
     }
     Fmap->val();

     // empty the Lifo stack
     while (ipLifo2->getMax() > 0) ipLifo2->pop();

     // receive the convergence flag
     communicator_.recv(iter, "continueLoop", continueLoop);
     cout << "POW: Value of continueLoop : "<< continueLoop << " at iteration " << iter << endl;
   }
   cout << "POW: close the Calcium communicator" << endl ;
   communicator_.terminate();
   cout.precision(10);
   cout << "POW: converged K-effective=" << Keff_conv << endl ;
}
