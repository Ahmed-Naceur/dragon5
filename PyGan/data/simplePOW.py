#
# simplePOW: a simple multiphysics example with THM: module
#
import lifo
import lcm
import cle2000
from assertS import *
import numpy as np

# construct the Lifo stack for IniPowCompo
ipLifo1=lifo.new()
ipLifo1.pushEmpty("Fmap", "LCM")
ipLifo1.pushEmpty("Matex", "LCM")
ipLifo1.pushEmpty("Cpo", "LCM")
ipLifo1.pushEmpty("Track", "LCM")

# call IniPowCompo Cle-2000 procedure
IniPowCompo = cle2000.new('IniPowCompo',ipLifo1,1)
IniPowCompo.exec()
print("IniPowCompo execution completed")

# recover the output LCM objects
Fmap = ipLifo1.node("Fmap")
Matex = ipLifo1.node("Matex")
Cpo = ipLifo1.node("Cpo")
Track = ipLifo1.node("Track")
stateVector = Fmap["STATE-VECTOR"]
mylength = stateVector[0]*stateVector[1]
npar = stateVector[7]

# empty the Lifo stack
while ipLifo1.getMax() > 0:
  ipLifo1.pop();

# iteration loop
iter = 0
continueLoop = 1
powi = 17.3 # Reference at 17.3 MW
densB = 2000.0
ipLifo2 = lifo.new()
ipLifo3 = lifo.new()
PowComponent = cle2000.new('PowComponent',ipLifo2,1)
ThmComponent = cle2000.new('ThmComponent',ipLifo3,1)
conv = False
while not conv:
  iter += 1
  if iter > 5:
    raise Exception("simplePOW: maximum number of iterations is reached")
    
  print("POW: ITERATION NUMBER:", iter)

  # construct the Lifo stack for PowComponent
  ipLifo2.push(Fmap);
  ipLifo2.push(Matex);
  if iter == 1:
    Flux = ipLifo2.pushEmpty("Flux", "LCM")
  else:
    ipLifo2.push(Flux)

  ipLifo2.push(Cpo)
  ipLifo2.push(Track)
  ipLifo2.push(iter)
  ipLifo2.push(powi)
  ipLifo2.push(densB)
   
  # call PowComponent Cle-2000 procedure
  print("call PowComponent procedure")
  PowComponent.exec()
  print("PowComponent execution completed")
  Flux = ipLifo2.node("Flux")
  Keff_conv = Flux["K-EFFECTIVE"][0]
  print("POW: iter=", iter, " ------------- Keffective=", Keff_conv)

  # construct the Lifo stack for ThmComponent
  ipLifo3.push(Fmap);
  if iter == 1:
    Thm =ipLifo3.pushEmpty("Thm", "LCM");
  else:
    ipLifo3.push(Thm);

  ipLifo3.push(iter)
  ipLifo3.push(densB)
  ipLifo3.pushEmpty("CONV", "B")

  # call ThmComponent Cle-2000 procedure
  print("call ThmComponent procedure")
  ThmComponent.exec()
  conv = ipLifo3.node("CONV")
      
  print("ThmComponent execution completed. conv=", conv)

  # recover thermo-hydraulics information
  Thm = ipLifo3.node("Thm")
  Jpmap = Fmap["PARAM"];
  myIntPtr = np.array([2,], dtype='i')
  for ipar in range(0, npar):
    Kpmap = Jpmap[ipar]
    pname = Kpmap["P-NAME"]
    if pname == "T-FUEL":
      continue
    ptype = Kpmap["P-TYPE"]
    myArray = Kpmap["P-VALUE"]
    if pname == "T-FUEL":
      Kpmap.put("P-VALUE", myArray, mylength);
      Kpmap.put("P-TYPE", myIntPtr, 1);
    elif pname == "D-COOL":
      Kpmap.put("P-VALUE", myArray, mylength);
      Kpmap.put("P-TYPE", myIntPtr, 1);
    elif pname == "T-COOL":
      Kpmap.put("P-VALUE", myArray, mylength);
      Kpmap.put("P-TYPE", myIntPtr, 1);

  Fmap.val()

  # empty the ipLifo2 Lifo stack
  while ipLifo2.getMax() > 0:
    ipLifo2.pop();

  # empty the ipLifo3 Lifo stack
  while ipLifo3.getMax() > 0:
    ipLifo3.pop();

print("POW: converged K-effective=", Keff_conv)
assertS(Flux,'K-EFFECTIVE',0,1.011134)
print("test simplePOW completed")
