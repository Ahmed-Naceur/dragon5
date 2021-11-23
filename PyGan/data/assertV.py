#
# Assert procedure for non-regression testing
# Recover a value from a list of real arrays
# Author: A. Hebert
#
def assertV(lcmobj,key,iset,ipos,refvalue):
  lcmobj.lib()
  val=lcmobj[key][iset-1][ipos-1]
  delta=abs((val-refvalue)/refvalue)
  if delta < 1.0e-4:
    print("Test successful; delta=",delta)
  else:
    print("Reference=",refvalue,"Calculated=",val)
    print("------------")
    print("TEST FAILURE")
    print("------------")
    raise Exception("abort in assertV")
