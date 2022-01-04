#
# Assert procedure for non-regression testing
# Recover a value from a real array
# Author: A. Hebert
#
def assertS(lcmobj,key,ipos,refvalue):
  val=lcmobj[key][ipos-1]
  delta=abs((val-refvalue)/refvalue)
  if delta < 1.0e-4:
    print("Test successful; delta=",delta)
  else:
    print("Reference=",refvalue,"Calculated=",val)
    print("------------")
    print("TEST FAILURE")
    print("------------")
    raise Exception("abort in assertS")
