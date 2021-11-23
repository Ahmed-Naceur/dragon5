#
# fact: non regression testing for lifo and pygan classes
#
import lifo
import cle2000
my_lifo=lifo.new()
my_lifo.push(5)
my_lifo.push(int)
my_lifo.lib()
my_fact=cle2000.new('fact',my_lifo,1)
my_fact.exec()
print("factorial(5)=", my_lifo.node(1))
print("test fact completed")
