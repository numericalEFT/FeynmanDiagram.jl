using FeynmanDiagram

#a,b,c,d,e = set_variables(['a','b','c','d','e'],order = 2)
a, b, c, d, e = set_variables("a b c d e", order=2)
print(a + b)