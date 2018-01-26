q0_list = [0.,0.,0.]
q1_list = [1.,2.,3.]
q2_list = [4.,5.,6.]
q3_list = [7.,8.,9.]


d1 = {}

d1.update({"q0" : q0_list})
d1.update({"q1" : q1_list})
d1.update({"q2" : q2_list})
d1.update({"q3" : q3_list})

ransl = []

for i in range(4):
    ransl.append("q"+str(i))
	
eh = []
	
for i in range(2):
    ii = float(i)
    q0_list = [ii*0.,0.,0.]
    q1_list = [ii*1.,2.,3.]
    q2_list = [ii*4.,5.,6.]
    q3_list = [ii*7.,8.,9.]

    d1 = {}

    d1.update({"q0" : q0_list})
    d1.update({"q1" : q1_list})
    d1.update({"q2" : q2_list})
    d1.update({"q3" : q3_list})

    field = [[data for data in d1[s]] for s in ransl]
    eh.append(field)


#print(d1['q2'])
#plt.plot(eh[:][nt-1][2])
#print(eh)
print(eh[:][1][1])