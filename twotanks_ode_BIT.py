import numpy as np
import matplotlib.pyplot as plt
import sys
#from numpy import dot, array

# GIT TEST
# two-tanks mixing problem
# http://www.ees.nmt.edu/outside/courses/hyd510/PDFs/Lecture%20notes/Lectures%20Part%202.7%20SimultODE.pdf

class twotanks:

    def __init__(self,a11,a12,a21,a22,dt,tend):
	# initialize two tanks problem with model parameters
	    self.a11 = a11
	    self.a12 = a12
	    self.a21 = a21
	    self.a22 = a22
	    self.dt = dt
	    self.tend = tend

    def gauss_elimin(self,a, b):
      (rows, cols) = a.shape
      # elimination phase
      for row in range(0, rows-1): # pivot equation/row
        for i in range(row+1, rows):
          if a[i, row] != 0.0:
            factor = a [i, row]/a[row, row]
            a[i, row+1:rows] = a[i, row+1:rows] - factor*a[row, row+1:rows]
            b[i] = b[i] - factor*b[row]
      # back substitution 
      for k in range(rows-1,-1,-1):
        b[k] = (b[k] - np.dot(a[k, k+1:rows],b[k+1:rows]))/a[k, k]
      return b
		
		
dt = 1.
tend = 100.
prob_1 = twotanks(0.02,-0.02,-0.02,0.02,dt,tend)

# initial conditions

y1_sol_0 = 150.
y2_sol_0 = 0.

amatrix = np.array([[1.+0.02*dt,-0.02*dt],[-0.02*dt,1.+0.02*dt]])
bvector = np.array([y1_sol_0,y2_sol_0]).T

t = 0.

time = []
y1_sol = []
y2_sol = []

# integration

while t < tend:   
#    solution = prob_1.gauss_elimin(amatrix,bvector)
    solution = np.linalg.solve(amatrix,bvector)
    bvector = np.array([solution[0],solution[1]]).T
    y1_sol.append(solution[0])
    y2_sol.append(solution[1])		
    t = t + dt
    time.append(t)
	
# calculate analytic solution
	
y1analytic = 75.+75.*np.exp(-0.04*np.asarray(time))
y2analytic = 75.-75.*np.exp(-0.04*np.asarray(time))
	
# plot results
	
plt.plot(time,y1_sol,color='r')
plt.plot(time,y2_sol,color='g')
plt.plot(time,y1analytic,color='k',linewidth=0.5,linestyle='--')
plt.plot(time,y2analytic,color='k',linewidth=0.5,linestyle='--')

plt.hlines(y=75.,xmin=np.min(time),xmax=np.max(time),color='k',linewidth=1)

plt.xlabel('t (s)')
plt.ylabel('y')
plt.legend(loc=1,prop={'size':14})

# display

#tst1 = np.array([[3,2],[-6,6]])
#tst2 = np.array([7,6])

tst1 = np.array([[3.,1.], [1.,2.]],dtype = np.float64)
tst2 = np.array([[9.,8.]],dtype = np.float64).T

print(sys.version)

print(tst1)
print(tst2)
print("---")

sola =prob_1.gauss_elimin(tst1,tst2)
solb = np.linalg.solve(tst1,tst2)

print sola
print solb

#print sola[0],sola[1]
#print solb[0],solb[1]

#print(3.*sola[0]+2.*sola[1],-6.*sola[0]+6.*sola[1])
#print(3.*solb[0]+2.*solb[1],-6.*solb[0]+6.*solb[1])



#plt.show()



