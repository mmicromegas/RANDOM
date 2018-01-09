import numpy as np
import matplotlib.pyplot as plt

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
		
    def LHS(self,y1_lhs_sol,y2_lhs_sol):
	    return np.array(y1_lhs_sol,y2_lhs_sol)
		
    def RHS(self,y1_rhs,y2_rhs):
        I = np.array([[1.,0.],[0.,1.]])
        A = np.array([[self.a11,self.a12],[self.a21,self.a22]])
        dt = self.dt
        Y = np.array([y1_rhs,y2_rhs]).T
        return np.dot((I - dt*A),Y)
		
dt = 0.1
tend = 100.
prob_1 = twotanks(0.02,-0.02,-0.02,0.02,dt,tend)

# initial conditions

y1_lhs_sol_0 = 150.
y2_lhs_sol_0 = 0.

t = 0.

time = []
y1_sol = []
y2_sol = []

while t < tend:   
    solution = prob_1.RHS(y1_lhs_sol_0,y2_lhs_sol_0)
    y1_lhs_sol_0 = solution[0]
    y2_lhs_sol_0 = solution[1]
    y1_sol.append(solution[0])
    y2_sol.append(solution[1])		
    t = t + dt
    time.append(t)
	
y1analytic = 75.+75.*np.exp(-0.04*np.asarray(time))
y2analytic = 75.-75.*np.exp(-0.04*np.asarray(time))
	
plt.plot(time,y1_sol,color='r')
plt.plot(time,y2_sol,color='g')
plt.plot(time,y1analytic,color='k',linewidth=0.5,linestyle='--')
plt.plot(time,y2analytic,color='k',linewidth=0.5,linestyle='--')

plt.hlines(y=75.,xmin=np.min(time),xmax=np.max(time),color='k',linewidth=1)

plt.xlabel('t (s)')
plt.ylabel('y')
plt.legend(loc=1,prop={'size':14})




plt.show()



