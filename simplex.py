import numpy as np
np.warnings.filterwarnings('ignore')
#ignore divide by zero warnings as it returns np.inf which is fine


#compare two rows lexicographically
def lexLess(a, b):
    if (a == b).all():
        return 0
    idx = np.where( (a>b) != (a<b) )[0][0]
    if a[idx] < b[idx]: 
        return 1
    return 0

class Simplex():

    def __init__(self, A, b, c, epsilon=1e-5):
        self.A = A.astype(float)
        self.b = b.astype(float)
        self.c = c.astype(float)

        self.epsilon = epsilon

        self.soln = None
        self.minObj = np.inf

    #methods

    def solve(self):
        self.soln, self.minObj = self.solveLP(self.A, self.b, self.c)
        return self.getRes()

    def isFeasible(self):
        return self.hasFeasibleSoln(self.A, self.b, self.c)

    #helper functions:

    #main helper function for solving LP's
    def solveLP(self, A, b, c, onlySoln = True):

        n, m = A.shape
        solnVars = set()
        
        s_which_col = np.ones(n).astype(int) * -1
        cost = np.array(0.0)
        
        badSoln = 0
        
        #finds the initial basic variables
        #also calculate initial cost row and cost of tableau
        for i in range(m):
            if np.sum(np.abs(A[:,i])) == 1:
                for j in range(n):
                    if s_which_col[j] != -1:
                        continue
                    if abs(A[j][i] - 1) < self.epsilon:
                        solnVars.add(i)
                        s_which_col[j] = i
                        
                        cost -= b[j] * c[i]
                        c -= A[j] * c[i] 
                        
                        if b[j] < 0:
                            badSoln = 1
                 
        
        #if our starting soln has negative vars,
        #we need to find an actual feasible soln
        if badSoln:
            hasFeas = self.hasFeasibleSoln(A, b, c, cost)
            if hasFeas == 0:
                #no feasible soln
                self.soln = "no feasible solution", np.inf
                return "no feasible solution"
            else:
                #found feasible soln
                #start solving LP with new starting tableau
                return self.solveLP(A, b, c)
        
        #starting simplex method from feasible solution
        
        while(1):
            done = 1
            for i in c:
                if i > self.epsilon:
                    done = 0
            if done:
                #all elements of cost <= 0, so we found optimal solution
                break
                
            #pivot to element with largest value in cost row
            #picks first element if there are multiple
            piv = np.argmax(c)
                
            bound = self.doPivot(A, b, c, solnVars, cost, piv)
            if bound == 0:
                #LP is unbounded
                return "unbounded", -np.inf
                
        soln = np.zeros(m)
        for i in solnVars:
            corRow = np.argmax(A[:, i])
            soln[i] = b[corRow]

        self.soln = soln
        self.minObj = -cost
            
        if onlySoln:
            #returning -cost because objective function is equal to
            #-cost of tableau
            return soln, -cost
            
        #returning soln and full tableau
        return soln, solnVars, A, b, c, cost
   

    #find a starting feasible soln for LP
    #by solving a sub-LP
    #set given tableau to new feasible point
    def hasFeasibleSoln(self, A0, b0, c0, cost = 0):
        #starting sub-LP to find feasible solution to LP
        A = np.copy(A0)
        b = np.copy(b0)
        c = np.copy(c0)


        n, m = A.shape
        
        #create sub-LP tableau
        c = np.zeros(m + n)
        c[0:n] = -np.ones(n)
        
        for i in range(n):
            if b[i] < 0:
                A[i] *= -1
                b[i] *= -1
                
        newA = np.zeros([n, m+n])
        for i in range(n):
            curRow = np.zeros(n)
            curRow[i] = 1
            newA[i] = np.concatenate((curRow, A[i]), axis=0)
        
        #solve sub-LP using our solve function
        #is guaranteed not to infinitely recurse
        #as sub-LP will have good starting soln
        _, soln0, resultA, resultB, resultC, cost = self.solveLP(newA, b, c, onlySoln = False)
        
        if cost < self.epsilon:
            #original LP is feasible
            #start it from found feasible soln
            A0[:] = resultA[:, n:]
            b0[:] = resultB[:]
            return 1
        else:
            #original LP is infeasible
            return 0
        
    #pivot to a given variable
    #appropriately modify tableau and basic variables
    def doPivot(self, A, b, c, solnVars, cost, piv):
        n, m = A.shape

        #get indicator values for each row
        ind = b / A[:, piv]
        for i in range(n):
            if A[i, piv] < 0:
                ind[i] = np.inf
            if A[i, piv] == 0:
                ind[i] = np.inf
        
        if np.min(ind) == np.inf:
            #LP is unbounded
            return 0
        
        lowInd = np.argmin(ind)

        #get all rows with the lowest ind value
        lowInds = []
        for i in range(n):
            if abs(ind[i] - ind[lowInd]) < self.epsilon:
                lowInds.append(i)

        #use lexicographic method for breaking ties
        for i in lowInds:
            if lexLess(A[i], A[lowInd]):
                lowInd = i
        
        #change basic variables according to our pivoting
        curSoln = list(solnVars)
        for i in curSoln:
            if abs(A[lowInd][i] - 1) < self.epsilon:
                solnVars.remove(i)
                break
        solnVars.add(piv)
                    
        #change tableau
        div = A[lowInd][piv]
        A[lowInd, :] /= div
        b[lowInd] /= div
        cost -= b[lowInd] * c[piv]
        c -= A[lowInd] * c[piv]
        for i in range(n):
            if i == lowInd:
                continue
            mult = A[i][piv]
            A[i] -= A[lowInd] * A[i][piv]
            b[i] -= b[lowInd] * mult

        return 1
        
    def getRes(self):
        return self.soln, self.minObj
