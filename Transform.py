from __future__ import division
import math
import random
import time

sqrt = math.sqrt
exp = math.exp
log = math.log
abs = math.fabs
fact = math.factorial

class sourceBox:

    def __init__(self, index, Nside, S, p):
        self.index = index # box index
        self.left = [i/Nside for i in index] # left endpoints in each dimension of box
        self.right = [(i+1)/Nside for i in index] # right endpoints
        self.center = [(i+0.5)/Nside for i in index] # center of box
        self.sources = self.getSources(S)
        self.weights = [1]*len(self.sources) # weight/charge for each source, automatically use 1
        self.Hermite = [[0 for col in range(p)] for row in range(p)] # Hermitian coefficients for box

    def getSources(self, S):
        sources = []
        ind = 0

        '''
            getSources iterates through the list of sources, finding sources that belong in the box
            and removing them from the list of sources so that after all boxes the source list
            will be empty
        '''
        while ind < len(S):
            flag = True

            for i, x in enumerate(S[ind]):
                if not (self.left[i] <= x < self.right[i]):
                    flag = False
                    break

            if flag:
                sources.append(S.pop(ind))
            else:
                ind += 1

        return sources

    def printSources(self):
        print 'Sources in Box', self.index
        print self.sources

class targetBox:

    def __init__(self, index, Nside, T, p):
        self.index = index
        self.left = [i/Nside for i in index]
        self.right = [(i+1)/Nside for i in index]
        self.center = [(i+0.5)/Nside for i in index]
        self.targets = self.getTargets(T)
        self.potentials = [0]*len(self.targets) # potential at  target
        self.Taylor = [[0 for col in range(p)] for row in range(p)] # Taylor coefficients for box

    # identical to getSources
    def getTargets(self, T):
        targets = []
        ind = 0

        while ind < len(T):
            flag = True

            for i, x in enumerate(T[ind]):
                if not (self.left[i] <= x < self.right[i]):
                    flag = False
                    break

            if flag:
                targets.append(T.pop(ind))
            else:
                ind += 1

        return targets

    def printTargets(self):
        print 'Targets in Box', self.index
        print self.targets

def norm(x, y):
    if len(x) == len(y):
        normsq = 0

        for i in range(len(x)):
            normsq += (x[i] - y[i])**2

        return sqrt(normsq)

    else:
        print 'vectors not same length'

def hermite(x, y, delta, p, d):
    # computes Hermite functions using recursion
    h = [[0 for col in range(p)] for row in range(d)]

    for k in range(d):
        h[k][0] = exp(-((x[k] - y[k])/sqrt(delta))**2)
        h[k][1] = 2*(x[k] - y[k])/sqrt(delta)*h[k][0]
        
        for l in range(2,p):
            h[k][l] = 2*(x[k] - y[k])/sqrt(delta)*h[k][l-1] - 2*(l-1)*h[k][l-2]

    return h


def FastGauss(N, M, d, delta, epsilon, r):
    """
        Implements the Fast Gauss Transform.

        Inputs:
            N - number of sources
            M - number of targets
            d - dimension (NOTE: hard-coded in some spots for 2 dimensions)
            delta - variance
            epsilon - tolerance (NOTE: results only yield about half the digits that tolerance specifies)

    """
    
    Q = N
#    r = 0.5
    Nside = math.ceil(1/(r*sqrt(2*delta)))
    r = 1/(Nside*sqrt(2*delta))
    n = int(math.ceil(sqrt(-log(epsilon)/(2*r**2))))
    K = 1.09**d
    weights = [1]*N


    p = 1
    while (K*(1/fact(p))**(d/2)*(r**(p+1)/(1-r))**d > epsilon):
        p += 1

#    p=10

    # cutoff values, use direct evaluation or plug in for boxes with sources (NF) or targets (ML) below number
    NF = p**2
    ML = p**2

    S = [[random.random() for col in range(d)] for row in range(N)]
    T = [[random.random() for col in range(d)] for row in range(M)]

    # calculating potentials directly to check errors
    S1 = S[:]
    T1 = T[:]
    checkPot =[]
    
    print 'Computing potentials directly for N =', N
    # start timing direct solve
    start = time.clock()
    for target in T1:
        potential = 0

        # for each source, add to potential at current target
        for index, source in enumerate(S1):
            potential += weights[index]*exp(-norm(target, source)**2/delta)

        checkPot.append(target+[potential])

    elapsedDirect = time.clock() - start

    checkPot.sort()

    print 'Finished direct. Computing potentials using Fast Gauss Transform.'
    # start timing Fast Gauss
    start = time.clock()

    BS = [sourceBox([i,j], Nside, S, p) for j in range(int(Nside)) for i in range(int(Nside))]
    BT = [targetBox([i,j], Nside, T, p) for j in range(int(Nside)) for i in range(int(Nside))]

    # counts keep track of how interactions are computed, numbers refer to methods in paper
    Total = 0
    by1 = 0
    by2 = 0
    by3 = 0
    by4 = 0

    for B in BS:
        NB = len(B.sources)

        # use Gaussian for each source in B
        if NB <= NF:
            for C in BT:
                # if C is in the interaction list of B. NOTE: we hard code this for 2 dimensions. Can change later
                if (abs(B.index[0] - C.index[0]) <= n) and (abs(B.index[1] - C.index[1]) <= n):
                    Total += 1
                    MC = len(C.targets)
                    
                    # directly compute each source's effect on potential of each target
                    if MC <= ML:
                        #print 'Computing B', B.index, 'on  C', C.index, 'by 1 (direct evaluation)'
                        by1 += 1

                        for targetInd, target in enumerate(C.targets):
                            for sourceInd, source in enumerate(B.sources):
                                C.potentials[targetInd] += B.weights[sourceInd]*exp(-norm(target,source)**2/delta)

                    # for each source add to the Taylor expansion for all targets in C
                    else:
                        #print 'Computing B', B.index, 'on  C', C.index, 'by 2 (add to Taylor series for each source)'
                        by2 += 1

                        # compute and add to Taylor series coefficients for terms indexed by i, j, stored in C.Taylor
                        for i in range(p):
                            for j in range(p):
                                for sourceInd, source in enumerate(B.sources):
                                    # for each source create Hermite functions using recursion
                                    h = hermite(C.center, source, delta, p, d)
                                    C.Taylor[i][j] += B.weights[sourceInd]*(-1)**(i+j)/(fact(i)*fact(j))*h[0][i]*h[1][j]

        # NB > NF
        else:
            # form the Hermite coefficients
            for i in range(p):
                for j in range(p):
                    sum = 0
                    for sourceInd, source in enumerate(B.sources):
                        sum += B.weights[sourceInd]*((source[0] - B.center[0])/sqrt(delta))**i*((source[1] - B.center[1])/sqrt(delta))**j

                    B.Hermite[i][j] = sum/(fact(i)*fact(j))

#####################################################################
#            print 'testing Hermite coefficients,', B.index
#            directsum = 0
#            for sourceInd, source in enumerate(B.sources):
#                 directsum += B.weights[sourceInd]*exp(-norm([1,0], source)**2/delta)
#            print 'directsum =', directsum
#
#            h = hermite([1,0], B.center, delta, p, d)
#            hermitesum = 0
#            for i in range(p):
#                for j in range(p):
#                    hermitesum += B.Hermite[i][j]*h[0][i]*h[1][j]
#
#            print 'hermitesum =', hermitesum
#
#####################################################################

            for C in BT:
                # if C is in interaction list of B
                if (abs(B.index[0] - C.index[0]) <= n) and (abs(B.index[1] - C.index[1]) <= n):
                    Total += 1
                    MC = len(C.targets)

                    # evaluate Hermite expansion at each target
                    if MC <= ML:
                        #print 'Computing B', B.index, 'on C', C.index, 'by 3 (Hermite for B, evaluate at target)'
                        by3 += 1

                        for targetInd, target in enumerate(C.targets):
                            # for each target create Hermite functions using recursion
                            h = hermite(target, B.center, delta, p, d)

                            # for each target add to the potential the effect from the box, indexed over terms in expansion
                            for i in range(p):
                                for j in range(p):
                                    C.potentials[targetInd] += B.Hermite[i][j]*h[0][i]*h[1][j]

                    # Hermite expansion in B and Taylor series about C
                    else:
                        #print 'Computing B', B.index, 'on C', C.index, 'by 4 (Hermite for B, Taylor about C)'
                        by4 += 1

                        # create the Hermite functions between B and C
                        h = hermite(C.center, B.center, delta, 2*p, d)
                        
                        # i, j index the Taylor series coefficient we add to
                        for i in range(p):
                            for j in range(p):
                                sum = 0

                                # k, l index the Hermite expansion coefficients we computed
                                for k in range(p):
                                    for l in range(p):
                                        sum += B.Hermite[k][l]*h[0][i+k]*h[1][j+l]

                                C.Taylor[i][j] += (-1)**(i+j)/(fact(i)*fact(j))*sum

    # Sum up Taylor coefficients
    for C in BT:
        if len(C.targets) > ML:
            for targetInd, target in enumerate(C.targets):
                sum = 0

                for i in range(p):
                    for j in range(p):
                        sum += C.Taylor[i][j]*((target[0] - C.center[0])/sqrt(delta))**i*((target[1] - C.center[1])/sqrt(delta))**j

                C.potentials[targetInd] += sum

    # stop timing
    elapsedFG = time.clock() - start

    # gather potentials to check error
    Pot = []
    for C in BT:
        for targetInd, target in enumerate(C.targets):
            Pot.append(target+[C.potentials[targetInd]])

    Pot.sort()

    for index, line in enumerate(checkPot):
        line.append(Pot[index][2])

    # compute relative error
    sum = 0
    for line in checkPot:
        sum += abs(line[2]-line[3])/line[2]
    relError = sum/len(checkPot)

    print 'r =', r
    print 'Nside =', Nside
    print 'n =', n
    print 'p =', p
    print 'Relative Error =', relError
    print Total, 'interactions total'
    print 'solved', by1, 'interactions by 1: direct evaluation'
    print 'solved', by2, 'interactions by 2: Taylor expansion'
    print 'solved', by3, 'interactions by 3: Hermite expansion and evaluate'
    print 'solved', by4, 'interactions by 4: Hermite then Taylor'

    print 'direct solve seconds:', elapsedDirect
    print 'FG solve seconds:', elapsedFG
    print '-----------------------------------------------------------------------------------------'

#    print 'Printing Fast Gauss targets and potentials'
#    for line in Pot:
#        print line
    
    return [elapsedDirect, elapsedFG, relError]

if __name__ == '__main__':

    runs = 7
    delta = 1
    epsilon = 10**-8
    table = []

    for i in range(runs):
        [d, fg, e] = FastGauss(100*2**i, 100*2**i, 2, delta, epsilon, 0.5)
        table.append([100*2**i]+[d, fg, e])

    print 'On 2 dimensional grid, with delta =', 1, 'and epsilon =', epsilon
    print '[N=M, Time(Direct), Time(FastGauss), Relative Error]'
    for line in table:
        print line



