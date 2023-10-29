import matplotlib as plt
from matplotlib.animation import FuncAnimation
import random
import math

folded = [[3,3],[2,3],[1,3],[0,3],[0,2],[1,2],[2,2],[3,2],[3,1],[3,0],[2,0],[2,1],[1,1],[1,0],[0,0],[0,1]] # Using x,y coordinates on a lattice, list is ordered

def dist(res1, res2):
    return ((res1[0]-res2[0])**2 + (res1[1]-res2[1])**2)**(0.5) # Euclidean distance between two points


def metro(e0, e1): # Defining metropolis criterion, returns boolean
    deltaE = e0 - e1
    w = math.exp(deltaE)
    r = random.uniform(0,1)
    if (w>1):
        return True
    elif (w>r):
        return True
    else:
        return False


# energy = num of non-covalent interactions * interaction energy

def avgX(grid):
    x = 0
    for coord in grid:
        x+=coord[0]
    return x/len(grid)

def avgY(grid):
    y = 0
    for coord in grid:
        y+=coord[1]
    return y/len(grid)

def avgCoord(grid):
    return [avgX(grid), avgY(grid)]

def edgeMove(res, grid, idx):
    
    resAfter = grid[idx+1]
    resBefore = grid[idx-1]
    # residues before and after must be at distance of 1 unit
    # can only move diagonally

    moves = [[-1,1],[1,1],[1,-1],[-1,-1]]

    validMoves = []
    
    for i in moves:
        tryMove = i
        newCoord = [tryMove[0]+res[0],tryMove[1]+res[1]] # where the residue would go after making the move
        if (newCoord not in grid): # Checking if lattice point is already occupied
            if ((dist(newCoord,resBefore)==1) and (dist(newCoord,resAfter)==1)): # Checking if chain remains intact
                validMoves.append(tryMove)
    
    if (len(validMoves)!=0): # If the residue has any possible allowed move
        return random.choice(validMoves)

    # elif (idx!=1):
    #    return edgeMove(grid[idx-1], grid, idx-1)
        
    else:
        # return cornerMove(grid[0], grid, 0)
        return [0,0]
    

def cornerMove(res, grid, idx):

    
    if (idx==0):
        checkDistRes = grid[idx+1] # First residue, check next residue for structure
    else:
        checkDistRes = grid[idx-1] # Last residue, check previous residue for structure

    moves = [[-1,1],[1,1],[1,-1],[-1,-1]]

    validMoves = []
    
    for i in moves: # Same as edge move conditions, except we only check one nearest neighbour of the residue
        tryMove = i
        newCoord = [tryMove[0]+res[0],tryMove[1]+res[1]]
        if (newCoord not in grid):
            if ((dist(newCoord,checkDistRes)==1)):
                validMoves.append(tryMove)
    
    if (len(validMoves)!=0):
        return random.choice(validMoves)

    else:
        return [0,0]
    
nativeinteractions = [[0,7],[1,6],[2,5],[5,12],[8,11],[10,13],[12,15]]

def numNative(grid):
    s=0
    for pair in nativeinteractions:
        if (dist(grid[pair[0]], grid[pair[1]])==1):
            s+=1
    return s

def ratioNative(grid):
    s=numNative(grid)
    return s/9

numSteps = 100000

grid = folded.copy() # Initializing to completely folded state

interaction_e = -0.25
# nc = totalNC(folded)

avgE = 0 # Average number of steps needed for reaching unfolded state

microstates = []

steps = []

for i in range(50):
    grid = folded.copy()
    for step in range(numSteps):
        res = random.choice(grid) # Pick a random residue
        idx = grid.index(res)

        
        testGrid = grid.copy()

        

        e0 = numNative(grid)*interaction_e # Calculate energy before the move is made (Ei)

        

        if ((idx==0)or(idx==(len(testGrid)-1))): # Seeing whether the randomly picked residue is an end bead or not, moving accordingly
            move = cornerMove(res, testGrid, idx)
        else:
            move = edgeMove(res, testGrid, idx)

        newCoord = [move[0]+res[0], move[1]+res[1]] # Trying out the move

        testGrid[idx] = newCoord # Make move on a replica of our lattice (not actually made yet, just to check energy values)

        e1 = numNative(testGrid)*interaction_e # Energy after making the move

        approve = metro(e0, e1) # Check metropolis criterion

        if (approve==True): # If approved by metropolis criterion
            grid = testGrid.copy()

        currentE = numNative(grid)*interaction_e # Update the energy value

        ratio = ratioNative(grid)

        rc = (dist(grid[0], grid[-1]))*(dist(grid[0], avgCoord(grid)))*(dist(grid[-1], avgCoord(grid)))

        if (rc not in microstates):
            microstates.append(rc)

        print(f"Exp {i}: Iteration {step+1}/{numSteps}, Energy = {currentE}/{9*interaction_e}, Reaction Coordinate = {rc}")

        if ((currentE == 0)or(step==99999)):
            steps.append(step+1)
            break

totalSteps = 0
for i in steps: totalSteps += i

avgSteps = totalSteps/len(steps)

steps.sort()
mid = len(steps) // 2
median = (steps[mid] + steps[~mid]) / 2

mode = max(set(steps), key=steps.count)


print(f"FOR ENERGY = {interaction_e}, STEPS DISTRIBUTION: MEAN = {avgSteps}, MEDIAN = {median}, MODE = {mode} ")
print(steps)

