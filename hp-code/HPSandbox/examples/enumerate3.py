#! /usr/bin/python

usage = """
enumerate3.py <configfile>

Try:  enumerate3.py enumerate.conf

This program will read in an HP chain specified in the configure file,
and perform a full enumeration of conformational space.

The problem tablulates:

    1) the density of states (in energies/contacts)
    
    2) the number density of unique contact states, i.e. disjoint collections
       of microscopic conformations all sharing a unique set of interresidue contacts. 

These values are printed as output.

"""


import sys
sys.path.append("/home/charles/repositories/HPSandbox3/hp-code/HPSandbox/")

from Config3 import *
from Chain3 import *
from Monty3 import *
from Replica3 import *
from Trajectory3 import *
import copy

import random
import string
import math
import os

g = random.Random(randseed)


if len(sys.argv) < 2 or '-h' in sys.argv:
    print('Usage:  enumerate3.py <configfile> [options]')
    print('Make sure the config file is the first argument.')
    print('Possible options: ')
    print('-mpl: Matplotlib, plot native structure (or example with max number of contacts) in matplotlib')
    print('-nc: No contacts, do not show contacts explicitly in plots')
    print('-full: Full mode, print all unique conformations found (I would not try this for huge sequences, say, much bigger than a 7-mer)')
    print('-h: Help mode, print this help message')
    print('-q: Quiet mode, do not print the config file at the beginning')
    print('-s: Silent mode, only print out the final density table and the plots')
    print('-nn: No native, do not print an example with the max number of contacts')
    print('-pe: Print examples, print an example with each number of contacts')
    print('-t: Thick, used with -mpl flag to plot thick covalent bonds')
    print('-nb: No boundary, used with -mpl flag to remove figure boundary box')
    sys.exit(1)

# default options
VERBOSE = 1
mplflag = 0
fullflag = 0
contactflag = 1
silentflag = 0
nativeflag = 1
printeveryflag = 0
thickflag = 0
nbflag = 0

if '-q' in sys.argv:
    VERBOSE = 0

if '-s' in sys.argv:
    VERBOSE = 0
    silentflag = 1

if '-mpl' in sys.argv:
    mplflag = 1

if '-full' in sys.argv:
    fullflag = 1

if '-nc' in sys.argv:
    contactflag = 0

if '-nn' in sys.argv:
    nativeflag = 0

if '-pe' in sys.argv:
    nativeflag = 0
    printeveryflag = 1

if '-t' in sys.argv:
    thickflag = 1

if '-nb' in sys.argv:
    nbflag = 1

# define the function that takes in a vector list and spits out a structure string

def vec2stxstring(maxvec,sequence,xycoordsoption=False):
    """Turns a vector for a sequence (telling where the bonds go, 0 = north, 1 = west, 2 = south, 3 = east) into a structure that we can print in the terminal. xycoordsoption denotes whether we want to return a dictionary with xy coords of each monomer that we could easily plot in matplotlib as balls and sticks."""
    # to do this, we will make a matrix of strings, then write out the sequence inside that matrix in a way that connects the whole chain
    # first, find the total number of spaces that we need
    # the total number of rows we need is 2*(#0s) + 2*(#2s) + 1 (each 0 adds a ball and a stick north, each 2 adds a ball and a stick south, plus the origin position)
    # the total number of columns we need is 2*(#1s) + 2*(#3s) + 1 (again, ball and stick)
    numrows = 2 * sum([1 for kk in maxvec if kk == 0 or kk == 2]) + 1
    numcolumns = 2 * sum([1 for kk in maxvec if kk == 1 or kk == 3]) + 1
    structuremat = [[' ' for jj in range(numcolumns)] for kk in range(numrows)]
    # need to find the origin's index now
    # row index is 2*(#0s) (this automatically accounts for the +1, e.g. if we have one 0 then the origin is in row 2)
    originrow = 2 * sum([1 for kk in maxvec if kk == 0])
    # column index is 2*(#1s)
    origincolumn = 2 * sum([1 for kk in maxvec if kk == 1])
    # put first monomer at origin
    structuremat[originrow][origincolumn] = sequence[0]
    # set current position at origin
    currentposition = [originrow,origincolumn]
    # get storage for xy coords
    xycoords = {0:(0,0)} # string index:(x,y)
    sequenceindx = 1 # for tracking where we are in the string
    # now loop through the moves and place their corresponding monomers
    for (mv,toplace) in zip(maxvec,sequence[1:]):
        currentrow = currentposition[0]
        currentcolumn = currentposition[1]
        if mv == 0:
            # moving to north, which are the two previous rows
            structuremat[currentrow-1][currentcolumn] = '|'
            structuremat[currentrow-2][currentcolumn] = toplace
            currentposition = [currentrow-2,currentcolumn]

        if mv == 2:
            # moving south, which are the two next rows
            structuremat[currentrow+1][currentcolumn] = '|'
            structuremat[currentrow+2][currentcolumn] = toplace
            currentposition = [currentrow+2,currentcolumn]

        if mv == 1:
            # moving to west, which are the two previous columns
            structuremat[currentrow][currentcolumn-1] = '-'
            structuremat[currentrow][currentcolumn-2] = toplace
            currentposition = [currentrow,currentcolumn-2]

        if mv == 3:
            # moving east, which are the two next rows
            structuremat[currentrow][currentcolumn+1] = '-'
            structuremat[currentrow][currentcolumn+2] = toplace
            currentposition = [currentrow,currentcolumn+2]

        # now we want to put in contacts.
        newrow = currentposition[0]
        newcolumn = currentposition[1]
        if toplace == 'H':
            # only do it if we have an H.
            # Now check each position +/- 2 (when structuremat actually exists in that direction)
            # only add the contact c if that direction wasn't the move
            if mv != 2 and newrow - 1 > 0 and structuremat[newrow-2][newcolumn] == 'H':
                # first check to the north and make sure the move wasn't south (2)
                structuremat[newrow-1][newcolumn] = 'c'

            if mv != 0 and newrow + 1 < numrows and structuremat[newrow+2][newcolumn] == 'H':
                # now check to the south and make sure the move wasn't north (0)
                structuremat[newrow+1][newcolumn] = 'c'

            if mv != 3 and newcolumn - 1 > 0 and structuremat[newrow][newcolumn-2] == 'H':
                # now check to the west and make sure the move wasn't east (3)
                structuremat[newrow][newcolumn-1] = 'c'

            if mv != 1 and newcolumn + 1 < numcolumns and structuremat[newrow][newcolumn+2] == 'H':
                # now check to the east and make sure the move wasn't west (1)
                structuremat[newrow][newcolumn+1] = 'c'

        # now store xycoords
        xycoords[sequenceindx] = (int(-1*(newcolumn - origincolumn)/2),int(-1*(newrow - originrow)/2))
        sequenceindx += 1

    # join the columns
    structurelist = [''.join(kk) for kk in structuremat]
    # cut empty rows
    emptyrow = ''.join([' ' for jj in range(numcolumns)])
    structurelist = [kk for kk in structurelist if kk != emptyrow]
    # join the rows with newlines to keep the structure
    structurestring = '\n'.join(structurelist)
    # optionally return xycoords
    if xycoordsoption:
        return (structurestring,xycoords)
    return structurestring
    

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(1)
	
    configfile = sys.argv[1]
    config = Config( filename=configfile)
    if VERBOSE:
        config.print_config()
        print()
        print("Executed command: {} ".format(" ".join(sys.argv)))
        print()
    
    # create a single Replica
    replicas = [ Replica(config,0) ]
    
    traj = Trajectory(replicas,config)	# a trajectory object to write out trajectories

    nconfs = 0
    contact_states = {}		# dictionary of {repr{contact state}: number of conformations}
    contacts = {}               # dictionary of {number of contacts: number of conformations}


    #################
    #
    # This is a useful subroutine for enumerating all conformations of an HP chain
    #
    # NOTE: in order for this to work correctly, the initial starting vector must be [0,0,0,....,0]
    # 
    done = 0
    maxcontacts = 0        # Track the max number of contacts and print the conf if we found a new one
    maxvectors = []    # keep a list with examples of each amount of contacts
    maxstatecontacts = [] # keep a list with the contact states as well
    fullstates = {} # keep dictionary of ncontacts:[all state vectors with this many contacts] which we will display if we have the fullflag. Also, only store them if we have the fullflag.
    while not(done):
	    	    
        if len(replicas[0].chain.vec) == replicas[0].chain.n-1:    
            if replicas[0].chain.viable:		
                if replicas[0].chain.nonsym():
		    
		    # tally the number of contacts
                    state = replicas[0].chain.contactstate()
                    ncontacts = len(state)
                    if (ncontacts in contacts) == False:
                        contacts[ncontacts] = 1
                    else:
                        contacts[ncontacts] = contacts[ncontacts] + 1

		    # tally the contact state
                    this_state_repr = repr(state)
                    if (this_state_repr in contact_states) == False:
                        contact_states[this_state_repr] = 1
                    else:
                        contact_states[this_state_repr] = contact_states[this_state_repr] + 1

		    # tally the number of conformations
                    nconfs = nconfs + 1

                    # store the chain.vec if we have the fullflag
                    if fullflag:
                        if ncontacts not in fullstates.keys():
                            # not there, so we need to make the list
                            fullstates[ncontacts] = list()
                            
                        fullstates[ncontacts].append(copy.deepcopy(replicas[0].chain.vec))

		    # write to trajectory
                    if (nconfs % config.TRJEVERY) == 0:
                        traj.queue_trj(replicas[0])
		    # print progress, do it if not silent
                    if (nconfs % config.PRINTEVERY) == 0 and not silentflag:
                        print('%-4d confs  %s'%(nconfs,replicas[0].chain.vec))
                    # print if new max contacts and not silent flag
                    if ncontacts > maxcontacts:
                        maxcontacts = ncontacts
                        if not silentflag:
                            print('New max, %d contacts, config  %s'%(ncontacts,replicas[0].chain.vec))
                            
                        maxvectors.append((maxcontacts,copy.deepcopy(replicas[0].chain.vec)))
                        maxstatecontacts.append(copy.deepcopy(state))
    
                done = replicas[0].chain.shift()
		    
            else:
                done = replicas[0].chain.shift()

        else:
            if replicas[0].chain.viable:
                replicas[0].chain.grow()
            else:
                done = replicas[0].chain.shift()

        if replicas[0].chain.vec[0] == 1:    # skip the other symmetries
            break	
    #
    #
    #################
        
    
    # write the last of the trj and ene buffers
    # and close all the open trajectory file handles
    traj.cleanup(replicas)
    
    # print out the density of contact states
    print()
    print('DENSITY of CONTACT STATES:')
    print('%-40s %s'%('contact state','number of conformations'))
    for state in list(contact_states.keys()):
        print('%-40s %d'%(state, contact_states[state]))
    
    # print out the density of states (energies)
    print() 
    print('DENSITY of STATES (in energies/contacts):')
    print('%-20s %-20s %s'%('number of contacts','energy (kT)','number of conformations'))
    for c in list(contacts.keys()):
        print('%-20d %-20d %d'%(c,config.eps*c,contacts[c]))

    print()
    print('at T = %4.1f K'%config.T)
    print()
    # now we want to print the native structure (if it uniquely exists; if not, this just prints an example of the lowest energy conformations)
    maxvec = maxvectors[-1][1] # get the lowest energy conformation we stored
    sequence = config.HPSTRING # get the sequence
    nativecontacts = maxstatecontacts[-1] # get the contacts
    (structurestring,xycoords) = vec2stxstring(maxvec,sequence,xycoordsoption=True)
    if nativeflag: # print the example native state if available
        print('Example conformation vector with %d contacts: %s'%(maxcontacts,maxvectors[-1][1]))
        print("NOTE on conformation vectors: 0 = north, 1 = west, 2 = south, 3 = east.")
        print()
        print("Here is a structure with this conformation vector.")
        # display with contacts if desired
        if contactflag:
            print("c denotes a contact.")
            print()
            print("Coordinates for sequence " + sequence + ":")
            print([xycoords[kk] for kk in range(len(xycoords.keys()))])
            print()
            print(structurestring)
            print()
        else: # no contacts
            print()
            print("Coordinates for sequence " + sequence + ":")
            print([xycoords[kk] for kk in range(len(xycoords.keys()))])
            print()
            print(structurestring.replace('c',' '))
            print()

    # now, if we have the printeveryflag, display each maxvec in maxvectors
    if printeveryflag:
        print('Using -pe option, now printing examples at each number of contacts')
        print("NOTE on conformation vectors: 0 = north, 1 = west, 2 = south, 3 = east.")
        if contactflag:
            print("c denotes a contact.")
        print()
        for (contax,mvec) in maxvectors:
            (structurestring,xycoords) = vec2stxstring(mvec,sequence,xycoordsoption=True)
            print('Example conformation vector with %d contacts: %s'%(contax,mvec))
            # display with contacts if desired
            if contactflag:
                print("Coordinates for sequence " + sequence + ":")
                print([xycoords[kk] for kk in range(len(xycoords.keys()))])
                print()
                print(structurestring)
                print()
            else: # no contacts
                print("Coordinates for sequence " + sequence + ":")
                print([xycoords[kk] for kk in range(len(xycoords.keys()))])
                print()
                print(structurestring.replace('c',' '))
                print()
            
    # now, if we have the fullflag, display all unique conformations at each energy level
    if fullflag:
        # say what we are printing
        print('Using -full option, now printing all unique conformations')
        print()
        # loop over all numbers of contacts
        for ke in fullstates.keys():
            # say that we are now printing a certain energy level
            print('Now printing states with number of contacts = {}'.format(ke))
            print()
            counter = 0
            # loop over all vectors with that many contacts
            for chainvec in fullstates[ke]:
                # get the string that prints to the structure
                (structurestring,xycoords) = vec2stxstring(chainvec,sequence,xycoordsoption=True)
                # print the structure with or without contacts
                counter += 1 # just for labeling the printouts for reference
                if contactflag:
                    print("({}) Coordinates:".format(counter))
                    print([xycoords[kk] for kk in range(len(xycoords.keys()))])
                    print()
                    print(structurestring)
                    print()
                else: # no contacts
                    print("Coordinates:")
                    print([xycoords[kk] for kk in range(len(xycoords.keys()))])
                    print()
                    print(structurestring.replace('c',' '))
                    print()

    # now, if we have the mpl flag, do the plotting
    # we also leave the import statements for matplotlib down here so that it is not a dependency for the full code, only this section
    if mplflag:
        (structurestring,xycoords) = vec2stxstring(maxvec,sequence,xycoordsoption=True)
        import matplotlib.pyplot as plt
        # loop over the sequence and plot the appropriately colored circle at the correct location
        # make the figure square and size it so the circles are radius 0.05 in
        xcoords = [xy[0] for xy in xycoords.values()]
        ycoords = [xy[1] for xy in xycoords.values()]
        boxmax = max([max(xcoords),max(ycoords)])+1
        boxmin = min([min(xcoords),min(ycoords)])-1
        # boxsize needs to be such that 2/3 = 1 in
        boxlength = boxmax - boxmin
        boxsize = boxlength * 1 / (2/3)
        plt.figure(1,(boxsize,boxsize))
        plt.xlim([boxmin,boxmax])
        plt.ylim([boxmin,boxmax])
        ax = plt.gca()
        for kk in range(len(sequence)):
            # blue is #2c38ff
            # red is #f60000
            # starting red pink is #f1cccc
            if sequence[kk] == 'H':
                # red circle
                if kk == 0:
                    # plot pink if first
                    circ = plt.Circle(xycoords[kk],0.33,facecolor='#f1cccc',edgecolor='black',linewidth=4,zorder=10)
                else:
                    circ = plt.Circle(xycoords[kk],0.33,facecolor='#f60000',edgecolor='black',linewidth=4,zorder=10)
            if sequence[kk] == 'P':
                # blue circle
                circ = plt.Circle(xycoords[kk],0.33,facecolor='#2c38ff',edgecolor='black',linewidth=4,zorder=10)

            ax.add_patch(circ)

            # now draw a line to the next point if it exists, and make sure it is beneath the circles by giving circles zorder 10
            if kk < len(sequence) - 1:
                # black line, use correct thickness
                linewidthchoice = 6
                if thickflag:
                    linewidthchoice = 24
                    
                plt.plot([xycoords[kk][0],xycoords[kk+1][0]],[xycoords[kk][1],xycoords[kk+1][1]],'k-',linewidth=linewidthchoice)

        # now display the contacts as dashed red lines, if desired
        if contactflag:
            for pair in nativecontacts:
                plt.plot([xycoords[pair[0]][0],xycoords[pair[1]][0]],[xycoords[pair[0]][1],xycoords[pair[1]][1]],'r--',linewidth=4)

        # remove boundary box, if desired
        if nbflag:
            ax.axis("off")
            
        plt.show()
