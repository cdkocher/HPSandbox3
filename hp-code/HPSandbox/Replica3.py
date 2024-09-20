#! /usr/bin/python
#
# Replica.py    

from __future__ import absolute_import
from builtins import str
from builtins import object
from Config3 import *
import Chain3 as Chain
import Monty3 as Monty
from Trajectory3 import *

import random
import string
import math
import sys
import os

class Replica(object):
    """A container object, to hold the Chain() and Monty() objects"""

    def __init__(self,config,repnum):
        """Initialize the Replica() object."""
	
        temp = config.REPLICATEMPS[repnum]
        self.repnum = repnum
        #self.repname = 'rep' + string.zfill(str(repnum),2)
        self.repname = 'rep' + str(repnum).zfill(2)
        self.repdir = config.DATADIR + self.repname
        self.chain = Chain.Chain(config)
        self.mc = Monty.Monty(config,temp,self.chain)    
    
