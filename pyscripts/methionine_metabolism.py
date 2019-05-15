# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 14:58:48 2018
Human genome scale model of methionine metabolism 
@author: scampit
"""

import cobrapy
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import seaborn as sns

from cobra import Model, Reaction, Metabolite

model = Model('Methionine Metabolism')

Reaction = Reaction("")
