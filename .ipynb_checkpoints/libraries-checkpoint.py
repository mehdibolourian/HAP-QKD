#!/usr/bin/env python3.12

## General libraries
import time
import random
import threading
import concurrent.futures
import copy
import math
import pickle
import os
import sys
print(sys.executable)
import cdsapi
import pprint
import torch
import multiprocessing       as     mp
import xarray                as     xr
import numpy                 as     np
import gurobipy              as     gp
import matplotlib            as     mpl
import matplotlib.pyplot     as     plt
import matplotlib.patches    as     mpatches
import plotly.express        as     px
import networkx              as     nx
import xml.etree.ElementTree as     ET
import pandas                as     pd
import matplotlib.ticker     as     ticker
import torch.nn              as     nn
import torch.optim           as     optim
import torch.nn.functional   as     F
from   torch.distributions   import Categorical
import gymnasium             as     gym
from   gurobipy              import GRB
from   collections           import Counter
from   collections           import deque
from   matplotlib.ticker     import ScalarFormatter
from   tabulate              import tabulate
from   itertools             import product, combinations
from   scipy.interpolate     import interp1d, CubicSpline
from   pyproj                import Transformer
from   mpl_toolkits.mplot3d  import Axes3D
from   matplotlib.animation  import FuncAnimation
from   matplotlib.lines      import Line2D
from   concurrent.futures    import ThreadPoolExecutor,ProcessPoolExecutor, as_completed
from   tqdm                  import tqdm
from   sklearn.cluster       import KMeans
from   matplotlib.lines      import Line2D

from   mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.append("./balloon_qnet")
import transmittance_simulation as ts

## Project library files
from data   import *
from plot   import *
from helper import *
from setup  import *
from problems.offline  import *
from problems.online   import *
from problems.planning import *