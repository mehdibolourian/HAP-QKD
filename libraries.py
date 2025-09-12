#!/usr/bin/env python3.11

## General libraries
import time
import random
import threading
import multiprocessing
import concurrent.futures
import copy
import math
import pickle
import os
import cdsapi
import xarray                as     xr
import numpy                 as     np
import gurobipy              as     gp
import matplotlib            as     mpl
import matplotlib.pyplot     as     plt
import matplotlib.patches    as     mpatches
import networkx              as     nx
import xml.etree.ElementTree as     ET
import pandas                as     pd
import matplotlib.ticker     as     ticker
from   gurobipy              import GRB
from   collections           import Counter
from   matplotlib.ticker     import ScalarFormatter
from   tabulate              import tabulate
from   itertools             import product
from   scipy.interpolate     import interp1d, CubicSpline
from   pyproj                import Transformer
from   mpl_toolkits.mplot3d  import Axes3D
from   matplotlib.animation  import FuncAnimation

## Project library files
from data  import *
from setup import *
from plot  import *
from problems.offline_wind import *
from problems.planning     import *