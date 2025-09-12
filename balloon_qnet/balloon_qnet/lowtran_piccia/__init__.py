"""
Python wrapper of the venerable LOWTRAN7 atmospheric absorption and solar transmission
model circa 1992.

Note: specified Lowtran7 model limitations include
wlcminvstep >= 20 cm^-1
0 <= wlcminv <= 50000

historical note:
developed on CDC CYBER, currently runs on 32-bit single float--this can cause loss of numerical
precision, future would like to ensure full LOWTRAN7 code can run at 64-bit double float.

user manual:
https://apps.dtic.mil/sti/pdfs/ADA206773.pdf
"""

from .base import *
