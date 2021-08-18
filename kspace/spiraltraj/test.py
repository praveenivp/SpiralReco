#!/usr/bin/env python

import spiraltraj
import matplotlib.pyplot as plt

print('basic test: spiral out trajectory calculation with default settings:')
res = spiraltraj.calc_traj(spiraltype=1)
print('shape of returned list: [', len(res), ', ', len(res[0]), ']')

plt.plot(res)