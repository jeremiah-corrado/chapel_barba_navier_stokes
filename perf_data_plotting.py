import matplotlib.pyplot as plt
import csv
import sys
import re
from collections import OrderedDict

chpl_perf = {}
py_perf = {}

with open(sys.argv[1]) as csv_data:
    r = csv.reader(csv_data, delimiter=',')
    header = True
    for row in r:
        if header:
            header = False
        else:
            problem_size = int(re.search(r"nx=([0-9]+)", row[1]).group(1)) * int(re.search(r"ny=([0-9]+)", row[1]).group(1))
            if row[0] == 'chpl':
                chpl_perf[problem_size] = float(row[-1].strip())
            elif row[0] == 'py':
                py_perf[problem_size] = float(row[-1].strip())

chpl_perf = OrderedDict(sorted(chpl_perf.items()))
py_perf = OrderedDict(sorted(py_perf.items()))

f1 = plt.figure(1)
ax1 = f1.add_subplot(1, 1, 1)
ax1.plot(chpl_perf.keys(), chpl_perf.values(), label='chpl')
ax1.plot(py_perf.keys(), py_perf.values(), label='py')
ax1.legend()
ax1.set_xlabel("Problem Size")
ax1.set_ylabel("Walltime (sec)")
ax1.set_title("Performance Comparison (Chapel vs. Python) - {deets}".format(deets = sys.argv[2]))
f1.savefig("perf_plots/" + sys.argv[3] + ".png")
f1.show()

f2 = plt.figure(2)
ax2 = f2.add_subplot(1, 1, 1)
ax2.semilogy(chpl_perf.keys(), chpl_perf.values(), label='chpl')
ax2.semilogy(py_perf.keys(), py_perf.values(), label='py')
ax2.legend()
ax2.set_xlabel("Problem Size")
ax2.set_ylabel("Walltime (sec) - Log Scale")
ax2.set_title("Performance Comparison (Chapel vs. Python) - {deets}".format(deets = sys.argv[2]))
f2.savefig("perf_plots/" + sys.argv[3] + "_log.png")
f2.show()
