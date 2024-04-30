__author__ = 'vetrot'

import sys
import os

table = sys.argv[1]
best_function = sys.argv[2]
function_name = sys.argv[3]
identity_var = sys.argv[4]
#k141_1000000_1_442_-	3.30E-22	[K00074]
linked_tab = sys.argv[5]

#read functions...
fun = {}
for n, line in enumerate(open(best_function)):
    vals = line.rstrip().split('\t')
    if len(vals) == 3:
        fun[vals[0]] = vals[1] +'\t' + vals[2]

print 'Functions processed... '+str(len(fun))

#link it...
nfun = 0
fp = open(linked_tab, 'w')
for n, line in enumerate(open(table)):
    line = line.rstrip()
    line_new = ''
    if n == 0:
        #header
        line_new = line +'\t' + identity_var + '\t' +function_name
    else:
        vals = line.rstrip().split('\t')
        fun_line = 'NaN'+ '\t' +'-'
        if fun.has_key(vals[0]):
            nfun = nfun + 1
            fun_line = fun[vals[0]]
        line_new = line + '\t'+fun_line
    fp.write(line_new + "\n")
fp.close()

print 'Done... used functions: '+str(nfun)+'/'+str(len(fun))