#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys

"""Run the Model Predictive Control problem on encrypted data in a Two-Server architecture. 
	The two servers are emulated by different threads.
	See https://arxiv.org/pdf/1803.09891.pdf for more details"""
 
if __name__ == '__main__':
	proc1 = Popen(['python3','server1.py'])
	proc2 = Popen(['python3','server2.py'],stdout=PIPE)
	# stdout_value = proc2.communicate()[0]
	stdout_value = proc1.communicate()[0]

