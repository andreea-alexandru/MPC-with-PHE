#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys

"""Run the Model Predictive Control problem on encrypted data in a Client-Server architecture. 
	The Client and Server are emulated by different threads.
	See https://arxiv.org/pdf/1803.09891.pdf for more details"""
 
if __name__ == '__main__':
	proc1 = Popen(['python3','server.py'],stdout=PIPE)
	proc2 = Popen(['python3','client.py'])
	stdout_value = proc2.communicate()[0]
	# stdout_value = proc1.communicate()[0]


