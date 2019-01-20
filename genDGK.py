#!/usr/bin/env python3
import random
from phe.util import invert, powmod, getprimeover
import numpy
import time
from gmpy2 import mpz

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False

DEFAULT_KEYSIZE = 512
DEFAULT_MSGSIZE = 16
DEFAULT_SECURITYSIZE = 160

def isprime(x):
	if (n <= 1): return 0
	if (n <= 3):  return 1
	if (n%2 == 0 or n%3 == 0): return 0;

	i = 5
	while i*i<=n:
		if (n%i == 0 or n%(i+2) == 0): return 0
		i = i+6 
	return 1

def keysDGK(n_length=DEFAULT_KEYSIZE,l=DEFAULT_MSGSIZE,t=DEFAULT_SECURITYSIZE):
	u = 2**l
	len_vp = t+1
	while len_vp > t:
		vp = getprimeover(t)
		len_vp = vp.bit_length()
	len_vq = t+1
	while len_vq > t:
		vq = getprimeover(t)
		len_vq = vq.bit_length()
	n_len = 0
	prime = 0
	# print(u,vp,vq)
	while n_len != n_length:
		fp = getprimeover((n_length//2) - l - t)
		p = u*vp*fp + 1
		# print(p.bit_length())
		while (not(gmpy2.is_prime(p))):
			fp = getprimeover((n_length//2) - l - t)
			p = u*vp*fp + 1
		fq = getprimeover((n_length//2) - l - t + 1)
		q = u*vq*fq + 1
		# print(q.bit_length())
		while (not(gmpy2.is_prime(q))):
			fq = getprimeover((n_length//2) - l - t + 1)
			q = u*vq*fq + 1		
		n = p*q
		n_len = n.bit_length()
		# print(n_len)
	
	n = p*q
	g,h = find_gens(p,q,u,vp,vq,fp,fq,n)
	return p,q,u,vp,vq,fp,fq,g,h

def loadkey(file):
	with open(file, 'r') as fin:
		data=[line.split() for line in fin]
	# data = numpy.loadtxt(file, delimiter=' ')
	p = data[0][0]
	q = data[1][0]
	u = data[2][0]
	vp = data[3][0]
	vq = data[4][0]
	fp = data[5][0]
	fq = data[6][0]	
	g = data[7][0]
	h = data[8][0]		
	p = mpz(p)
	q = mpz(q)
	u = mpz(u)
	vp = mpz(vp)
	vq = mpz(vq)
	fp = mpz(fp)
	fq = mpz(fq)	
	g = mpz(g)
	h = mpz(h)					
	return p,q,u,vp,vq,fp,fq,g,h


def find_gens(p,q,u,vp,vq,fp,fq,n,n_length=DEFAULT_KEYSIZE,l=DEFAULT_MSGSIZE,t=DEFAULT_SECURITYSIZE):
	found = False
	# compute g from CRT with g mod p = xp, g mod q = xq, g = g^(fp*fq) mod n
	# xp is a prime element from Z^*_p, that does not have order a divisor of fp*u*vp
	exp1 = mpz(fp*vp*u/2)
	exp2 = mpz(vp*u)
	exp3 = mpz(fp*vp)
	p_len = p.bit_length()
	while not(found):
		xp = getprimeover(p_len)
		if (xp < p):
			found = True
		if (found):
			order = gmpy2.powmod(xp,exp1,p)
			if (order==1):
				found = False
			else:
				order = gmpy2.powmod(xp,exp2,p)
				if (order==1):
					found = False
				else:
					order = gmpy2.powmod(xp,exp3,p)
					if (order==1):
						found = False
	# xq is a prime element from Z^*_q, that does not have order a divisor of fq*u*vq
	found = False
	exp1 = mpz(fq*vq*u/2)
	exp2 = mpz(vq*u)
	exp3 = mpz(fq*vq)
	q_len = q.bit_length()
	while not(found):
		xq = getprimeover(q_len)
		if (xq < q):
			found = True
		if (found):
			order = gmpy2.powmod(xq,exp1,q)
			if (order==1):
				found = False
			else:
				order = gmpy2.powmod(xq,exp2,q)
				if (order==1):
					found = False
				else:
					order = gmpy2.powmod(xq,exp3,q)
					if (order==1):
						found = False
	# CRT: g = xp*q*(q^{-1}mod p) + xq*p*(p^{-1}mod q) mod n, g = g^(fp*fq) mod n
	inv_q = gmpy2.invert(q, p)
	inv_p = gmpy2.invert(p, q)
	tmp = xp*q*inv_q % n
	tmp2 = xq*p*inv_p % n
	g = (tmp + tmp2) % n
	g = gmpy2.powmod(g,fp*fq,n)
	# print(g)

	# compute h as xh^(fp*fq*u) mod n
	# xh is a prime element from Z^*_n
	found = False
	while not(found):
		xh = getprimeover(n_length)
		if(xh < n):
			found = True
	h = gmpy2.powmod(xh,fp*fq*u,n)
	# print(h)
	return g,h

def inverses_mod(a,b):
	inv = []
	for d in range(1, b):
		r = (d * a) % b
		if r == 1:
			inv = inv.append(d)
		# else:
		# 	raise ValueError('%d has no inverse mod %d' % (a, b))
	return inv

def main():

	file = 'DGK_keys.txt'
	p,q,u,vp,vq,fp,fq,g,h = keysDGK()
	n = p*q
	v = vp*vq
	with open(file, 'w') as f:
		f.write("%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d" % (p, q, u, vp, vq, fp, fq, g, h))

	# p,q,u,vp,vq,fp,fq,g,h = loadkey(file)


main()