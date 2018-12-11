#!/usr/bin/env python3
import random
from phe.util import invert, powmod, getprimeover, isqrt
from gmpy2 import mpz


DEFAULT_KEYSIZE = 512
DEFAULT_SECURITYSIZE = 160

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False


class DGKpubkey:
	def __init__(self,n,g,h,u,t=DEFAULT_SECURITYSIZE):
		self.n = n
		self.g = g;
		self.h = h;
		self.u = u;
		self.t = t;

	def raw_encrypt(self, plaintext, r_value=None):
		"""DGK encryption of a positive integer plaintext .
		You probably should be using :meth:`encrypt` instead, because it handles positive and negative ints and floats.

		Args:
			plaintext (int): a positive integer < :attr:`n` to be DGK
			encrypted. Typically this is an encoding of the actual 
			number you want to encrypt.
			r_value (int): obfuscator for the ciphertext; by default (i.e.
			r_value is None), a random value of 2t bits is used.
		Returns:
			int: DGK encryption of plaintext.

		Raises:
			TypeError: if plaintext is not an int or mpz.
		"""
		if not isinstance(plaintext, int) and not isinstance(plaintext, type(mpz(1))):
			raise TypeError('Expected int type plaintext but got: %s' %
                            type(plaintext))

		nude_ciphertext = powmod(self.g, plaintext, self.n)
		r = r_value or powmod(self.h, self.get_random_lt_2t(), self.n) # Pass the precomputed obfuscator
		obfuscator = r		

		return (nude_ciphertext * obfuscator) % self.n

	def get_random_lt_2t(self):
		"""Return a cryptographically random number less than :attr:`n`"""
		t2 = 2*DEFAULT_SECURITYSIZE
		return random.SystemRandom().randrange(1, 2**t2)

class DGKprivkey:
	def __init__(self,p,q,v,pubkey):
		self.p = p
		self.q = q
		self.v = v
		self.pubkey = pubkey

	def raw_decrypt0(self, ciphertext):
		"""Decrypt raw ciphertext and return raw plaintext.

		Args:
			ciphertext (int): (usually from :meth:`EncryptedNumber.ciphertext()`)
				that is to be DGK decrypted.

		Returns:
			int: DGK decryption of ciphertext. 0 if the plaintext is 0 by 
			checking cyphertext^v mod p == 1, 1 otherwise
		"""

		c = powmod(ciphertext, self.v, self.p)
		if c==1:
			return 0
		else: return 1



def loadkey(file):
	# To generate the keys, run genDGK.py
	with open(file, 'r') as fin:
		data=[line.split() for line in fin]
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

def add_encrypted(a,b,pubkey):
	return gmpy2.t_mod(gmpy2.mul(a,b),pubkey.n)

def diff_encrypted(a,b,pubkey):
	return add_encrypted(a, gmpy2.invert(b,pubkey.n), pubkey)

def mul_sc_encrypted(a,b,pubkey):
	return gmpy2.powmod(a,b,pubkey.n)
