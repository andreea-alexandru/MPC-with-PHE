#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
import numpy as np
import time
import DGK
from pathlib import Path
import os


DEFAULT_KEYSIZE = 512						# set here the default number of bits of the RSA modulus
DEFAULT_MSGSIZE = 64 						# set here the default number of bits the plaintext can have
DEFAULT_SECURITYSIZE = 100					# set here the default number of bits for the one time pads
DEFAULT_PRECISION = int(DEFAULT_MSGSIZE/2)	# set here the default number of fractional bits
DEFAULT_DGK = 160							# set here the default security size of DGK
# The message size of DGK has to be greater than 2*log2(DEFAULT_MSGSIZE), check u in DGK_pubkey
NETWORK_DELAY = 0 							# set here the default network delay

seed = 42	# pick a seed for the random generator

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False


def encrypt_vector(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.encrypt(y) for y in x]
	else: 
		return [pubkey.encrypt(y,coins.pop()) for y in x]

def encrypt_matrix(pubkey, x, coins=None):
	if (coins==None):
		return [[pubkey.encrypt(y) for y in z] for z in x]
	else: return [[pubkey.encrypt(y,coins.pop()) for y in z] for z in x]

def decrypt_vector(privkey, x):
    return np.array([privkey.decrypt(i) for i in x])

def sum_encrypted_vectors(x, y):
	return [x[i] + y[i] for i in range(np.size(x))]

def diff_encrypted_vectors(x, y):
	return [x[i] - y[i] for i in range(len(x))] 

def mul_sc_encrypted_vectors(x, y): # x is encrypted, y is plaintext
    return [y[i]*x[i] for i in range(len(x))]    

def dot_sc_encrypted_vectors(x, y): # x is encrypted, y is plaintext
    return sum(mul_sc_encrypted_vectors(x,y))

def dot_m_encrypted_vectors(x, A):
    return [dot_sc_encrypted_vectors(x,vec) for vec in A]

def encrypt_vector_DGK(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.raw_encrypt(y) for y in x]
	else: return [pubkey.raw_encrypt(y,coins.pop()) for y in x]

def decrypt_vector_DGK(privkey, x):
    return np.array([privkey.raw_decrypt0(i) for i in x])

"""We take the convention that a number x < N/3 is positive, and that a number x > 2N/3 is negative. 
	The range N/3 < x < 2N/3 allows for overflow detection.""" 

def Q_s(scalar,prec=DEFAULT_PRECISION):
	return int(scalar*(2**prec))/(2**prec)

def Q_vector(vec,prec=DEFAULT_PRECISION):
	if np.size(vec)>1:
		return [Q_s(x,prec) for x in vec]
	else:
		return Q_s(vec,prec)

def Q_matrix(mat,prec=DEFAULT_PRECISION):
	return [Q_vector(x,prec) for x in mat]

def fp(scalar,prec=DEFAULT_PRECISION):
	if prec < 0:
		return gmpy2.t_div_2exp(mpz(scalar),-prec)
	else: return mpz(gmpy2.mul(scalar,2**prec))

def fp_vector(vec,prec=DEFAULT_PRECISION):
	if np.size(vec)>1:
		return [fp(x,prec) for x in vec]
	else:
		return fp(vec,prec)

def fp_matrix(mat,prec=DEFAULT_PRECISION):
	return [fp_vector(x,prec) for x in mat]

def retrieve_fp(scalar,prec=DEFAULT_PRECISION):
	return scalar/(2**prec)

def retrieve_fp_vector(vec,prec=DEFAULT_PRECISION):
	return [retrieve_fp(x,prec) for x in vec]

def retrieve_fp_matrix(mat,prec=DEFAULT_PRECISION):
	return [retrieve_fp_vector(x,prec) for x in mat]

class Client:
	def __init__(self, l=DEFAULT_MSGSIZE):
		"""This would generate the keys on the spot"""
		# keypair = paillier.generate_paillier_keypair(n_length=DEFAULT_KEYSIZE)
		# self.pubkey, self.privkey = keypair
		# file = 'Keys/pubkey'+str(DEFAULT_KEYSIZE)+".txt"
		# with open(file, 'w') as f:
		# 	f.write("%d" % (self.pubkey.n))
		# file = 'Keys/privkey'+str(DEFAULT_KEYSIZE)+".txt"			
		# with open(file, 'w') as f:
		# 	f.write("%d\n%d" % (self.privkey.p,self.privkey.q))
		filepub = "Keys/pubkey"+str(DEFAULT_KEYSIZE)+".txt"
		with open(filepub, 'r') as fin:
			data=[line.split() for line in fin]
		Np = mpz(data[0][0])
		pubkey = paillier.PaillierPublicKey(n=Np)
		self.pubkey = pubkey

		filepriv = "Keys/privkey"+str(DEFAULT_KEYSIZE)+".txt"
		with open(filepriv, 'r') as fin:
			data=[line.split() for line in fin]
		p = mpz(data[0][0])
		q = mpz(data[1][0])
		self.privkey = paillier.PaillierPrivateKey(pubkey, p, q)

	def load_data(self,n,m,N):
		fileparam = "Data/x0"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		x0 = np.loadtxt(fileparam, delimiter='\n')
		self.x0 = x0;
		self.enc_x0 = encrypt_vector(self.pubkey,fp_vector(x0))
		filew0 = "Data/w0"+str(n)+"_"+str(m)+"_"+str(N)+".txt"		
		w0 = np.loadtxt(filew0, delimiter=',')
		hu = np.concatenate([w0[2*i*m:(2*i+1)*m] for i in range(0,N)])
		lu = np.concatenate([-w0[(2*i+1)*m:2*(i+1)*m] for i in range(0,N)])
		self.hu = hu; self.lu = lu
		fileA = "Data/A"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		A = np.loadtxt(fileA, delimiter=',')
		self.A = A
		fileB = "Data/B"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		B = np.loadtxt(fileB, delimiter=',')
		self.B = B

	def closed_loop(self,u):
		u = retrieve_fp_vector(decrypt_vector(self.privkey,u))
		print("Last input: ", ["%.8f"% i for i in u])
		with np.errstate(invalid='ignore'): self.x0 = np.dot(self.A,self.x0) + np.dot(self.B,u)
		print("Next state: ", ["%.8f"% i for i in self.x0])
		self.enc_x0 = encrypt_vector(self.pubkey,fp_vector(self.x0))


class Server1:
	def __init__(self,n,m,N,T,l=DEFAULT_MSGSIZE,sigma=DEFAULT_SECURITYSIZE):
		self.l = l
		self.sigma = sigma
		filepub = "Keys/pubkey"+str(DEFAULT_KEYSIZE)+".txt"
		with open(filepub, 'r') as fin:
			data=[line.split() for line in fin]
		Np = mpz(data[0][0])
		self.Np = Np
		pubkey = paillier.PaillierPublicKey(n=Np)
		self.pubkey = pubkey	
		self.N_len = Np.bit_length()		
		fileH = "Data/H"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		H = np.loadtxt(fileH, delimiter=',')
		fileF = "Data/F"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		F = np.loadtxt(fileF, delimiter=',')		
		fileG0 = "Data/G0"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		G0 = np.loadtxt(fileG0, delimiter=',')
		fileK = "Data/K"+str(n)+"_"+str(m)+"_"+str(N)+".txt"		
		K = np.loadtxt(fileK, delimiter=',')	
		Kc = K[0];	Kw = K[1]
		self.Kc = int(Kc);	self.Kw = int(Kw); self.T = T
		self.m = m
		nc = m*N
		self.nc = nc
		Hq = Q_matrix(H)
		eigs = np.linalg.eigvals(Hq)
		L = np.real(max(eigs))
		mu = np.real(min(eigs))
		cond = Q_s(L/mu)
		eta = Q_s((np.sqrt(cond)-1)/(np.sqrt(cond)+1))
		Hf = Q_matrix([[h/Q_s(L) for h in hv] for hv in Hq])
		Ft = F.transpose()
		Ff = Q_matrix([[Q_s(h)/Q_s(L) for h in hv] for hv in Ft])
		self.eta = eta
		self.Hf = Hf
		mFf = np.negative(Ff)
		self.mFft = fp_matrix(mFf,2*DEFAULT_PRECISION)

		coeff_z = np.eye(nc) - Hf;
		self.coeff_z = fp_matrix(coeff_z)

	def gen_rands(self,DGK_pubkey): ### CHECK SIZES
		self.DGK_pubkey = DGK_pubkey
		T = self.T
		nc = self.nc
		m = self.m
		l = self.l
		lf = DEFAULT_PRECISION
		sigma = self.sigma
		Kc = self.Kc
		Kw = self.Kw
		random_state = gmpy2.random_state(seed)
		filePath = Path('Randomness/'+str(l + sigma)+'.txt')
		if filePath.is_file():		
			with open(filePath) as file:
				# Noise for updating the iterate
				rn1 = [[[int(next(file)), int(next(file))] for x in range(0,2*nc)] for y in range(0,Kc+(T-1)*Kw)]
				# Noise for comparison
				rn2 = [[int(next(file)) for x in range(0,nc)] for y in range(0,2*Kc+2*(T-1)*Kw)]
		else:
			rn1 = [[[gmpy2.mpz_urandomb(random_state,l + sigma),gmpy2.mpz_urandomb(random_state,l + sigma)] for i in range(0,2*nc)] for k in range(0,Kc+(T-1)*Kw)]
			rn2 = [[gmpy2.mpz_urandomb(random_state,l + sigma) for i in range(0,nc)] for k in range(0,2*Kc+2*(T-1)*Kw)]
		self.obfuscations = rn1
		self.rn = rn2
		# Noise for Paillier encryption
		filePath = Path('Randomness/'+str(self.N_len)+'.txt')
		if filePath.is_file():		
			with open(filePath) as file:
				coinsP = [int(next(file)) for x in range(0,4*(T-1)*nc*Kw+ 4*nc*Kc)] 
		else:
			coinsP = [gmpy2.mpz_urandomb(random_state,self.N_len-1) for i in range(0,4*(T-1)*nc*Kw+ 4*nc*Kc)]	
		coinsP = [gmpy2.powmod(x, self.Np, self.pubkey.nsquare) for x in coinsP]
		# Noise for DGK encryption
		filePath = Path('Randomness/'+str(2*DEFAULT_DGK)+'.txt')
		if filePath.is_file():		
			with open(filePath) as file:
				coinsDGK = [int(next(file)) for x in range(0,3*(l+1)*nc*Kc + 3*(l+1)*nc*Kw*(T-1))]
		else:
			coinsDGK = [gmpy2.mpz_urandomb(random_state,2*DEFAULT_DGK) for i in range(0,3*(l+1)*nc*Kc + 3*(l+1)*nc*Kw*(T-1))]
		coinsDGK = [gmpy2.powmod(self.DGK_pubkey.h, x, self.DGK_pubkey.n) for x in coinsDGK]
		self.coinsDGK = coinsDGK
		# Noise for truncation
		filePath = Path('Randomness/'+str(l+2*lf+sigma)+'.txt')
		if filePath.is_file():		
			with open(filePath) as file:
				rn = [int(next(file)) for x in range(0,nc*Kc + nc*Kw*(T-1))]
		else:
			rn = [gmpy2.mpz_urandomb(random_state,l+2*lf+sigma) for i in range(0,nc*Kc + nc*Kw*(T-1))]
		self.fixedNoise = encrypt_vector(self.pubkey, rn) # ,coinsP[-2*nc*K:])
		er = [-fp(x,-2*lf) for x in rn]
		er = encrypt_vector(self.pubkey,er) # ,coinsP[-2*nc*K:-nc*K])
		self.er = er
		# coinsP = coinsP[:-3*nc*K]
		self.coinsP = coinsP


	def compute_coeff(self,x0):
		coeff_0 = np.dot(self.mFft,x0)
		self.coeff_0 = coeff_0

	def t_iterate(self,z):
		return sum_encrypted_vectors(np.dot(self.coeff_z,z),self.coeff_0)

	def z_iterate(self,new_U,U):
		new_z = [fp(1+self.eta)*v for v in new_U]
		z = [fp(-self.eta)*v for v in U]
		return sum_encrypted_vectors(new_z,z)

	def temporary_prec_t(self):
		nc = self.nc
		pubkey = self.pubkey
		r = [self.fixedNoise.pop() for i in range(0,nc)]
		temp_t = sum_encrypted_vectors(self.t,r)	
		return temp_t

	def randomize(self,limit):
		nc = self.nc
		a = [0]*nc
		b = [0]*nc
		for i in range(0,nc):
			a[i],b[i] = np.random.permutation([limit[i]+self.pubkey.encrypt(0),self.t[i]])
		self.a = a
		self.b = b
		# a and b have to be numbers of l bits
		return self.a,self.b

	def init_comparison_s1(self,limit):
		nc = self.nc
		l = self.l
		pubkey = self.pubkey
		r = self.r 		
		a,b = self.randomize(limit)
		z = diff_encrypted_vectors(b,a)
		z = sum_encrypted_vectors(z,encrypt_vector(pubkey,r,self.coinsP[-nc:]))
		z = sum_encrypted_vectors(z,encrypt_vector(pubkey,[2**l]*nc,self.coinsP[-2*nc:-nc]))
		self.coinsP = self.coinsP[:-2*nc]
		alpha = [gmpy2.t_mod_2exp(x,l) for x in r]
		alpha = [x.digits(2) for x in alpha]
		for i in range(0,nc):
			if (len(alpha[i]) < l):
				alpha[i] = "".join(['0'*(l-len(alpha[i])),alpha[i]])
		self.alpha = alpha
		return z

	def obfuscate(self):
		nc = self.nc
		self.a2 = [0]*nc
		self.b2 = [0]*nc
		for i in range(0,nc):
			r = self.obfuscation[i]
			self.a2[i] = self.a[i]+self.pubkey.encrypt(r[0])
			self.b2[i] = self.b[i]+self.pubkey.encrypt(r[1])
		return self.a2, self.b2

	def update_max(self,v):
		new_U = [0]*self.nc
		for i in range(0,self.nc):
			r = self.obfuscation[i]
			new_U[i] = v[i] + (self.t_comp[i]-1)*r[0] + self.t_comp[i]*(-r[1]) 
		return new_U

	def update_min(self,v):
		t = [0]*self.nc
		for i in range(0,self.nc):
			r = self.obfuscation[i]
			t[i] = v[i] + (self.t_comp[i]-1)*r[1] + self.t_comp[i]*(-r[0]) 
		return t

	def DGK_s1(self,b): 
		l = self.l
		nc = self.nc
		self.delta_A = [0]*nc
		c_all = [[0]*l]*nc
		for k in range(0,nc):
			beta = b[k]
			alpha = self.alpha[k]
			DGK_pubkey = self.DGK_pubkey
			delta_A = np.random.randint(0,2)
			self.delta_A[k] = delta_A
			prod = [0]*l
			c = [DGK_pubkey.raw_encrypt(0)]*l
			# index 0 is the MSB
			for i in range(0,l):
				if (int(alpha[i]) == 0):
					prod[i] = beta[i]
				else: prod[i] = DGK.diff_encrypted(DGK_pubkey.raw_encrypt(1,self.coinsDGK.pop()),beta[i],DGK_pubkey)
				if (int(delta_A)==int(alpha[i])):
					if i==0: c[i] = DGK_pubkey.raw_encrypt(0,self.coinsDGK.pop())
					else: 
						for iter in range(0,i):
							c[i] = DGK.add_encrypted(c[i],prod[iter],DGK_pubkey)
					if (int(delta_A) == 0):
						diff = DGK.diff_encrypted(DGK_pubkey.raw_encrypt(1,self.coinsDGK.pop()),beta[i],DGK_pubkey)
						c[i] = DGK.add_encrypted(c[i],diff,DGK_pubkey)
					else: c[i] = DGK.add_encrypted(c[i],beta[i],DGK_pubkey)
			for i in range(0,l):
				if (int(delta_A)==int(alpha[i])):
					r = gmpy2.mpz_urandomb(gmpy2.random_state(),self.sigma+self.sigma)
					c[i] = DGK.mul_sc_encrypted(c[i],r,DGK_pubkey)
				else: 
					c[i] = DGK_pubkey.raw_encrypt(gmpy2.mpz_urandomb(gmpy2.random_state(),self.sigma+self.sigma),self.coinsDGK.pop()) 
			c_all[k] = np.random.permutation(c)
		return c_all

	def compute_tDGK(self,delta_B,zdivl):
		t_comp = [0]*self.nc
		for i in range(0,self.nc):
			if (self.delta_A[i] == 1):
				t_comp[i] = delta_B[i]
			else: t_comp[i] = self.pubkey.encrypt(1) - delta_B[i]
			t_comp[i] = zdivl[i] - self.pubkey.encrypt(mpz(gmpy2.t_div_2exp(self.r[i],self.l))) - t_comp[i]
		self.t_comp = t_comp
		return t_comp

def key(serialised):
	received_dict = json.loads(serialised)
	pk = received_dict['public_key_DGK']
	n = mpz(pk['n']); g = mpz(pk['g']); h = mpz(pk['h']); u = mpz(pk['u']);
	DGK_pubkey = DGK.DGKpubkey(n=n,g=g,h=h,u=u)
	return DGK_pubkey

def send_encr_data(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	enc_with_one_pub_key = {}
	enc_with_one_pub_key = [str(x.ciphertext()) for x in encrypted_number_list]
	return json.dumps(enc_with_one_pub_key)

def send_plain_data(data):
	time.sleep(NETWORK_DELAY)
	return json.dumps([str(x) for x in data])

def recv_size(the_socket):
	#data length is packed into 4 bytes
	total_len=0;total_data=[];size=sys.maxsize
	size_data=sock_data=bytes([]);recv_size= 4096
	while total_len<size:
		sock_data=the_socket.recv(recv_size)
		if not total_data:
			if len(sock_data)>4:
				size=struct.unpack('>i', sock_data[:4])[0]
				recv_size=size
				if recv_size>262144:recv_size=262144
				total_data.append(sock_data[4:])
			else:
				size_data+=sock_data

		else:
			total_data.append(sock_data)
		total_len=sum([len(i) for i in total_data ])
	return b''.join(total_data)

def get_enc_data(received_dict,pubkey):
	return [paillier.EncryptedNumber(pubkey, int(x)) for x in received_dict]

def send_DGK_data(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	encrypted = {}
	encrypted = [str(x) for x in encrypted_number_list]
	return json.dumps(encrypted)

def send_DGK_matrix(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	encrypted = {}
	encrypted = [[str(y) for y in x] for x in encrypted_number_list]
	return json.dumps(encrypted)

def get_comp_data(received_dict):
	return [mpz(x) for x in received_dict]

def get_comp_matrix(received_dict):
	return [[mpz(y) for y in x] for x in received_dict]

def main():
	# Make sure the default parameters are the same as in server2.py
	lf = DEFAULT_PRECISION
	n = 5	# set the number of states
	m = 5	# set the number of control inputs
	N = 7	# set the horizon length
	T = 1 	# set the number of time steps
	s1 = Server1(n,m,N,T)
	s1.Kc = 50; s1.Kw = 20	# set the number of cold start iterations and warm start iterations
	Kc = s1.Kc; Kw = s1.Kw
	nc = s1.nc
	pubkey = s1.pubkey
	U = [0]*nc

	client = Client()
	client.n = n; client.m = m; client.N = N; client.Kc = Kc; client.Kw = Kw; client.T = T
	client.nc = nc
	client.load_data(n,m,N)

	s1.hu = encrypt_vector(client.pubkey,fp_vector(client.hu))
	s1.lu = encrypt_vector(client.pubkey,fp_vector(client.lu))

	# Create a TCP/IP socket
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	port = 10000

	# Connect the socket to the port where the server2 is listening
	localhost = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2] if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]]) if l][0][0]
	server_address = (localhost, port)
	print('Server1: Connecting to {} port {}'.format(*server_address))
	sock.connect(server_address)

	cont = 1		
	try:
		while cont:
			# Send n,m,N,Kc,Kw,T
			data = send_plain_data([n,m,N,Kc,Kw,T])
			sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
			U = encrypt_vector(pubkey,fp_vector(U))
			z = [uz*(2**lf) for uz in U]
			K = Kc
			# Get DGK_pubkey
			data = recv_size(sock)
			DGK_pubkey = key(data)
			s1.gen_rands(DGK_pubkey)
			sec = [0]*T
			time_s1 = [0]*K
			time_s2 = [0]*K

			start = time.time()
			# Time steps
			for i in range(0,T):
				# print("i = ", i)
				x0 = client.enc_x0
				s1.compute_coeff(x0)	
				# Optimization steps
				for k in range(0,K):
					# print("k = ", k)
					start_s1 = time.time()
					s1.t = s1.t_iterate(z)
					s1.obfuscation = s1.obfuscations[k]
					s1.r = s1.rn[k]
					temp_t = s1.temporary_prec_t()
					# Send temp_t to the target
					data = send_encr_data(temp_t)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					start_s2 = time.time()
					# Receive [(temp_t + r)*2^{-2lf}]
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					temp_tr = get_enc_data(data,pubkey)
					s1.t = sum_encrypted_vectors(temp_tr,[s1.er.pop() for i in range(0,nc)]) # t = int(t*2**16)

					# Projection on hu
					# Send z_DGK
					z_DGK = s1.init_comparison_s1(s1.hu)
					data = send_encr_data(z_DGK)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					start_s2 = time.time()
					# Receive b
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					b = get_comp_matrix(data)
					c = s1.DGK_s1(b)
					# Send c
					serialized_data = send_DGK_matrix(c)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					start_s2 = time.time()
					# Receive delta_B, zvdil
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					merged = get_enc_data(data,pubkey)
					delta_B = merged[:nc];zdivl = merged[nc:]
					t_comp = s1.compute_tDGK(delta_B,zdivl)
					# Send t_comp,a2,b2
					a2,b2 = s1.obfuscate()
					data = send_encr_data(t_comp+a2+b2)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					start_s2 = time.time()
					# Receive v
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					v = get_enc_data(data,pubkey)
					s1.t = s1.update_min(v)			

					# Projection on lu
					# Send z_DGK
					z_DGK = s1.init_comparison_s1(s1.lu)
					data = send_encr_data(z_DGK)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					start_s2 = time.time()
					# Receive b
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					b = get_comp_matrix(data)
					c = s1.DGK_s1(b)
					# Send c
					serialized_data = send_DGK_matrix(c)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					start_s2 = time.time()
					# Receive delta_B, zvdil
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					merged = get_enc_data(data,pubkey)
					delta_B = merged[:nc];zdivl = merged[nc:]
					t_comp = s1.compute_tDGK(delta_B,zdivl)
					# Send t,a2,b2
					a2,b2 = s1.obfuscate()
					data = send_encr_data(t_comp+a2+b2)
					time_s1[k] += time.time() - start_s1
					sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					start_s2 = time.time()
					# Receive v
					data = json.loads(recv_size(sock))
					time_s2[k] += time.time() - start_s2
					start_s1 = time.time()
					v = get_enc_data(data,pubkey)
					new_U = s1.update_max(v)			# [[U_{k+1}]]

					#  New [[U_{k+1}]]
					z = s1.z_iterate(new_U,U)
					U = new_U
					time_s1[k] += time.time() - start_s1
				u = U[:m]
				client.closed_loop(u);
				U = list(U[m:]) + list([pubkey.encrypt(0)]*m)
				z = [el*2**lf for el in U]
				K = Kw
				sec[i] = time.time() - start
				start = time.time()
			print('total time', sec)
			with open(os.path.abspath(str(DEFAULT_KEYSIZE)+'_'+str(DEFAULT_PRECISION)+'_results_SS'+'.txt'),'a+') as f: 
				f.write("%d, %d, %d, %d, %d, %d: " % (n,m,N,Kc,Kw,T));
				for item in sec:
  					f.write("total time %.2f " % item)
				f.write("\n")
				f.write("avg. time FGM iteration for S1: %.3f\n" % np.mean(time_s1))
				f.write("avg. time FGM iteration for S2: %.3f\n" % np.mean(time_s2))

			cont = 0				

	finally:
	# Clean up the sock
		print('Server1: Closing sock')
		sock.close()
# main()
if __name__ == '__main__':
	main()