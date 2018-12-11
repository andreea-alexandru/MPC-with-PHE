#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
from pathlib import Path
import numpy as np
import time
import DGK
import os


DEFAULT_KEYSIZE = 512						# set here the default number of bits of the RSA modulus
DEFAULT_MSGSIZE = 64 						# set here the default number of bits the plaintext can have
DEFAULT_SECURITYSIZE = 100					# set here the default number of bits for the one time pads
DEFAULT_PRECISION = int(DEFAULT_MSGSIZE/2)	# set here the default number of fractional bits
DEFAULT_DGK = 160							# set here the default security size of DGK
# The message size of DGK has to be greater than 2*log2(DEFAULT_MSGSIZE), check u in DGK_pubkey
NETWORK_DELAY = 0 		

seed = 43	# pick a seed for the random generator


try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False


def encrypt_vector(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.encrypt(y) for y in x]
	else: return [pubkey.encrypt(y,coins.pop()) for y in x]

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

class Server2:
	def __init__(self, l=DEFAULT_MSGSIZE,t_DGK=DEFAULT_DGK,sigma=DEFAULT_SECURITYSIZE):
		filepub = "Keys/pubkey"+str(DEFAULT_KEYSIZE)+".txt"
		with open(filepub, 'r') as fin:
			data=[line.split() for line in fin]
		Np = mpz(data[0][0])
		self.N_len = Np.bit_length()		
		pubkey = paillier.PaillierPublicKey(n=Np)
		self.pubkey = pubkey

		filepriv = "Keys/privkey"+str(DEFAULT_KEYSIZE)+".txt"
		with open(filepriv, 'r') as fin:
			data=[line.split() for line in fin]
		p = mpz(data[0][0])
		q = mpz(data[1][0])
		self.privkey = paillier.PaillierPrivateKey(pubkey, p, q)
		self.l = l
		self.sigma = sigma
		self.t_DGK = t_DGK
		self.generate_DGK()


	def params(self,n,m,N,Kc,Kw,T):	### CHECK SIZES FOR COINS
		self.Kc = Kc
		self.Kw = Kw
		nc = m*N
		self.nc = nc
		t2 = 2*self.t_DGK
		N_len = self.N_len
		random_state = gmpy2.random_state(seed)
		# Noise for Paillier encryption
		filePath = Path('Randomness/'+str(N_len)+'.txt')
		if filePath.is_file():
			with open(filePath) as file:
				coinsP = [int(next(file)) for x in range(0,7*(T-1)*nc*Kw + 7*nc*Kc)]
		else: 
			coinsP = [gmpy2.mpz_urandomb(random_state,N_len-1) for i in range(0,7*(T-1)*nc*Kw + 7*nc*Kc)]
		coinsP = [gmpy2.powmod(x, self.pubkey.n, self.pubkey.nsquare) for x in coinsP]
		self.coinsP = coinsP
		filePath = Path('Randomness/'+str(t2)+'.txt')
		if filePath.is_file():		
			with open('Randomness/'+str(t2)+'.txt') as file:
				coinsDGK = [int(next(file)) for x in range(0,2*(self.l+1)*nc*Kc + 2*(self.l+1)*nc*Kw*(T-1))]
		else:
			coinsDGK = [gmpy2.mpz_urandomb(random_state,t2) for i in range(0,2*(self.l+1)*nc*Kc + 2*(self.l+1)*nc*Kw*(T-1))]
		coinsDGK = [gmpy2.powmod(self.DGK_pubkey.h, x, self.DGK_pubkey.n) for x in coinsDGK]
		self.coinsDGK = coinsDGK
		# self.delta_B = [0]*nc

	def init_comparison_s2(self,msg):
		l = self.l
		z = decrypt_vector(self.privkey,msg)
		z = [mpz(x) for x in z]
		self.z = z
		beta = [gmpy2.t_mod_2exp(x,l) for x in z]
		beta = [x.digits(2) for x in beta]
		for i in range(0,self.nc):
			if (len(beta[i]) < l):
				beta[i] = "".join(['0'*(l-len(beta[i])),beta[i]])
		self.beta = beta


	def generate_DGK(self):
		file = 'Keys/DGK_keys.txt'
		p,q,u,vp,vq,fp,fq,g,h = DGK.loadkey(file)
		n = p*q
		self.DGK_pubkey = DGK.DGKpubkey(n,g,h,u)
		self.DGK_privkey = DGK.DGKprivkey(p,q,vp,self.DGK_pubkey)


	def DGK_s2(self,c_all):
		l = self.l
		nc = self.nc
		for i in range(0,nc):
			c = c_all[i]
			self.delta_B[i] = 0
			for j in range(0,l):
				if (int(self.DGK_privkey.raw_decrypt0(c[j])) == 0):
					self.delta_B[i] = 1
					break
		db = encrypt_vector(self.pubkey,self.delta_B,self.coinsP[-nc:]); z = encrypt_vector(self.pubkey,[mpz(gmpy2.t_div_2exp(self.z[i],l)) for i in range(0,nc)],self.coinsP[-2*nc:-nc])
		self.coinsP = self.coinsP[:-2*nc]
		return db,z


	def choose_max(self,a,b):
		nc = self.nc
		v = [0]*nc
		for i in range(0,nc):
			if int(self.t_comp[i])==0: 
				v[i] = a[i] + self.pubkey.encrypt(0,self.coinsP.pop())
			else: v[i] = b[i] + self.pubkey.encrypt(0,self.coinsP.pop())
		return v

	def choose_min(self,a,b):
		nc = self.nc
		v = [0]*nc
		for i in range(0,nc):
			if int(self.t_comp[i])==1: 
				v[i] = a[i] + self.pubkey.encrypt(0,self.coinsP.pop())
			else: v[i] = b[i] + self.pubkey.encrypt(0,self.coinsP.pop())
		return v

def keys(DGK_pubkey):
	pubkeys = {}
	pubkeys['public_key_DGK'] = {'n': int(DGK_pubkey.n), 'g':int(DGK_pubkey.g),'h':int(DGK_pubkey.h), 'u':int(DGK_pubkey.u)}
	serialized_pubkeys = json.dumps(pubkeys)
	return serialized_pubkeys

def get_enc_data(received_dict,pubkey):
	return [paillier.EncryptedNumber(pubkey, int(x)) for x in received_dict]

def get_plain_data(data):
	return [int(x) for x in data]

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

def send_encr_data(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	encrypted = {}
	encrypted = [str(x.ciphertext()) for x in encrypted_number_list]
	return json.dumps(encrypted)

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

def get_DGK_data(received_dict):
	return [mpz(x) for x in received_dict]

def get_DGK_matrix(received_dict):
	return [[mpz(y) for y in x] for x in received_dict]

def main():
	# Make sure the default parameters are the same as in server1.py
	lf = DEFAULT_PRECISION
	s2 = Server2()
	l = s2.l
	pubkey = s2.pubkey
	privkey = s2.privkey
	DGK_pubkey = s2.DGK_pubkey
	serialized_pubkey = keys(DGK_pubkey)

	# Create a TCP/IP socket
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	print('Server2: Socket successfully created')
	port = 10000
	# Bind the socket to the port
	localhost = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2] if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]]) if l][0][0]
	server_address = (localhost, port)
	print('Server2: Starting up on {} port {}'.format(*server_address))
	sock.bind(server_address)

	# Listen for incoming connections
	sock.listen(1)      
	print('Server2: Socket is listening')
	connection, client_address = sock.accept()	
	try:
		print('Server2: Connection from', client_address)
		# data = recv_size(connection)
		data = json.loads(recv_size(connection))
		if data:
			n,m,N,Kc,Kw,T = get_plain_data(data)
			s2.params(n,m,N,Kc,Kw,T)
			nc = m*N
			K = Kc
			# Send DGK public key
			connection.sendall(struct.pack('>i', len(serialized_pubkey))+serialized_pubkey.encode('utf-8'))		
			for i in range(0,T):
				for k in range(0,K):
					# Receive temp_t + r
					data = json.loads(recv_size(connection))
					temp_tr = get_enc_data(data,pubkey)
					temp_tr = decrypt_vector(privkey,temp_tr)
					temp_tr = fp_vector(temp_tr,-2*lf)
					temp_tr = encrypt_vector(s2.pubkey,temp_tr,s2.coinsP[-nc:])
					s2.coinsP = s2.coinsP[:-nc]
					# Send temp_tr
					serialized_data = send_encr_data(temp_tr)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))

					# Projection on hu
					s2.delta_B = [0]*nc
					# Receive z_DGK
					data = json.loads(recv_size(connection))
					z_DGK = get_enc_data(data,pubkey)
					s2.init_comparison_s2(z_DGK)
					s2.coinsDGK = s2.coinsDGK[:-nc]
					b = [[0]*l]*nc
					b = [encrypt_vector_DGK(DGK_pubkey,[int(s2.beta[i][j]) for j in range(0,l)],s2.coinsDGK[-(i+1)*l:-i*l] or s2.coinsDGK[-l:]) for i in range(0,nc)]
					s2.coinsDGK = s2.coinsDGK[:-l*nc]
					# Send b = bits of beta
					serialized_data = send_DGK_matrix(b)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# Receive c
					data = json.loads(recv_size(connection))
					c = get_DGK_matrix(data)
					delta_B, zdivl = s2.DGK_s2(c)
					# Send delta_B, zdivl
					serialized_data = send_encr_data(delta_B+zdivl)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# Receive t,a2,bs
					data = json.loads(recv_size(connection))
					merged = get_enc_data(data,pubkey)
					t_comp = merged[:nc]; a2 = merged[nc:2*nc]; b2 = merged[2*nc:]
					s2.t_comp = decrypt_vector(s2.privkey,t_comp)
					v = s2.choose_min(a2,b2)
					# Send v
					serialized_data = send_encr_data(v)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))

					# Projection on lu
					s2.delta_B = [0]*nc
					# Receive z_DGK
					data = json.loads(recv_size(connection))
					z_DGK = get_enc_data(data,pubkey)
					s2.init_comparison_s2(z_DGK)
					b = [[0]*l]*nc
					b = [encrypt_vector_DGK(DGK_pubkey,[int(s2.beta[i][j]) for j in range(0,l)],s2.coinsDGK[-(i+1)*l:-i*l] or s2.coinsDGK[-l:]) for i in range(0,nc)]
					s2.coinsDGK = s2.coinsDGK[:-l*nc]
					# Send b
					serialized_data = send_DGK_matrix(b)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# Receive c
					data = json.loads(recv_size(connection))
					c = get_DGK_matrix(data)
					delta_B, zdivl = s2.DGK_s2(c)
					# Send delta_B, zdivl
					serialized_data = send_encr_data(delta_B+zdivl)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# Receive t,a2,bs
					data = json.loads(recv_size(connection))
					merged = get_enc_data(data,pubkey)
					t_comp = merged[:nc]; a2 = merged[nc:2*nc]; b2 = merged[2*nc:]
					s2.t_comp = decrypt_vector(s2.privkey,t_comp)
					v = s2.choose_max(a2,b2)
					# Send v
					serialized_data = send_encr_data(v)
					connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))

					K = Kw				


	finally:
		print('Server2: Closing socket')
		connection.close()				


# main()
if __name__ == '__main__':
	main()