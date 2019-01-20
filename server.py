#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
import numpy as np
import time
import random
try:
	import gmpy2
	HAVE_GMP = True
except ImportError:
	HAVE_GMP = False

DEFAULT_KEYSIZE = 512						# set here the default number of bits of the RSA modulus
DEFAULT_MSGSIZE = 64 						# set here the default number of bits the plaintext can have
DEFAULT_SECURITYSIZE = 100					# set here the default number of bits for the one time pads
DEFAULT_PRECISION = int(DEFAULT_MSGSIZE/2)	# set here the default number of fractional bits
NETWORK_DELAY = 0 							# set here the default network delay

seed = 42	# pick a seed for the random generator

def encrypt_vector(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.encrypt(y) for y in x]
	else: return [pubkey.encrypt(y,coins.pop()) for y in x]

def sum_encrypted_vectors(x, y):
	return [x[i] + y[i] for i in range(np.size(x))]

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
	return mpz(scalar*(2**prec))

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

def decrypt_vector(privkey, x):
    return np.array([privkey.decrypt(i) for i in x])

class Server:
	def __init__(self, n, m, N, l=DEFAULT_MSGSIZE, sigma = DEFAULT_KEYSIZE):
		filepub = "Keys/pubkey"+str(DEFAULT_KEYSIZE)+".txt"
		with open(filepub, 'r') as fin:
			data=[line.split() for line in fin]
		Np = mpz(data[0][0])
		pubkey = paillier.PaillierPublicKey(n=Np)
		self.pubkey = pubkey		
		fileH = "Data/H"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		H = np.loadtxt(fileH, delimiter=',')
		fileF = "Data/F"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		F = np.loadtxt(fileF, delimiter=',')		
		fileG0 = "Data/G0"+str(n)+"_"+str(m)+"_"+str(N)+".txt"
		G0 = np.loadtxt(fileG0, delimiter=',')
		fileK = "Data/K"+str(n)+"_"+str(m)+"_"+str(N)+".txt"		
		K = np.loadtxt(fileK, delimiter=',')	
		Kc = K[0];	Kw = K[1]
		self.Kc = int(Kc);	self.Kw = int(Kw)

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

	def compute_coeff(self,x0):
		coeff_0 = np.dot(self.mFft,x0)
		self.coeff_0 = coeff_0

	def t_iterate(self,z):
		return sum_encrypted_vectors(np.dot(self.coeff_z,z),self.coeff_0)

	def z_iterate(self,new_U,U):
		new_z = [fp(1+self.eta)*v for v in new_U]
		z = [fp(-self.eta)*v for v in U]
		return sum_encrypted_vectors(new_z,z)


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
	size_data=sock_data=bytes([]);recv_size=4096
	while total_len<size:
		sock_data=the_socket.recv(recv_size)
		if not total_data:
			if len(sock_data)>4:
				size=struct.unpack('>i', sock_data[:4])[0]
				recv_size=size
				if recv_size>4096:recv_size=4096
				total_data.append(sock_data[4:])
			else:
				size_data+=sock_data

		else:
			total_data.append(sock_data)
		total_len=sum([len(i) for i in total_data ])
	return b''.join(total_data)

def get_enc_data(received_dict,pubkey):
	return [paillier.EncryptedNumber(pubkey, int(x)) for x in received_dict]

def main():
	# Make sure the default parameters are the same as in client.py
	lf = DEFAULT_PRECISION
	n = 5	# set the number of states
	m = 5	# set the number of control inputs
	N = 7	# set the horizon length
	T = 1 	# set the number of time steps
	server = Server(n,m,N)
	server.Kc = 50; server.Kw = 20	# set the number of cold start iterations and warm start iterations
	Kc = server.Kc; Kw = server.Kw
	nc = server.nc
	server.m = m
	pubkey = server.pubkey
	U = [0]*nc

	# Create a TCP/IP socket
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	port = 10000

	# Connect the socket to the port where the server is listening
	localhost = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2] if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]]) if l][0][0]
	server_address = (localhost, port)
	print('Server: Connecting to {} port {}'.format(*server_address))
	sock.connect(server_address)	
	cont = 1		
	start = time.time()
	try:		
		while cont:
			# Send n,m,N,Kc,Kw,T
			data = send_plain_data([n,m,N,Kc,Kw,T])
			sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
			U = encrypt_vector(pubkey,fp_vector(U))
			z = [u*(2**lf) for u in U]
			K = Kc
			for i in range(0,T):
				# Receive [[x0]]
				data = json.loads(recv_size(sock))
				x0 = get_enc_data(data,pubkey)
				server.compute_coeff(x0)				
				for k in range(0,K):
					print(k)
					t = server.t_iterate(z)
					# Send [[t_k]]
					data = send_encr_data(t)
					sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					# Receive new [[U_{k+1}]]
					data = json.loads(recv_size(sock))
					new_U = get_enc_data(data,pubkey)
					z = server.z_iterate(new_U,U)
					U = new_U
				U = list(U[m:]) + list([pubkey.encrypt(0)]*m)
				z = [el*2**lf for el in U] 
				K = Kw
			cont = 0
		print(time.time() - start)
	finally:
		print('Server: Closing socket')
		sock.close()


if __name__ == '__main__':
	main()