#!/usr/bin/env python3
# Adapted from pyphe. Transform it to be a fixed-point library, no encoding
#
"""Paillier encryption library for partially homomorphic encryption."""
import random
import hashlib
import math
import sys
import numpy
try:
    from collections.abc import Mapping
except ImportError:
    Mapping = dict

from util import invert, powmod, getprimeover, isqrt
from gmpy2 import mpz

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False

DEFAULT_KEYSIZE = 1024


def generate_paillier_keypair(private_keyring=None, n_length=DEFAULT_KEYSIZE):
    """Return a new :class:`PaillierPublicKey` and :class:`PaillierPrivateKey`.

    Add the private key to *private_keyring* if given.

    Args:
      private_keyring (PaillierPrivateKeyring): a
        :class:`PaillierPrivateKeyring` on which to store the private
        key.
      n_length: key size in bits.

    Returns:
      tuple: The generated :class:`PaillierPublicKey` and
      :class:`PaillierPrivateKey`
    """
    p = q = n = None
    n_len = 0
    while n_len != n_length:
        p = getprimeover(n_length // 2)
        q = p
        while q == p:
            q = getprimeover(n_length // 2)
        n = p * q
        n_len = n.bit_length()

    public_key = PaillierPublicKey(n)
    private_key = PaillierPrivateKey(public_key, p, q)

    if private_keyring is not None:
        private_keyring.add(private_key)

    return public_key, private_key

class PaillierPublicKey(object):
    """Contains a public key and associated encryption methods.

    Args:

      n (int): the modulus of the public key - see Paillier's paper.

    Attributes:
      g (int): part of the public key - see Paillier's paper.
      n (int): part of the public key - see Paillier's paper.
      nsquare (int): :attr:`n` ** 2, stored for frequent use.
      max_int (int): Maximum int that may safely be stored. This can be
        increased, if you are happy to redefine "safely" and lower the
        chance of detecting an integer overflow.
    """
    def __init__(self, n):
        self.g = n + 1
        self.n = n
        self.nsquare = n * n
        self.max_int = n // 3 - 1

    def __repr__(self):
        nsquare = self.nsquare.to_bytes(1024, 'big')
        g = self.g.to_bytes(1024, 'big')
        publicKeyHash = hashlib.sha1(nsquare + g).hexdigest()
        return "<PaillierPublicKey {}>".format(publicKeyHash[:10])

    def __eq__(self, other):
        return self.n == other.n

    def __hash__(self):
        return hash(self.n)

    def raw_encrypt(self, plaintext, r_value=None):
        """Paillier encryption of a positive integer plaintext < :attr:`n`.

        You probably should be using :meth:`encrypt` instead, because it
        handles positive and negative ints and floats.

        Args:
          plaintext (int): a positive integer < :attr:`n` to be Paillier
            encrypted. Typically this is an encoding of the actual
            number you want to encrypt.
          r_value (int): obfuscator for the ciphertext; by default (i.e.
            r_value is None), a random value is used.

        Returns:
          int: Paillier encryption of plaintext.

        Raises:
          TypeError: if plaintext is not an int or mpz.
        """
        if not isinstance(plaintext, int) and not isinstance(plaintext, type(mpz(1))) and not isinstance(plaintext, numpy.int64):
            raise TypeError('Expected int type plaintext but got: %s' %
                            type(plaintext))

        if self.n - self.max_int <= plaintext < self.n:
            # Very large plaintext, take a sneaky shortcut using inverses
            neg_plaintext = self.n - plaintext  # = abs(plaintext - nsquare)
            neg_ciphertext = (self.n * neg_plaintext + 1) % self.nsquare
            nude_ciphertext = invert(neg_ciphertext, self.nsquare)
        else:
            # we chose g = n + 1, so that we can exploit the fact that 
            # (n+1)^plaintext = n*plaintext + 1 mod n^2
            nude_ciphertext = (self.n * plaintext + 1) % self.nsquare

        # r = r_value or self.get_random_lt_n()
        # obfuscator = powmod(r, self.n, self.nsquare)
        r = r_value or powmod(self.get_random_lt_n(), self.n, self.nsquare) # Pass the precomputed obfuscator
        obfuscator = r

        return (nude_ciphertext * obfuscator) % self.nsquare

    def get_random_lt_n(self):
        """Return a cryptographically random number less than :attr:`n`"""
        return random.SystemRandom().randrange(1, self.n)

    def encrypt(self, value, r_value=None): ### Do raw_encrypt
        """Encode and Paillier encrypt a real number *value*.

        Args:
          value: an int or mpz to be encrypted.
            If int, it must satisfy abs(*value*) < :attr:`n`/3.
          r_value (int): obfuscator for the ciphertext; by default (i.e.
            if *r_value* is None), a random value is used.

        Returns:
          EncryptedNumber: An encryption of *value*.

        Raises:
          ValueError: if *value* is out of range or *precision* is so
            high that *value* is rounded to zero.
        """

        obfuscator = r_value or 1
        if value < 0:
            value = value + self.n
        ciphertext = self.raw_encrypt(value, r_value=obfuscator)
        encrypted_number = EncryptedNumber(self, ciphertext)
        if r_value is None:
            encrypted_number.obfuscate()
        return encrypted_number

        return self.encrypt_encoded(encoding, r_value)


class PaillierPrivateKey(object):
    """Contains a private key and associated decryption method.

    Args:
      public_key (:class:`PaillierPublicKey`): The corresponding public
        key.
      p (int): private secret - see Paillier's paper.
      q (int): private secret - see Paillier's paper.

    Attributes:
      public_key (PaillierPublicKey): The corresponding public
        key.
      p (int): private secret - see Paillier's paper.
      q (int): private secret - see Paillier's paper.
      psquare (int): p^2
      qsquare (int): q^2
      p_inverse (int): p^-1 mod q
      hp (int): h(p) - see Paillier's paper.
      hq (int): h(q) - see Paillier's paper.
    """
    def __init__(self, public_key, p, q):
        if not p*q == public_key.n:
            raise ValueError('given public key does not match the given p and q.')
        if p == q: #check that p and q are different, otherwise we can't compute p^-1 mod q
            raise ValueError('p and q have to be different')
        self.public_key = public_key
        if q < p: #ensure that p < q. 
            self.p = q
            self.q = p
        else:
            self.p = p
            self.q = q
        self.psquare = self.p * self.p
        
        self.qsquare = self.q * self.q
        self.p_inverse = invert(self.p, self.q)
        self.hp = self.h_function(self.p, self.psquare);
        self.hq = self.h_function(self.q, self.qsquare);
        self.n = public_key.n

    @staticmethod
    def from_totient(public_key, totient):
        """given the totient, one can factorize the modulus
        
        The totient is defined as totient = (p - 1) * (q - 1),
        and the modulus is defined as modulus = p * q
        
        Args:
          public_key (PaillierPublicKey): The corresponding public
            key
          totient (int): the totient of the modulus
          
        Returns:
          the :class:`PaillierPrivateKey` that corresponds to the inputs
          
        Raises:
          ValueError: if the given totient is not the totient of the modulus
            of the given public key
        """
        p_plus_q = public_key.n - totient + 1
        p_minus_q = isqrt(p_plus_q * p_plus_q - public_key.n * 4)
        q = (p_plus_q - p_minus_q) // 2
        p = p_plus_q - q
        if not p*q == public_key.n:
            raise ValueError('given public key and totient do not match.')
        return PaillierPrivateKey(public_key, p, q)
        
    def __repr__(self):
        pub_repr = repr(self.public_key)
        return "<PaillierPrivateKey for {}>".format(pub_repr)

    def decrypt(self, encrypted_number):
        """Return the decrypted & decoded plaintext of *encrypted_number*.

        Args:
          encrypted_number (EncryptedNumber): an
            :class:`EncryptedNumber` with a public key that matches this
            private key.

        Returns:
          the int or float that `EncryptedNumber` was holding. N.B. if
            the number returned is an integer, it will not be of type
            float.

        Raises:
          TypeError: If *encrypted_number* is not an
            :class:`EncryptedNumber`.
          ValueError: If *encrypted_number* was encrypted against a
            different key.
        """
        if not isinstance(encrypted_number, EncryptedNumber):
            raise TypeError('Expected encrypted_number to be an EncryptedNumber'
                            ' not: %s' % type(encrypted_number))

        if self.public_key != encrypted_number.public_key:
            raise ValueError('encrypted_number was encrypted against a '
                             'different key!')

        return self.raw_decrypt(encrypted_number.ciphertext(be_secure=False))
        # encoded = Encoding(self.public_key, encoded,
        #                      encrypted_number.exponent)
        # return encoded.decode()

    def raw_decrypt(self, ciphertext):
        """Decrypt raw ciphertext and return raw plaintext.

        Args:
          ciphertext (int): (usually from :meth:`EncryptedNumber.ciphertext()`)
            that is to be Paillier decrypted.

        Returns:
          int: Paillier decryption of ciphertext. This is a positive
          integer < :attr:`public_key.n`.

        Raises:
          TypeError: if ciphertext is not an int.
        """
        if not isinstance(ciphertext, int) and not isinstance(ciphertext, type(mpz(1))) and not isinstance(scalar, numpy.int64):
            raise TypeError('Expected ciphertext to be an int, not: %s' %
                type(ciphertext))

        decrypt_to_p = self.l_function(powmod(ciphertext, self.p-1, self.psquare), self.p) * self.hp % self.p
        decrypt_to_q = self.l_function(powmod(ciphertext, self.q-1, self.qsquare), self.q) * self.hq % self.q
        value = self.crt(decrypt_to_p, decrypt_to_q)
        if value < self.n/3:
            return value
        else:
            return value - self.n
    
    def h_function(self, x, xsquare):
        """Computes the h-function as defined in Paillier's paper page 12, 
        'Decryption using Chinese-remaindering'.
        """
        return invert(self.l_function(powmod(self.public_key.g, x - 1, xsquare),x), x)
            
    
    def l_function(self, x, p):
        """Computes the L function as defined in Paillier's paper. That is: L(x,p) = (x-1)/p"""
        return (x - 1) // p
    
    def crt(self, mp, mq):
        """The Chinese Remainder Theorem as needed for decryption. Returns the solution modulo n=pq.
        
        Args:
           mp(int): the solution modulo p.
           mq(int): the solution modulo q.
       """
        u = (mq - mp) * self.p_inverse % self.q
        return mp + (u * self.p)

    def __eq__(self, other):
        return (self.p == other.p and self.q == other.q)

    def __hash__(self):
        return hash((self.p, self.q))

class PaillierPrivateKeyring(Mapping):
    """Holds several private keys and can decrypt using any of them.

    Acts like a dict, supports :func:`del`, and indexing with **[]**,
    but adding keys is done using :meth:`add`.

    Args:
      private_keys (list of PaillierPrivateKey): an optional starting
        list of :class:`PaillierPrivateKey` instances.
    """
    def __init__(self, private_keys=None):
        if private_keys is None:
            private_keys = []
        public_keys = [k.public_key for k in private_keys]
        self.__keyring = dict(zip(public_keys, private_keys))

    def __getitem__(self, key):
        return self.__keyring[key]

    def __len__(self):
        return len(self.__keyring)

    def __iter__(self):
        return iter(self.__keyring)

    def __delitem__(self, public_key):
        del self.__keyring[public_key]

    def add(self, private_key):
        """Add a key to the keyring.

        Args:
          private_key (PaillierPrivateKey): a key to add to this keyring.
        """
        if not isinstance(private_key, PaillierPrivateKey):
            raise TypeError("private_key should be of type PaillierPrivateKey, "
                            "not %s" % type(private_key))
        self.__keyring[private_key.public_key] = private_key

    def decrypt(self, encrypted_number):
        """Return the decrypted & decoded plaintext of *encrypted_number*.

        Args:
          encrypted_number (EncryptedNumber): encrypted against a known public
            key, i.e., one for which the private key is on this keyring.

        Returns:
          the int or float that *encrypted_number* was holding. N.B. if
          the number returned is an integer, it will not be of type
          float.

        Raises:
          KeyError: If the keyring does not hold the private key that
            decrypts *encrypted_number*.
        """
        relevant_private_key = self.__keyring[encrypted_number.public_key]
        return relevant_private_key.decrypt(encrypted_number)


class EncryptedNumber(object):
    """Represents the Paillier encryption of a float or int.

    Typically, an `EncryptedNumber` is created by
    :meth:`PaillierPublicKey.encrypt`. You would only instantiate an
    `EncryptedNumber` manually if you are de-serializing a number
    someone else encrypted.


    Paillier encryption is only defined for non-negative integers less
    than :attr:`PaillierPublicKey.n`. :class:`EncodedNumber` provides
    an encoding scheme for floating point and signed integers that is
    compatible with the partially homomorphic properties of the Paillier
    cryptosystem:

    1. D(E(a) * E(b)) = a + b
    2. D(E(a)**b)     = a * b

    where `a` and `b` are ints or floats, `E` represents encoding then
    encryption, and `D` represents decryption then decoding.

    Args:
      public_key (PaillierPublicKey): the :class:`PaillierPublicKey`
        against which the number was encrypted.
      ciphertext (int): encrypted representation of the encoded number.
      exponent (int): used by :class:`EncodedNumber` to keep track of
        fixed precision. Usually negative.

    Attributes:
      public_key (PaillierPublicKey): the :class:`PaillierPublicKey`
        against which the number was encrypted.
      exponent (int): used by :class:`EncodedNumber` to keep track of
        fixed precision. Usually negative.

    Raises:
      TypeError: if *ciphertext* is not an int, or if *public_key* is
        not a :class:`PaillierPublicKey`.
    """
    def __init__(self, public_key, ciphertext):
        self.public_key = public_key
        self.__ciphertext = ciphertext
        self.__is_obfuscated = False
        if isinstance(self.ciphertext, EncryptedNumber):
            raise TypeError('ciphertext should be an integer')
        if not isinstance(self.public_key, PaillierPublicKey):
            raise TypeError('public_key should be a PaillierPublicKey')

    def __add__(self, other):
        """Add an int, float, `EncryptedNumber` or `EncodedNumber`."""
        if isinstance(other, EncryptedNumber):
            return self._add_encrypted(other)
        else:
            return self._add_scalar(other)

    def __radd__(self, other):
        """Called when Python evaluates `34 + <EncryptedNumber>`
        Required for builtin `sum` to work.
        """
        return self.__add__(other)

    def __mul__(self, other):
        """Multiply by an int."""
        if isinstance(other, EncryptedNumber):
            raise NotImplementedError('Good luck with that...')
        if other < 0:
            other = other + self.public_key.n
        product = self._raw_mul(other)

        return EncryptedNumber(self.public_key, product)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __sub__(self, other):
        return self + (other * -1)

    def __rsub__(self, other):
        return other + (self * -1)

    def __truediv__(self, scalar):
        return self.__mul__(1 / scalar)

    def ciphertext(self, be_secure=True):
        """Return the ciphertext of the EncryptedNumber.

        Choosing a random number is slow. Therefore, methods like
        :meth:`__add__` and :meth:`__mul__` take a shortcut and do not
        follow Paillier encryption fully - every encrypted sum or
        product should be multiplied by r **
        :attr:`~PaillierPublicKey.n` for random r < n (i.e., the result
        is obfuscated). Not obfuscating provides a big speed up in,
        e.g., an encrypted dot product: each of the product terms need
        not be obfuscated, since only the final sum is shared with
        others - only this final sum needs to be obfuscated.

        Not obfuscating is OK for internal use, where you are happy for
        your own computer to know the scalars you've been adding and
        multiplying to the original ciphertext. But this is *not* OK if
        you're going to be sharing the new ciphertext with anyone else.

        So, by default, this method returns an obfuscated ciphertext -
        obfuscating it if necessary. If instead you set `be_secure=False`
        then the ciphertext will be returned, regardless of whether it
        has already been obfuscated. We thought that this approach,
        while a little awkward, yields a safe default while preserving
        the option for high performance.

        Args:
          be_secure (bool): If any untrusted parties will see the
            returned ciphertext, then this should be True.

        Returns:
          an int, the ciphertext. If `be_secure=False` then it might be
            possible for attackers to deduce numbers involved in
            calculating this ciphertext.
        """
        if be_secure and not self.__is_obfuscated:
            self.obfuscate()

        return self.__ciphertext

    def obfuscate(self):
        """Disguise ciphertext by multiplying by r ** n with random r.

        This operation must be performed for every `EncryptedNumber`
        that is sent to an untrusted party, otherwise eavesdroppers
        might deduce relationships between this and an antecedent
        `EncryptedNumber`.

        For example::

            enc = public_key.encrypt(1337)
            send_to_nsa(enc)       # NSA can't decrypt (we hope!)
            product = enc * 3.14
            send_to_nsa(product)   # NSA can deduce 3.14 by bruteforce attack
            product2 = enc * 2.718
            product2.obfuscate()
            send_to_nsa(product)   # NSA can't deduce 2.718 by bruteforce attack
        """
        r = self.public_key.get_random_lt_n()
        r_pow_n = powmod(r, self.public_key.n, self.public_key.nsquare)
        self.__ciphertext = self.__ciphertext * r_pow_n % self.public_key.nsquare
        self.__is_obfuscated = True

    def _add_scalar(self, scalar):
        """Returns E(a + b), given self=E(a) and b.

        Args:
          scalar: an int or float b, to be added to `self`.

        Returns:
          EncryptedNumber: E(a + b), calculated by encrypting b and
            taking the product of E(a) and E(b) modulo
            :attr:`~PaillierPublicKey.n` ** 2.
        """

        a, b = self, scalar

        # Don't bother to salt/obfuscate in a basic operation, do it
        # just before leaving the computer.
        encrypted_scalar = a.public_key.raw_encrypt(b, 1)

        sum_ciphertext = a._raw_add(a.ciphertext(False), encrypted_scalar)
        return EncryptedNumber(a.public_key, sum_ciphertext)

    def _add_encrypted(self, other):
        """Returns E(a + b) given E(a) and E(b).

        Args:
          other (EncryptedNumber): an `EncryptedNumber` to add to self.

        Returns:
          EncryptedNumber: E(a + b), calculated by taking the product
            of E(a) and E(b) modulo :attr:`~PaillierPublicKey.n` ** 2.

        Raises:
          ValueError: if numbers were encrypted against different keys.
        """
        if self.public_key != other.public_key:
            raise ValueError("Attempted to add numbers encrypted against "
                             "different public keys!")

        a, b = self, other

        sum_ciphertext = a._raw_add(a.ciphertext(False), b.ciphertext(False))
        return EncryptedNumber(a.public_key, sum_ciphertext)

    def _raw_add(self, e_a, e_b):
        """Returns the integer E(a + b) given ints E(a) and E(b).

        N.B. this returns an int, not an `EncryptedNumber`, and ignores
        :attr:`ciphertext`

        Args:
          e_a (int): E(a), first term
          e_b (int): E(b), second term

        Returns:
          int: E(a + b), calculated by taking the product of E(a) and
            E(b) modulo :attr:`~PaillierPublicKey.n` ** 2.
        """
        return e_a * e_b % self.public_key.nsquare

    def _raw_mul(self, plaintext):
        """Returns the integer E(a * plaintext), where E(a) = ciphertext

        Args:
          plaintext (int): number by which to multiply the
            `EncryptedNumber`. *plaintext* is typically an encoding.
            0 <= *plaintext* < :attr:`~PaillierPublicKey.n`

        Returns:
          int: Encryption of the product of `self` and the scalar
            encoded in *plaintext*.

        Raises:
          TypeError: if *plaintext* is not an int.
          ValueError: if *plaintext* is not between 0 and
            :attr:`PaillierPublicKey.n`.
        """
        if not isinstance(plaintext, int) and not isinstance(plaintext, type(mpz(1))) and not isinstance(plaintext, numpy.int64):
            raise TypeError('Expected ciphertext to be int, not %s' %
                type(plaintext))

        if plaintext < 0 or plaintext >= self.public_key.n:
            raise ValueError('Scalar out of bounds: %i' % plaintext)

        if self.public_key.n - self.public_key.max_int <= plaintext:
            # Very large plaintext, play a sneaky trick using inverses
            neg_c = invert(self.ciphertext(False), self.public_key.nsquare)
            neg_scalar = self.public_key.n - plaintext
            return powmod(neg_c, neg_scalar, self.public_key.nsquare)
        else:
            return powmod(self.ciphertext(False), plaintext, self.public_key.nsquare)
