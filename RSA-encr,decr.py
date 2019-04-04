import random


'''
Euclid's algorithm for determining the greatest common divisor
Use iteration to make it faster for larger integers
'''
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a


'''
Euclid's extended algorithm for finding the multiplicative inverse of two numbers
'''
def multiplicative_inverse(e, phi):
    d = 0
    x1 = 0
    x2 = 1
    y1 = 1
    temp_phi = phi
    while e > 0:
        temp1 = temp_phi / e
        temp2 = temp_phi - temp1 * e
        temp_phi = e
        e = temp2
        x = x2 - temp1 * x1
        y = d - temp1 * y1
        x2 = x1
        x1 = x
        d = y1
        y1 = y
    if temp_phi == 1:
        return d + phi


'''
Tests to see if a number is prime.
'''
def is_prime(num):
    if num == 2:
        return True
    if num < 2 or num % 2 == 0:
        return False
    for n in xrange(3, int(num ** 0.5) + 2, 2):
        if num % n == 0:
            return False
    return True


def generate_keypair(p, q):
    '''if not (is_prime(p) and is_prime(q)):
        raise ValueError('Both numbers must be prime.')
    elif p == q:
        raise ValueError('p and q cannot be equal')'''
    # n = pq
    n = p * q

    # Phi is the totient of n
    phi = (p - 1) * (q - 1)

    # Choose an integer e such that e and phi(n) are coprime
    e = random.randrange(1, phi)

    # Use Euclid's Algorithm to verify that e and phi(n) are comprime
    g = gcd(e, phi)
    while g != 1:
        e = random.randrange(1, phi)
        g = gcd(e, phi)

    # Use Extended Euclid's Algorithm to generate the private key
    d = multiplicative_inverse(e, phi)

    # Return public and private keypair
    # Public key is (e, n) and private key is (d, n)
    return ((e, n), (d, n))


def encrypt(pk, plaintext):
    # Unpack the key into it's components
    key, n = pk
    # Generate the ciphertext based on the plaintext and key using a^b mod m
    cipher = [pow(plaintext,key,n)]
    return cipher


def decrypt(pk, ciphertext):
    # Unpack the key into its components
    key, n = pk
    # Generate the plaintext based on the ciphertext and key using a^b mod m
    plain = pow(ciphertext, key, n)
    plaintext=''
    while(plain):
        plaintext = plaintext+chr((plain%1000))
        plain = plain/1000
    plaintext = plaintext[::-1]
    return plaintext


if __name__ == '__main__':
    print ("RSA - Encryption,Decryption\n")

    p = int(input("Enter a prime number (p): "))
    q = int(input("Enter another prime number (q): "))

    print ("\nGenerating your public and private keypairs now...")
    public, private = generate_keypair(p, q)
    print "Public key {e,n}: ", public, "\nPrivate key {d,n}: ", private

    e,n=public

    while(1):
        message = raw_input("\nEnter the message to encrypt with your public key: ")
        mess = ''
        for i in message:
            if (ord(i) < 100):
                mess = mess + '0' + str(ord(i))
            else:
                mess = mess + str(ord(i))

        print mess

        if (n <= int(mess)):
            print "Size of m should be less than the size of n"
            continue

        encrypted_msg = encrypt(public, int(mess))

        print ("\nEncrypting message with public key...")
        print "Encrypted message: ", ''.join(map(lambda x: str(x), encrypted_msg))

        print ("\nDecrypting message with private key...")
        print "Decrypted message: ", decrypt(private, int(''.join(map(lambda x: str(x), encrypted_msg))))