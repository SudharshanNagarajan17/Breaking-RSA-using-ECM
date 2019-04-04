def encrypt(pk, plaintext):
    # Unpack the key into it's components
    key, n = pk
    # Generate the ciphertext based on the plaintext and key using a^b mod m
    cipher = [pow(plaintext,key,n)]
    return cipher

if __name__ == '__main__':
    print "\nEnter the public key {e,n} for encryption:\n"
    e = int(input("Enter e: "))
    n = int(input("Enter n: "))
    public = (e, n)
    while(1):
        message = raw_input("\nEnter the message (m): ")
        mess = ''
        for i in message:
            if (ord(i) < 100):
                mess = mess + '0' + str(ord(i))
            else:
                mess = mess + str(ord(i))

        if (n <= int(mess)):
            print "Size of m should be less than the size of n"
            continue

        encrypted_msg = encrypt(public, int(mess))
        print "Encrypted message: ", ''.join(map(lambda x: str(x), encrypted_msg))