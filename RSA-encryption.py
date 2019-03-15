def encrypt(pk, plaintext):
    # Unpack the key into it's components
    key, n = pk
    # Generate the ciphertext based on the plaintext and key using a^b mod m
    cipher = pow(plaintext,key,n)
    return cipher

if __name__ == '__main__':
    print "\nEnter the public key {e,n} for encryption:\n"
    e = int(input("Enter e: "))
    n = int(input("Enter n: "))
    #e=58240451197
    #n=819708416353
    nlength = len(str(n))
    public = (e, n)
    while(1):
        message = raw_input("\nEnter the message (m): ")
        
        mess = ''
        for i in message:
            if( ord(i) < 100 ):
                mess = mess + '0'+str(ord(i))
            else:
                mess=mess+str(ord(i))
        print "ASCII value for encryption: "+mess
        
        if( n <= int(mess)):
            print "Size of m should be less than the size of n"
            continue
        
        encrypted_msg = encrypt(public, int(mess))
        print "Encrypted Message: "+str(encrypted_msg)
        