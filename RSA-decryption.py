def decrypt(pk, ciphertext):
    # Unpack the key into its components
    key, n = pk
    # Generate the plaintext based on the ciphertext and key using a^b mod m
    plain = pow(ciphertext, key, n)
    plaintext=''
    while(plain):
        plaintext = plaintext+chr((plain%100))
        plain = plain/100
    plaintext = plaintext[::-1]
    return plaintext

if __name__ == '__main__':
    print "\nEnter the private key {d,n} for decryption:\n"
    d = int(input("Enter d: "))
    n = int(input("Enter n: "))
    while(1):
        encrypted_msg = int(input("\nEnter the cipher text (c): "))
        private = (d, n)
        print "Decrypted message: ", decrypt(private, encrypted_msg)