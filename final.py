from ECCPrimeFactorization import ecc


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


if __name__ == '__main__':
    print "\nEnter the public key {e,n}:\n"
    e = int(input("Enter e: "))
    n = int(input("Enter n: "))

    num = str(n)
    e1 = ecc()
    p,q=e1.compute(num)

    if(p==-1 and q==-1):
        print "Invalid n : n should be a product of two primes"
        exit()
    print "\np = ",p
    print "q = ",q

    phi = (p-1)*(q-1)
    print "\nPhi = ",phi

    d=multiplicative_inverse(e,phi)
    print "\nd = ",d

    print "\nPublic Key {e,n} = {",e,",",n,"}"
    print "Private Key {d,n} = {",d,",",n,"}"
