from ECMPrimeFactorization import ecc
import multiprocessing

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
    else:
        return -1


if __name__ == '__main__':
    print "\nEnter the public key {e,n}:\n"
    e = int(input("Enter e: "))
    n = int(input("Enter n: "))

    num = str(n)
    e1 = ecc()

    q1=multiprocessing.Queue()

    p1 = multiprocessing.Process(target=e1.compute, args=(num,1.0,q1))
    p2 = multiprocessing.Process(target=e1.compute, args=(num,2.0,q1))
    p3 = multiprocessing.Process(target=e1.compute, args=(num,3.0,q1))

    p1.start()
    p2.start()
    p3.start()

    while (p1.is_alive() and p2.is_alive() and p3.is_alive()):
        a = 1

    p1.terminate()
    p2.terminate()
    p3.terminate()

    p=q1.get()
    q=q1.get()

    if(p==-1 and q==-1):
        print "Invalid n : n should be a product of two primes"
        exit()
    print "\np = ",p
    print "q = ",q

    phi = (p-1)*(q-1)
    print "\nPhi = ",phi

    d=multiplicative_inverse(e,phi)
    if(d==-1):
        print "Invalid e: e should be relatively prime to phi(n)"
        exit()
    print "\nd = ",d

    print "\nPublic Key {e,n} = {",e,",",n,"}"
    print "Private Key {d,n} = {",d,",",n,"}"