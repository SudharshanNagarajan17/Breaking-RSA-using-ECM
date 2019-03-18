import threading

#print len(str(245841236512478852752909734912575581815967630033049838269083))
def fun(q,p):
    i = q
    k = pow(10,7)
    while i < p:
        i=i+1
        if(i%k==0):
            print i

if __name__ == "__main__":
    '''fun(0,pow(10,8))
    '''
    t1 = threading.Thread(target=fun, args=(0,pow(10,8)/2), name='t1')
    t2 = threading.Thread(target=fun, args=(pow(10,8)/2,pow(10,8)),name='t2')

    t1.start()
    t2.start()