import multiprocessing

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
    p1 = multiprocessing.Process(target=fun, args=(0,pow(10,8)/3))
    p2 = multiprocessing.Process(target=fun, args=(pow(10,8)/3,pow(10,8)*2/3))
    p3 = multiprocessing.Process(target=fun, args=(pow(10,8)*2/3,pow(10,8)))

    p1.start()
    p2.start()
    p3.start()