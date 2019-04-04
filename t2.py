import multiprocessing

def fun(q,p,que):
    i = q
    k = pow(10,7)
    while i < p:
        i=i+1
        '''if(i%k==0):
            print i'''
    que.put(q)
    que.put(100)

if __name__ == "__main__":
    '''fun(0,pow(10,8))
    '''
    q1=multiprocessing.Queue()

    p1 = multiprocessing.Process(target=fun, args=(0,pow(10,8)/3,q1))
    p2 = multiprocessing.Process(target=fun, args=(pow(10,8)/3,pow(10,8)*2/3,q1))
    p3 = multiprocessing.Process(target=fun, args=(pow(10,8)*2/3,pow(10,8),q1))

    p1.start()
    p2.start()
    p3.start()

    while(p1.is_alive() and p2.is_alive() and p3.is_alive()):
        a=1


    p1.terminate()
    p2.terminate()
    p3.terminate()

    print q1.get()
    print q1.get()