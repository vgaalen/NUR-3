import numpy as np

def n(x, a, b, c, Nsat, A):
    return A * Nsat * (x/b)**(a-3) * np.exp(-(x/b)**c)

def N(x, a, b, c, Nsat, A):
    return 4*np.pi*x**2*n(x)

def x2n(x, a, b, c, Nsat, A):
    #print("x2n")
    return A * Nsat * x**(a-1) * b**(3-a) * np.exp(-(x/b)**c)

def Ndx(x ,xmin,xmax,Nsat,A=None):
    a,b,c = x[0],x[1],x[2]
    if A is None:
        A = 1/Romberg_integration(lambda x: x2n(x, a, b, c, 1, 1),0,2.5)
    return 4*np.pi*Romberg_integration(lambda x: x2n(x, a, b, c, Nsat, A),xmin,xmax)

def Romberg_integration(f,a,b,*args,m=5,h=1):
    """
    Implementation of Romberg Integration.
    Numerical integration of function f(x), from x=a until x=b
    """
    if a>b:
        a,b=b,a
    
    h=h*(b-a)
    r = np.zeros(m)
    r[0] = 0.5*h*(f(a,*args)+f(b,*args))

    N_p = 1
    for i in range(1,m):
        delta = h
        h = 0.5*h
        x = a + h
        for j in range(N_p):
            r[i] += f(x,*args)
            x += delta
        r[i] = 0.5*(r[i-1]+delta*r[i])
        N_p *= 2
    
    N_p = 1
    for i in range(1,m):
        N_p *= 4
        for j in range(m-i):
            r[j] = (N_p*r[j+1]-r[j])/(N_p-1)
    return r[0]

def golden_section_search(f, x0, x1, x2, target_accuracy=None, max_itt=1000):

    # define the golden ratio
    w=1/(1+((1+np.sqrt(5))*0.5))

    if target_accuracy is None:
        target_accuracy = np.abs(1e-8*x1)

    # make sure the x-values are sorted
    if x0>x1 or x0>x2:
        (x0, x1, x2) = np.sort([x0,x1,x2])
    elif x2<x0 or x2<x1:
        (x0, x1, x2) = np.sort([x0,x1,x2])

    # find the initial y-values
    y0 = f(x0)
    y1 = f(x1)
    y2 = f(x2)

    # find the smallest interval (we will tighten this one)
    if np.abs(x0-x1) > np.abs(x2-x1):
        smallest = 1
    else:
        smallest = 2
    
    for itt in range(max_itt):
        if smallest==1:
            # define the point at the golden ratio
            x3 = x1 + (x0-x1)*w
            # if we tighten this interval the other will now be the smallest
            smallest = 2
            # update the edge-point
            y3 = f(x3)
            if y3<y1:
                x2,x1 = x1,x3
                y2,y1 = y1,y3
            else:
                x0 = x3
                y0 = y3
        else: # same approach if the second interval is the smallest
            x3 = x1 + (x2-x1)*w
            smallest = 1
            y3 = f(x3)
            if y3<y1:
                x0,x1 = x1,x3
                y0,y1 = y1,y3
            else:
                x2 = x3
                y2 = y3

        # if the target accuracy is reached, return the x-value associated with the smallest y-value
        if np.abs(x2-x0)<=target_accuracy:
            if y3 < y1:
                print('a')
                return x3
            else:
                print('b')
                return x1
    if y3 < y1:
        print('c')
        return x3
    else:
        print('d')
        return x1

def multid_quick_sort(f,x):
    """
    Quick sort algorithm

    Parameters
    ----------
    f : function
        The function generating the values to sort on
    x : list
        The list that should be sorted based on f(x)

    Returns
    -------
    list
        The sorted list
    """

    if x.shape[0] <= 1:
        return x
    else:
        set = []
        for i in range(np.array(x).shape[0]):
            set.append(f(x[i,:]))

        if set[len(set)//2]<set[0]:
            if set[-1]<=set[len(set)//2]:
                set[0], set[len(set)//2], set[-1] = set[-1], set[len(set)//2], set[0]
                x[0,:], x[len(set)//2,:], x[-1,:] = x[-1,:], x[len(set)//2,:], x[0,:]
            elif set[-1]<set[0]:
                set[0], set[len(set)//2], set[-1] = set[len(set)//2], set[-1], set[0]
                x[0,:], x[len(set)//2,:], x[-1,:] = x[len(set)//2,:], x[-1,:], x[0,:]
            else:
                set[0], set[len(set)//2], set[-1] = set[len(set)//2], set[0], set[-1]
                x[0,:], x[len(set)//2,:], x[-1,:] = x[len(set)//2,:], x[0,:], x[-1,:]
        else:
            if set[-1]<=set[0]:
                set[0], set[len(set)//2], set[-1] = set[-1], set[0], set[len(set)//2]
                x[0,:], x[len(set)//2,:], x[-1,:] = x[-1,:], x[0,:], x[len(set)//2,:]
            elif set[-1]<set[len(set)//2]:
                set[0], set[len(set)//2], set[-1] = set[0], set[-1], set[len(set)//2]
                x[0,:], x[len(set)//2,:], x[-1,:] = x[0,:], x[-1,:], x[len(set)//2,:]
        
        i=0
        j=len(set)-1
        pivot = len(set)//2
        while i<j:
            if set[i]>=set[pivot]:
                if set[j]<=set[pivot]:
                    set[i], set[j] = set[j], set[i]
                    x[i,:], x[j,:] = x[j,:], x[i,:]
                    if i==pivot:
                        pivot=j
                    elif j==pivot:
                        pivot=i
                    i+=1
                    j-=1
                else:
                    j-=1
            else:
                i+=1
                if set[j]<=set[pivot]:
                    j-=1
        return np.concatenate((multid_quick_sort(f,x[:pivot,:]), [x[pivot,:]], multid_quick_sort(f,x[pivot+1:,:])), axis=0)

def DownhillSimplex(f,x,target_accuracy, num_itt=1000):
    if x.shape[0]!=x.shape[1]+1:
        raise ValueError("The input array x should have shape (n+1,n)")
    for itt in range(num_itt):
        # 1. order the points and calculate the mean (excluding the worst point)
        x = multid_quick_sort(f,x)
        mean = np.mean(x[:-1],axis=0)

        # 2. check if the fractional range in f(x) (no in x!), this is |f(x_N)-f(x_0)|/[0.5*|f(x_N)+f(x_0)|], is within target accuracy and if so terminate
        if np.abs(f(x[-1])-f(x[0])) / (0.5*np.abs(f(x[-1])+f(x[0]))) < target_accuracy:
            return x[0]
        
        # 3. propose a new point by reflecting x_N:x_try=2mean-x_N
        x_try = 2*mean - x[-1]
        if f(x[0])<=f(x_try)<f(x[-1]):
            x[-1] = x_try
        elif f(x_try)<f(x[0]):
            x_exp = 2*x_try-mean
            if f(x_exp)<f(x_try):
                x[-1] = x_exp
            else:
                x[-1] = x_try
        elif f(x_try)>=f(x[-1]):
            x_try = 0.5*(mean+x[-1])
            if f(x_try)<f(x[-1]):
                x[-1] = x_try
        else:
            x[1:] = 0.5*(x[0]+x[1:])
    return x[0]

def chi2(x, data, xmin, xmax, Nsat):
    chi2 = 0
    for i in range(data.size):
        val = Ndx(x,xmin,xmax,Nsat)
        if np.isnan(x).any():
            return np.inf
        chi2 += (data[i]-val)**2
    return chi2



if __name__=='__main__':
    ########
    ## 1a ##
    ########

    a = 2.4
    b = 0.25
    c = 1.6
    x_max = 5
    Nsat = 100
    A = 256/(5*np.pi**(3/2))

    # Find the maximum of N(x)
    #out = downhill_simplex(lambda x: -1*N(x), np.array([0,2.5,5]),0.001)
    out = golden_section_search(lambda x: -1*N(x, a, b, c, Nsat), 1e-32, 2.5, 5, target_accuracy=0.001)
    print(out, N(out))
    with open('output/ex1a.txt','w') as f:
        f.write(f'The maximum of N(x) is at x = {out} with N(x) = {N(out)}')
    #out = downhill_simplex(lambda x: -1*N(x), np.array([[0,2.5,5],[0.1,3,1],[2,3,4],[4,3,2]]),0.001)
    #print(out, N(out))

    ########
    ## 1b ##
    ########
    def readfile(filename, num_bins): # Courtesy of Marcel van Daalen
        f = open(filename, 'r')
        data = f.readlines()[3:] #Skip first 3 lines 
        nhalo = int(data[0]) #number of halos
        radius = []
        for i,line in enumerate(data[1:]):
            if line[:-1]!='#':
                radius.append(float(line.split()[0]))

        radius = np.array(radius, dtype=float)    
        f.close()
        return radius, nhalo #Return the virial radius for all the satellites in the file, and the number of halos

    num_bins = 75 #why?
    x_max = 2.5
    bin_edges = np.linspace(0,x_max,num_bins+1)
    import matplotlib.pyplot as plt

    mean_num_sat = np.zeros((5,num_bins))
    std_num_sat = np.zeros((5,num_bins))

    for i,file_name in enumerate(['satgals_m11.txt','satgals_m12.txt','satgals_m13.txt','satgals_m14.txt','satgals_m15.txt']):
        radius, nhalo = readfile('data/'+file_name, num_bins)
        hist, bins = np.histogram(radius, bins=num_bins)
        bin_centers = 0.5*(bins[1:]+bins[:-1])
        
        mean_num_sat[i] = hist/nhalo #hist/nhalo
        #std_num_sat[i] = np.std(satellites_per_halo, axis=0)

    x_init = np.array([[2.4,0.25,1.6],
                    [2.4,0.25,1.6],
                    [2.4,0.25,1.6],
                    [2.4,0.25,1.6]])

    best_val = []
    for i in range(5):
        x_init2 = x_init + 0.5*np.random.rand(4,3)
        out = DownhillSimplex(lambda x: chi2(x,mean_num_sat[i],0,2.5,np.sum(mean_num_sat[i])),x_init2,1e-30,num_itt=100)
        print("[a,b,c] = ",out, Ndx(out,0,2.5,100), chi2(out,mean_num_sat[i],0,2.5,np.sum(mean_num_sat[i])))
        best_val.append(out)
        with open('output/ex1b.txt','w') as f:
            f.write(f'The best fit parameters for file {i} are: [a,b,c] = {out} with Nsat = {np.sum(mean_num_sat[i])}')

        plt.bar(bin_centers,mean_num_sat[i],width=2.5/75,alpha=0.5)
        print(mean_num_sat[i].shape)
        print(bin_centers.shape)
        vals = []
        for j in range(num_bins):
            vals.append(Ndx(out,bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])))
        plt.plot(bin_centers,vals,color='red')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('plots/1b_'+str(i)+'.png')
        plt.close()


    ########
    ## 1c ##
    ########
    def log_likelihood(f,x,y):
        return np.sum( y*np.log(f(x)) - f(x) )

    best_val2 = []
    for i in range(5):
        x_init2 = x_init + 0.5*np.random.rand(4,3)
        out = DownhillSimplex(lambda x: -1*log_likelihood(lambda x: Ndx(x,0,2.5,np.sum(mean_num_sat[i])),x,mean_num_sat[i]),x_init2,1e-30,num_itt=100)
        print("[a,b,c] = ",out, Ndx(out,0,2.5,100), log_likelihood(lambda x: Ndx(x,0,2.5,np.sum(mean_num_sat[i])),out,mean_num_sat[i]))
        best_val2.append(out)
        with open('output/ex1c.txt','w') as f:
            f.write(f'The best fit parameters for file {i} are: [a,b,c] = {out} with Nsat = {np.sum(mean_num_sat[i])}')

        plt.bar(bin_centers,mean_num_sat[i],width=2.5/75,alpha=0.5)
        print(mean_num_sat[i].shape)
        print(bin_centers.shape)
        vals = []
        for j in range(num_bins):
            vals.append(Ndx(out,bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])))
        plt.plot(bin_centers,vals,color='red')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('plots/1c_'+str(i)+'.png')
        plt.close()

    ########
    ## 1d ##
    ########
    def Gtest(obs,exp):
        mask = obs>0
        #print(obs,exp)
        return 2*np.sum(obs[mask]*np.log(obs[mask]/exp[mask]))

    with open('output/ex1d.txt','w') as f:
        for i in range(5):
            print("G-test for file ",i)
            vals = []
            vals2= []
            for j in range(num_bins):
                vals.append(Ndx(best_val[i],bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])))
                vals2.append(Ndx(best_val2[i],bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])))
            print(Gtest(mean_num_sat[i],np.array(vals)))
            f.write(f'G-test for file {i} with best fit parameters from b: {Gtest(mean_num_sat[i],np.array(vals))}\n')
            print(Gtest(mean_num_sat[i],np.array(vals2)))
            f.write(f'G-test for file {i} with best fit parameters from c: {Gtest(mean_num_sat[i],np.array(vals2))}\n')