import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import warnings

def n(x, a, b, c, Nsat, A):
    return A * Nsat * (x/b)**(a-3) * np.exp(-(x/b)**c)

def N(x, a, b, c, Nsat, A):
    return 4*np.pi*x**2*n(x, a, b, c, Nsat, A)

def x2n(x, a, b, c, Nsat, A):
    if a<1:
        raise ValueError("A value for parameter a below 1 is not allowed")
    if b<0:
        raise ValueError("Negative values for b are not supported")
    if c<0:
        raise ValueError("Negative values for c are not supported")
    return A * Nsat * x**(a-1) * b**(3-a) * np.exp(-(x/b)**c)

def Ndx(parameters, x, A=None):
    a,b,c = parameters[0], parameters[1], parameters[2]
    xmin, xmax, Nsat = x
    if A is None:
        norm = 4*np.pi*Romberg_integration(lambda x: x2n(x, a, b, c, Nsat, 1),0,2.5)
        if norm!=0:
            A = Nsat / norm
        else:
            A = 1
    return 4*np.pi*Romberg_integration(lambda x: x2n(x, a, b, c, Nsat, A),xmin,xmax)

def Romberg_integration(f,a,b,*args,m=5,h=None):
    """
    Implementation of Romberg Integration.
    Numerical integration of function f(x), from x=a until x=b
    """
    if a>b:
        a,b=b,a
    
    if h is None:
        h=(b-a)
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
    """
    Implementation of Golden Section Search.
    This function finds a minimum in f using the bracket [x0,x1,x2]
    """

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
                return x3
            else:
                return x1
    if y3 < y1:
        return x3
    else:
        return x1

def multid_quick_sort(y,x):
    """
    Quick sort algorithm

    Parameters
    ----------
    y : list
        The values to sort with
    x : 2d numpy array
        The list that should be sorted based on the values in y
    """

    if x.shape[0] <= 1:
        return y, x
    else:
        # sort the first, middle, and last values to make sure the pivot is the maximum or minimum
        if y[len(y)//2]<y[0]:
            if y[-1]<=y[len(y)//2]:
                y[0], y[-1] = y[-1], y[0]
                x[0,:], x[-1,:] = x[-1,:].copy(), x[0,:].copy()
            elif y[-1]<y[0]:
                y[0], y[len(y)//2], y[-1] = y[len(y)//2], y[-1], y[0]
                x[0,:], x[len(y)//2,:], x[-1,:] = x[len(y)//2,:].copy(), x[-1,:].copy(), x[0,:].copy()
            else:
                y[0], y[len(y)//2] = y[len(y)//2], y[0]
                x[0,:], x[len(y)//2,:] = x[len(y)//2,:].copy(), x[0,:].copy()
        else:
            if y[-1]<=y[0]:
                y[0], y[len(y)//2], y[-1] = y[-1], y[0], y[len(y)//2]
                x[0,:], x[len(y)//2,:], x[-1,:] = x[-1,:].copy(), x[0,:].copy(), x[len(y)//2,:].copy()
            elif y[-1]<y[len(y)//2]:
                y[len(y)//2], y[-1] = y[-1], y[len(y)//2]
                x[len(y)//2,:], x[-1,:] = x[-1,:], x[len(y)//2,:].copy()

        i=0
        j=len(y)-1
        pivot = len(y)//2
        # sort the arrays relative to the pivot
        while i<j:
            if y[i]>=y[pivot]:
                if y[j]<=y[pivot]:
                    y[i], y[j] = y[j], y[i]
                    x[i,:], x[j,:] = x[j,:], x[i,:].copy()
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
                if y[j]<=y[pivot]:
                    j-=1
        
        # the pivot is in the right location now, sort the rest of the array
        start_y, start_x = multid_quick_sort(y[:pivot],x[:pivot,:])
        end_y, end_x = multid_quick_sort(y[pivot+1:],x[pivot+1:,:])
        return np.concatenate((start_y, [y[pivot]], end_y), axis=0), np.concatenate((start_x, [x[pivot,:]], end_x), axis=0)

def DownhillSimplex(f,x,target_accuracy, num_itt=1000):
    """
    Implementation of the Downhill Simplex minimalisation algorithm.
    """
    if x.shape[0]!=x.shape[1]+1:
        raise ValueError("The input array x should have shape (n+1,n)")
    
    y = [0]*x.shape[0]
    for i in range(len(y)):
        y[i] = f(x[i])
    
    for itt in range(num_itt):
        # 1. order the points and calculate the mean (excluding the worst point)
        y,x = multid_quick_sort(y,x)
        mean = np.mean(x[:-1],axis=0)

        # 2. check if the fractional range in f(x) (no in x!), this is |f(x_N)-f(x_0)|/[0.5*|f(x_N)+f(x_0)|], is within target accuracy and if so terminate
        if np.isnan(y[-1]) is False and (y[0] == y[1] or np.abs(y[-1] - y[0]) / (0.5*np.abs(y[-1] + y[0]))) < target_accuracy:
            return x[0]
        
        # 3. propose a new point by reflecting x_N:x_try=2mean-x_N
        x_try = 2*mean - x[-1]
        if f(x[0])<=f(x_try)<f(x[-1]):
            x[-1] = x_try
            y[-1] = f(x_try)
        elif f(x_try)<f(x[0]):
            x_exp = 2*x_try-mean
            if f(x_exp)<f(x_try):
                x[-1] = x_exp
                y[-1] = f(x_exp)
            else:
                x[-1] = x_try
                y[-1] = f(x_try)
        else:
            x_try = 0.5*(mean+x[-1])
            if f(x_try)<f(x[-1]):
                x[-1] = x_try
                y[-1] = f(x_try)
            else:
                x[1:] = 0.5*(x[0]+x[1:])
                for i in range(1,len(y)):
                    y[i] = f(x[i])
    return x[0]

class RNG():
    def __init__(self, seed):
        self.seed = seed
        self.x = np.uint64(seed)
        self.a1 = np.uint64(21)
        self.a2 = np.uint64(35)
        self.a3 = np.uint64(4)

        self.a = np.uint64(1664525)
        self.c = np.uint64(1013904223)
        
    
    def __call__(self, shape=None):
        """
        Generate random numbers using XOR_shift and LCG between 0 and 1
        """
        if shape is None:
            return self._LCG(self._XOR_shift())/18446744073709551615
        else:
            return np.array([self._LCG(self._XOR_shift())/18446744073709551615 for _ in range(np.prod(shape))]).reshape(shape)
    
    def get_state(self):
        return {"seed": self.seed, "a1": self.a1, "a2": self.a2, "a3": self.a3, "a": self.a, "c": self.c}

    def _XOR_shift(self):
        """
        Implementation of XOR_shift Random Number Generator
        """
        x = self.x
        x = x ^ ( x >> self.a1 )
        x = x ^ ( x << self.a2 )
        x = x ^ ( x >> self.a3 )
        self.x = x
        return x
    
    def _LCG(self, x):
        """
        Implementation of Linear Congruential Generator
        """
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="overflow encountered in scalar multiply")
            return (self.a*x + self.c)

def TestDegeneracy(x):
    """Test whether a set of points is degenerate
    returns False if not degenerate"""
    for i in range(x.shape[0]):
        x_i = np.concatenate([x[:i], x[i+1:]], axis=0) - x[i]

        angles = np.zeros((x_i.shape[0],x.shape[1]-1))
        for j in range(x.shape[1]-1):
            for k in range(x_i.shape[0]):
                angles[k,j] = np.arctan(x_i[k,j]/x_i[k,j+1])
                for l in range(k):
                    if (angles[l]==angles[k]).all(): # 2 points lay on a line as seen from a 3rd point, so the set is degenerate
                        return True
    return False

def chi2(x, model, x_data, y_data, error=None):
    if x[0]<1 or x[1]<0 or x[2]<0: # check if the model is defined for these input values, otherwise return inf
        return np.inf

    if len(x_data.shape)==1:
        y_model = model(x, x_data)
    else:
        y_model = np.zeros(y_data.shape)
        for i in range(y_model.shape[0]):
            y_model[i] = model(x, x_data[i])
    if error is None:
        error = np.sqrt(y_data)
    mask = (error!=0) # mask locations where the error is zero, as these are (probably) unintentional and lead to division by zero.
    return np.sum(np.power((y_data[mask]-y_model[mask]),2)/np.power(error[mask],2))



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
    out = golden_section_search(lambda x: -1*N(x, a, b, c, Nsat, A), 1e-32, 2.5, 5, target_accuracy=0.001)
    with open('output/ex1a.txt','w') as f:
        f.write(f'The maximum of N(x) is at x = {out} with N(x) = {N(out, a, b, c, Nsat, A)}')

    ########
    ## 1b ##
    ########
    def readfile(filename, num_bins): # This function was supplied in the exercise
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

    num_bins = 20
    x_max = 2.5
    bin_edges = np.logspace(-1,np.log10(x_max),num_bins+1)

    mean_num_sat = np.zeros((5,num_bins))
    std_num_sat = np.zeros((5,num_bins))

    file_names = ['satgals_m11.txt','satgals_m12.txt','satgals_m13.txt','satgals_m14.txt','satgals_m15.txt']
    names = ['Mhalo = 1e11 Msun','Mhalo = 1e12 Msun','Mhalo = 1e13 Msun','Mhalo = 1e14 Msun','Mhalo = 1e15 Msun']
    for i,file in enumerate(file_names):
        radius, nhalo = readfile('data/'+file, num_bins)
        hist, bins = np.histogram(radius, bins=bin_edges)
        bin_centers = 0.5*(bins[1:]+bins[:-1])
        mean_num_sat[i] = hist/nhalo

    # Plot the dependance of the distribution on parameters a, b, and c in order to determine good initial values for the Downhill Simplex
    for a in np.linspace(1,10,10):
        vals = []
        for i in range(len(bin_centers)):
            vals.append(Ndx([a,0.25,1.6],np.array([bin_edges[i],bin_edges[i+1],100]).T))
        plt.plot(bin_centers, vals,label=a)
    plt.legend()
    plt.title('relation between Ndx and a')
    plt.savefig('plots/Ndx_a.pdf')
    plt.close()

    for b in np.linspace(0.15,1,10):
        vals = []
        for i in range(len(bin_centers)):
            vals.append(Ndx([2.4,b,1.6],np.array([bin_edges[i],bin_edges[i+1],100]).T))
        plt.plot(bin_centers, vals, label=b)
    plt.legend()
    plt.title('relation between Ndx and b')
    plt.savefig('plots/Ndx_b.pdf')
    plt.close()

    for c in np.linspace(0.5,5,10):
        vals = []
        for i in range(len(bin_centers)):
            vals.append(Ndx([2.4,0.25,c],np.array([bin_edges[i],bin_edges[i+1],100]).T))
        plt.plot(bin_centers, vals,label=c)
    plt.legend()
    plt.title('relation between Ndx and c')
    plt.savefig('plots/Ndx_c.pdf')
    plt.close()

    # Use random initial values within the range found using these plots
    rng = RNG(2390478)
    x_init = np.array(( (10-1)*rng(4)+1,
                        (0.75-0.15)*rng(4)+0.15,
                        (5-1)*rng(4)+1 )).T

    while TestDegeneracy(x_init):
        x_init = np.array(( (10-1)*rng(4)+1,
                        (0.75-0.15)*rng(4)+0.15,
                        (5-1)*rng(4)+1 )).T
    
    best_val = []
    for i in range(5):
        x_init2 = x_init.copy() # to make sure these values are not overwritten during the downhill simplex itterations
        out = DownhillSimplex(lambda x: chi2(x,Ndx,np.array((bin_edges[:-1],bin_edges[1:],[np.sum(mean_num_sat[i])]*len(bin_centers))).T,mean_num_sat[i]),x_init2,1e-20,num_itt=100)
        #print("[a,b,c] = ",out, chi2(out,Ndx,np.array((bin_edges[:-1],bin_edges[1:],[np.sum(mean_num_sat[i])]*len(bin_centers))).T,mean_num_sat[i]))
        best_val.append(out)
        with open('output/ex1b.txt','a') as f:
            f.write(f'The best fit parameters for {names[i]} are: [a,b,c] = {out} with Nsat = {np.sum(mean_num_sat[i])}\n')

        plt.bar(np.log10(bin_centers),mean_num_sat[i],width=(bin_edges[-1]-bin_edges[0])/(num_bins+1)-0.05,zorder=0)
        vals = []
        for j in range(num_bins):
            vals.append(Ndx(out,[bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])]))
        plt.plot(np.log10(bin_centers),vals,color='red',zorder=1)
        plt.xlabel('log10(x)')
        plt.ylabel(r'$\widetilde{N}_i$')
        plt.yscale('log')
        plt.savefig('plots/1b_'+str(i)+'.png')
        plt.close()

    ########
    ## 1c ##
    ########
    def log_likelihood(x, model, x_data, y_data, error=None):
        if x[0]<1 or x[1]<0 or x[2]<0:
            return np.inf
        if len(x_data.shape)==1:
            y_model = model(x, x_data)
        else:
            y_model = np.zeros(y_data.shape)
            for i in range(y_model.shape[0]):
                y_model[i] = model(x, x_data[i])
        mask = (y_model!=0)
        return -1*np.sum( y_data[mask]*np.log(y_model[mask]) - y_model[mask] )

    best_val2 = []
    for i in range(5):
        x_init2 = x_init.copy()
        out = DownhillSimplex(lambda x: log_likelihood(x,Ndx,np.array((bin_edges[:-1],bin_edges[1:],[np.sum(mean_num_sat[i])]*len(bin_centers))).T,mean_num_sat[i]),x_init2,1e-20,num_itt=100)
        #print("[a,b,c] = ",out, -1*log_likelihood(out,Ndx,np.array((bin_edges[:-1],bin_edges[1:],[np.sum(mean_num_sat[i])]*len(bin_centers))).T,mean_num_sat[i]))
        best_val2.append(out)
        with open('output/ex1c.txt','a') as f:
            f.write(f'The best fit parameters for {names[i]} are: [a,b,c] = {out} with Nsat = {np.sum(mean_num_sat[i])}\n')

        plt.bar(np.log10(bin_centers),mean_num_sat[i],width=(bin_edges[-1]-bin_edges[0])/(num_bins+1)-0.05,zorder=0)
        vals = []
        for j in range(num_bins):
            vals.append(Ndx(out,[bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])]))
        plt.plot(np.log10(bin_centers),vals,color='red',zorder=1)
        plt.xlabel('log10(x)')
        plt.ylabel(r'$\widetilde{N}_i$')
        plt.yscale('log')
        plt.savefig('plots/1c_'+str(i)+'.png')
        plt.close()
    
    ########
    ## 1d ##
    ########
    from scipy.special import gammainc

    def Gtest(obs,exp):
        mask = obs>0
        return 2*np.sum(obs[mask]*np.log(obs[mask]/exp[mask]))
    
    def Qvalue(x,k):
        return gammainc(k/2,x/2) # scipy.special's implementation of the incomplete gamma-function is regularized (has a factor 1/\Gamma)

    # here k=3, as we are fitting 3 parameters
    k=3

    gtest = np.zeros((2,5))
    for i in range(5):
        vals = []
        vals2 = []
        for j in range(num_bins):
            vals.append(Ndx(best_val[i],[bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])]))
            vals2.append(Ndx(best_val2[i],[bin_edges[j],bin_edges[j+1],np.sum(mean_num_sat[i])]))
        gtest[:,i] = [Gtest(mean_num_sat[i],np.array(vals)), Gtest(mean_num_sat[i],np.array(vals2))]
    
    qval = Qvalue(gtest, k)

    with open('output/ex1d.txt','w') as f:
        f.write('file, G Gaussian, G Poisson, Q-value Gaussian, Q-value Poisson\n')
        for i in range(5):
            f.write(f'{file_names[i]}, {gtest[0,i]}, {gtest[1,i]}, {qval[0,i]}, {qval[1,i]}\n')
