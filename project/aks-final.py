
from sympy import perfect_power
from sympy.ntheory import isprime
import math
import mpmath
import sympy
from mpmath import mp 
import numpy as np
test = np.zeros((10000, 10000), dtype=np.int)
import time
import copy
#from sympy import fft, ifft
from sympy import convolution
from sympy import re
from sympy.discrete.convolutions import convolution_fft
from sympy.discrete.convolutions import convolution_ntt
import copy
import cProfile
import re
from mpmath import sin, cos
from random import randint
from math import floor,log,sqrt,ceil,factorial
# f
mp.dps = 100
def fft( b,  invert=False) :
    a=[mp.mpc(0)]*len(b);
    a[0:len(b)]=b;
    n=len(a);
    if (n == 1):
        return a;
    a0=[mp.mpc(0)]*(n//2);
    a1=[mp.mpc(0)]*(n//2);
    for i in range(n//2) :
        a0[i] = mp.mpc(a[2*i]);
        a1[i] = mp.mpc(a[2*i+1]);
    a0=fft(a0, invert);
    a1=fft(a1, invert);
    factor=[-1,1];
    ang = (mp.mpf(2) * mp.pi / n) * factor[int(invert)];
    w=mp.mpc(1);
    wn=mp.exp(mp.mpc(real=0,imag=ang));
    
    for i in range(n//2) :
        a[i] = a0[i] + w * a1[i];
        a[i + n//2] = a0[i] - w * a1[i];
        if (invert) :
            a[i] /= mp.mpc(2);
            a[i + n//2] /= mp.mpc(2);
        w *= wn;

    
    return a;
    
def C(n,r):
    if test[n][r]!=0:
        return test[n][r]
    if n == r : 
        return 1
    elif r==1:
        return n
    else: 
        test[n][r]=C(n-1, r) + C(n-1, r-1)
        return test[n][r]
def C1(n,r):
    return factorial(n)//(factorial(r)*factorial(n-r))
fac=np.zeros(100000000)
def C2(n,r):
    if fac[n]==0:
        fac[n]=factorial(n)
    if fac[r]==0:
        fac[r]=factorial(r)
    if fac[n-r]==0:
        fac[n-r]=factorial(n-r)
    return fac[n]//(fac[n-r]*fac[r])
def checkRCandidate(r,n):
    len=n.bit_length()+1
    for i in range(1,len):
        if pow(n,i,r)==1 :
            return False
    return True
def findR(n):# step 2
    r=2
    while r<ceil(sqrt(n)) :
        if n%r==0 :
            return False
        else:
            q=r-1
            rand=min(5,q)
            while rand>0:
                s=randint(1,q)
                if (C1(q+s-1,s)>=pow(n,floor(2*sqrt(r))))&((pow(n,(r-1)/q)%r)>1):
                    break
                rand=rand-1
        r=r+1
    return r

def isPerfectPower(n):
    if perfect_power(n)==False :#返回整数(b，e)的元组，使得b^e==n。
        return False
    return True

def polyMul(p1,p2,n):
    r=len(p1);
    result = [0]*r;
    for i in range(r):
        for j in range(r):
            result[(i+j)%r]+=p1[i]*p2[j];
            result[(i+j)%r]%=n;
    return result;



def polyMulFast(p1,p2,n):
    r=len(p1);
    m=2**(r.bit_length());
    p1F=[0]*m;
    p2F=[0]*m;
    p1F[0:len(p1)]=[mp.mpc(x) for x in p1];
    p2F[0:len(p1)]=[mp.mpc(x) for x in p2];
    p1F=fft(p1F);
    p2F=fft(p2F);
    result=[0]*m;
    for i in range(m):
        result[i]=p1F[i]*p2F[i];
    result=fft(result,True);
    res=[0]*r;
    for i in range(len(result)):
        res[i%r]+=int(round(result[i].real));
       # res[i%r]%=n;
    #res[i%len(p1)]%=n;
    return res;



def polypow(a,n,r,m,fast=False): # calculates (x+a)**n %(n,x**r-1)
    x = [0]*r;
    result=[0]*r;
    x[0]=a;    
    x[1]=1;
    result[0]=1;
    while n > 0 :
        if n%2== 1 :
            if fast:
                result=polyMulFast(result,x,m);
            else :
                result=polyMul(result,x,m);
        n = n//2;
        if fast:  
            x = polyMulFast(x,x,m);
        else:
            x = polyMul(x,x,m);
    return result; 

def checkWitness(a,n,r,fast=False): 
    LHS=polypow(a,n,r,n,fast);
    RHS = [0]*r;
    RHS[0]=a;
    RHS[n%r]=1;
    for i in range(r):
        if LHS[i]!=RHS[i] :
             return True;
    return False;


def aks(n,fast=False):

    if n%2==0|n%3==0|n%5==0|n==1:#简单的情况下，直接判断
        return False
    if( isPerfectPower(n)) :# step 1再判断是否有能整除的
        return False

    r=findR(n)
    if r==False :
        return False
    if r>=ceil(sqrt(n)) :# step 2
        return True
    #print("n=",n,"r=",r)
    end=min(math.ceil(math.sqrt(r))*(n.bit_length()+1),n)
    rand=5
    while rand>0:
        a=randint(1,end-1)
        if math.gcd(a,n)>1 :
            return False
        if checkWitness(a,n,r,fast) :
            return False
        rand=rand-1
    return True



n=10000

end=int(input("输入上限："))
start_time = time.time()
for i in range(n,n+end):
    #if aks(i)!=isprime(i):
    #    print("False")
    #    break
    if aks(i)==True:
        print(i)

print("--- %s seconds ---" % (time.time() - start_time))

