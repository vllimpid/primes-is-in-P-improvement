from sympy import perfect_power
from sympy.ntheory import isprime
import math
import mpmath
import sympy
from mpmath import mp 

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
from sys import argv
from random import randint
from math import floor,log

def compute_power(a,u,p):
    ans=1
    u=int(u)
    while u>0:
        if u%2==1 :
            ans=(ans*a)%p
        a=(a*a)%p
        u=u//2
    return ans 

def miller_rabin_witness(a, p):
    if p == 1:
        return False
    if p == 2:
        return True
 
    n = p - 1
    t = int(floor(log(n, 2)))
    u = 1
    while t > 0:
        u = n / 2**t
        if n % 2**t == 0 and u % 2 == 1:
            break
        t = t - 1
 
    b1 = b2 = compute_power(a, u, p)
    for i in range(1, t + 1):
        b2 = b1**2 % p
        if b2 == 1 and b1 != 1 and b1 != (p - 1):
            return False
        b1 = b2
    if b1 != 1:
        return False
 
    return True
 
def prime_test_miller_rabin(p, k):
    while k > 0:
        a = randint(1, p - 1)
        if not miller_rabin_witness(a, p):
            return False
        k = k - 1
    return True

base=10000000
n=1
#100000000 +38
end=5000
if __name__ == "__main__":
    start_time = time.time()
    for i in range(base+1,base+end):
        if prime_test_miller_rabin(i,5)==True:
           print(i) 
    print("--- %s seconds ---" % (time.time() - start_time))