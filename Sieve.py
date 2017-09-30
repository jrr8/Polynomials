from math import sqrt
from time import time, sleep
import datetime
#import comp182

def sieve(n):
	"""
	Creates a sieve of Erasthoces from 1 to n.
	Returns an ordered list of primes less than
	or equal to n.
	"""
	if n < 2 or n%1:
		return []

	sieve = {}
	for i in range(2, int(sqrt(n)) + 1):
		if sieve.get(i, True):
			for j in range(2, n/i + 1):
				sieve[i*j] = False

	primes = [i for i in range(2,n+1) if sieve.get(i, True)]
	return primes

def distribution(sieve):
	"""
	'sieve' is object returned by function sieve
	"""
	s = time()
	dist = {}
	for i in range(sieve[1]-1):
		gap = sieve[2][i+1]-sieve[2][i]
		dist[gap] = dist.get(gap,0) + 1
	return time()-s, dist

# a = sieve(50000000)  ## Checked at https://primes.utm.edu/lists/small/100000.txt
# dist_a = distribution(a)  ## 4.6 seconds for n = 5,000,000. 57.5 seconds for n = 50,000,000.
# print "Elapsed time:", a[0]+dist_a[0]
# print "Primes less than 50,000,000:", a[1]
# comp182.plot_dist_linear(dist_a[1], "Distribution of Prime Gaps", "Size of Gap", "Frequency", 'prime_dist_50mil')

def prime_factor(n):
	"""
	Returns the prime facorization of
	n as a list of tuples (p, r), where
	p is prime and p^r divides n, r being
	maximal
	This algorithm is inefficient and expects small inputs
	"""
	factors = {}
	while n-1:
		primes = sieve(int(sqrt(n)))
		_ = 1

		for prime in primes:

			if n % prime == 0:
				_ = 0

				j = 1
				n = n / prime
				while n % prime == 0:
					j += 1
					n = n / prime

				factors[prime] = j
				break

		if _:
			factors[n] = 1
			break


	factors = factors.items()
	factors.sort(key = lambda x: x[0])
	return factors

# def num_factors(n):
# 	"""
# 	Helper function to expedite counting number
# 	of prime factors for each number from 1 to n
# 	"""
# 	num_factors = {}
# 	checked = {}
# 	for i in range(2,n+1):

# 		if checked.get(i, False) == False:
# 			factors = prime_factor(i)
# 			num_factors[len(factors)] = num_factors.get(len(factors), 0) + 1
			
# 			if i > n/2:
# 				continue

# 			new_factors = range(2, n/i + 1)
# 			for factor in factors:
# 				for j in range()
# 				new_factors.remove(factor[0])

# 			for factor in new_factors:




# comp182.plot_dist_linear(dist_a[1], "Distribution of Prime Gaps", "Size of Gap", "Frequency", 'prime_dist_50mil')


a = b if b is not None else 1

