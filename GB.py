"""
This code sample is just a hobby project I've been working on
since learning about the computational tools afforded by 
Groebner Bases that we learned about in my Abstract Algebra II
class this fall
"""

from math import sqrt

## Helper functions for Poly class
## These functions don't check for proper inputs, but instead
## these checks are done in the Poly class methods

def xgcd(b, n):
	"""
	Extended Euclidean algorithm.
	Returns g, x, y such that g = gcd(b,n) and
	bx + ny = g
	"""
	x0, x1, y0, y1 = 1, 0, 0, 1
	while n != 0:
		q, b, n = b // n, n, b % n
		x0, x1 = x1, x0 - q * x1
		y0, y1 = y1, y0 - q * y1
	return  b, x0, y0

def mult_inv(b, n):
	"""
	Returns the multiplicative inverse of b
	mod n (assuming b is a unit mod n)
	"""
	if n is None:
		return b

	return xgcd(b,n)[1] % n

def trim(poly):
	"""
	Removes leading zeros from an iterable
	"""
	for i in range(len(poly)):
		if poly[i]:
			return list(poly[i:])

	return [0]

def is_prime(n):
	"""
	Primality test for n based on the sieve of Erasthoces.
	"""
	sieve = {2:True}
	for i in range(2, int(sqrt(n))+1):
		if sieve.get(i, True):
			for j in range(2, n/i+1):
				sieve[i*j] = False

	return sieve.get(n, True)

##

class Poly:
	"""
	Class for working with polynomials. Input is a list or tuple
	representing coefficients, with leading term at index 0 and
	constant term at end of iterable.

	Currently, errors are handled by returning the null element
	of whatever type the method returns, and an error message is printed.
	Consider modifying to have more robust error codes.
	"""
	
	def __init__(self, poly, mod=None):
		if not isinstance(poly,  (list,tuple)):
			print "Improper input. Unknown type " + str(type(poly))
			return

		try:
			mod = abs(int(mod))
			if not is_prime(mod):
				print "Error: Modulus must be prime integer or None"
				return
		except:
			if mod is not None:
				print "Error: Modulus must be prime integer or None"
				return

		self.poly = trim(poly)
		self.update(mod)

	def __str__(self):
		"""
		This has lots of lines of code mainly to conform to my unnecessaryily complicated
		ideas of how polynomials should be printed. Particularly:
			1) There should be spaces between monomials
			2) +- 1 is to be ommitted as a coefficient for x^n (n > 0)
			3) The leading term should be preceded by no spaces
		"""

		## Constant term
		if self.degree < 1:
			return str(self.poly[0])

		else:
			end_str = ''
			if self.poly[-1] < 0:
				end_str += ' - ' + str(abs(self.poly[-1]))
			if self.poly[-1] > 0:
				end_str += ' + ' + str(self.poly[-1])

		## Linear term
		mid_str = ''
		if self.poly[-2] < 0:

			if self.poly[-2] == -1:
				if self.degree == 1:
					return '-x' + end_str
				mid_str += ' - x'

			else:
				if self.degree == 1:
					return str(self.poly[0]) + 'x' + end_str
				mid_str += ' - ' + str(abs(self.poly[-2])) + 'x'

		if self.poly[-2] > 0:

			if self.poly[-2] == 1:
				if self.degree == 1:
					return 'x' + end_str
				mid_str += ' + x'

			else:
				if self.degree == 1:
					return str(self.poly[0]) + 'x' + end_str
				mid_str += ' + ' + str(abs(self.poly[-2])) + 'x'


		## Higher order terms
		beg_str = ''
		if self.poly[0] == 1:
			beg_str += 'x^'+str(self.degree)

		elif self.poly[0] == -1:
			beg_str += '-x^'+str(self.degree)

		else:
			beg_str += str(self.poly[0])+"x^"+str(self.degree)

		for i in range(1, self.degree - 1):
			if self.poly[i] == 0:
				continue

			if self.poly[i] == 1:
				beg_str += " + x^"+str(self.degree-i)

			elif self.poly[i] < 0:
				if self.poly[i] == -1:
					beg_str += " - x^"+str(self.degree-i)
				else:
					beg_str += " - "+str(abs(self.poly[i]))+'x^'+str(self.degree-i)

			else:
				beg_str += " + "+str(self.poly[i])+'x^'+str(self.degree-i)

		return beg_str + mid_str + end_str

	def update(self, mod=None):
		try:
			mod = abs(int(mod))
			if not is_prime(mod):
				print "Error: Modulus must be prime integer or None"
				return
		except:
			if mod is not None:
				print "Error: Modulus must be prime integer or None"
				return

		if mod is not None:
			self.poly = [i % mod for i in self.poly]

		self.poly = trim(self.poly)
		self.set_degree()

	def set_degree(self):
		"""
		Sets degree attribute apporpriately. This method
		is called by __init__ and update
		"""
		if len(self.poly) > 1:
			self.degree = len(self.poly) - 1

		elif self.poly[0] == 0:
			self.degree = float('-inf')

		else:
			self.degree = 0

	def get_monomial(self, deg):
		"""
		Returns the deg degree term of self, 
		as a Poly object.
		"""
		try:
			deg = int(deg)
		except:
			if deg != float('-inf'):
				print "Degree must be non-negative integer"
				return Poly([])

		if 0 <= deg <= self.degree:
			return Poly([self.poly[self.degree - deg]] + [0]*deg)

		return Poly([])

	def leading_term(self):
		"""
		Wrapper for self.get_monomial(self.degree)
		"""
		return self.get_monomial(self.degree)

	def is_monic(self, mod=None):
		"""
		Checks if leading coefficient is a unit
		"""
		try:
			mod = abs(int(mod))
			if not is_prime(mod):
				print "Error: Modulus must be prime integer or None"
				return
		except:
			if mod is not None:
				print "Error: Modulus must be prime integer or None"
				return

		if mod is not None: ## Check if leading term is a unit for given modulus
			return xgcd(self.poly[0], mod)[0] == 1

		return abs(self.poly[0]) == 1

	def const_factor(self):
		"""
		Tries to factor out a constant
		if self is not monic. This is solely
		used as a bonus feature of div to try to
		divide by polynomials that are constant
		multiples of a monic polynomial.
		"""
		if not any([i % self.poly[0] for i in self.poly]):
			return self.poly[0], Poly([i / self.poly[0] for i in self.poly])

		return 1, self

	def add(self, poly2, mod=None):
		"""
		Modifies self by adding poly2, using modulus mod if supplied.
		Returns none.
		"""
		## Check that poly2 is a Poly object
		if not isinstance(poly2, Poly):
			print "Must pass Polynomial object"
			return Poly([])

		## Make a copy of poly2 to leave poly2 unchanged
		poly3 = Poly(poly2.poly, mod)
		self.update(mod)

		## If adding zero, do nothing
		if poly3.degree < 0:
			return

		## If adding to zero, the result is poly2
		if self.degree < 0:
			self.poly = poly3.poly
			self.update(mod)
			return

		diff = self.degree - poly3.degree
		a = [0]*(max(self.degree, poly3.degree)+1)

		if diff < 0:
			diff = abs(diff)
			a[:diff] = poly3.poly[:diff]
			for i in range(self.degree + 1):
				a[i + diff] = poly3.poly[i + diff] + self.poly[i]

		elif diff > 0:
			a[:diff] = self.poly[:diff]
			for i in range(poly3.degree + 1):
				a[i + diff] = self.poly[i + diff] + poly3.poly[i]

		else:
			for i in range(self.degree + 1):
				a[i] = self.poly[i] + poly3.poly[i]

		self.poly = a
		self.update(mod)
		return

	def subtract(self, poly2, mod=None):
		"""
		A wrapper for addition
		"""
		self.add(Poly([-1*i for i in poly2.poly]), mod)
		return

	def mult(self, poly2, mod=None):
		"""
		Updates self to be product of self and poly2,
		using modulus mod if supplied
		"""
		## Check that poly2 is a Poly object
		if not isinstance(poly2, Poly):
			print "Must pass Polynomial object"
			return

		## Make a copy of poly2 to leave poly2 unchanged
		poly3 = Poly(poly2.poly, mod)
		self.update(mod)

		## If either factor is the zero polynomial, multiplication is easy
		if min(poly3.degree, self.degree) < 0:
			self.poly = []
			self.update(mod)
			return

		## Multiply by each monomial of poly2 then add them all together
		a = []
		for i in range(poly3.degree + 1):
			if poly3.poly[i] == 0:
				continue

			b = [poly3.poly[i]*j for j in self.poly] + [0]*(poly3.degree - i)
			a.append(Poly(b))

		## Now add up all of the monomial partial sums
		self.poly = a.pop().poly
		self.update(mod)
		for i in a:
			self.add(i, mod)

		return

	def div(self, poly2, mod=None):
		"""
		Return quotient and remainder of dividing by poly2.
		Does not modify self
		"""
		if not isinstance(poly2, Poly):
			print "Must pass Polynomial object"
			return Poly([]), Poly([])

		## Copying self and poly2 to leave them unchanged
		self2, poly3 = Poly(self.poly, mod), Poly(poly2.poly, mod)

		if poly3.degree < 0:
			print "Cannot divide by 0 polynomial"
			return Poly([]), Poly([])

		## Try to factor out any constants from the divisor
		const, poly3 = poly3.const_factor()
		poly3.update(mod)

		if not poly3.is_monic(mod):
			print "Divisor must be monic"
			return Poly([]), Poly([])

		remainder = Poly(self2.poly)
		quotient = Poly([])

		while remainder.degree >= poly3.degree:
			lead_coeff = remainder.poly[0]*mult_inv(poly3.poly[0], mod)

			## Check to see that factoring by constant yields quotient in Z[x]
			if lead_coeff % const:
				print "Divisor must be monic"
				return Poly([]), Poly([])

			new_term = Poly([lead_coeff] + [0]*(remainder.degree-poly3.degree))
			quotient.add(new_term, mod)
			new_term.mult(poly3, mod)
			remainder.subtract(new_term, mod)

		## Factoring constant factor of divisor back into quotient
		quotient.poly = [i / const for i in quotient.poly]
		quotient.update(mod)

		return quotient, remainder

	def quotient(self, poly2, mod=None):
		"""
		Updates self to the quotient of division by poly2
		"""
		self.poly = self.div(poly2, mod)[0]
		self.update(mod)
		return

	def remainder(self, poly2, mod=None):
		"""
		Updates self to the remainder of division by poly2
		"""
		self.poly = self.div(poly2, mod)[1]
		self.update(mod)
		return



## For Groebner bases, will need to generalize to multivariate polynomials. This will
## require implementing a monomial ordering and general polynomial division, among others
			



