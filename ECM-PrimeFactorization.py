import math
import sys
import random

try:
   import psyco
   psyco.full()
   PSYCO_EXISTS = True
except ImportError:
   PSYCO_EXISTS = False

try:
   from gmpy2 import isqrt as sqrt
   from gmpy2 import iroot as root
   from gmpy2 import gcd, invert, mpz, next_prime
   import gmpy2
   GMPY_EXISTS = True
except ImportError:
   try:
      from gmpy import gcd, invert, mpz, next_prime, sqrt, root
      GMPY_EXISTS = True
   except ImportError:
      GMPY_EXISTS = False

# Constants declarations
INV_C = 1.4
LOG_2 = math.log(2)
LOG_4 = math.log(4)
LOG_3_MINUS_LOG_LOG_2 = math.log(3) - math.log(LOG_2)
LOG_4_OVER_9 = LOG_4 / 9
_3_OVER_LOG_2 = 3 / LOG_2
_5_LOG_10 = 5 * math.log(10)
_7_OVER_LOG_2 = 7 / LOG_2
BIG = 2.0**512
BILLION = 10**9
MULT = math.log(3) / LOG_2
ONE = mpz(1)
SMALL = 2.0**(-30)
SMALLEST_COUNTEREXAMPLE_FASTPRIME = 2047
T = (type(mpz(1)), type(1), type(1L))
DUMMY = 'dummy'
_12_LOG_2_OVER_49 = 12 * math.log(2) / 49
RECORD = 1162795072109807846655696105569042240239

class ts:

   def __init__(self, degree, acc, p):
      self.acc = acc
      self.coefficients = p[:degree + 1]
      while len(self.coefficients) <= degree:
         self.coefficients.append(0)

   def add(self, a, b):
      '''Adds a and b'''
      b_ = b.coefficients[:]
      a_ = a.coefficients[:]
      self.coefficients = []

      while len(b_) > len(a_):
         a_.append(0)
      while len(b_) < len(a_):
         b_.append(0)

      for i in xrange(len(a_)):
         self.coefficients.append(a_[i] + b_[i])

      self.acc = a.acc

   def ev(self, x):
      '''Returns a(x)'''
      answer = 0
      for i in xrange(len(self.coefficients) - 1, -1, -1):
         answer *= x
         answer += self.coefficients[i]
      return answer

   def evh(self):
      '''Returns a(1/2)'''
      answer = 0
      for i in xrange(len(self.coefficients) - 1, -1, -1):
         answer >>= 1
         answer += self.coefficients[i]
      return answer

   def evmh(self):
      '''Returns a(-1/2)'''
      answer = 0
      for i in xrange(len(self.coefficients) - 1, -1, -1):
         answer = - answer >> 1
         answer += self.coefficients[i]
      return answer

   def int(self):
      '''Replaces a by an integral of a'''
      self.coefficients = [0] + self.coefficients
      for i in xrange(1, len(self.coefficients)):
         self.coefficients[i] /= i

   def lindiv(self, a):
      '''a.lindiv(k) -- sets a/(x-k/2) for integer k'''
      for i in xrange(len(self.coefficients) - 1):
         self.coefficients[i] <<= 1
         self.coefficients[i] /= a
         self.coefficients[i + 1] -= self.coefficients[i]
      self.coefficients[-1] <<= 1
      self.coefficients[-1] /= a

   def neg(self):
      '''Sets a to -a'''
      for i in xrange(len(self.coefficients)):
         self.coefficients[i] = - self.coefficients[i]

   def set(self, a):
      '''a.set(b) sets a to b'''
      self.coefficients = a.coefficients[:]
      self.acc = a.acc

   def simp(self):
      '''Turns a into a type of Taylor series that can be fed into ev, but cannot be computed with further.'''
      for i in xrange(len(self.coefficients)):
         shift = max(0, int(math.log(abs(self.coefficients[i]) + 1) / LOG_2) - 1000)
         self.coefficients[i] = float(self.coefficients[i] >> shift)
         shift = self.acc - shift
         for _ in xrange(shift >> 9):
            self.coefficients[i] /= BIG
         self.coefficients[i] /= 2.0**(shift & 511)
         if abs(self.coefficients[i] / self.coefficients[0]) <= SMALL:
            self.coefficients = self.coefficients[:i]
            break

def add(p1, p2,  n):
   inv = range(len(p1))

   for i in xrange(len(p1)):
      inv[i] = p1[i][0] - p2[i][0]

   inv = parallel_invert(inv, n)

   if not isinstance(inv, list):
      return inv

   for i in xrange(len(p1)):
      m = ((p1[i][1] - p2[i][1]) * inv[i]) % n
      p2[i][0] = (m * m - p1[i][0] - p2[i][0]) % n
      p2[i][1] = (m * (p1[i][0] - p2[i][0]) - p1[i][1]) % n

   return p2

def add_sub_x_only(p1, p2,  n):
   sums = range(len(p1))
   difs = range(len(p1))

   for i in xrange(len(p1)):
      sums[i] = p2[i][0] - p1[i][0]

   sums = parallel_invert(sums, n)

   if not isinstance(sums, list):
      return (sums, None)

   for i in xrange(len(p1)):
      ms = ((p2[i][1] - p1[i][1]) * sums[i]) % n
      md = ((p2[i][1] + p1[i][1]) * sums[i]) % n
      sums[i] = (ms * ms - p1[i][0] - p2[i][0]) % n
      difs[i] = (md * md - p1[i][0] - p2[i][0]) % n

   sums = tuple(sums)
   difs = tuple(difs)

   return (sums, difs)

def atdn(a, d, n):
   x = 1
   pos = int(math.log(d) / LOG_2)

   while pos >= 0:
      x = (x * x) % n
      if (d >> pos) & 1:
         x *= a
      pos -= 1

   return x % n

def copy(p):
   answer = []
   for i in p:
      answer.append(i[:])

   return answer

def could_be_prime(n):
   if n < 2:
      return False
   if n == 2:
      return True
   if not n & 1:
      return False

   product = ONE
   log_n = int(math.log(n)) + 1
   bound = int(math.log(n) / (LOG_2 * math.log(math.log(n))**2)) + 1
   if bound * log_n >= n:
      bound = 1
      log_n = int(sqrt(n))
   prime_bound = 0
   prime = 3

   for _ in xrange(bound):
      p = []
      prime_bound += log_n
      while prime <= prime_bound:
         p.append(prime)
         prime = next_prime(prime)
      if p != []:
         p = prod(p)
         product = (product * p) % n

   return gcd(n, product) == 1

def double(p, n):
   inv = range(len(p))

   for i in xrange(len(p)):
      inv[i] = p[i][1] << 1

   inv = parallel_invert(inv, n)

   if not isinstance(inv, list):
      return inv

   for i in xrange(len(p)):
      x = p[i][0]
      m = (x * x) % n
      m = ((m + m + m + p[i][2]) * inv[i]) % n
      p[i][0] = (m * m - x - x) % n
      p[i][1] = (m * (x - p[i][0]) - p[i][1]) % n

   return p

def fastprime(n):
   if not could_be_prime(n):
      return False
   if n == 2:
      return True

   j = 1
   d = n >> 1

   while not d & 1:
      d >>= 1
      j += 1

   p = 1
   pos = int(math.log(d) / LOG_2)

   while pos >= 0:
      p = (p * p) % n
      p <<= (d >> pos) & 1
      pos -= 1

   if p in (n - 1, n + 1):
      return True

   for _ in xrange(j):
      p = (p * p) % n

      if p == 1:
         return False
      elif p == n - 1:
         return True

   return False

def greatest_n(phi_max):
   phi_product = 1
   product = 1
   prime = 1
   while phi_product <= phi_max:
      prime = next_prime(prime)
      phi_product *= prime - 1
      product *= prime

   n_max = (phi_max * product) / phi_product

   phi_values = range(n_max)

   prime = 2
   while prime <= n_max:
      for i in xrange(0, n_max, prime):
         phi_values[i] -= phi_values[i] / prime

      prime = next_prime(prime)

   for i in xrange(n_max - 1, 0, -1):
      if phi_values[i] <= phi_max:
         return i

def inv_const(n):
   return int(INV_C * math.log(n)**0.42)

def naf(d):
   g = 0L
   while d:
      g <<= 2
      g ^= ((d & 2) & (d << 1)) ^ (d & 1)
      d += (d & 2) >> 1
      d >>= 1
   return g

def parallel_invert(l, n):
   l_ = l[:]
   for i in xrange(len(l)-1):
      l[i+1] = (l[i] * l[i+1]) % n

   try:
      inv = invert(l[-1], n)
   except ZeroDivisionError:
      inv = 0
   if inv == 0:
      return gcd(l[-1], n)

   for i in xrange(len(l)-1, 0, -1):
      l[i] = (inv * l[i-1]) % n
      inv = (inv * l_[i]) % n
   l[0] = inv

   return l

def prod(p):
   jump = 1
   while jump < len(p):
      for i in xrange(0, len(p) - jump, jump << 1):
         p[i] *= p[i + jump]
         p[i + jump] = None
      jump <<= 1
   return p[0]

def rho_ev(x, ts):
   return ts[int(x)].ev(x - int(x) - 0.5)

def rho_ts(n):

   f = ts(10, 10, [])
   answer = [ts(10, 10, [1])]
   for _ in xrange(n):
      answer.append(ts(10, 10, [1]))
   deg = 5
   acc = 50 + n * int(1 + math.log(1 + n) + math.log(math.log(3 + n)))
   r = 1
   rho_series = ts(1, 10, [0])
   while r != rho_series.coefficients[0]:
      deg = (deg + (deg << 2)) / 3
      r = rho_series.coefficients[0]
      rho_series = ts(deg, acc, [(1L) << acc])
      center = 0.5
      for i in xrange(1, n+1):
         f.set(rho_series)
         center += 1
         f.lindiv(int(2*center))
         f.int()
         f.neg()
         d = ts(deg, acc, [rho_series.evh() - f.evmh()])
         f.add(f, d)
         rho_series.set(f)
         f.simp()
         answer[i].set(f)
      rho_series.simp()

   return answer

def sub_sub_sure_factors(f, u, curve_parameter):

   while not (f & 1):
      yield 2
      f >>= 1

   while not (f % 3):
      yield 3
      f /= 3

   if isprime(f):
      yield f
      return

   log_u = math.log(u)
   u2 = int(_7_OVER_LOG_2 * u * log_u / math.log(log_u))
   primes = []
   still_a_chance = True
   log_mo = math.log(f + 1 + sqrt(f << 2))

   g = gcd(curve_parameter, f)
   if g not in (1, f):
      for factor in sub_sub_sure_factors(g, u, curve_parameter):
         yield factor
      for factor in sub_sub_sure_factors(f/g, u, curve_parameter):
         yield factor
      return

   g2 = gcd(curve_parameter**2 - 5, f)
   if g2 not in (1, f):
      for factor in sub_sub_sure_factors(g2, u, curve_parameter):
         yield factor
      for factor in sub_sub_sure_factors(f / g2, u, curve_parameter):
         yield factor
      return

   if f in (g, g2):
      yield f

   while still_a_chance:
      p1 = get_points([curve_parameter], f)
      for prime in primes:
         p1 = multiply(p1, prime, f)
         if not isinstance(p1, list):
            if p1 != f:
               for factor in sub_sub_sure_factors(p1, u, curve_parameter):
                  yield factor
               for factor in sub_sub_sure_factors(f/p1, u, curve_parameter):
                  yield factor
               return
            else:
               still_a_chance = False
               break

      if not still_a_chance:
         break

      prime = 1
      still_a_chance = False
      while prime < u2:
         prime = next_prime(prime)
         should_break = False
         for _ in xrange(int(log_mo / math.log(prime))):
            p1 = multiply(p1, prime, f)
            if not isinstance(p1, list):
               if p1 != f:
                  for factor in sub_sub_sure_factors(p1, u, curve_parameter):
                     yield factor
                  for factor in sub_sub_sure_factors(f/p1, u, curve_parameter):
                     yield factor
                  return

               else:
                  still_a_chance = True
                  primes.append(prime)
                  should_break = True
                  break
         if should_break:
            break

   for i in xrange(2, int(math.log(f) / LOG_2) + 2):
      r = root(f, i)
      if r[1]:
         for factor in sub_sub_sure_factors(r[0], u, curve_parameter):
            for _ in xrange(i):
               yield factor
         return

   a = 1 + sqrt(f)
   bsq = a * a - f
   iter = 0

   while bsq != sqrt(bsq)**2 and iter < 3:
      a += 1
      iter += 1
      bsq += a + a - 1

   if bsq == sqrt(bsq)**2:
      b = sqrt(bsq)
      for factor in sub_sub_sure_factors(a - b, u, curve_parameter):
         yield factor
      for factor in sub_sub_sure_factors(a + b, u, curve_parameter):
         yield factor
      return

   yield f
   return

def sub_sure_factors(f, u, curve_params):

   if len(curve_params) == 1:
      for factor in sub_sub_sure_factors(f, u, curve_params[0]):
         yield factor
      return

   c1 = curve_params[:len(curve_params) >> 1]
   c2 = curve_params[len(curve_params) >> 1:]

   if mainloop(f, u, c1) == 1:
      for factor in sub_sure_factors(f, u, c2):
         yield factor
      return

   if mainloop(f, u, c2) == 1:
      for factor in sub_sure_factors(f, u, c1):
         yield factor
      return

   for factor in sub_sure_factors(f, u, c1):
      if isprime(factor):
         yield factor
      else:
         for factor_of_factor in sub_sure_factors(factor, u, c2):
            yield factor_of_factor

   return

def subtract(p1, p2,  n):

   inv = range(len(p1))

   for i in xrange(len(p1)):
      inv[i] = p2[i][0] - p1[i][0]

   inv = parallel_invert(inv, n)

   if not isinstance(inv, list):
      return inv

   for i in xrange(len(p1)):
      m = ((p1[i][1] + p2[i][1]) * inv[i]) % n
      p2[i][0] = (m * m - p1[i][0] - p2[i][0]) % n
      p2[i][1] = (m * (p1[i][0] - p2[i][0]) + p1[i][1]) % n

   return p2

def sure_factors(n, u, curve_params, veb, ra, ov, tdb, pr):

   f = mainloop(n, u, curve_params)

   if f == 1:
      return

   if isprime(f):
      yield f
      n /= f
      if isprime(n):
         yield n
      return

   for factor in sub_sure_factors(f, u, curve_params):
      if isprime(factor):
         yield factor
      else:
         for factor_of_factor in ecm(factor, True, ov, veb, tdb, pr):
            yield factor_of_factor
      n /= factor

   if isprime(n):
      yield n

   return

def to_tuple(p):
   answer = []
   for i in p:
      answer.append((i[0], i[1]))

   return tuple(answer)

def mainloop(n, u, p1):
   k = inv_const(n)
   log_u = math.log(u)
   log_log_u = math.log(log_u)
   log_n = math.log(n)
   u2 = int(_7_OVER_LOG_2 * u * log_u / log_log_u)
   ncurves = len(p1)
   w = int(math.sqrt(_3_OVER_LOG_2 * ncurves / k) - 0.5)
   number_of_primes = int((ncurves << w) * math.sqrt(LOG_4_OVER_9 * log_n / k) / log_u) # Lagrange multipliers!
   number_of_primes = min(number_of_primes, int((log_n / math.log(log_n))**2 * ncurves / log_u), int(u / log_u))
   number_of_primes = max(number_of_primes, 1)
   m = math.log(number_of_primes) + log_log_u
   w = min(w, int((m - 2 * math.log(m) + LOG_3_MINUS_LOG_LOG_2) / LOG_2))
   w = max(w, 1)
   max_order = n + sqrt(n << 2) + 1 # By Hasse's theorem.
   det_bound = ((1 << w) - 1 + ((w & 1) << 1)) / 3
   log_mo = math.log(max_order)
   p = range(number_of_primes)
   prime = mpz(2)

   p1 = get_points(p1, n)
   if not isinstance(p1, list):
      return p1

   for _ in xrange(int(log_mo / LOG_2)):
      p1 = double(p1, n)
      if not isinstance(p1, list):
         return p1

   for i in xrange(1, det_bound):
      prime  = (i << 1) + 1
      if isprime(prime):
         for _ in xrange(int(log_mo / math.log(prime))):
            p1 = multiply(p1, prime, n)
            if not isinstance(p1, list):
               return p1

   while prime < sqrt(u) and isinstance(p1, list):
      for i in xrange(number_of_primes):
         prime = next_prime(prime)
         p[i] = prime ** max(1, int(log_u / math.log(prime)))
      p1 = fast_multiply(p1, prod(p),  n, w)

   if not isinstance(p1, list):
      return p1

   while prime < u and isinstance(p1, list):
      for i in xrange(number_of_primes):
         prime = next_prime(prime)
         p[i] = prime
      p1 = fast_multiply(p1, prod(p),  n, w)

   if not isinstance(p1, list):
      return p1

   del p

   small_jump = int(greatest_n((1 << (w + 2)) / 3))
   small_jump = max(120, small_jump)
   big_jump = 1 + (int(sqrt((5 << w) / 21)) << 1)
   total_jump = small_jump * big_jump
   big_multiple = max(total_jump << 1, ((int(next_prime(prime)) - (total_jump >> 1)) / total_jump) * total_jump)
   big_jump_2 = big_jump >> 1
   small_jump_2 = small_jump >> 1
   product = ONE

   psmall_jump = multiply(p1, small_jump, n)
   if not isinstance(psmall_jump, list):
      return psmall_jump

   ptotal_jump = multiply(psmall_jump, big_jump, n)
   if not isinstance(ptotal_jump, list):
      return ptotal_jump

   pgiant_step = multiply(p1, big_multiple, n)
   if not isinstance(pgiant_step, list):
      return pgiant_step

   small_multiples = [None]
   for i in xrange(1, small_jump >> 1):
      if gcd(i, small_jump) == 1:
         tmp = multiply(p1, i, n)
         if not isinstance(tmp, list):
            return tmp
         for i in xrange(len(tmp)):
            tmp[i] = tmp[i][0]
         small_multiples.append(tuple(tmp))
      else:
         small_multiples.append(None)
   small_multiples = tuple(small_multiples)

   big_multiples = [None]
   for i in xrange(1, (big_jump + 1) >> 1):
      tmp = multiply(psmall_jump, i, n)
      if not isinstance(tmp, list):
         return tmp
      big_multiples.append(to_tuple(tmp))
   big_multiples = tuple(big_multiples)

   psmall_jump = to_tuple(psmall_jump)
   ptotal_jump = to_tuple(ptotal_jump)

   while big_multiple < u2:
      big_multiple += total_jump
      center_up = big_multiple
      center_down = big_multiple
      pgiant_step = add(ptotal_jump, pgiant_step, n)
      if not isinstance(pgiant_step, list):
         return pgiant_step

      prime_up = next_prime(big_multiple - small_jump_2)
      while prime_up < big_multiple + small_jump_2:
         s = small_multiples[abs(int(prime_up) - big_multiple)]
         for j in xrange(ncurves):
            product *= pgiant_step[j][0] - s[j]
            product %= n
         prime_up = next_prime(prime_up)

      for i in xrange(1, big_jump_2 + 1):
         center_up += small_jump
         center_down -= small_jump

         pmed_step_up, pmed_step_down = add_sub_x_only(big_multiples[i], pgiant_step, n)
         if pmed_step_down == None:
            return pmed_step_up

         while prime_up < center_up + small_jump_2:
            s = small_multiples[abs(int(prime_up) - center_up)]
            for j in xrange(ncurves):
               product *= pmed_step_up[j] - s[j]
               product %= n
            prime_up = next_prime(prime_up)

         prime_down = next_prime(center_down - small_jump_2)
         while prime_down < center_down + small_jump_2:
            s = small_multiples[abs(int(prime_down) - center_down)]
            for j in xrange(ncurves):
               product *= pmed_step_down[j] - s[j]
               product %= n
            prime_down = next_prime(prime_down)

   if gcd(product, n) != 1:
      return gcd(product, n)

   return 1

def fast_multiply(p, d, n, w):

   mask = (1 << (w << 1)) - 1
   flop = mask / 3
   g = naf(d) >> 4
   precomp = {}
   m = copy(p)
   p = double(p, n)

   for i in xrange((flop >> w) + (w & 1)):
      key = naf((i << 1) + 1)
      precomp[key] = to_tuple(m)
      precomp[((key & flop) << 1) ^ key] = precomp[key]
      m = add(p, m, n)

   while g > 0:
      if g & 1:
         t = g & mask
         sh = 1 + int(math.log(t) / LOG_4)
         for _ in xrange(sh):
            p = double(p, n)

         if g & 2:
            p = subtract(precomp[t], p, n)
         else:
            p = add(precomp[t], p,  n)

         g >>= (sh << 1)
         if not isinstance(p, list):
            return p
      else:
         p = double(p, n)
         g >>= 2

   return p

def get_points(p1, n):
   p1 = list(p1)
   invs = p1[:]
   ncurves = len(p1)

   for j in xrange(ncurves):
      sigma = mpz(p1[j])
      u = (sigma**2 - 5) % n
      v = sigma << 2
      i = (((u * u) % n) * ((v * u << 2) % n)) % n
      p1[j] = [u, v, i]
      invs[j] = (i * v) % n

   invs = parallel_invert(invs, n)
   if not isinstance(invs, list):
      return invs

   for j in xrange(ncurves):
      u, v, i = p1[j]
      inv = invs[j]

      a = (((((((v - u)**3 % n) * v) % n) * (u + u + u + v)) % n) * inv - 2) % n # <-- This line is a thing of beauty
      x_0 = (((((u * i) % n) * inv) % n) ** 3) % n # And this one gets second place
      b = ((((x_0 + a) * x_0 + 1) % n) * x_0) % n
      x_0 = (b * x_0) % n
      y_0 = (b**2) % n

      while a % 3:
         a += n

      x_0 = (x_0 + a * b / 3) % n
      c = (y_0 * ((1 - a**2 / 3) % n)) % n

      p1[j] = [x_0, y_0, c]

   return p1

def isprime(n):
   if not fastprime(n):
      return False
   elif n < SMALLEST_COUNTEREXAMPLE_FASTPRIME:
      return True

   do_loop = False
   j = 1
   d = n >> 1
   a = 2
   bound = int(0.75 * math.log(math.log(n)) * math.log(n)) + 1

   while not d & 1:
      d >>= 1
      j += 1

   while a < bound:
      a = next_prime(a)
      p = atdn(a, d, n)

      if p == 1 or p == n - 1:
         continue

      for _ in xrange(j):
         p = (p * p) % n

         if p == 1:
            return False
         elif p == n - 1:
            do_loop = True
            break

      if do_loop:
         do_loop = False
         continue

      return False

   return True

def multiply(p1, d, n):
   pos = int(math.log(d) / LOG_2) - 1
   p = copy(p1)

   while pos >= 0:
      p = double(p, n)
      if not isinstance(p, list):
         return p
      if (d >> pos) & 1:
         p = add(p1, p,  n)
         if not isinstance(p, list):
            return p
      pos -= 1

   return p

def ecm(n, ra, ov, veb, tdb, pr):
   if veb:
      looking_for = 0
   k = inv_const(n)

   if ra:
      sigma = 6 + random.randrange(BILLION)
   else:
      sigma = 6

   for factor in sure_factors(n, k, range(sigma, sigma + k), veb, ra, ov, tdb, pr):
      yield factor
      n /= factor

   if n == 1:
      return

   if ra:
      sigma += k + random.randrange(BILLION)
   else:
      sigma += k

   x_max = 0.5 * math.log(n) / math.log(k)
   t = rho_ts(int(x_max))
   prime_probs = []
   nc = 1 + int(_12_LOG_2_OVER_49 * ov * ov * k)
   eff_nc = nc / pr

   for i in xrange(1 + (int(math.log(n)) >> 1)):
      if i < math.log(tdb):
         prime_probs.append(0)
      else:
         prime_probs.append(1.0/i)

   for i in xrange(len(prime_probs)):
      p_success = rho_ev((i - 2.65) / math.log(k), t)
      p_fail = max(0, (1 - p_success * math.log(math.log(k)))) ** (k / pr)
      prime_probs[i] = p_fail * prime_probs[i] / (p_fail * prime_probs[i] + 1 - prime_probs[i])

   while n != 1:
      low = int(k)
      high = n
      while high > low + 1:
         u = (high + low) >> 1
         sum = 0
         log_u = math.log(u)
         for i in xrange(len(prime_probs)):
            log_p = i - 2.65
            log_u = math.log(u)
            quot = log_p / log_u
            sum += prime_probs[i] * (rho_ev(quot - 1, t) - rho_ev(quot, t) * log_u)
         if sum < 0:
            high = u
         else:
            low = u

      if ra:
         sigma += nc + random.randrange(BILLION)
      else:
         sigma += nc

      for factor in sure_factors(n, u, range(sigma, sigma + nc), veb, ra, ov, tdb, pr):
         yield factor
         n /= factor

      for i in xrange(len(prime_probs)):
         p_success = rho_ev((i - 2.65) / math.log(u), t)
         p_fail = max(0, (1 - p_success * math.log(math.log(u)))) ** eff_nc
         prime_probs[i] = p_fail * prime_probs[i] / (p_fail * prime_probs[i] + 1 - prime_probs[i])
      prime_probs = prime_probs[:1 + (int(math.log(n)) >> 1)]

      if veb and n != 1:
         m = max(prime_probs)
         for i in xrange(len(prime_probs)):
            if prime_probs[i] == m:
               break

         new_looking_for = (int(i / _5_LOG_10) + 1)
         new_looking_for += new_looking_for << 2
         if new_looking_for != looking_for:
            looking_for = new_looking_for

   return

def factors(n, veb, ra, ov, pr):

   while not n & 1:
      n >>= 1
      yield 2

   n = mpz(n)
   k = inv_const(n)
   prime = 2
   trial_division_bound = max(10 * k**2, 100)

   while prime < trial_division_bound:
      prime = next_prime(prime)
      while not n % prime:
         n /= prime
         yield prime

   if isprime(n):
      yield n
      return

   if n == 1:
      return

   for factor in ecm(n, ra, ov, veb, trial_division_bound, pr):
      yield factor

def interactive(veb, ra, ov, pr):
   print "ECC - Prime Factorization"
   response = raw_input("\nEnter number to factorize (n): ")
   n = eval(response)
   int(n)
   print '\nFactoring number %d:' % n
   if n == 0:
      print '0 does not have a well-defined factorization.'
   else:
      if n < 0:
         print -1
         n = -n
      elif n == 1:
         print 1
      if ov == DUMMY:
         ov = 2 * math.log(math.log(n))
      nf = 0
      for factor in factors(n, veb, ra, ov, pr):
         nf = nf + 1
         print factor
      print "\nNumber of factors = ", nf

def main():
   ra = veb = False
   pr = 1.0
   ov = DUMMY
   interactive(veb, ra, ov, pr)

if __name__ == '__main__':
   try:
      main()
   except (EOFError, KeyboardInterrupt):
      sys.exit()