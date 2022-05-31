# Bell operators
# Functions that build expected values of Bell operators
# Requires numpy, math, cmath, scipy.integrate and scipy.special

import numpy as np
from scipy import special
from scipy import integrate
import math
import cmath

# Auxiliary functions

def Power(a,b): # Useful for analytic expressions coming from Mathematica
    return a**b

# GKMR Bell operators
## Exact functions of the covariance matrix
## Functions of the covariance matrix as a dictionary
## Works for both Minkowski and de Sitter spacetime

def Szz_GKMR(cov): # <SzSz>
	dgam = cov["dgam"] # extract required variable
	return 1/(4*np.sqrt(dgam))

def Sxx_GKMR(cov): # <SxSx>
	# extract required variables
	ax = cov["ax"]
	ay = cov["ay"]
	axy = cov["axy"]
	return (-2/np.pi*np.arctan(axy*(ax*ay-axy**2)**(-1/2)))

## Approximated expressions

def asymptote(delta): # asymptote on top of which the approximation is built
	return 2 * (6*Power(np.pi,2))/(5.*np.sqrt((Power(4*Power(delta,2) + 4*Power(delta,3) + Power(delta,4) + 
           2*Power(delta,2)*(3 + 3*delta + Power(delta,2))*np.log(delta) - 
           4*delta*(4 + 6*delta + 4*Power(delta,2) + Power(delta,3))*np.log(2*(1 + delta)) - 4*np.log(4*(1 + delta)) + 
           8*np.log(2 + delta) + 16*delta*np.log(2 + delta) + 18*Power(delta,2)*np.log(2 + delta) + 
           10*Power(delta,3)*np.log(2 + delta) + 2*Power(delta,4)*np.log(2 + delta),2)*
         Power(112*Power(delta,2) + 224*Power(delta,3) + 196*Power(delta,4) + 84*Power(delta,5) + 
           14*Power(delta,6) + Power(delta,4)*(15 + 15*delta + 4*Power(delta,2))*np.log(delta) - 
           8*delta*(6 + 15*delta + 20*Power(delta,2) + 15*Power(delta,3) + 6*Power(delta,4) + Power(delta,5))*
            np.log(2*(1 + delta)) - 8*np.log(4*(1 + delta)) + 16*np.log(2 + delta) + 48*delta*np.log(2 + delta) + 
           120*Power(delta,2)*np.log(2 + delta) + 160*Power(delta,3)*np.log(2 + delta) + 
           105*Power(delta,4)*np.log(2 + delta) + 33*Power(delta,5)*np.log(2 + delta) + 4*Power(delta,6)*np.log(2 + delta)
           ,2))/(Power(delta,8)*Power(10 + 10*delta + 5*Power(delta,2) + Power(delta,3),4))))

def Bell_GKMR_largealpha(alpha,delta): # Large alpha approximation for GKMR Bell
	asym = asymptote(delta)
	return 2*np.sqrt((asym/2)**2+(8*(1+delta)/(9*np.pi*alpha**2))**2)

def Bell_GKMR_smalldelta(alpha,logdelta): # Small delta approximation for GKMR Bell
	return (8*np.pi**2/(9*np.abs((1-2*logdelta+2*np.log(2))))*(1+2/alpha**4*(4/81+((1-2*logdelta+2*np.log(2))/np.pi**3)**2)))

# Larsson Bell operators
## Functions of the parameter l and the covariance matrix
## Computed as a series with cut-off given by max

def Szz_L(l,max,cov): # <SzSz> is a series with a cut-off

	# extract required variables
	ax = cov["ax"]
	ay = cov["ay"]
	axy = cov["axy"]
	dgam = cov["dgam"]
	invdgampi = cov["invdgampi"]
	
	def J(n,m,l): # Auxiliary to be summed
		def temp(x):
			return math.exp(-1/2 * x**2 * (ax - axy**2/ay)) * (special.erf(1/math.sqrt(2*ay) * (axy * x + ay * (m+1)*l )) - special.erf(1/math.sqrt(2*ay) * (axy * x + ay * m*l )))
		return math.sqrt(np.pi/(2*ay)) * integrate.quad(lambda x : temp(x), n*l, l*(n+1))[0]

	# return the sum
	return 1/(2*np.pi*math.sqrt(dgam*invdgampi)) * np.sum([[(-1)**(n+m)*J(n,m,l) for n in range(-max,max+1)] for m in range(-max,max+1)])

def Sxx_L(l,max,cov): # <SxSx> is a series with a cut-off

	# extract required variables
	ax = cov["ax"]
	ay = cov["ay"]
	axy = cov["axy"]
	dgam = cov["dgam"]
	invdgampi = cov["invdgampi"]
	invgam = cov["invgam"]
	invgampi = cov["invgampi"]
	invgampiinv = cov["invgampiinv"]

	# some auxiliaries
	def atx(eps):
		return 1/2 * (invgam[0,1] * invgampiinv[0,0] * eps[0] + invgam[0,1] * invgampiinv[0,1] * eps[1] + 
			invgam[0,3] * invgampiinv[1,0] * eps[0] + invgam[0,3] * invgampiinv[1,1] * eps[1])

	def aty(eps):
		return 1/2 * (invgam[2,1] * invgampiinv[0,0] * eps[0] + invgam[2,1] * invgampiinv[0,1] * eps[1] + 
			invgam[2,3] * invgampiinv[1,0] * eps[0] + invgam[2,3] * invgampiinv[1,1] * eps[1])

	def K(n,m,l):
		def temp(x,eps):
			return cmath.exp(-1/2 * x**2 * (ax - axy**2/ay) + x*(1j*l/ay*(axy*aty(eps)-ay*atx(eps))) - (l*aty(eps))**2/(2*ay) -1/2*l**2 * (eps@invgampi@eps) ) * (special.erf(1/math.sqrt(2*ay) * (axy * x + l*aty(eps)*1j + ay * (2*m+3/2)*l )) - special.erf(1/math.sqrt(2*ay) * (axy * x + l*aty(eps)*1j + ay * (2*m+1/2)*l )))
		def temp2(x):
			return np.real(temp(x,[1,1]) + temp(x,[-1,1]) + temp(x,[1,-1]) + temp(x,[-1,-1]))
		return math.sqrt(np.pi/(2*ay)) * integrate.quad(lambda x : temp2(x), l/2+(2*n)*l, l/2+l*(2*n+1))[0]

	# return the sum
	return 1/(2*np.pi*math.sqrt(dgam*invdgampi)) * np.sum([[K(n,m,l) for n in range(-max,max+1)] for m in range(-max,max+1)])

## Approximated expressions for small l

def Sxx_L_app(l,cov):
	invgampi = cov["invgampi"] # extract
	total = 0
	for eps in np.array([[1,1],[1,-1],[-1,1],[-1,-1]]): # sum
		total += np.exp(-1/2*l**2*eps@invgampi@eps)
	return 1/4*total # return sum

def Szz_L_app(l,cov):

	#extract
	ax = cov["ax"]
	ay = cov["ay"]
	axy = cov["axy"]
	
	# return exp
	return np.exp(-np.pi**2*(ax+ay-2*axy)/(2*l**2*(ax*ay-axy**2)))
























