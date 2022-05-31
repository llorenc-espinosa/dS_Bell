# Quantum discord
# Functions that compute quantum discord and related quantities
# Requires numpy

import numpy as np

# Auxiliary functions

def purity(cov):
	dgam = 16 * cov["dgam"]
	# extract covariance matrix
	# note that the factor 2**4 is required due to different conventions
	return 1/np.sqrt(dgam)

def symplectic(cov): # input is covariance matrix dictionary
	
	gam = 2 * cov["gam"]
	# extract covariance matrix
	# note that the factor 2 is required due to different conventions

	# compute symplectic eigenvalues
	sigma12 = np.sqrt(gam[0,0]*gam[1,1]-gam[1,2]**2)
	sigmaplus = np.sqrt((gam[0,0]+gam[0,2])*(gam[1,1]+gam[1,3])-(gam[0,1]+gam[0,3])**2)
	sigmaminus = np.sqrt((gam[0,0]-gam[0,2])*(gam[1,1]-gam[1,3])-(gam[0,1]-gam[0,3])**2)
	
	# return symplectic eigenvalues
	return sigma12, sigmaplus, sigmaminus

def f(x): # needed later
	return (x + 1)/2 *np.log2((x + 1)/2) - (x - 1)/2 * np.log2((x - 1)/2)


# Measures of mutual information
## Iinfo and Jinfo are alternative quantum counterparts to classical mutual information
## Both take covariance matrix dictionary and return value of mutual information

def Iinfo(cov):
	sigma12, sigmaplus, sigmaminus = symplectic(cov)
	return 2*f(sigma12) - f(sigmaplus) - f(sigmaminus)

def Jinfo(cov):
	
	gam = 2 * cov["gam"] # extract covariance matrix

	def gamA(gam):
		return np.array([[gam[0,0],gam[0,1]],[gam[0,1],gam[1,1]]])
	
	def gamC(gam):
		return np.array([[gam[0,2],gam[0,3]],[gam[1,2],gam[1,3]]])
	
	sigma12, sigmaplus, sigmaminus = symplectic(cov)
	
	sigma3 = np.array([[1,0],[0,-1]])
	
	def gamCtilde(gam):
		product = np.linalg.matrix_power(gamA(gam),-1/2)@gamC(gam)@np.linalg.matrix_power(gamA(gam),-1/2)
		if np.linal.det(gamC(gam)) < 0:
			return sigma3@product
		else:
			return product
	
	def frake(gam):
		np.sqrt(np.min(np.linalg.eigvals(gamCtilde(gam)@np.linalg.transpose(gamC(tilde(gam))))))
	
	def frakf(gam):
		np.sqrt(np.max(np.linalg.eigvals(gamCtilde(gam)@np.linalg.transpose(gamC(tilde(gam))))))
		
	def fraks(gam):
		return ((sigma12**2 + 1)*np.linalg.det(gamC(gam))**2*(sigma12**2+sigmaplus**2*sigmaminus**2)
				- (sigma12**4 - sigmaplus**2 * sigmaminus**2)**2
				)
	
	def eplus(gam):
		return ((2 * np.linalg.det(gamC(gam)**2) + (sigma12**2-1)*(sigmaplus**2 * sigmaminus**2 - sigma12**2)
		 + 2*np.abs(np.linalg.det(gamC(gam)))*np.sqrt(np.linalg.det(gamC(gam))**4
		+ (sigma12**2-1) * (sigmaplus**2 * ssigmaminus**2 - sigma12**2)))/(sigma12**2-1)**2
				)
	
	def eminus(gam):
		return ( (sigma12**4 - np.linalg.det(gamC(gam))**2 + sigmaplus**2*sigmaminus**2 - 
				np.sqrt(np.linalg.det(gamC(gam))**4 + (sigma12**4 - sigmaplus**2*sigmaminus**2)**2 - 
				2 * np.linalg.det(gamC(gam))**2 * (sigma12**4+sigmaplus**2*sigmaminus**2)))/(2*sigma12**2)
					   )
	
	def emin(gam):
		if fraks(gam) > 0:
			return eplus(gam)
		else:
			return eminus(gam)
	
	return f(sigma12) - f(np.sqrt(emin(gam)))

# Quantum discord is simply the difference

def discord(cov):
	return Iinfo(cov)-Jinfo(cov)
