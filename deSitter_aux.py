# Auxiliary functions
# Description of 2-point correlation functions in deSitter spacetime
# Requires numpy 

import numpy as np 

# Window and related functions

def F(delta):
	return 1/4 * (delta + 2)*(delta**2 + 2*delta + 2)

def G(delta):
	return 8*(delta**3 + 5*delta**2 + 10*delta + 10)/(5*(delta+2)**2*(delta**2+2*delta+2)**2)

# Caligraphic functions
## The functions Ka, valid for small beta and delta

def Ka(beta,mu,delta):
    if mu == -3:
        return (
            1/(2*beta**2) + 1/200*(40*(delta+1)*np.log(2*beta)-73*delta+40*np.euler_gamma*(delta +1) -93)
            - 3/350*beta**2*(2*delta+1)+beta**4/4725*(1+3*delta)
        )
    if mu == -1:
        return (
            -np.log(2*beta)-np.euler_gamma+7/4-delta/2-7*delta**2/48+beta**2/10*(1+delta+5*delta**2/6)
        )
    if mu == 1:
        return (
            9/4*(1-delta)-beta**2/2
        )
    if mu == 3:
        return (
            9/4*(1-2*np.log(delta/2)) + 9*delta*np.log(9/2)-(9/4+np.log(512))*delta-beta**4/4
        )

# The function L, valid for small beta and delta

## Notice that the difference between small and large rho is placed at rho = 1,
## in theory they are limits, but the transition works and is sufficiently smooth

def L(beta,mu,delta,rho):
    if np.abs(rho) < 1:
        if mu == -3:
            return (
                1/(2*beta**2)+(np.log(2*beta)+np.euler_gamma)/5*(1+delta)-(93+73*delta)/200-3/350*(1+2*delta)*beta**2
                + ((np.log(2*beta)+np.euler_gamma)/6-7/24+delta/12)*rho**2 - beta**2*rho**2/60*(1+delta)
            )
        if mu == -1:
            return (
                7/4 - np.log(2*beta)-np.euler_gamma-delta/2 + beta**2/10*(1+delta) + (beta**2/12 - 3/8 +3/8*delta)*rho**2
                - 3/160*(delta + (2-4*delta)*np.log(delta/2)-1)*rho**4
            )
    else:
        if mu == -3:
            return (
                1/(2*beta**2) + ((np.log(beta*np.abs(rho))+np.euler_gamma)/6-11/36)*rho**2 + (np.log(beta*np.abs(rho)) + np.euler_gamma - 1)/5*(1+delta)
                +3/175*(1+2*delta)/rho**2 - rho**4*beta**2/240 - rho**2*beta**2/60*(1+delta)-3/350*beta**2*(1+2*delta)
            )
        if mu == -1:
            return (
                1 - np.euler_gamma - np.log(np.abs(rho)*beta)-(1+delta)/(5*rho**2)+rho**2*beta**2/12-6/(175*rho**4)*(1+2*delta)
            )
        if mu == 1:
            return (
                1/rho**2 + 4*(delta**4+4*delta**3+7*delta**2+6*delta+3)/(15*(delta**2+2*delta+2)*rho**4)-beta**2/2+72*(1+2*delta)/(175*rho**6)
            )
        if mu == 3:
            return (
                -2/rho**4 - 16*(delta**4+4*delta**3+7*delta**2+6*delta+3)/(5*(delta**2+2*delta+2)*rho**6)-beta**4/4+rho**2*beta**6/36
            )
        
# The function M, valid for small beta and delta

## Notice that the difference between small and large rho is placed at rho = 1,
## in theory they are limits, but the transition works and is sufficiently smooth

def M(beta,mu,delta,rho):
    if np.abs(rho) < 1:
        if mu == -3:
            return (
                1/(2*beta**2) + 1/200*(40*(delta+1)*np.log(2*beta)-73*delta+40*np.euler_gamma*(delta+1)-93)
                +(4*np.log(2*beta)+4*np.euler_gamma-7+2*delta)/8*rho**2 - 3/350*beta**2*(2*delta+1)+3/32*(1-delta)*rho**4
                -(1+delta)/20*beta**2*rho**2+beta**4*(3*delta+1)/4725+1/320*(delta+(2-4*delta)*np.log(delta/2)-1)*rho**6
            )
        if mu == -1:
            return (
                7/4 - np.log(2*beta) - np.euler_gamma -delta/2+beta**2/10*(1+delta)+(beta**2/4 - 9/8 +9*delta/8)*rho**2
                +3/32*(1-delta-2*(1-2*delta)*np.log(delta/2))*rho**4
            )
    else:
        if mu == -3:
            return (
                1/(2*beta**2)+rho**2/4*(-3+2*np.euler_gamma + 2*np.log(beta*np.abs(rho)))+(np.euler_gamma + np.log(beta*np.abs(rho)))/5*(1+delta) - 3/175*(1+2*delta)/rho**2
                - beta**2*rho**4/48 - (1+delta)/20*beta**2*rho**2 - 3/350*(1+2*delta)*beta**2
            )
        if mu == -1:
            return (
                - np.log(beta*np.abs(rho)) - np.euler_gamma + (1+delta)/(5*rho**2) + rho**2*beta**2/4 + (1+delta)*beta**2/10
            )

# Function that builds the covariance matrix
## Takes arguments (alpha, beta, HR, delta) characterizing the quantum state
## Returns a dictionary with the matrix itself and useful derived quantities

def compute_gammas_exact(al, beta, HR, delta):
	gam11 = 2/(3*np.pi*G(delta))*HR**2*(Ka(beta,-1,delta) + Ka(beta,1,delta)/HR**2)
	gam12 = -2/(3*np.pi*G(delta))*HR*Ka(beta,1,delta)
	gam21 = gam12
	gam22 = 2/(3*np.pi*G(delta))*Ka(beta,3,delta)
	gam13 = 2/(3*np.pi*G(delta))*HR**2*(L(beta,-1,delta,al) + L(beta,1,delta,al)/HR**2)
	gam14 = -2/(3*np.pi*G(delta))*HR*L(beta,1,delta,al)
	gam24 = 2/(3*np.pi*G(delta))*L(beta,3,delta,al)
	gam23 = gam14
	gam31 = gam13
	gam32 = gam23
	gam33 = gam11
	gam34 = gam12
	gam41 = gam14
	gam42 = gam24
	gam43 = gam34
	gam44 = gam22

	gam = 1/2*np.array([[gam11, gam12, gam13, gam14],[gam21, gam22, gam23, gam24],[gam31, gam32, gam33, gam34], [gam41, gam42, gam43, gam44]],dtype=np.float64)
	gam1 = 1/2*np.array([[gam11, gam12],[gam21, gam22]],dtype=np.float64)
	gam2 = 1/2*np.array([[gam33, gam34],[gam43, gam44]],dtype=np.float64)
	gamoff12 = 1/2*np.array([[gam13, gam14],[gam23, gam24]])
	gamoff21 = np.transpose(gamoff12)
	gamphi = 1/2*np.array([[gam11,gam13],[gam31,gam33]])
	gampi = 1/2*np.array([[gam22,gam24],[gam42,gam44]])

	invgam = np.linalg.inv(gam)
	invgamphi = np.array([[invgam[0,0],invgam[0,2]],[invgam[2,0],invgam[2,2]]])
	invgampi = np.array([[invgam[1,1],invgam[1,3]],[invgam[3,1],invgam[3,3]]])
	invgamphipi = np.array([[invgam[0,1],invgam[0,3]],[invgam[2,1],invgam[2,3]]])
	invgampiphi = np.array(([invgam[1,0],invgam[1,2]],[invgam[3,0],invgam[3,2]]))
	invgamphiinv = np.linalg.inv(invgamphi)
	invgampiinv = np.linalg.inv(invgampi)

	a = invgamphi - invgamphipi@invgampiinv@invgampiphi

	ax = a[0,0]
	ay = a[1,1]
	axy = a[0,1]

	dgam = np.linalg.det(gam)
	invdgampi = np.linalg.det(invgampi)

	cov = {} # Initiate the dictionary to be returned
	cov["gam"] = gam
	cov["gam1"] = gam1
	cov["gam2"] = gam2
	cov["gamoff12"] = gamoff12
	cov["gamoff21"] = gamoff21
	cov["gamphi"] = gamphi
	cov["gampi"] = gampi
	cov["invgam"] = invgam
	cov["invgamphi"] = invgamphi
	cov["invgampi"] = invgampi
	cov["invgamphiinv"] = invgamphiinv
	cov["invgampiinv"] = invgampiinv
	cov["ax"] = ax
	cov["ay"] = ay
	cov["axy"] = axy
	cov["dgam"] = dgam
	cov["invdgampi"] = invdgampi

	return cov # Returns the full dictionary