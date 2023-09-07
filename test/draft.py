# Input parameters 
N       = 384 # number of Chebyshev polynomials

R  = 180         # Reynolds number
kx = 1           # streamwise wavenumber
kz = 0           # spanwise wavenumber
Ri = 0.0        # Richardson number
Pr = 0.71        # Prantl number
Ra = 8*R*R/Pr*Ri # Rayleigh number
print('Rayleigh number is '+ str(Ra))

bf = 2   # set = 1 for Couette, 2 for Poiseuille, 3 for quiescent