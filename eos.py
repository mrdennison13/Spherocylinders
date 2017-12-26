import numpy as np
import pandas as pd


def free_energy(alpha, rho, poly, order=2):
    """Calculate the free energy for a given density and alpha value"""
    if alpha >150:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / alpha)
        sig = np.log(alpha) - 1.0
    elif alpha < 1e-10:
        x = 0.0
        sig = 0.0
    else:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / (alpha*np.tanh(alpha)))
        sig = np.log(alpha/np.tanh(alpha)) + np.arctan(np.sinh(alpha))/np.sinh(alpha) - 1.0
    fe = sig + np.log(rho) - 1.0
    for i in range (2,order+1):
        Bi = 'B'+str(i)
        B = func(x,poly[Bi][0],poly[Bi][1],poly[Bi][2],poly[Bi][3],poly[Bi][4],poly[Bi][5],poly[Bi][6])
        fe += B*rho**(i-1)/(i-1)
    return fe




def chemical_potential(alpha, rho, poly, order=2):
    """Calculate the chemical potential for a given density and alpha value"""
    if alpha >150:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / alpha)
        sig = np.log(alpha) - 1.0
    elif alpha < 1e-10:
        x = 0.0
        sig = 0.0
    else:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / (alpha*np.tanh(alpha)))
        sig = np.log(alpha/np.tanh(alpha)) + np.arctan(np.sinh(alpha))/np.sinh(alpha) - 1.0
    mu = sig + np.log(rho)
    for i in range (2,order+1):
        Bi = 'B'+str(i)
        B = func(x,poly[Bi][0],poly[Bi][1],poly[Bi][2],poly[Bi][3],poly[Bi][4],poly[Bi][5],poly[Bi][6])
        mu += B*i*rho**(i-1)/(i-1)
    return mu



def pressure(alpha, rho, poly, order=2):
    """Calculate the pressure for a given density and alpha value"""
    P = rho
    if alpha >150:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / alpha)
    elif alpha < 1e-10:
        x = 0.0
    else:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / (alpha*np.tanh(alpha)))
    for i in range (2,order+1):
        Bi = 'B'+str(i)
        B = func(x,poly[Bi][0],poly[Bi][1],poly[Bi][2],poly[Bi][3],poly[Bi][4],poly[Bi][5],poly[Bi][6])
        P += B*rho**i
    return P




def nematic_order(alpha):
    """Calculate the nematic order parameter for a given alpha value"""
    if alpha > 150:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / alpha)
    elif alpha < 1e-4:
        x = 0.0
    else:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / (alpha*np.tanh(alpha)))
    if x < 0.0:
        x = -x
    return x






def calc_EOS(poly, rho_min=0.0001, rho_max=0.1, n_points=500, order=2):
    from scipy.optimize import minimize_scalar
    rho = np.linspace(rho_min, rho_max, n_points)
    EOS = pd.DataFrame({"rho":rho})
    alpha = np.empty(n_points)
    S2 = np.empty(n_points)
    Pn = np.empty(n_points)
    Pi = np.empty(n_points)
    fen = np.empty(n_points)
    fei = np.empty(n_points)
    mun = np.empty(n_points)
    mui = np.empty(n_points)
    for j in range(len(rho)):
        alpha[j] = (minimize_scalar(free_energy, args = (rho[j], poly, order), bounds = (-1e-10,2e3), method='bounded')).x
        S2[j] = nematic_order(alpha[j])
        Pn[j] = pressure(alpha[j], rho[j], poly, order)
        Pi[j] = pressure(0.0, rho[j], poly, order)
        fen[j] = free_energy(alpha[j], rho[j], poly, order)
        fei[j] = free_energy(0.0, rho[j], poly, order)
        mun[j] = chemical_potential(alpha[j], rho[j], poly, order)
        mui[j] = chemical_potential(0.0, rho[j], poly, order)
    EOS["alpha"] = alpha
    EOS["S2"] = S2
    EOS["Pr_i"] = Pi
    EOS["Pr_n"] = Pn
    EOS["mu_i"] = mui
    EOS["mu_n"] = mun
    EOS["fe_i"] = fei
    EOS["fe_n"] = fen
    return EOS





def d_mu(alpha, rho, poly, order=2):
    """Calculate the derivative of chemical potential for a given density and alpha value"""
    if alpha >150:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / alpha)
        sig = np.log(alpha) - 1.0
    elif alpha < 1e-10:
        x = 0.0
        sig = 0.0
    else:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / (alpha*np.tanh(alpha)))
        sig = np.log(alpha/np.tanh(alpha)) + np.arctan(np.sinh(alpha))/np.sinh(alpha) - 1.0
    dmu = 1.0/rho
    for i in range (2,order+1):
        Bi = 'B'+str(i)
        B = func(x,poly[Bi][0],poly[Bi][1],poly[Bi][2],poly[Bi][3],poly[Bi][4],poly[Bi][5],poly[Bi][6])
        dmu += B*i*rho**(i-2)
    return dmu



def d_pr(alpha, rho, poly, order=2):
    """Calculate the derivative of pressure for a given density and alpha value"""
    P = 1.0
    if alpha >150:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / alpha)
    elif alpha < 1e-10:
        x = 0.0
    else:
        x = 1.0 + (3.0 / (alpha**2)) - (3.0 / (alpha*np.tanh(alpha)))
    for i in range (2,order+1):
        Bi = 'B'+str(i)
        B = func(x,poly[Bi][0],poly[Bi][1],poly[Bi][2],poly[Bi][3],poly[Bi][4],poly[Bi][5],poly[Bi][6])
        P += i*B*rho**(i-1)
    return P






def func(x,p0,p1,p2,p3,p4,p5,p6):
    return (p0+p1*x+p2*x**2+p3*x**3)/(1+p4*x+p5*x**2+p6*x**3)


def dfunc(x,p0,p1,p2,p3,p4,p5,p6):
    A = (p0+p1*x+p2*x**2+p3*x**3)
    B = (1+p4*x+p5*x**2+p6*x**3)
    dA = p1 + 2.0*p2*x + 3.0*p3*x**2
    dB = p4 + 2.0*p5*x + 3.0*p6*x**2
    deriv = (dA*B - A*dB)/(B*B)
    return deriv
