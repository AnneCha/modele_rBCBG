#!/usr/bin/python
# -*- coding: utf-8 -*-

# Computation of the different parameters for the rBCBG model
# author : Anne Chadoeuf
# 9th March 2009

import csv
from math import sqrt, cosh, exp, pi

solutions = csv.DictReader(open("solutions_simple.csv"))

# fixed parameters
A_GABA=0.25
A_AMPA=1
A_NMDA=0.025
D_GABA=5
D_AMPA=5
D_NMDA=100
Ri=200.E-2
Rm=20000.E-4

N_GPe = 251000.0
N_STN =  77000.0
N_FSI =  532000.0
N_MSN = 26448000.0
N_GPi = 143000.0

NUCLEI=['MSN','FSI','STN','GPe','GPi']
dx={'MSN':1.E-6,'FSI':1.5E-6,'STN':1.5E-6,'GPe':1.7E-6,'GPi':1.2E-6}
lx={'MSN':619E-6,'FSI':961E-6,'STN':750E-6,'GPe':865E-6,'GPi':1132E-6}
LX={}
for n in NUCLEI:
    LX[n]=lx[n]*sqrt((4*Ri)/(dx[n]*Rm))

P_MSN_GPe = 0.82
P_STN_GPe = 0.83
P_STN_GPi = 0.72
P_STN_MSN = 0.17
P_STN_FSI = 0.17
P_GPe_GPe = 0.84
P_GPe_GPi = 0.84
P_GPe_MSN = 0.16
P_GPe_FSI = 0.16
P_MSN_GPi = 0.82

solutions2=[]
for row in solutions:
    solutions2.append(row)

line=3

alpha_GPe_STN = float(solutions2[line]['NU_GPe_STN'])
nu_GPe_STN =  N_GPe / N_STN * alpha_GPe_STN

alpha_GPe_GPe = float(solutions2[line]['NU_GPe_GPe'])
nu_GPe_GPe = P_GPe_GPe * alpha_GPe_GPe

alpha_GPe_FSI = float(solutions2[line]['NU_GPe_FSI'])
nu_GPe_FSI = P_GPe_FSI * N_GPe / N_FSI * alpha_GPe_FSI

alpha_GPe_MSN = float(solutions2[line]['NU_GPe_MSN'])
nu_GPe_MSN = P_GPe_MSN * N_GPe / N_MSN * alpha_GPe_MSN

alpha_MSN_GPe = float(solutions2[line]['NU_MSN_GPe'])
nu_MSN_GPe = P_MSN_GPe * (N_MSN / N_GPe) * alpha_MSN_GPe

alpha_MSN_GPi = float(solutions2[line]['NU_MSN_GPi'])
nu_MSN_GPi =  (N_MSN / N_GPi) * alpha_MSN_GPi

alpha_STN_GPi = float(solutions2[line]['NU_STN_GPi'])
nu_STN_GPi = P_STN_GPi * (N_STN / N_GPi) * alpha_STN_GPi

alpha_STN_GPe = float(solutions2[line]['NU_STN_GPe'])
nu_STN_GPe = P_STN_GPe * (N_STN / N_GPe) * alpha_STN_GPe

alpha_STN_FSI = float(solutions2[line]['NU_STN_FSI'])
nu_STN_FSI = P_STN_FSI * (N_STN / N_FSI) * alpha_STN_FSI

alpha_STN_MSN = float(solutions2[line]['NU_STN_MSN'])
nu_STN_MSN = P_STN_MSN * (N_STN / N_MSN) * alpha_STN_MSN

alpha_GPe_GPi = float(solutions2[line]['NU_GPe_GPi'])
nu_GPe_GPi = P_GPe_GPi * (N_GPe / N_GPi) * alpha_GPe_GPi

alpha_FSI_MSN = float(solutions2[line]['NU_FSI_MSN'])
nu_FSI_MSN = (N_FSI / N_MSN) * alpha_FSI_MSN

alpha_FSI_FSI = float(solutions2[line]['NU_FSI_FSI'])
nu_FSI_FSI =  alpha_FSI_FSI

alpha_MSN_MSN = float(solutions2[line]['NU_MSN_MSN'])
nu_MSN_MSN =  alpha_MSN_MSN


# computation of the weight of the different connexions

W_STN_GPi=nu_STN_GPi*(cosh(LX['GPi']-float(solutions2[line]['DIST_STN_GPi'])*LX['GPi'])/cosh(LX['GPi']))*(A_AMPA*D_AMPA+A_NMDA*D_NMDA)
W_STN_GPe=nu_STN_GPe*(cosh(LX['GPe']-float(solutions2[line]['DIST_STN_GPe'])*LX['GPe'])/cosh(LX['GPe']))*(A_AMPA*D_AMPA+A_NMDA*D_NMDA)
W_STN_MSN=nu_STN_MSN*(cosh(LX['MSN']-float(solutions2[line]['DIST_STN_MSN'])*LX['MSN'])/cosh(LX['MSN']))*(A_AMPA*D_AMPA+A_NMDA*D_NMDA)
W_STN_FSI=nu_STN_FSI*(cosh(LX['FSI']-float(solutions2[line]['DIST_STN_FSI'])*LX['FSI'])/cosh(LX['FSI']))*(A_AMPA*D_AMPA+A_NMDA*D_NMDA)
W_MSN_GPe=nu_MSN_GPe*(cosh(LX['GPe']-float(solutions2[line]['DIST_MSN_GPe'])*LX['GPe'])/cosh(LX['GPe']))*A_GABA*D_GABA
W_MSN_GPi=nu_MSN_GPi*(cosh(LX['GPi']-float(solutions2[line]['DIST_MSN_GPi'])*LX['GPi'])/cosh(LX['GPi']))*A_GABA*D_GABA
W_GPe_GPe=nu_GPe_GPe*(cosh(LX['GPe']-float(solutions2[line]['DIST_GPe_GPe'])*LX['GPe'])/cosh(LX['GPe']))*A_GABA*D_GABA
W_GPe_GPi=nu_GPe_GPi*(cosh(LX['GPi']-float(solutions2[line]['DIST_GPe_GPi'])*LX['GPi'])/cosh(LX['GPi']))*A_GABA*D_GABA
W_FSI_FSI=nu_FSI_FSI*(cosh(LX['FSI']-float(solutions2[line]['DIST_FSI_FSI'])*LX['FSI'])/cosh(LX['FSI']))*A_GABA*D_GABA
W_GPe_STN=nu_GPe_STN*(cosh(LX['STN']-float(solutions2[line]['DIST_GPe_STN'])*LX['STN'])/cosh(LX['STN']))*A_GABA*D_GABA
W_GPe_MSN=nu_GPe_MSN*(cosh(LX['MSN']-float(solutions2[line]['DIST_GPe_MSN'])*LX['MSN'])/cosh(LX['MSN']))*A_GABA*D_GABA
W_GPe_FSI=nu_GPe_FSI*(cosh(LX['MSN']-float(solutions2[line]['DIST_GPe_FSI'])*LX['MSN'])/cosh(LX['MSN']))*A_GABA*D_GABA
W_FSI_MSN=nu_FSI_MSN*(cosh(LX['MSN']-float(solutions2[line]['DIST_FSI_MSN'])*LX['MSN'])/cosh(LX['MSN']))*A_GABA*D_GABA
W_MSN_MSN=nu_MSN_MSN*(cosh(LX['MSN']-float(solutions2[line]['DIST_MSN_MSN'])*LX['MSN'])/cosh(LX['MSN']))*A_GABA*D_GABA


print 'w_msn_gpe = {}' .format(W_MSN_GPe)
print 'w_msn_gpi = {}' .format(W_MSN_GPi)
print 'w_stn_gpe = {}' .format(W_STN_GPe)
print 'w_stn_gpi = {}' .format(W_STN_GPi)
print 'w_stn_msn = {}' .format(W_STN_MSN)
print 'w_stn_fsi = {}' .format(W_STN_FSI)
print 'w_gpe_stn = {}' .format(W_GPe_STN)
print 'w_gpe_gpi = {}' .format(W_GPe_GPi)
print 'w_gpe_msn = {}' .format(W_GPe_MSN)
print 'w_gpe_fsi = {}' .format(W_GPe_FSI)
print 'w_gpe_gpe = {}' .format(W_GPe_GPe)
print 'w_fsi_msn = {}' .format(W_FSI_MSN)
print 'w_fsi_fsi = {}' .format(W_FSI_FSI)
print 'w_msn_msn = {}' .format(W_MSN_MSN)
print 'nu gpe stn = {}' .format(nu_GPe_STN)

#Firing rate at rest:
Phi_CSN = 2		#Hz (s-1)
Phi_PTN = 15		#Hz (s-1)
Phi_CMPf = 4		#Hz (s-1)

#Constantes pour les Medium Spiny Neurons (MSN):

Smax_MSN = 300		#Hz
Smax_STN = 300		#Hz
Smax_GPe = 400		#Hz
Smax_GPi = 400		#Hz

sigma = 7*1e3			#microV
sigma_prime = (sigma * sqrt(3) / pi) 	#microV

TAU_CTX_STR = 7#00
TAU_CTX_STN = 3#00
TAU_STR_GPe = 7#00
TAU_STR_GPi = 11#00
TAU_STN_GPe = 3#00
TAU_STN_GPi = 3#00
TAU_GPe_STN = 10#00
TAU_GPe_GPi = 3#00
TAU_STN_MSN = 3#00
TAU_STN_FSI = 3#00
TAU_GPe_MSN = 3#00
TAU_GPe_FSI = 3#00

theta_MSN=float(solutions2[line]['THETA_MSN'])
theta_FS=float(solutions2[line]['THETA_FSI'])
theta_STN=float(solutions2[line]['THETA_STN'])
theta_GPe=float(solutions2[line]['THETA_GPe'])
theta_GPi=float(solutions2[line]['THETA_GPi'])
