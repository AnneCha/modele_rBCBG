#!/usr/bin/python
# -*- coding: utf-8 -*-

# Rewriting in python of the basal ganglia CBG model (Girard et al., 2008)
# author : Benoît Girard <benoit.girard@isir.fr>
# modified by Jean Bellot & Anne CHadoeuf for the rBCBG
# 9th March 2009

import math
import numpy

from connexweights import *

#-------------------------------------------
class BasalGanglia:
  #-----------------------------------------
  def __init__(self,NbChannels,line,invtau = 200.,model='rBCBG'):
    #NbChannels : number of channels of the model
    #line : line of the chosen solutions
    #invtau : neurons' time constant

    self.invTau = invtau

    self.NbChannels = NbChannels # number of channels in competition

    self.model = model           # model type, always rBCBG

    print self.model
    self.paramInit(line)   # parameter initialisation (connection weights, neuron biases)

    self.stateReset()            # reset to 0 of all the internal variables

    self.f=open('log/BG_' + str(self.model),'w')		# log file where the internal state will be stored if logAll function is used
    self.f.write(str(self.model) + '\n' + str(self.NbChannels)+ '\n')

    self.f.write('FS;MSN;STN;GPe;GPi;\n')

  #-----------------------------------------
  def __del__(self):
    self.f.close()

  #-----------------------------------------
  def stateReset(self):
      self.FS = numpy.zeros((self.NbChannels))
      self.old_FS = numpy.zeros((self.NbChannels))
      self.olds_FS = [numpy.zeros((self.NbChannels))]*15

      self.olds_MSN  = [numpy.zeros((self.NbChannels))]*15
      self.olds_MSN_D1  = [numpy.zeros((self.NbChannels))]*15
      self.olds_MSN_D2  = [numpy.zeros((self.NbChannels))]*15

      self.olds_STN = [numpy.zeros((self.NbChannels))]*15
      self.olds_GPe = [numpy.zeros((self.NbChannels))]*15
      self.olds_GPi = [numpy.zeros((self.NbChannels))]*1000

      self.old_MSN  = numpy.zeros((self.NbChannels))
      self.old_MSN_D1  = numpy.zeros((self.NbChannels))
      self.old_MSN_D2  = numpy.zeros((self.NbChannels))

      self.old_STN = numpy.zeros((self.NbChannels))
      self.old_GPe = numpy.zeros((self.NbChannels))
      self.old_GPi = numpy.zeros((self.NbChannels))

      self.MSN  = numpy.zeros((self.NbChannels))	# medium spiny neurons of the striatum
      self.MSN_D1  = numpy.zeros((self.NbChannels))	#with D1 dopamine receptors
      self.MSN_D2  = numpy.zeros((self.NbChannels))	#with D2 dopamine receptors

      self.STN = numpy.zeros((self.NbChannels))   # Sub-Thalamic Nucleus
      self.GPe = numpy.zeros((self.NbChannels))   # external Globus Pallidus
      self.GPi = numpy.zeros((self.NbChannels))   # internal Globus Pallidus & reticular Substantia Nigra

  #-----------------------------------------
  def paramInit(self,line):

    # invTau are 1/tau, tau being the neurons' time constants

    # W_A_B is the projection weight from neuron A to neuron B

    # I_A is the bias applied to neuron A
    solutions = numpy.loadtxt(open("compact_weights.csv","rb"),delimiter=";",skiprows=1)

    self.W_salience_MSN = solutions[line][0]*numpy.ones((self.NbChannels))
    self.W_salience_MSN_D1 = solutions[line][0]*numpy.ones((self.NbChannels))
    self.W_salience_MSN_D2 = solutions[line][0]*numpy.ones((self.NbChannels))
    self.W_CSN_MSN	  = solutions[line][0]*numpy.ones((self.NbChannels))
    self.W_CSN_MSN_D1 = solutions[line][0]*numpy.ones((self.NbChannels))
    self.W_CSN_MSN_D2 = solutions[line][0]*numpy.ones((self.NbChannels))
    self.W_CSN_FS = solutions[line][1]
    self.W_PTN_STN = solutions[line][2]
    self.W_MSN_GPe = solutions[line][3]
    self.W_MSN_GPi = solutions[line][4]
    self.W_STN_GPe = solutions[line][5]/self.NbChannels
    self.W_STN_GPi = solutions[line][6]/self.NbChannels
    self.W_STN_MSN = solutions[line][7]/self.NbChannels
    self.W_STN_FS = solutions[line][8]/self.NbChannels
    self.W_GPe_STN = solutions[line][9]
    self.W_GPe_GPi = solutions[line][10]/self.NbChannels
    self.W_GPe_MSN = solutions[line][11]/self.NbChannels
    self.W_GPe_FS = solutions[line][12]/self.NbChannels
    self.W_GPe_GPe = solutions[line][13]/self.NbChannels
    self.W_FS_MSN = solutions[line][14]/self.NbChannels
    self.W_FS_FS = solutions[line][15]/self.NbChannels
    self.W_MSN_MSN = solutions[line][16]/self.NbChannels
    self.W_CMPf_MSN = solutions[line][17]/self.NbChannels
    self.W_CMPf_FS = solutions[line][18]/self.NbChannels
    self.W_CMPf_STN = solutions[line][19]/self.NbChannels
    self.W_CMPf_GPe = solutions[line][20]/self.NbChannels
    self.W_CMPf_GPi = solutions[line][21]/self.NbChannels
    self.W_PTN_MSN = solutions[line][22]
    self.W_PTN_FS = solutions[line][23]

    self.theta_MSN = solutions[line][24]*1e3
    self.theta_FS = solutions[line][25]*1e3
    self.theta_STN = solutions[line][26]*1e3
    self.theta_GPe = solutions[line][27]*1e3
    self.theta_GPi = solutions[line][28]*1e3
    self.Smax_FS = solutions[line][29]

    print '#----------------------------------------------------'
    print  'Parametres de la simulation :'
    print 'W_CSN_MSN = %.2f' %self.W_CSN_MSN[0]
    print 'W_CSN_FS = %.2f' %self.W_CSN_FS
    print 'W_PTN_STN = %.2f' %self.W_PTN_STN
    print 'W_MSN_GPe = %.2f' %self.W_MSN_GPe
    print 'W_MSN_GPi = %.2f' %self.W_MSN_GPi
    print 'W_STN_GPe = %.2f' %(self.W_STN_GPe*self.NbChannels)
    print 'W_STN_GPi = %.2f' %(self.W_STN_GPi*self.NbChannels)
    print 'W_STN_MSN = %.2f' %(self.W_STN_MSN*self.NbChannels)
    print 'W_STN_FS = %.2f' %(self.W_STN_FS*self.NbChannels)
    print 'W_GPe_STN = %.2f' %self.W_GPe_STN
    print 'W_GPe_GPi = %.2f' %(self.W_GPe_GPi*self.NbChannels)
    print 'W_GPe_MSN = %.2f' %(self.W_GPe_MSN*self.NbChannels)
    print 'W_GPe_FS = %.2f' %(self.W_GPe_FS*self.NbChannels)
    print 'W_GPe_GPe = %.2f' %(self.W_GPe_GPe*self.NbChannels)
    print 'W_FS_MSN = %.2f' %(self.W_FS_MSN*self.NbChannels)
    print 'W_FS_FS = %.2f' %(self.W_FS_FS*self.NbChannels)
    print 'W_MSN_MSN = %.2f' %(self.W_MSN_MSN*self.NbChannels)
    print 'W_CMPf_MSN = %.2f' %(self.W_CMPf_MSN*self.NbChannels)
    print 'W_CMPf_FS = %.2f' %(self.W_CMPf_FS*self.NbChannels)
    print 'W_CMPf_STN = %.2f' %(self.W_CMPf_STN*self.NbChannels)
    print 'W_CMPf_GPe = %.2f' %(self.W_CMPf_GPe*self.NbChannels)
    print 'W_CMPf_GPi = %.2f' %(self.W_CMPf_GPi*self.NbChannels)
    print 'W_PTN_MSN = %.2f' %self.W_PTN_MSN
    print 'W_PTN_FS = %.2f' %self.W_PTN_FS

    print 'theta_MSN = %.2f' %self.theta_MSN
    print 'theta_STN = %.2f' %self.theta_STN
    print 'theta_GPe = %.2f' %self.theta_GPe
    print 'theta_GPi = %.2f' %self.theta_GPi
    print 'Smax_FS = %.2f' %self.Smax_FS
    print 'sigma_prime = %.2f' %sigma_prime
    print '------------------------------------------------------'


  #-----------------------------------------
  # updates the model state, integrating over timestep "dt" and salience input "salience",
  # using the (very) basic Euler method.
  # "FC_Input" : excitatory input from the frontal cortex
  # the update for the CBG and CBGcustom is based on lPDS neurons
  # the update for the GPR is based on leaky-integrator neurons
  def stepCompute(self,dt,saliences,olds_CSN,olds_PTN,olds_CMPf):


        ##########################################
        ## 	Calcul des Delta_V_y(t)		##
        ##########################################
        #Version simplifiée
        Delta_V_MSN =   self.W_PTN_MSN * olds_PTN[TAU_CTX_STR - 1] \
                       + self.W_STN_MSN * self.olds_STN[TAU_STN_MSN-1].sum() \
                       - self.W_GPe_MSN * self.olds_GPe[TAU_GPe_MSN-1].sum() \
                       - self.W_MSN_MSN  * self.olds_MSN[0].sum() \
                       - self.W_FS_MSN * self.olds_FS[0].sum() \
                       + self.W_CMPf_MSN * olds_CMPf[0].sum() \
                       + self.W_CSN_MSN * olds_CSN[TAU_CTX_STR - 1]


        #Fast Spiking Interneurons:
        Delta_V_FS = self.W_STN_FS * self.olds_STN[TAU_STN_FSI-1].sum() \
        - self.W_GPe_FS * self.olds_GPe[TAU_GPe_FSI-1].sum() \
        - self.W_FS_FS * self.olds_FS[0].sum() \
        + self.W_CSN_FS * olds_CSN[TAU_CTX_STR-1] \
        + self.W_PTN_FS * olds_PTN[TAU_CTX_STR-1] \
        + self.W_CMPf_FS * olds_CMPf[0].sum()



        #SubThalamic Nucleus:
        Delta_V_STN = - self.W_GPe_STN * self.olds_GPe[TAU_GPe_STN - 1] \
        + self.W_PTN_STN * olds_PTN[TAU_CTX_STN - 1] \
        + self.W_CMPf_STN * olds_CMPf[0].sum()
        + self.W_PTN_STN * olds_CSN[TAU_CTX_STN - 1] \

        #Globus Pallidus externe:
        Delta_V_GPe = self.W_STN_GPe * self.olds_STN[TAU_STN_GPe - 1].sum() \
        - self.W_MSN_GPe * self.olds_MSN[TAU_STR_GPe - 1] \
        - self.W_GPe_GPe * self.olds_GPe[0].sum() \
        + self.W_CMPf_GPe * olds_CMPf[0].sum()

        #Aucune asymétrie
        #Globus Pallidus interne:
        Delta_V_GPi = self.W_STN_GPi * self.olds_STN[TAU_STN_GPi - 1].sum() \
        - 0.82*self.W_MSN_GPi * self.olds_MSN[TAU_STR_GPi - 1] \
        - self.W_GPe_GPi * self.olds_GPe[TAU_GPe_GPi - 1].sum() \
        + self.W_CMPf_GPi * olds_CMPf[0].sum()


        # Computation of a(t + dt)
        self.MSN = self.MSN + self.invTau * dt * ( Delta_V_MSN - self.MSN )
        #self.MSN_D1 = self.MSN_D1 + self.invTau * dt * ( Delta_V_MSN_D1 - self.MSN_D1 )
        #self.MSN_D2 = self.MSN_D2 + self.invTau * dt * ( Delta_V_MSN_D2 - self.MSN_D2 )
        self.FS  = self.FS  + self.invTau * dt * ( Delta_V_FS  - self.FS  )
        self.STN = self.STN + self.invTau * dt * ( Delta_V_STN - self.STN )
        self.GPe = self.GPe + self.invTau * dt * ( Delta_V_GPe - self.GPe )
        self.GPi = self.GPi + self.invTau * dt * ( Delta_V_GPi - self.GPi )



        # Computation of y=f(a)
        self.old_MSN = Smax_MSN     / (1 + numpy.exp((self.theta_MSN - self.MSN) / sigma_prime))
        self.old_FS  = self.Smax_FS / (1 + numpy.exp((self.theta_FS  - self.FS ) / sigma_prime))
        self.old_STN = Smax_STN / (1 + numpy.exp((self.theta_STN - self.STN) / sigma_prime))
        self.old_GPe = Smax_GPe / (1 + numpy.exp((self.theta_GPe - self.GPe) / sigma_prime))
        self.old_GPi = Smax_GPi / (1 + numpy.exp((self.theta_GPi - self.GPi) / sigma_prime))

        #----------------------------
        #On gère la pile permettant de garder l'activité sur 11ms pour tous les noyaux
        #(pas besoin de plus dans un premier temps).
        #Permet de modéliser la latence entre les différents noyaux.
        #----------------------------

        self.olds_MSN.pop()
        self.olds_FS.pop()
        self.olds_STN.pop()
        self.olds_GPe.pop()
        self.olds_GPi.pop()

        self.olds_MSN = [self.old_MSN] + self.olds_MSN
        self.olds_FS  = [self.old_FS]  + self.olds_FS
        self.olds_STN = [self.old_STN] + self.olds_STN
        self.olds_GPe = [self.old_GPe] + self.olds_GPe
        self.olds_GPi = [self.old_GPi] + self.olds_GPi



  #-----------------------------------------
  def readGPi(self):
    return self.old_GPi
  #-----------------------------------------
  def readGPe(self):
    return self.old_GPe
  #-----------------------------------------
  def readSTN(self):
    return self.old_STN

  #-----------------------------------------
  # logs the internal state of the module
  # easily visualized with gnuplot : splot 'log/BG' matriw with lines
  def logAll(self):
      self.f.writelines(' '.join([str(self.old_FS[i])  for i in range(self.NbChannels)])+' ')
      self.f.writelines(' '.join([str(self.old_MSN[i]) for i in range(self.NbChannels)])+' ')
      self.f.writelines(' '.join([str(self.old_STN[i]) for i in range(self.NbChannels)])+' ')
      self.f.writelines(' '.join([str(self.old_GPe[i]) for i in range(self.NbChannels)])+' ')
      self.f.writelines(' '.join([str(self.old_GPi[i]) for i in range(self.NbChannels)]) +'\n')

  def restActivity(self,lineNb,writer,olds_CSN,olds_PTN,olds_CMPf):
      saliences = numpy.zeros((6))
      saliences[0]=0.4

      BG = BasalGanglia(6,line = lineNb)
      print lineNb
      dt=0.001
      for t in range(500):
          BG.stepCompute(dt,saliences,olds_CSN,olds_PTN,olds_CMPf)
          BG.logAll()
      print '============================================='
      print 'Rest activity of the different nuclei:'
      print 'MSN activity: %.2f' %BG.old_MSN[0]
      print 'FSI activity: %.2f' %BG.old_FS[0]
      print 'STN activity: %.2f' %BG.old_STN[0]
      print 'GPe activity: %.2f' %BG.old_GPe[0]
      print 'GPi activity: %.2f' %BG.old_GPi[0]
      print '============================================='
      writer.writerow([BG.old_MSN[0],BG.old_FS[0],BG.old_STN[0],BG.old_GPe[0],BG.old_GPi[0]])

#---------------------------

def main():

    model = 'rBCBG'
    dt = 0.001
    saliences = numpy.zeros((6))
    BG = BasalGanglia(10,line = 6)

    olds_CSN = [2*numpy.ones(10)]*15
    olds_PTN=[15*numpy.ones(10)]*15
    olds_CMPf=[4*numpy.ones(10)]*15

    activity_levels = open("activity_levels.csv","wb")
    writer_al = csv.writer(activity_levels)
    writer_al.writerow(["MSN","FSI","STN","GPe","GPi"])



#---------------------------

if __name__ == '__main__':
  # Import Psyco if available
  try:
    import psyco
    psyco.log()
    psyco.profile()
    psyco.full()
  except ImportError:
    print 'Psyco not available.'
  main()
