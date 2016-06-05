#!/usr/bin/python
# -*- coding: utf-8 -*-

# Rewriting in python of the basal ganglia CBG model (Girard et al., 2008, Neural Netwk.)
# This is the thalamus-frontal cortex module
# author : Benoît Girard <benoit.girard@isir.fr>
# modified by Jean Bellot & Anne CHadoeuf for the rBCBG
# 9th March 2009

import math
import numpy
from random import gauss

from connexweights import *

#-------------------------------------------
class thalamusFC:
  #-----------------------------------------
  def __init__(self,NbChannels,model='rBCBG'):
    self.NbChannels = NbChannels # number of channels in competition

    self.model = model           # model type, can be :
                                 # rBCBG only

    self.paramInit()   # parameter initialisation (connection weights, neuron biases)

    self.stateReset()            # reset to 0 of all the internal variables

    self.f=open('log/ThFC_'+model,'w')  # log file where the internal state will be stored if logAll function is used
    self.f.write(str(self.model) + '\n' + str(self.NbChannels)+ '\n')
    self.f.write('CMPf;CSN;PTN;\n')

  #-----------------------------------------
  def __del__(self):
    self.f.close()

  #-----------------------------------------
  def stateReset(self):

    Phi_PTN=numpy.array([15.,15.,15.,15.,15.,15.])

    self.Th = numpy.zeros((self.NbChannels)) # Neurons of the thalamic nucleus implied in the considered loop
    self.FC = numpy.zeros((self.NbChannels)) # Neurons of the frontal cortex area implied in the considered loop

    # Variables named "variable_old" are buffers used to store the previous output of the considered neurons
    self.old_Th = numpy.zeros((self.NbChannels))
    self.old_FC = numpy.zeros((self.NbChannels))

    # TRN: Thalmic Reticular Nucleus
    # The TRN is made of NbChannels neurons in the GPR,
    # while it is one neuron only in the CBG:
    if self.model == 'rBCBG':
      self.PTN = Phi_PTN*numpy.zeros((self.NbChannels))
      self.CSN = Phi_CSN*numpy.zeros((self.NbChannels))
      self.CMPf = Phi_CMPf*numpy.zeros((self.NbChannels))

      self.old_PTN = Phi_PTN*numpy.ones((self.NbChannels))
      self.old_CSN = Phi_CSN*numpy.ones((self.NbChannels))
      self.old_CMPf = Phi_CMPf*numpy.ones((self.NbChannels))

      self.olds_PTN  = [Phi_PTN*numpy.zeros((self.NbChannels))]*15
      self.olds_CSN  = [Phi_CSN*numpy.zeros((self.NbChannels))]*15
      self.olds_CMPf = [Phi_CMPf*numpy.zeros((self.NbChannels))]*15

      #Ajout du TRN:
      self.TRN = numpy.zeros((self.NbChannels))
      self.old_TRN = numpy.zeros((self.NbChannels))
      self.olds_TRN = [numpy.zeros((self.NbChannels))]*15

    else:
      print 'stateReset(): ', self.model, ' model type unknown'


  #-----------------------------------------
  def paramInit(self):

    if self.model == 'rBCBG':

	  self.PTN = Phi_PTN*numpy.ones((self.NbChannels))
	  self.CSN = Phi_CSN*numpy.ones((self.NbChannels))
	  self.CMPf = Phi_CMPf*numpy.ones((self.NbChannels))

	  self.old_PTN = Phi_PTN*numpy.ones((self.NbChannels))
	  self.old_CSN = Phi_CSN*numpy.ones((self.NbChannels))
	  self.old_CMPf = Phi_CMPf*numpy.ones((self.NbChannels))

	  self.olds_PTN  = [Phi_PTN*numpy.ones((self.NbChannels))]*15
	  self.olds_CSN  = [Phi_CSN*numpy.ones((self.NbChannels))]*15
	  self.olds_CMPf = [Phi_CMPf*numpy.ones((self.NbChannels))]*15

    else:
      print 'paramInit(): ', self.model, ' model type unknown'
      exit()

  #-----------------------------------------
  # updates the model state, integrating over timestep "dt" and salience input "salience",
  # using the (very) basic Euler method.
  # "BG_Input" : inhibitory input from the BG (from the GPi/SNr)


  def stepCompute(self,dt,saliences,BG_Input):

    if self.model == 'rBCBG':


      #Ici  seul  l'activité  de  CSN  change  en  fonction  de  la  salience
      self.PTN  =  Phi_PTN*numpy.ones(self.NbChannels) #+ 15*saliences    #+  numpy.random.normal(0,1,self.NbChannels)
      self.CSN  =  0.*Phi_CSN*numpy.ones(self.NbChannels) + 4*saliences  #+  Phi_CSN  #+  numpy.random.normal(0,1,self.NbChannels)
      self.CMPf  =  Phi_CMPf*numpy.ones(self.NbChannels)  #+  numpy.random.normal(0,1,self.NbChannels)
      self.TRN = 2

      self.old_PTN  = self.PTN
      self.old_CSN  = self.CSN
      self.old_CMPf = self.CMPf
      self.old_TRN  = self.TRN

          #Gestion  des  activités  plus  anciennes:
      self.olds_PTN.pop()
      self.olds_CSN.pop()
      self.olds_CMPf.pop()
      self.olds_TRN.pop()

      self.olds_PTN  =  [self.old_PTN]  +  self.olds_PTN
      self.olds_CSN  =  [self.old_CSN]  +  self.olds_CSN
      self.olds_CMPf =  [self.old_CMPf] +  self.olds_CMPf
      self.olds_TRN  =  [self.old_TRN]  +  self.olds_TRN


    #-----------------------------
    else:
      print 'paramInit(): ', self.model, ' model type unknown'

  #-----------------------------------------
  def readFC(self):
    return self.FC

  def readTHFC(self):
    return [self.olds_CSN, self.olds_PTN, self.olds_CMPf]

  #-----------------------------------------
  def logAll(self):
  # logs the internal state of the module
  # easily visualized with gnuplot : splot 'log/BG' matriw with lines
    self.f.writelines(' '.join([str(self.CMPf[i]) for i in range(self.NbChannels)])+' ')
    self.f.writelines(' '.join([str(self.CSN[i]) for i in range(self.NbChannels)])+' ')
    self.f.writelines(' '.join([str(self.PTN[i]) for i in range(self.NbChannels)])+'\n')

#---------------------------

def main():
  dt = 0.001
  THFC = thalamusFC(6)
  saliences = numpy.zeros((6))
  saliences[0] = 0.4
  BG_Input = numpy.zeros((6))

  for t in range(200):
    THFC.stepCompute(dt,saliences,BG_Input)
    THFC.logAll()

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
