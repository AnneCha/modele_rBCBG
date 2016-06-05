#!/usr/bin/python
# -*- coding: utf-8 -*-

# Simulation of a whole cortico-baso-thalamo-cortical loop (CBGTC)
# author : Benoît Girard <benoit.girard@isir.fr>
# modified by Jean Bellot & Anne Chadoeuf
# 25th September 2009

import math
import random
import numpy
import numpy.random
import basalganglia_connex
import thalamusFC
import pickle
import matplotlib.pyplot as pl
from matplotlib import rcParams
from plot_log import *

pl.rc("figure",facecolor="#ffffff") # make the background white

#-------------------------------------------
class CBGTC:
  #-----------------------------------------
  def __init__(self,NbChannels=6,model='rBCBG'):
    self.model = model
    self.NbChannels = NbChannels
    # creates the BG and Th-FC modules :
    self.BG = basalganglia_connex.BasalGanglia(NbChannels,line=1,model=model)
    self.THFC = thalamusFC.thalamusFC(NbChannels)
    # computes the inhibition at rest : any inhibition below this
    # level is considered as partial selection
    self.restInhibition = self.getInhibRest(0.001)

    print "============="
    print model+' model created'
    print self.NbChannels,'channels'
    print 'Inhibition at rest',self.restInhibition
    print "============="

  #-----------------------------------------
  # updates the model state, integrating over timestep "dt" and salience inumpyut "salience",
  # using the (very) basic Euler method.


  def __del__(self):
    del self.BG
    del self.THFC

  #-----------------------------------------
  # updates the model state, integrating over timestep "dt" and salience inumpyut "salience",
  # using the (very) basic Euler method.

  def stepCompute(self,dt,saliences):

    inhibs = self.BG.readGPi()

    [CSN,PTN,CMPf] = self.THFC.readTHFC()

    self.BG.stepCompute(dt,saliences,CSN,PTN,CMPf)
    self.THFC.stepCompute(dt,saliences,inhibs)
    return self.BG.readGPi()

  #---------------------------
  # simulates the CBGTC for a given number of steps (NbSteps)
  # logs the state of the model at each timestep if verbosity[0]=='v'
  # returns inhibition levels

  def nbStepsCompute(self,dt,NbSteps,saliences,verbosity='v'):

    for t in range(NbSteps):
      self.stepCompute(dt,saliences)
      if verbosity[0] == 'v':
        self.logAll()
    return self.BG.readGPi()

  #---------------------------
  # simulates the CBGTC loop until convergence of all channels
  # i.e. until |GPi(t+dt)-GPi(t)| < threshold
  # stops before convergence if t>3s
  # logs the state of the model at each timestep if verbosity[0]=='v'
  # returns time to convergence and inhibition levels
  def CvgCompute(self,dt,threshold,saliences,verbosity='v'):

    t = dt
    self.stepCompute(dt,saliences)
    if verbosity[0] == 'v':
        self.logAll()
    cvg = False

    while ((cvg == False) or (t<0.1)) and (t<3.0):
      t+=dt
      inhibs = self.BG.readGPi()

      self.stepCompute(dt,saliences)
      if verbosity[0] == 'v':
        self.logAll()

      new_inhibs = self.BG.readGPi()

      cvg = True
      for i in range(len(inhibs)) :
        if abs(inhibs[i]-new_inhibs[i]) >= threshold :
          cvg = False
          break

    #print t,new_inhibs
    return t,new_inhibs

  #---------------------------
  # returns the level of inhibition at rest in the GPi
  def getInhibRest(self,dt):
    saliences = numpy.zeros((self.NbChannels))
    inh = self.nbStepsCompute(dt,1000,saliences,'stfu')
    return inh[0]

  #-----------------------------------------
  # logs the internal state of the loop
  # easily visualized with gnuplot : splot 'log/moduleName' matrix with lines
  def logAll(self):
    self.BG.logAll()
    self.THFC.logAll()

  #--------------------------------------------------------
  # simulates the selection test from the (Gurney et al, 2001b) paper
  # returns a score between 0 and 1, depending on the completion of the success criteria
  # verbosity 'v' logs internal state
  # verbosity 'vv' prints step results on the terminal
  def simpleTest(self,dt,verbosity='vv'):
    score = 0

    # STEP 1
    #--------
    saliences = numpy.zeros((self.NbChannels))
    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)

    print inhibs

    if inhibs[0] > 0.01:
      score += 0.2
      if verbosity=='vv':
        print 'step 1 : inhibitory output at rest %.2f' %inhibs[0]
    else :
      if verbosity=='vv':
        print 'step 1 : no inhibitory output at rest'

    # STEP 2
    #--------
    saliences[0] = 0.4

    print self.BG.readGPi()
    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)

    print 'Step 2'
    print 'canal 1:' + str(inhibs[0])
    print 'canal 2:' + str(inhibs[1])
    print 'canal 3:' + str(inhibs[2])

    if (inhibs[0] < self.restInhibition) and (inhibs[1] >= self.restInhibition):
      score += 0.2


      if verbosity=='vv':
        print 'step 2 : channel 1 selected'
    else :
      if verbosity=='vv':
        print 'step 2 : channel 1 not selected'

    print self.BG.readGPi()

    # STEP 3
    #--------
    saliences[1] = 0.6
    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)
    print 'Step 3'
    print 'canal 1:' + str(inhibs[0])
    print 'canal 2:' + str(inhibs[1])
    print 'canal 3:' + str(inhibs[2])
    if (inhibs[0] > inhibs[1]) and  (inhibs[1] < self.restInhibition) :
      score+=0.1
      if inhibs[0] >= self.restInhibition:
	score+=0.1
        if verbosity=='vv':
          print 'step 3 : Channel 2 selected alone'
      else:
        if verbosity=='vv':
          print 'step 3 : Channel 2 more selected than channel 1'
    else:
      if verbosity=='vv':
        print 'step 3 : Channel 2 not selected, or channel 1 more selected than channel 2'


    # STEP 4
    #--------
    saliences[0] = 0.6
    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)
    print 'Step 4'
    print 'canal 1:' + str(inhibs[0])
    print 'canal 2:' + str(inhibs[1])
    print 'canal 3:' + str(inhibs[2])
    if (inhibs[0] < self.restInhibition) and (inhibs[1] < self.restInhibition):
      score+=0.1
      if (inhibs[0]-inhibs[1]<0.005):
        score+=0.1
        if verbosity=='vv':
          print 'step 4 : Channels 1 and 2 similarly selected'
      else:
        if verbosity=='vv':
          print 'step 4 : Channels 1 or 2 not similarly selected'
    else:
      if verbosity=='vv':
        print 'step 4 : Channels 1 or 2 not selected'

    # STEP 5
    #--------
    saliences[0] = 0.4
    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)
    print 'Step '
    print 'canal 1:' + str(inhibs[0])
    print 'canal 2:' + str(inhibs[1])
    print 'canal 3:' + str(inhibs[2])
    if (inhibs[0] > inhibs[1]) and (inhibs[1] < self.restInhibition) :
      score+=0.1
      if inhibs[0] >= self.restInhibition:
        score+=0.1
        if verbosity=='vv':
          print 'step 5 : Channel 2 selected alone'
      else:
        if verbosity=='vv':
          print 'step 5 : Channel 2 more selected than channel 1'
    else:
      if verbosity=='vv':
        print 'step 5 : Channel 2 not selected, or channel 1 more selected than channel 2'

    return score


  def testSaliences(self,dt = 0.001,verbosity = 'stfu',show = True, plot = True, saliences = numpy.array([0.1,0.14,0.2,0.3,0.23,0.15,0.135,0.12,0.1,0.07])):

    #saliences = numpy.zeros((self.NbChannels))
    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)

    if self.model == 'rBCBG':
      saliences = saliences * 5

    inhibs = self.nbStepsCompute(dt,2000,saliences,verbosity)
    m = 0
    if plot:
        if self.BG.DA < 0.4:
            s = 'low'
        elif self.BG.DA < 0.8:
            s = 'middle'
        else:
            s = 'high'

        #pl.subplot(311)
        #pl.bar(range(len(saliences)),saliences)
        #pl.title('Salience')
        #pl.ylabel('CNS input')
        #m = max(pl.ylim()[1],m)

        pl.subplot(311)
        pl.plot(range(len(inhibs)),inhibs,label = s )
        pl.title('Selection')
        pl.ylabel('GPi output')
        pl.xlabel('Channels')
        m = max(pl.ylim()[1],m)
        pl.legend()

        pl.subplot(312)
        pl.plot(range(len(inhibs)),proba(inhibs, model = self.model),label = s )
        pl.title('Selection')
        pl.ylabel('P(Ai|GPi)')
        pl.xlabel('Channels')

    if show and plot:
        pl.show()
    return inhibs

  def testSelection(self):
	dt = 0.001
	nbSteps = 51

	e1 = numpy.zeros((nbSteps,nbSteps))
	e2 = numpy.zeros((nbSteps,nbSteps))
	e3 = numpy.zeros((nbSteps,nbSteps))
	canal1 = numpy.zeros((nbSteps,nbSteps))
	canal2 = numpy.zeros((nbSteps,nbSteps))

	epsilon = self.restInhibition / 100.

	s = []
	saliences = numpy.zeros((self.NbChannels))
	for c1 in range(nbSteps):
		saliences[0] = 5*(c1/float(nbSteps - 1))

		print saliences[0]
		s += [saliences[0]]

		for c2 in range(nbSteps):
			saliences[1] = 5*(c2/float(nbSteps - 1))

			tcvg, inhibs = self.CvgCompute(dt,1e-5,saliences,'vv')

			if inhibs[1] < inhibs[0] - epsilon:
				e2[c1,c2] = 1
			elif inhibs[0] < inhibs[1] - epsilon:
				e2[c1,c2] = -1

			e1[c1,c2] = (min(inhibs))/ self.restInhibition
			e3[c1,c2] = abs(inhibs[0] - inhibs[1]) / self.restInhibition
			canal1[c1,c2] = abs(inhibs[0]) / self.restInhibition
			canal2[c1,c2] = abs(inhibs[1]) / self.restInhibition

	print numpy.arange(1,20,(1/float(nbSteps - 1))*25)[::-1]


	pickle.dump((e1,e2,e3,nbSteps), file=open('test_'+str(self.model)+'_selection.dat','w'))


  #-------------------------------------------------------
  # Computes the multiple successive vectors test
  #
  # * score[0] evaluates the capacity of the system of selecting the
  # channel with the highest inumpyut
  # * score[1] evaluates the capacity of the system of separating the
  # channel with the highest inumpyut from the channel with the second
  # highest inumpyut
  # * score[2] evaluates the amplification of the salience signal in the
  # winning FC channel
  # * score[3] evaluates the contrast of amplification between highest and second highest inumpyut
  # * score[4] is the average time of convergence
  # * score[5] is an histogram of the time of convergence (values longer than 1s are grouped in the last bin)

  def TwoHundredMotelsTest(self,dt, steps, verbosity='stfu'):

    score = [0.,0.,0.,0.,0.,numpy.zeros((100))]
    numpy.random.seed(17) # you may change the seed at your convenience
    for i in range(steps):
      saliences = numpy.random.random_sample((self.NbChannels))
      tcvg, inhibs =  self.CvgCompute(dt,1e-5,saliences,'stfu')
      score[4] += tcvg
      score[5][min(int(tcvg*100.),99)]+=1

      #-----------------------------------------
      max1 = 0. # maximum salience
      max2 = 0. # second maximum
      i1 = []   # list of the indexes of the salience maximum in the salience vector
      i2 = []   # the same for the second maximum
      for j in range(len(saliences)):
        if saliences[j]>max1:
          max2 = max1
          i2 = i1
          max1 = saliences[j]
          i1 = [j]
        elif saliences[j] == max1:
          i1.append(j)
        elif saliences[j]>max2:
          max2 = saliences[j]
          i2 = [j]
        elif saliences[j] == max2:
          i2.append(j)

      if verbosity=='vv':
        print '---------------------------'
        print 'Step :',i
        print 'Saliences   :',saliences
        print 'Inhibitions :',inhibs
        print 'FC          :',self.THFC.readFC()
        print 'Amplification Contrast :', ((float(self.THFC.readFC()[i1[0]]-max1) / max1) - (float(self.THFC.readFC()[i2[0]]-max2) / max2))/ (float(self.THFC.readFC()[i1[0]]-max1) / max1)

      #-----------------------------------------
      if (saliences.max() < self.restInhibition) :
        score[0] += 1.
        score[1] += 1.
      else:
        for m1 in i1:
          if (inhibs[m1]<self.restInhibition) and (inhibs[saliences.argmin()]>inhibs[m1]):
            score[0] += 1. / len(i1)
          for m2 in i2:
            #print inhibs, min(max(0.,(inhibs[m2]-inhibs[m1])/(self.restInhibition-inhibs[m1])),1) / (len(i1)*len(i2))
            score[1] += min(max(0.,(inhibs[m2]-inhibs[m1])/(self.restInhibition-inhibs[m1])),1) / (len(i1)*len(i2))
            if (max1>0.) and (max2>0.) and (score[2]>0.) :
              score[3] += (   (float(self.THFC.readFC()[m1]-max1) / max1)
                              - (float(self.THFC.readFC()[m2]-max2) / max2)
                              ) \
                              / (len(i1)*len(i2)) \
                              / (float(self.THFC.readFC()[m1]-max1) / max1)
            if max1>0. :
              score[2] += float(self.THFC.readFC()[m1]-max1) / max1 / len(i1)

    if verbosity[0]=='v':
      print '=============================='
      print 'Selection of the max inumpyut:',score[0]/steps
      print 'Selection contrast:        ',score[1]/steps
      print 'Amplification of the max:  ',score[2]/steps
      print 'Amplification Contrast:    ',score[3]/steps
      print 'T cvg:                     ',score[4]*1000./steps
    return score[0]/steps, score[1]/steps, score[2]/steps, score[3]/steps, score[4]/steps, score[5]/steps

  #---------------------------
  # computes selection efficiency as in the test defined in (Prescott et al 2006 Neural Netw)

  def evaluate2ChannelsCompetition(self,dt):

    nbsteps = 21
    e1=numpy.zeros((nbsteps,nbsteps))
    e2=numpy.zeros((nbsteps,nbsteps))

    saliences = numpy.zeros((self.NbChannels))
    for c1 in range(0,nbsteps):
      print 'column',c1
      for c2 in range(0,nbsteps):
        saliences[0]= c1/float(nbsteps-1)
        saliences[1]= c2/float(nbsteps-1)
        tcvg, inhibs = self.CvgCompute(dt,1e-5,saliences,'stfu')
        #inhibs = self.nbStepsCompute(dt,2000,saliences,'stfu')
        e1[c1,c2] = min(1,max(1 - inhibs[0]/self.restInhibition,0))
        e2[c1,c2] = min(1,max(1 - inhibs[1]/self.restInhibition,0))

    f1 = open('log/e1_'+self.model,'w')
    f1.writelines(' '.join([str(e1[i,j]) for i in range(0,nbsteps)]) + '\n' for j in range(0,nbsteps))
    f1.close()

    f2 = open('log/e2_'+self.model,'w')
    f2.writelines(' '.join([str(e2[i,j]) for i in range(0,nbsteps)]) + '\n' for j in range(0,nbsteps))
    f2.close()

def plot_test_selection(model):

  e1,e2,e3,nbSteps = pickle.load(open('test_rBCBG_selection.dat','r')) #'test_'+model+'_selection.dat'


  params = {
      'axes.labelsize': 8,
      'text.fontsize': 8,
      'legend.fontsize': 10,
      'xtick.labelsize': 13,
      'ytick.labelsize': 13,
      'text.usetex': False,
      'figure.figsize': [5, 5],
      'text.latex.unicode': True
      }

  rcParams.update(params)

  ax = pl.subplot(131)
  pl.title('Canal selectionne')
  im=pl.imshow(e2[::-1,:], interpolation = 'nearest',cmap = 'binary')
  pl.setp(ax, xticks = numpy.arange(0,nbSteps,5), xticklabels = numpy.arange(0,11,(1/float(nbSteps - 1))*50), yticks = numpy.arange(0,nbSteps,5),yticklabels = numpy.arange(0,11,(1/float(nbSteps - 1))*50)[::-1])
  pl.xlabel('salience du canal 1')
  pl.ylabel('salience du canal 2')
  pl.colorbar(im, fraction=0.046, pad=0.04)

  ax = pl.subplot(132)
  pl.title('Activite du canal le plus inhibe')
  im = pl.imshow(e1[::-1,:], interpolation = 'nearest', cmap = 'binary')
  pl.setp(ax, xticks = numpy.arange(0,nbSteps,5), xticklabels = numpy.arange(0,11,(1/float(nbSteps - 1))*50), yticks = numpy.arange(0,nbSteps,5),yticklabels = numpy.arange(0,11,(1/float(nbSteps - 1))*50)[::-1])
  pl.xlabel('salience du canal 1')
  pl.ylabel('salience du canal 2')
  pl.colorbar(im, fraction=0.046, pad=0.04)

  ax = pl.subplot(133)
  pl.title('Contraste')
  im = pl.imshow(e3[::-1,:], interpolation = 'nearest',cmap = 'binary')
  pl.setp(ax, xticks = numpy.arange(0,nbSteps,5), xticklabels = numpy.arange(0,11,(1/float(nbSteps - 1))*50), yticks = numpy.arange(0,nbSteps,5,),yticklabels = numpy.arange(0,11,(1/float(nbSteps - 1))*50)[::-1])
  pl.xlabel('salience du canal 1')
  pl.ylabel('salience du canal 2')
  pl.colorbar(im, fraction=0.046, pad=0.04)
  pl.show()


def plotSimpleTest(model,channels = [1,2,3]):

    evo = False

    myCBGTC = CBGTC()


    s = myCBGTC.simpleTest(1e-3,verbosity = 'vv')

    del myCBGTC

    nuclei = ['CSN','MSN','FS','STN','GPe','GPi']
    plot_log(nuclei,channels, model = model)


#=========================================
def main():
#=========================================

  dt = 0.001
  NbChannels = 6
  modeltype = 'rBCBG'

  saliences = numpy.zeros((6))
  myCBGTC = CBGTC()


  #myCBGTC.simpleTest(dt,'vv')
  #myCBGTC.TwoHundredMotelsTest(dt,200,'v')
  #myCBGTC.evaluate2ChannelsCompetition(dt) # can be pretty long

  #myCBGTC.testSelection()
  #plot_test_selection(modeltype)

  plotSimpleTest(modeltype)

  #myCBGTC.testSaliences()



  exit()

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
