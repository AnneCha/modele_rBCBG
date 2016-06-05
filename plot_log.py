#!/usr/bin/python
# -*- coding: utf-8 -*-

import pylab as pl
import numpy as np
import pickle

'''
Conseil d'utilisation:

Model rBCBG:
	plot_log(['MSN','FS','STN','GPe','GPi'], [1,2,3])


'''


#TODO: Ajouter ce qu'il faut pour la partie THFC
def plot_log(keys, channels, model = 'rBCBG', save = False, chemin = 'firing_rate'):
	f_BG = open('log/BG_'+str(model),'rb')
	f_BG.readline()
	NbChannels = int(f_BG.readline())

	nuclei_BG = f_BG.readline().split(';')[0:-1]

	log_BG = np.loadtxt(f_BG,delimiter=" ",skiprows=3)
	f_BG.close()


	f_THFC = open('log/ThFC_'+str(model),'rb')
	f_THFC.readline()
	NbChannels = int(f_THFC.readline())

	nuclei_THFC = f_THFC.readline().split(';')[0:-1]
	log_THFC = np.loadtxt(f_THFC,delimiter=" ",skiprows=3)
	f_THFC.close()


	for indk,k in enumerate(keys):
		if(k in nuclei_BG):
			log = log_BG
			nuclei = nuclei_BG
		elif(k in nuclei_THFC):
			log = log_THFC
			nuclei = nuclei_THFC
		if(not(k in nuclei_BG or k in nuclei_THFC)):
			print 'Le noyau ' + str(k) + 'est inconnu\n'
			break
		ind = nuclei.index(k)

		for indc,c in enumerate(channels):
			if c > NbChannels:
				print 'La chaine ' + str(c) + ' ne peux pas être afficher car le log ne comprend que ' + str(NbChannels) + ' en compétition.\n'
			else:
				pl.subplot(len(keys),len(channels), indk * len(channels) + indc + 1)
				if indc == 0:
					pl.ylabel(str(k) + ' (Hz)')
				if indk == 0:
					pl.title('Channels ' + str(c))

				#Dans CBG FS n'est représenter que par un neurones (indépendant du nombre de competing channels)
				#De même dans THFC le TRN n'est qu'un seul neurones
				if model == 'CBG' or model == 'CBGcustom':
					activity = log[:,ind * NbChannels + indc - (NbChannels - 1)]
				else:
					activity = log[:,ind * NbChannels + indc]
				pl.plot(range(len(activity)),activity)
				pl.ylim((0,pl.ylim()[1]))

	if save:
		pl.savefig(chemin + '.pdf')


	pl.show()



if __name__ == '__main__':

	plot_log(['MSN','FS','STN','GPe','GPi'], [1])
