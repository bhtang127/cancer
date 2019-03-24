# Two Carrying capacities one for symmetric division and one for asymetric division
# -*- coding: utf-8 -*-
#cancerMutations2.py
#Feedback Model
#Last Updated: 4/26/2015, Added feature to print out population of individual type 3 clones

    


import os
from os import path
import argparse
import numpy as np
import scipy as sp
import numpy.matlib as npmat
import time
import multiprocessing as mp
import warnings
from itertools import product
import sklearn.gaussian_process as gp
import sklearn.gaussian_process.kernels as kn
import pyDOE as pyl
import datetime
import jit_aux as aux
warnings.filterwarnings("ignore")


np.random.seed(seed=111)

class tissue:
    
    #parameters:
    #N: wild type population size
    #ngrp: number of mutation groups
    #pmute: probability of mutation
    #rho0: slope parameter for asymetric differentiation
    #delta: update time step for simulation
    
    #pgrp: conditional probability of groups given mutation. If
    #      length=1,all groups are equi-likely 
    #T: maximum lifetime (in weeks)
    #t3Threshold: number of mutations in last type for cancer to
    #             occur.
    #epsilon0: radius of lowest to highest rate intervals for normal
    #          population homeostasis
    #cang: angio-genesis effect (increased resource)
    #capp: apoptosis effect (death rate decrease)
    #cpro: proliferation effect (division increase)
    
    def __init__(self, N = 1, dcell = 1e9, ngrp=18, pmute =1e-6,pmuteLOH=1e-5,pmutesuperLOH=10**(-3),  
                       pgrp = [1], alpha=1., capp=1.5 ,cappInter=1., targ_size=10.,targ_size1h=100., delta=0.5, 
                       T=52*80,lastTime=40,deltaMute=5*1e-6,deltaOnc1=0.001,psym=0.0525,pasymdiff=0.9,psymdiff=0.1-0.0525,
                       deltapsym=1.7e-3,PowerRate=0.7,TargetDiv=0.07,mult=16,powersize=9.301):
        
        self.Tsind=np.append([0,2,3],np.arange(10,18))
        self.Oncind=np.append([1],np.arange(4,10)) 
        self.PathwaysCF=[{0,1,2},{3}]
        self.PathwaysCS=[{4,5,6},{7},{8,9},{10}]
        self.PathwaysGM=[{11},{12},{13},{14},{15},{16},{17}]   
        self.N = int(N)
        self.crypt_size = 1
        self.nb_crypt = self.N /self.crypt_size  
        
        
        self.dcell = int(dcell)
        
        self.ngrp = ngrp
        #self.pchange = pchange
        self.nmut = 0
        
        
        self.alpha = alpha
        self.capp = capp
        self.cappInter=cappInter
        self.sigma_app = float(targ_size)*((self.alpha+1)*targ_size + (self.alpha-1)*self.crypt_size)/((self.alpha-1)*targ_size + (self.alpha+1)*self.crypt_size)
        
        self.lastTime=lastTime
        
        
        self.delta = delta
        self.deltapsym=deltapsym
        if len(pgrp) == 1:
            self.pgrp = np.append(np.ones(ngrp-1)/(ngrp-1),0)
            self.pgrp[4] *=1
            self.pgrp =self.pgrp/self.pgrp.sum()
        else:
            self.pgrp = pgrp
            
        self.pmute = pmute
        self.pmuteLOH=pmuteLOH
        self.pmutesuperLOH=pmutesuperLOH
        self.deltaMute=deltaMute
        self.T = T
        self.PowerRate=PowerRate
        self.TargetDiv=TargetDiv
        self.groupsize = np.array([N] + [0]*self.ngrp, dtype=int)
        self.kchange = .1
        self.maxAge = 500
        self.deltaOnc1=deltaOnc1
        
        self.psym=psym
        self.pasymdiff=pasymdiff
        self.psymdiff=psymdiff
        self.age=0
        self.mult=mult
        self.powersize=powersize
        
        
       
        
        
   
    
    #basic division/death rate with total population feedback
    def hfun(self, n, a, s, ep):
        #a=np.log(1+((s-1)*float(n)/float(s)))/np.log(s)
        #u = float(7)/(3.5+a)
        if (float(10*n)/(np.sqrt(1)*10**self.powersize*s))<1 and (self.age<=20*52):
           #u =self.mult*np.exp(np.log(4/self.mult)/((52*20)**(1))*(20*52)**(1))*(3*np.exp(1)*0.45/3.5*np.exp(-1/(1-(float(10*n)/(np.sqrt(1)*10**self.powersize*s))**(self.PowerRate)))+self.TargetDiv)
           u=np.exp(1)*(1.75-self.TargetDiv)*np.exp(-1/(1-(float(10*n)/(np.sqrt(1)*10**self.powersize*s))**(self.PowerRate)))+self.TargetDiv
        elif (float(10*n)/(10**self.powersize*s))<1 and (self.age>20*52):
           #u=4*(3*np.exp(1)*0.45/3.5*np.exp(-1/(1-(float(10*n)/(np.sqrt(1)*10**self.powersize*s))**(self.PowerRate)))+self.TargetDiv)  
           u=np.exp(1)*(1.75-self.TargetDiv)*np.exp(-1/(1-(float(10*n)/(np.sqrt(1)*10**self.powersize*s))**(self.PowerRate)))+self.TargetDiv
        else:
           u=1*self.TargetDiv
        
        
            #a * self.crypt_size + b* float(n))/(self.crypt_size + float(n))
            # s is the target size
            # n is 
        return u
    
    
        
    

    #
    def CleanList(self,L):
       ind=0
       l=len(L)
       
       for i in range(l):
           if ((L[ind].sum())==L[ind][0,0,0]):
               L.pop(ind)
           else:
               ind +=1
       return L
    
    
    def getPreRateGeneric(self):
    
        Rates=np.zeros((3,3,2,4))
        AsymRates=np.zeros((3,3,2))
        deltaProb=0.01/2
        self.psymdiff=self.psym-deltaProb
        psym0=self.psym
        psymdiff0=self.psymdiff
        pasymdiff0=self.pasymdiff       
        tauDiv0=1
        tauDeathc=(psym0-psymdiff0)
        pmute0=self.pmute
        
        for i,j,k in product(range(3),range(3),range(2)):
            psym=psym0+i*self.deltapsym
            psymdiff=psymdiff0-(1*(i==1)+(2)*(i==2))*self.deltapsym
            pasymdiff=pasymdiff0
            tauDiv=tauDiv0*(self.capp**(1*(j==1)+(2)*(j==2)))
            pmute=pmute0+k*self.deltaMute            
            Rates[i,j,k,0]=tauDeathc
            Rates[i,j,k,1]=tauDiv*psymdiff    
            Rates[i,j,k,2]=tauDiv*psym*(1-pmute)
            Rates[i,j,k,3]=tauDiv*pasymdiff*float(pmute)/2+tauDiv*psym*float(pmute)
            AsymRates[i,j,k]=tauDiv*pasymdiff
            
           
        self.Rates=self.delta*Rates
        self.AsymRates=self.delta*AsymRates
        return self.delta*Rates
    
    

        
    
    def GenericReactions(self):
        GenericRates=np.zeros((3,3,2,5))    
        for i,j,k in product(range(3),range(3),range(2)):
            GenericRates[i,j,k,0]=self.Rates[i,j,k,0]+self.Rates[i,j,k,1]
            GenericRates[i,j,k,1]=self.Rates[i,j,k,2]
            GenericRates[i,j,k,2]=self.Rates[i,j,k,3]*(4*10**(-3))*(1/(1+i))
            GenericRates[i,j,k,3]=self.Rates[i,j,k,3]*(7*10**(-3))*(1/(1+j))       
            GenericRates[i,j,k,4]=self.Rates[i,j,k,3]*(7*10**(-3))                

        self.ReactionsRate=GenericRates
        return GenericRates
 
    def TypesIndices(self):
        Type1=[[],[],[]]
        Type2=[[],[],[]]
        Type3=[[],[],[]]
        for i,j,k in product(range(3),range(3),range(2)):
            if (i+j+k)==1:
                Type1[0].append(i)
                Type1[1].append(j)
                Type1[2].append(k)
            elif  (i+j+k)==2:
                Type2[0].append(i)
                Type2[1].append(j)
                Type2[2].append(k)
            elif (i+j+k)>=3:
                Type3[0].append(i)
                Type3[1].append(j)
                Type3[2].append(k)
        self.IndType1=Type1
        self.IndType2=Type2
        self.IndType3=Type3
            
                
            
        
    def Likelihoods(self,D,Y):

        return aux.Likelihoods(D,Y)
    
     
    def MaxLikelihood(self,K,Y,tol):
        
        return aux.MaxLikelihood(K,Y,tol)
    
    
    def EPApprox(self,K,Y,tol):
        
        return aux.EPApprox(K,Y,tol)
               

    def PosteriorMeanVariance(self,K,Y,tol,x,X,C,l):  
        
        return aux.PosteriorMeanVariance(K,Y,tol,x,X,C,l)
    

    def FitGradient(self,L_K,Y,tol,X,C,l,Const,LearningRate,iterations):
        
        return aux.FitGradient(L_K,Y,tol,X,C,l,Const,LearningRate,iterations)
            
    
    def ChooseBestParam(self,Y,tol,X,X_new,C,l,Const):
        
        return aux.ChooseBestParam(Y,tol,X,X_new,C,l,Const)
                


    # computes probability of group hit given mutation    
    # computes probability of group hit given mutation    
    
        
    
    def CancerTest(self,Tabmutations):
        
        Type1=np.any(Tabmutations[self.IndType1[0],self.IndType1[1],self.IndType1[2],:]>100,axis=0)
        Type2=np.any(Tabmutations[self.IndType2[0],self.IndType2[1],self.IndType2[2],:]>100,axis=0)
        Type3=np.any(Tabmutations[self.IndType3[0],self.IndType3[1],self.IndType3[2],:]>100,axis=0)
        
        #self.Plastic=np.logical_and(self.CrTot>3*10**8,self.CrTot/self.CrSizes<10)
        #WinType=np.logical_and(Type3,self.Plastic)      
        return np.any(Type3),Type3,Type2,Type1
        
    
    
    def TestCancerClonal(self):
        self.TypesIndices()
        Time=time.time()
        t=0
        self.age=20*52
        self.getPreRateGeneric()
        self.GenericReactions()
        self.age=0
        L=set()
        TabMutations=np.zeros((3,3,2,int(1.1*10**7)))
        TabDifferentiated=np.zeros((7,int(1.1*10**7)))
        TabHistory=5000*np.ones((3,int(1.1*10**7)))
        Tmutations=np.zeros((2,2,2,2))
        self.Diagnostic=False
        self.status=False
        self.CrTot45=np.array([10**11, 10**11])
        self.Nold=0
        self.N=0
        while t<self.T and not(self.status):
              self.Nold=self.N
              
              
              self.Diagnostic,self.Type3,self.Type2,self.Type1=self.CancerTest(TabMutations[:,:,:,list(L)])
              
              
              IndType1=np.where(self.Type1)[0]
              IndType2=np.where(self.Type2)[0]
              IndType3=np.where(self.Type3)[0]
              
              if len(IndType1)>0:
                 
                 TabHistory[0,np.array(list(L))[IndType1]]=np.minimum(TabHistory[0,np.array(list(L))[IndType1]],t)
                 
              if len(IndType2)>0:
                 TabHistory[1,np.array(list(L))[IndType2]]=np.minimum(TabHistory[1,np.array(list(L))[IndType2]],t)
                 
              if len(IndType3)>0:
                 TabHistory[2,np.array(list(L))[IndType3]]=np.minimum(TabHistory[2,np.array(list(L))[IndType3]],t)
                 
              if t<36:
                self.N=np.ceil((1+10*1388*2*t)/.1)
              elif (t>=36) and (t<(20*52+36)):
                self.N=np.ceil((1+10*1388*72+10*463*2*(t-36))/.1)
              else:
                self.N=np.ceil((1+10*1388*72+10*463*2*(52*20))/.1)
              
                
              
              self.nb_crypt=np.ceil(self.N/10)
              self.nb_cryptOld=np.ceil(self.Nold/10)
              #print(self.nb_crypt-self.nb_cryptOld)
              
              self.crypt_size=int(np.ceil(self.N/self.nb_crypt))
              

              #print(int(self.nb_crypt-self.nb_cryptOld))
              if t<=(20*52+36):
                 self.IndInheritence=np.intersect1d(np.random.randint(0,self.nb_cryptOld+1,int(self.nb_crypt-self.nb_cryptOld)),np.array(list(L)))
                 if len(self.IndInheritence)>0:
                    #print(len(self.IndInheritence))
                    TabMutations[:,:,:,np.arange(int(self.nb_cryptOld),int(self.nb_cryptOld)+len(self.IndInheritence))]=TabMutations[:,:,:,self.IndInheritence]
                    L=L.union(np.arange(int(self.nb_cryptOld),int(self.nb_cryptOld)+len(self.IndInheritence)))
                 
                          
              
              
              #print(t)
              if t<30*52:
              
              
                 N=np.random.poisson(self.N*self.hfun(10,1,10,1)*self.ReactionsRate[0,0,0,2]*np.array([1,7/4,1]))
              
              
                 A0=np.random.randint(0,self.nb_crypt+1,N[0])
                 A1=np.random.randint(0,self.nb_crypt+1,N[1])
                 A2=np.random.randint(0,self.nb_crypt+1,N[2])
              
              
              
                 unique0,count0=np.unique(A0,return_counts=True)
                 unique1,count1=np.unique(A1,return_counts=True)
                 unique2,count2=np.unique(A2,return_counts=True)
               
                 TabMutations[1,0,0,unique0] +=count0
                 TabMutations[0,1,0,unique1] +=count1 
                 TabMutations[0,0,1,unique2] +=count2
              
                 L =L.union(set(unique0),set(unique1),set(unique2))
              

              
              
              if len(L)>0:
                  
                  
                  
                  
                  
                  ind=list(L)
                  
                  Tmutations=np.copy(TabMutations[:,:,:,(ind)])
                  
                  TDifferentiated=np.copy(TabDifferentiated[:,ind])
                  
                  TDifferentiated=np.roll(TDifferentiated,1,axis=0)
                  
                  
                  
                  CrSizes=np.sum(Tmutations,axis=(0,1,2))
                  
                  self.CrSizes=CrSizes
                  
                  #print(CrSizes)
                  DivRates=np.zeros(len(L))
                  NumberDiff=np.dot(2**np.arange(7),TabDifferentiated[:,list(L)])
                  
                  for k in range(len(L)):
                      
                      DivRates[k]=self.hfun(CrSizes[k]+NumberDiff[k],1,10,1)   
                  
                  
                  
                  
                    
                  
                      
                  ReactionsMutants0=np.tile(np.reshape(self.ReactionsRate[:,:,:,2],(3,3,2,1)),(1,1,1,len(L)))*Tmutations
                                           
                  ReactionsMutants1=np.tile(np.reshape(self.ReactionsRate[:,:,:,3],(3,3,2,1)),(1,1,1,len(L)))*Tmutations
                                           
                  ReactionsMutants2=np.tile(np.reshape(self.ReactionsRate[:,:,:,4],(3,3,2,1)),(1,1,1,len(L)))*Tmutations
                                           
                  DivRatesTable=np.tile(np.reshape(DivRates,(1,1,1,len(L))),(3,3,2,1))
                  
                  ReactionsMutants0 *=DivRatesTable
                  
                  ReactionsMutants1 *=DivRatesTable
                  
                  ReactionsMutants2 *=DivRatesTable
                  
                  ReactionsDeath1=Tmutations*DivRatesTable*np.tile(np.reshape(self.Rates[:,:,:,0],(3,3,2,1)),(1,1,1,len(L)))
                  ReactionsDeath2=Tmutations*DivRatesTable*np.tile(np.reshape(self.Rates[:,:,:,1],(3,3,2,1)),(1,1,1,len(L)))
                  
                  ReactionsReplicate=Tmutations*DivRatesTable*np.tile(np.reshape(self.Rates[:,:,:,2],(3,3,2,1)),(1,1,1,len(L)))
                  
                  #print("d1 is: {0}".format((np.tile(np.reshape(self.Rates[:,:,:,0],(3,3,2,1)),(1,1,1,len(L)))+DivRatesTable*np.tile(np.reshape(self.Rates[:,:,:,1],(3,3,2,1)),(1,1,1,len(L)))).mean()))
                  
                  #print("d2 is: {0}".format((DivRatesTable*np.tile(np.reshape(self.Rates[:,:,:,2],(3,3,2,1)),(1,1,1,len(L)))).mean()))
                  
                  ReactionsDiff1=Tmutations*DivRatesTable*np.tile(np.reshape(self.AsymRates,(3,3,2,1)),(1,1,1,len(L)))
                  
                  #print(np.max(ReactionsDiff1))
                  
                  
                  PoissonReactionsMutants0=np.random.poisson(ReactionsMutants0)
                  
                  PoissonReactionsMutants1=np.random.poisson(ReactionsMutants1)
                  
                  PoissonReactionsMutants2=np.random.poisson(ReactionsMutants2)
                  
                  PoissonReactionsDeath1=np.random.poisson(ReactionsDeath1)
                  PoissonReactionsDeath2=np.random.poisson(ReactionsDeath2)
                  PoissonReactionsDeath=PoissonReactionsDeath1+PoissonReactionsDeath2
                  
                  
                  
                  
                  
                  PoissonReactionsReplicate=np.random.poisson(ReactionsReplicate)
                  
                  
                  
                  #print(PoissonReactionsReplicate.mean())
                  
                  PoissonReactionsDiff=np.random.poisson(ReactionsDiff1)+2*PoissonReactionsDeath2
                                                        
                  TDifferentiated[0,:]=np.sum(PoissonReactionsDiff,axis=(0,1,2))                                     
                  
                  
                  MutantsUpdate0=np.roll(PoissonReactionsMutants0,1,axis=0)
                  MutantsUpdate0[1,:,:,:] +=MutantsUpdate0[0,:,:,:]
                  MutantsUpdate0[0,:,:,:] = 0
                                
                  MutantsUpdate1=np.roll(PoissonReactionsMutants1,1,axis=1)
                  MutantsUpdate1[:,1,:,:] +=MutantsUpdate1[:,0,:,:]
                  MutantsUpdate1[:,0,:,:] = 0
                                
                  MutantsUpdate2=np.roll(PoissonReactionsMutants2,1,axis=2)
                  MutantsUpdate2[:,:,1,:] +=MutantsUpdate2[:,:,0,:]
                  MutantsUpdate2[:,:,0,:] = 0
                  
                  Tmutations +=(MutantsUpdate0+MutantsUpdate1+MutantsUpdate2+PoissonReactionsReplicate-PoissonReactionsDeath)
                  #print(Tmutations.shape)
                  Tmutations=np.maximum(Tmutations,np.zeros((3,3,2,len(L))))
                  #print(Tmutations.shape)
                  TabMutations[:,:,:,(ind)]=np.copy(Tmutations)
                  TabDifferentiated[:,ind]=np.copy(TDifferentiated)        
                  ind=np.array(ind)            
                  ind=np.copy(ind[np.sum(Tmutations,axis=(0,1,2))>0])
                  
                  
                  
                  L=set(ind)
                  
                  if t==(52*55):             
                     self.TabMutations45=TabMutations[:,:,:,list(L)]
                     self.CrSizes45=np.sum(self.TabMutations45,axis=(0,1,2))+self.crypt_size                                  
                     self.TabDifferentiated45=TabDifferentiated[:,list(L)]
                     self.DiffCrSizes45=np.dot(2**np.arange(7),self.TabDifferentiated45)+(self.crypt_size)*(2**7-1)
                     self.CrTot45=self.DiffCrSizes45+self.CrSizes45
                  
                  
                  
                  
                  
                  
              if self.Diagnostic:
                      self.TabMutations=TabMutations[:,:,:,list(L)]
                      self.CrSizes=np.sum(self.TabMutations,axis=(0,1,2))+self.crypt_size                                  
                      self.TabDifferentiated=TabDifferentiated[:,list(L)]
                      self.DiffCrSizes=np.dot(2**np.arange(7),self.TabDifferentiated)+(self.crypt_size)*(2**7-1)
                      self.CrTot=self.DiffCrSizes+self.CrSizes
                      self.TabHistory=TabHistory[:,list(L)]
                      
                      self.Plastic=np.logical_and(self.CrTot>10*10**8,self.CrTot/self.CrSizes<5000)
                      self.status=np.any(np.logical_and(self.CancerTest(self.TabMutations)[1],self.Plastic))
              t +=self.delta
              self.age=t
                       
                  
        self.Diagnostic,self.Type3,self.Type2,self.Type1=self.CancerTest(TabMutations[:,:,:,list(L)])
        
        
        self.TabMutations=TabMutations[:,:,:,list(L)]
        self.CrSizes=np.sum(self.TabMutations,axis=(0,1,2))+self.crypt_size                                  
        self.TabDifferentiated=TabDifferentiated[:,list(L)]
        self.DiffCrSizes=np.dot(2**np.arange(7),self.TabDifferentiated)+(self.crypt_size)*(2**7-1)
        self.CrTot=self.DiffCrSizes+self.CrSizes
        self.TabHistory=TabHistory[:,list(L)]                 
        self.L=L
        #print (np.min(self.CrTot/self.CrSizes))
              
        if self.Diagnostic:
                      self.TabMutations=TabMutations[:,:,:,list(L)]
                      self.CrSizes=np.sum(self.TabMutations,axis=(0,1,2))+self.crypt_size                                  
                      self.TabDifferentiated=TabDifferentiated[:,list(L)]
                      self.DiffCrSizes=np.dot(2**np.arange(7),self.TabDifferentiated)+(self.crypt_size)*(2**7-1)
                      self.CrTot=self.DiffCrSizes+self.CrSizes
                      self.TabHistory=TabHistory[:,list(L)]
                      
                      self.Plastic=np.logical_and(self.CrTot>10*10**8,self.CrTot/self.CrSizes<5000)
                      self.status=np.any(np.logical_and(self.CancerTest(self.TabMutations)[1],self.Plastic))
                  
        
        
#        self.Plastic=np.logical_and(self.CrTot>3*10**8,self.CrTot/self.CrSizes<10)
        self.ExitTime=t
        print (time.time()-Time) 



ExitTimes=list()
L=list() 
C=list() 
C1=list()
C45=list()
TH=list()
TM=list()
Ratio=list()
NPolyps1=list()
NPolyps145=list()
NPolyps2=list()
NPolyps245=list()
NPolyps3=list()
NPolyps345=list()


 

cases=0
TH_oneHit=list()
TH_twoHits=list()  
Label_oneHit=list()
Label_twoHits=list()

#
#
#

def mutationType(TM):
    types = ["F","S","M"]
    mtype = ""
    cur_count = 0
    cur_ind = np.zeros(3,dtype=np.int)
    for i,j,k in product(range(2),range(2),range(2)):
        if i + j + k == 1 and TM[i,j,k] > cur_count:
            cur_count = TM[i,j,k]
            cur_ind[0] = i; cur_ind[1] = j; cur_ind[2] = k
    if cur_count > 0:
        mtype += types[np.where(cur_ind > 0)[0][0]]
        cur_count = 0
        cur_ind_2 = np.zeros(3,dtype=np.int)
        for i,j,k in product(range(2),range(2),range(2)):
            if i + j + k == 1 and TM[i+cur_ind[0],j+cur_ind[1],k+cur_ind[2]] > cur_count:
                cur_count = TM[i+cur_ind[0],j+cur_ind[1],k+cur_ind[2]]
                cur_ind_2[0] = i; cur_ind_2[1] = j; cur_ind_2[2] = k
        if cur_count > 0:
            mtype += types[np.where(cur_ind_2 > 0)[0][0]]
            cur_count = 0
            cur_ind += cur_ind_2
            cur_ind_3 = np.zeros(3,dtype=np.int)
            for i,j,k in product(range(2),range(2),range(2)):
                if i + j + k == 1 and \
                   i+cur_ind[0] < 3 and \
                   j+cur_ind[1] < 3  and\
                   k+cur_ind[2] < 2 and \
                   TM[i+cur_ind[0],j+cur_ind[1],k+cur_ind[2]] > cur_count:
                    cur_count = TM[i+cur_ind[0],j+cur_ind[1],k+cur_ind[2]]
                    cur_ind_3[0] = i; cur_ind_3[1] = j; cur_ind_3[2] = k
            if cur_count > 0:
                mtype += types[np.where(cur_ind_3 > 0)[0][0]]
    return mtype

#
#
XSTAR=np.array([ 1.14794853 , 2.12     ,     5.0148502 ,  0.41173119])

Nindividuals=2500
print("Total iteration: ", Nindividuals)
cases=0
#

f = open('result.csv', 'w')

print("Individual, ExitTime, Minimum_ratio, Status, Size_largest_crypt, Size_largest_crypt_50, \
       Number_polyps1, Number_polyps1_50, Number_polyps2, Number_polyps2_50, \
       Number_polyps3, Number_polyps3_50, HitTime1, HitType1, HitTime2, HitType2, HitTime3, HitType3,\
       Case_Number", file=f)

for i in range(Nindividuals):
    print("Individual : ",i+1)
    print(i+1, end=", ", file=f)
    t1=tissue(N = 1, dcell = 1e9, ngrp=18, pmute =1e-6,pmuteLOH=1e-5,pmutesuperLOH=10**(-3),  
              pgrp = [1], alpha=1., cappInter=1., targ_size=10.,targ_size1h=100., delta=0.5,
              T=52*80,lastTime=40,deltaMute=5*1e-6,deltaOnc1=0.001,
              capp=XSTAR[0],psym=0.0525,pasymdiff=0.9,psymdiff=0.1-0.0525,deltapsym=XSTAR[1]*1e-3,
              PowerRate=XSTAR[3],TargetDiv=XSTAR[2]*0.03,mult=1,powersize=9.301)    
    t1.TestCancerClonal()
    print("{0:.1f}".format(t1.ExitTime), end=", ", file=f)
    ExitTimes.append(t1.ExitTime)

    print ("{0:.5f}".format((t1.CrTot[np.where(t1.CrSizes==max(t1.CrSizes))]/t1.CrSizes[np.where(t1.CrSizes==max(t1.CrSizes))])[0]), end=", ", file=f)
    Ratio.append(t1.CrTot[np.where(t1.CrSizes==max(t1.CrSizes))]/t1.CrSizes[np.where(t1.CrSizes==max(t1.CrSizes))])
    L.append(t1.status)
    print(t1.status)
    print(t1.status, end=", ", file=f)
    C.append((t1.CrTot).max())
    print("{0:.5f}".format(np.log((t1.CrTot).max())/np.log(10)), end=", ", file=f)
    C1.append(t1.CrSizes.max())
    print("{0:.5f}".format(np.log((t1.CrTot45).max())/np.log(10)), end=", ", file=f)
    C45.append(t1.CrTot45.max())
    
    
    NPolyps1.append(sum(t1.CrTot>10*10**8))
    print("{0:d}".format(sum(t1.CrTot>10*10**8)), end=", ", file=f)
    NPolyps145.append(sum(t1.CrTot45>10*10**8))
    print("{0:d}".format(sum(t1.CrTot45>10*10**8)), end=", ", file=f)
    NPolyps2.append(sum(t1.CrTot>10**9))
    print("{0:d}".format(sum(t1.CrTot>10**9)), end=", ", file=f)
    NPolyps245.append(sum(t1.CrTot45>10**9))
    print("{0:d}".format(sum(t1.CrTot45>10**9)), end=", ", file=f)
    NPolyps3.append(sum(t1.CrTot>3*10**9))
    print("{0:d}".format(sum(t1.CrTot>3*10**9)), end=", ", file=f)
    NPolyps345.append(sum(t1.CrTot45>3*10**9))
    print("{0:d}".format(sum(t1.CrTot45>3*10**9)), end=", ", file=f)
    
    
    HittingTimes1=np.where(np.logical_and(t1.Type1,t1.CrTot>10*10**8))[0]
    if HittingTimes1 !=[]:
       TH1 = t1.TabHistory[:,HittingTimes1[0]]
       TM1 = t1.TabMutations[:,:,:,HittingTimes1[0]]
       print("{0:.2f}".format(TH1[0]), end=", ", file=f)
       print("{0}".format(mutationType(TM1)[0]), end=", ", file=f)
    else:
        print(" ", end=", ", file=f)
        print(" ", end=", ", file=f)
       
       
    HittingTimes2=np.where(np.logical_and(t1.Type2,t1.CrTot>10*10**8))[0]
    if HittingTimes2 !=[]:
       TH2 = t1.TabHistory[:,HittingTimes2[0]]
       TM2 = t1.TabMutations[:,:,:,HittingTimes2[0]]
       print("{0:.2f}".format(TH2[1]), end=", ", file=f)
       print("{0}".format(mutationType(TM2)[1]), end=", ", file=f)
    else:
        print(" ", end=", ", file=f)
        print(" ", end=", ", file=f)   
    
    try:
        HittingTimes=np.where(np.logical_and(t1.Type3,t1.Plastic))[0]
        TH3 = t1.TabHistory[:,HittingTimes[0]]
        TM3 = t1.TabMutations[:,:,:,HittingTimes[0]] 
        print("{0:.2f}".format(TH3[2]), end=", ", file=f)
        print("{0}".format(mutationType(TM3)[2]), end=", ", file=f)
    except:
        print(" ", end=", ", file=f)
        print(" ", end=", ", file=f)

    if t1.status:
        print("{0:d}".format(cases+1), end=", ", file=f)
        cases +=1
    else:
        print(" ", end=", ", file=f)

    if i < Nindividuals-1:
        del t1
    print("", file=f)


f.close()

np.save("TM.npy", np.array(TM))



#Comments for Bohao: 
    

#The number of individuals simulates is in line 1198.

#Parameters aree already calibrated.

#The Times of cancers (in weeks) for each individual are stored in the list "ExitTimes". If the value is 4160, then there is no cancer. The person arrived to 80 years old without cancer.

# The size of the largest crypt at age 55 for each individual is stored in C45 (called 45 just because we started with 45 insstead of 55).

#The size of the larget crypt at age 80 is stored in the list C.

#The number of cancer are stored in the variables "cases".

#Given cancer, the timings of the hits (3 hits to get to cancer) are stored in the list TH.

#The types of Hits are stored in the list TM: TM[i] gives the types of hits for cancer number i