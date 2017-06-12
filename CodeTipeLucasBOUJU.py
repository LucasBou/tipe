#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lucas BOUJU
"""

import math 
import scipy
import numpy as np
from scipy.stats import norm

def normeinf(L):#norme infinie
    return max(list(map(abs,L)))

def mediane(L):#mediane
    S=sorted(L)
    l=len(L)
    if l%2==1:
        return S[int((l-1)/2)]
    else:
        return (S[int(l/2)]+S[int((l/2)-1)])/2
    
def moyenne(L):#moyenne
    return somme(L)/len(L)

def somme(L): #somme les elmts d'une liste
    s=0
    for k in L:
        s+=k
    return s

def calc1(n):#utilise dans frequencybloc
    return (n-0.5)**2

def mod(a,b):#utilise dans DFT
    return math.sqrt(a**2+b**2)

def cdf(x):#utilise dans cusum
    return norm.cdf(x)

def igamc(a,x):#utilise en test d'hypothese
    return scipy.special.gammaincc(a,x)

def chidist(chicarre,df):#utilise en test d'hypothese
    return scipy.special.gammaincc(df/2,(chicarre)/2)

def coeffFourierc(j,L):
    c=0
    n=len(L)
    for k in range (n):
        c+=L[k]*math.cos(2*k*j*math.pi/n)
    return c

def coeffFouriers(j,L):
    s=0
    n=len(L)
    for k in range (n):
        s+=L[k]*math.sin(2*k*j*math.pi/n)
    return s

def decompoFourierModule(L):#liste des n/2 premier coeffs de fourier (cf redondance)
    return [mod(coeffFouriers(j,L),coeffFourierc(j,L)) for j in range(int(len(L)/2))]

def variations(L):
    return [L[k+1]-L[k] for k in range(len(L)-1)]

def indicePremiereOccu(x,L):
    for k in range(len(L)):
        if L[k]==x:
            return k
    return "pas d'occurence de elmt dans la liste"

def listToBinaire(L):#utilise dans serial, bijection entre les séquences et leur emplacement de stockage de fréquence dans la liste
    s=0
    for k in range(len(L)):
        s+=L[k]*(2**k)
    return s

def longest(p):#trouve la taille de la plus grande chaine de 1 dans la liste p
    maxi=0
    compteur=0
    for k in p:
        if k==0:
            if compteur > maxi :
                maxi = compteur
            compteur = 0
        else:
            compteur+=1
    if compteur > maxi :
        return compteur
    return maxi


###transformation des donnees brutes en liste


def convertionbitmediane (L): #met des 1 pour toutes les valeurs au dessus de la mediane et 0 sinon, renvoie une nouvelle liste, marche bien pour des listes de nombres différents (cf si plusieurs pareils au niveau de la mediane on se retrouve pas avec autant de 0 que de 1) et de cardinal pair
    med=mediane(L)
    S=[]
    for k in L :
        if k>med :
            S.append(1)
        else :
            S.append(0)
    return S

def convertionbitmoyenne(L):#de meme mais pour la moyenne
    m=moyenne(L)
    for k in range(len(L)):
        if L[k]<m:
            L[k]=0
        else:
            L[k]=1
    return L

def listto2darray(L):#transforme une liste de 0 et 1 en array pret a etre transformé en image noir blanc
    a=int(math.sqrt(len(L)))+1
    A=[[2 for i in range(a)]for j in range (a)]
    s=0
    for i in range(a):#double boucle qui construit la matrice des o et 1
        for j in range (a):
            if s<len(L):
                A[i][j]=L[s]
                s+=1
    for i in range(len(A)):#double boucle qui la transforme en couleur
        for j in range(len(A)):
            if A[i][j]==0:
                A[i][j]=[0,0,0]
            else:
                A[i][j]=[255,255,255]
    return np.array(A)

def fichierToList(nomfichier):#format en fichierrandom.txt
    monfichier=open(nomfichier,"r")
    contenu=monfichier.read()
    monfichier.close()
    L=contenu.split("\n")
    L.pop()
    L.pop()
    for k in range(len(L)):
        L[k]=L[k].replace(',','.')
    X=list(map(float,L))
    return X

def randomOrgToList(nomfichier):
    monfichier=open(nomfichier,"r")
    contenu=monfichier.read()
    monfichier.close()
    contenu=contenu.replace(" ","")
    contenu=contenu.replace("\n","")
    X=[]
    for k in contenu:
        X.append(int(k))
    return X


###debut des tests


def frequencybloc (L,M):
    n=len(L)
    Lblocks=[]
    j=0
    for i in range(int(n/M)):#création de la liste des blocs
        L3=[]
        for k in range(M):
            L3.append(L[j+k])
        j+=M
        Lblocks.append(L3)
    listemoyennes=list(map(moyenne,Lblocks))
    print(listemoyennes)
    khicarre=4*M*somme(list(map(calc1,listemoyennes)))
    print(khicarre)
    df=int(n/M)-1
    print(df)
    pval=chidist(khicarre,df)
    return pval

def longestRunBlock(L):#pris pour M=128 donc K=5
    M=128
    n=len(L)
    Lblocks=[]
    j=0
    for i in range(int(n/M)):#création de la liste des blocs
        L3=[]
        for k in range(M):
            L3.append(L[j+k])
        j+=M
        Lblocks.append(L3)
    K=5
    PiAttendus=[0.1174,0.243,0.2493,0.1752,0.1027,0.1124]
    FreqObs=[0,0,0,0,0,0]
    for p in Lblocks:
        mu=longest(p)
        if mu <= 4:
            FreqObs[0]+=1
        elif mu == 5:
            FreqObs[1]+=1
        elif mu == 6:
            FreqObs[2]+=1
        elif mu == 7:
            FreqObs[3]+=1
        elif mu == 8:
            FreqObs[4]+=1
        else:
            FreqObs[5]+=1
    N=int(n/M)
    aux=[((FreqObs[i]-N*PiAttendus[i])**2)/(N*PiAttendus[i]) for i in range (K+1)]
    chicarre=somme(aux)
    pval=chidist(chicarre,K)
    return pval

def DFTtest(L):#utiliser sur une liste de cardinal pair
    n=len(L)
    X=list(map(lambda i:2*i-1,L))
    T=math.sqrt(2.995732274*n)
    N0=n*0.475
    N1=0
    modS=decompoFourierModule(X)
    for k in modS:
        if k<T:
            N1+=1
    d=(N1-N0)/(math.sqrt(n*0.95*0.05/2))
    pval=math.erfc(abs(d)/math.sqrt(2))
    return pval

def serialTest(L,m):#marche
    n=len(L)
    X=L+L[0:m-1]
    mum=[0 for k in range(2**m)]
    mum1=[0 for k in range(2**(m-1))]
    mum2=[0 for k in range (2**(m-2))]
    for k in range (n):#on gere les trois fenetres en meme temps
        mum[listToBinaire(X[k:m+k])]+=1  
        mum1[listToBinaire(X[k:m+k-1])]+=1
        mum2[listToBinaire(X[k:m+k-2])]+=1
    psicarrem=0
    psicarrem1=0
    psicarrem2=0
    for k in range (2**(m-2)):#de meme
        psicarrem+=(mum[k]-(n/(2**m)))**2
        psicarrem1+=(mum1[k]-(n/(2**(m-1))))**2
        psicarrem2+=(mum2[k]-(n/(2**(m-2))))**2
    for k in range (2**(m-2),2**(m-1)):
        psicarrem+=(mum[k]-(n/(2**m)))**2
        psicarrem1+=(mum1[k]-(n/(2**(m-1))))**2
    for k in range (2**(m-1),2**m):
        psicarrem+=(mum[k]-(n/(2**m)))**2
    psicarrem=psicarrem*(2**(m)/n)
    psicarrem1=psicarrem1*(2**(m-1)/n)
    psicarrem2=psicarrem2*(2**(m-2)/n)
    delta1=psicarrem-psicarrem1
    delta2=delta1-psicarrem1+psicarrem2
    return [igamc(2**(m-2),delta1/2),igamc(2**(m-3),delta2/2)]

def cusumTest(L,mode):
    X=list(map(lambda i:(2*i)-1,L))
    n=len(L)
    if mode==1:
        L.reverse()
    S=[X[0]]
    for k in range (1,n):
        S.append(S[k-1]+X[k])
    z=normeinf(S)
    somme1=0
    somme2=0
    for k in range(int(((-n/z)+1)/4),1+int(((n/z)-1)/4)):
        somme1+=cdf((z*(4*k+1))/math.sqrt(n))-cdf((z*(4*k-1))/math.sqrt(n))
    for k in range(int(((-n/z)-3)/4),1+int(((n/z)-1)/4)):
        somme2+=cdf((z*(4*k+3))/math.sqrt(n))-cdf((z*(4*k+1))/math.sqrt(n))
    pval=1-somme1+somme2
    return pval

def resume(L):
    print("Frequency Block", frequencybloc (L,500))#prendre 0,05n pour M
    print("Longest run of ones", longestRunBlock(L))
    print("Discrete Fourier Transform", DFTtest(L))
    print("Serial", serialTest(L,7))#m<log2(n)-2, pour 10 000, m=7 convient
    print("Cumulative Sums", cusumTest(L,0))
    return 0