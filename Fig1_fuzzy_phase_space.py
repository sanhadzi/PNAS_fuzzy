
"""
Created on Wed Aug 28 14:01:27 2019

@author: sanhadzi
"""
import time
from numpy import loadtxt, savetxt, array, transpose, hsplit, log, exp, diff, delete, matrix, arange, linspace, rint,arange, exp, array, log, zeros, insert, where, copy, column_stack,concatenate
from scipy.special import factorial
from sympy import Function, symbols, lambdify, diff, simplify, Matrix, expand, Mul
import matplotlib.pyplot as plt
from itertools import combinations 

start_time = time.time()
w, x, y= symbols("w x y")

def karkoli(vrednost, lista):
    a=False
    for z in lista:
        if z==vrednost:
            a=True
    return a

def VSOTA(sekvenca):
    vse_Z=0  ; st_odv_x=[] ; st_odv_y=[] ; v=0.048
    
    M3x = Matrix( [[w*x,v,0], [0,0,1], [v,v,1]]) 
    M3y = Matrix( [[w*y,v,0], [0,0,1], [v,v,1]]) 
    Q13= Matrix( [[0,0,1]] ) ;     COL_vektor3= Matrix( [[0], [1], [1]])
    M3 = Matrix( [[w,v,0], [0,0,1], [v,v,1]]) 
    
    produkt=1 ; x_y_xy=0
    for i in sekvenca:
        if i==0:
            M= M3
        if i==1:
            M= M3x ; st_odv_x.append(1)
        if i==2:
            M= M3y ; st_odv_y.append(1)
        produkt=M*produkt

    Z3= Q13*produkt*COL_vektor3 ; Z3=Z3[0]  ; NOx=sum(st_odv_x) ; NOy=sum(st_odv_y) 

    if NOx>0 and NOy>0:  
        dZ3_1=diff(Z3,x,NOx) ; dZ3_oba=diff(dZ3_1,y,NOy) ; Z3x_oba=(dZ3_oba*y**(NOy)*x**(NOx))/(factorial(NOy)*factorial(NOx)) 
        K3=expand(Z3x_oba) ; vse_Z=vse_Z+K3 ; x_y_xy=3
        
    if NOx>0 and NOy==0: 
        dZ3_1=diff(Z3,x,NOx)  ; Z3_x=(dZ3_1*x**NOx)/(factorial(NOx))
        K3=expand(Z3_x) ; vse_Z=vse_Z+K3 ; x_y_xy=1
        
    if NOx==0 and NOy>0: 
        dZ3_2=diff(Z3,y,NOy)  ; Z3_y=(dZ3_2*y**NOy)/(factorial(NOy))
        K3=expand(Z3_y) ; vse_Z=vse_Z+K3 ; x_y_xy=2
    
    return Z3, vse_Z, x_y_xy

def delckanje(Z, deli):
    j=0 ; suma1=0 ; suma2=0 ; suma3=0 ; suma4=0
    for i in (Z).atoms(Mul):
        if j<deli:
            suma1=suma1+i
        if j>=deli and j<(2*deli):
            suma2=suma2+i
        if j>=(2*deli) and j<(3*deli):
            suma3=suma3+i
        if j>=(3*deli):
            suma4=suma4+i
        j=j+1
    return suma1, suma2, suma3 ,suma4 


# SESTAVITI SET WT SEKVENC
N=20 ; sekvenca=zeros(N) ; sredina=int(rint(len(sekvenca)/2)) ; pipa=sekvenca[(sredina)]=1 ; no_hotspots=N-3
sekvence=[]; k=0 
for i in range(no_hotspots):
    if i % 2 == 0:
        k=k+1 
        sekvenca[(sredina-k)]=2
        if i>0:
            sekvenca[(sredina-(k-1))]=0
    else:
        sekvenca[(sredina+k)]=1
        sekvenca[(sredina+(k-1))]=0
    pipa=copy(sekvenca)
    sekvence.append(pipa)


ENTROPIJE=[] ; HELICNOSTI=[] ; SIGME=[] ; SEPARACIJE=[] ; POPULACIJE1=[] ; POPULACIJE2=[] ; CELI_dol=[] ; CELI_kr=[] ; POL=[] ; DGCH=[] ; AFINITETE=[]
ENTROPIJE_2=[] ; HELICNOSTI_2=[] ; SIGME_2=[] ; SEPARACIJE_2=[] ; POPULACIJE1_2=[] ; POPULACIJE2_2=[] ; CELI_dol_2=[] ; CELI_kr_2=[] ; POL_2=[] ; DG_INT=[] ; AFINITETE_2=[]

a=len(sekvence) ; b=a/2

print a, len(sekvence[(a/2):]), len(sekvence[:(a/2)])
for sekvenca in sekvence:
    print sekvenca 

sekvence1=sekvence[0:(b)]   # A
for sekvenca in sekvence1:
    print sekvenca 
print '>>>>>'
sekvence2=sekvence[(b):a]  # B
for sekvenca in sekvence2:
    print sekvenca


# ANALIZA WT SEKVENC
for sekvenca in sekvence2:
    N = float(len(sekvenca))
    enke=where(sekvenca == 1)[0] ; dvojke=where(sekvenca == 2)[0]; amin=(min(dvojke)) ; amax=(max(enke)) ; delz=float(amax-(amin-1)) ; separation=delz-2
    print 'nova', sekvenca, N, delz, separation
    
    Nhot=len(dvojke)+len(enke)
    vse_int=[] ; vse_int.append(enke) ; vse_int.append(dvojke) ; vse_int=concatenate(vse_int)
    sekvence_podsekvence=[] ; sekvenca0=zeros(int(N)) ; prilozna=copy(sekvenca0)

    for i in range(Nhot):
        i=Nhot-i ; comb = combinations(vse_int, i)
        for j in comb:
            for k in j:
                if karkoli(k, enke)==True:
                    prilozna[(k)]=1
                else:
                    prilozna[(k)]=2
                
                pipa=copy(prilozna)
            sekvence_podsekvence.append(pipa) ; prilozna=copy(sekvenca0)
    Z_vse=0 ; Z_ena=0 ; Z_dve=0
    
    for sekvenca in sekvence_podsekvence: 
        print sekvenca 
        Zu, Zint, x_y_xy =VSOTA(sekvenca)
        if x_y_xy==3:
            Z_dve=Z_dve+Zint
        else:
            Z_ena=Z_ena+Zint
        Z_vse=Z_vse+Zint
    

    K_Zu=expand(Zu) ;    K_Z_dve=expand(Z_dve) ; K_Z_vse=expand(Z_vse)
    LL=len(K_Zu.args) ; LL2=len(Z_ena.args)  ; LL3=len(Z_dve.args) ; LL4=len(Z_vse.args)
    print 'too', (LL), (LL2), (LL3), (LL4)
    #priprava 
    suma1, suma2, suma3, suma4= delckanje(K_Zu, int(LL/4)) ; lam_Zu1=lambdify((w, x, y), suma1, modules='numpy') ; lam_Zu2=lambdify((w, x, y), suma2, modules='numpy')  ; lam_Zu3=lambdify((w, x, y), suma3, modules='numpy') ; lam_Zu4=lambdify((w, x, y), suma4, modules='numpy')    
    suma1, suma2, suma3 ,suma4= delckanje(Z_dve, int(LL2/4)) ; lam_Zd1=lambdify((w, x, y), suma1, modules='numpy') ; lam_Zd2=lambdify((w, x, y), suma2, modules='numpy')  ; lam_Zd3=lambdify((w, x, y), suma3, modules='numpy') ; lam_Zd4=lambdify((w, x, y), suma4, modules='numpy')    
    suma1, suma2, suma3 ,suma4= delckanje(Z_ena, int(LL3/4)) ; lam_Ze1=lambdify((w, x, y), suma1, modules='numpy') ; lam_Ze2=lambdify((w, x, y), suma2, modules='numpy')  ; lam_Ze3=lambdify((w, x, y), suma3, modules='numpy') ; lam_Ze4=lambdify((w, x, y), suma4, modules='numpy')    
    suma1, suma2, suma3 ,suma4= delckanje(Z_vse, int(LL4/4)) ; lam_Zv1=lambdify((w, x, y), suma1, modules='numpy') ; lam_Zv2=lambdify((w, x, y), suma2, modules='numpy')  ; lam_Zv3=lambdify((w, x, y), suma3, modules='numpy') ; lam_Zv4=lambdify((w, x, y), suma4, modules='numpy')    

    #  >>>>>>>>>>>>>>>>>>>> skeniranje helix propensity <<<<<<<<<<<<<<<<<<

    #####    P A R A M E T R I ##########
    v = 0.048 ; T=273.15 ; R=1.987
    dGch_start=-50
    GCHwt=arange(dGch_start,+360, 10)
    dGx=-3000 ; dGy=-3000 ; xx=exp(dGx/(-R*T)) ; yy=exp(dGy/(-R*T))
    ####################################

    for dGch in GCHwt:
        wwt=float(exp(dGch/(-R*T)))
        #racunanje  
        ZuWT= (lam_Zu1(wwt, 1, 1)) +  (lam_Zu2(wwt, 1, 1))+  (lam_Zu3(wwt, 1, 1))+  (lam_Zu4(wwt, 1, 1)) 
        ZuWT=float(ZuWT+(K_Zu).args[0])
        Z_dve=(lam_Zd1(wwt, xx, yy))+(lam_Zd2(wwt, xx, yy))+(lam_Zd3(wwt, xx, yy))+(lam_Zd4(wwt, xx, yy))
        Z_ena=(lam_Ze1(wwt, xx, yy))+(lam_Ze2(wwt, xx, yy))+(lam_Ze3(wwt, xx, yy))+(lam_Ze4(wwt, xx, yy))
        Z_vse=(lam_Zv1(wwt, xx, yy))+(lam_Zv2(wwt, xx, yy))+(lam_Zv3(wwt, xx, yy))+(lam_Zv4(wwt, xx, yy))

        p_eni=Z_ena/Z_vse ; p_dva=Z_dve/Z_vse
    
        celiHEL_dol=[] ; celiHEL_kratki=[]; polHEL=[] 
        for i in ((K_Z_dve).atoms(Mul)):
            if i.args[0]<1e-8:
                break
            I=lambdify((w, x, y), i, modules='numpy') ; jj=(I(wwt,xx,yy)) ; p=jj/Z_vse    ; st_vezi=(i.args[3]).args[1]
            if st_vezi>(delz-1):
                if st_vezi>(N/2-1):
                    celiHEL_dol.append(p)
                else:
                    celiHEL_kratki.append(p)
            else:
                polHEL.append(p)
        
        hel=[]; S1=[] 
        for i in ((K_Z_vse).atoms(Mul)):
            II=diff(i, w) ; II=II*w ; nh=II/i
            I=lambdify((w, x, y), i, modules='numpy')
            jj=(I(wwt,xx,yy)) ; p=jj/Z_vse             
            Si=p*log(p) ; S1.append(Si)
            NH=p*nh ;     hel.append(NH)        
    
        abwtbound=float(sum(hel)/(N-2)) ; entr=-R*T*sum(S1)
        sig=[] 
        for i in ((K_Z_vse).atoms(Mul)):
            I=lambdify((w, x, y), i, modules='numpy')
            jj=(I(wwt,xx,yy)) ; p=jj/Z_vse  
            II=diff(i, w) ; II=II*w ; nh=II/i ; fH=nh/(N-2)
            sig.append(p*(fH-abwtbound)**2)
    
        ssigma=float(sum(sig)**0.5)    
        dGbindWT=-R*T*log(Z_vse/ZuWT) ; KbindWT=Z_vse/ZuWT

        pol=sum(polHEL) ;  celi_dol=sum(celiHEL_dol) ; celi_kratki=sum(celiHEL_kratki)

        ENTROPIJE.append(entr) ; CELI_dol.append(celi_dol) ; CELI_kr.append(celi_kratki) ; POL.append(pol) ; POPULACIJE1.append(p_eni) ; POPULACIJE2.append(p_dva) 
        SIGME.append(ssigma) ; HELICNOSTI.append(abwtbound) ; DGCH.append(dGch) ; SEPARACIJE.append(separation) ; AFINITETE.append(KbindWT)

    #  >>>>>>>>>>>>>>>>>>>> skeniranje interaction strength <<<<<<<<<<<<<<<<<<

    #####    P A R A M E T R I ##########
    dGch=50 ; wwt=float(exp(dGch/(-R*T))) 
    dGint_start=-4000    
    G_int=arange(dGint_start,-500, 100)
    ####################################

    for G_i in G_int:
        dGx=G_i ; dGy=G_i
        xx=exp(dGx/(-R*T)) ; yy=exp(dGy/(-R*T))
    
        #racunanje  
        ZuWT= (lam_Zu1(wwt, 1, 1)) +  (lam_Zu2(wwt, 1, 1))+  (lam_Zu3(wwt, 1, 1))+  (lam_Zu4(wwt, 1, 1)) 
        ZuWT=float(ZuWT+(K_Zu).args[0])
        Z_dve=(lam_Zd1(wwt, xx, yy))+(lam_Zd2(wwt, xx, yy))+(lam_Zd3(wwt, xx, yy))+(lam_Zd4(wwt, xx, yy))
        Z_ena=(lam_Ze1(wwt, xx, yy))+(lam_Ze2(wwt, xx, yy))+(lam_Ze3(wwt, xx, yy))+(lam_Ze4(wwt, xx, yy))
        Z_vse=(lam_Zv1(wwt, xx, yy))+(lam_Zv2(wwt, xx, yy))+(lam_Zv3(wwt, xx, yy))+(lam_Zv4(wwt, xx, yy))

        p_eni=Z_ena/Z_vse ; p_dva=Z_dve/Z_vse

        celiHEL_dol=[] ; celiHEL_kratki=[]; polHEL=[] 
        for i in ((K_Z_dve).atoms(Mul)):
            if i.args[0]<1e-8:
                break
            I=lambdify((w, x, y), i, modules='numpy') ; jj=(I(wwt,xx,yy)) ; p=jj/Z_vse    ; st_vezi=(i.args[3]).args[1]
            if st_vezi>(delz-1):
                if st_vezi>(N/2-1):
                    celiHEL_dol.append(p)
                else:
                    celiHEL_kratki.append(p)
            else:
                polHEL.append(p)
        
        hel=[]; S1=[] 
        for i in ((K_Z_vse).atoms(Mul)):
            II=diff(i, w) ; II=II*w ; nh=II/i
            I=lambdify((w, x, y), i, modules='numpy')
            jj=(I(wwt,xx,yy)) ; p=jj/Z_vse             
            Si=p*log(p) ; S1.append(Si)
            NH=p*nh ;     hel.append(NH)        
    
        abwtbound=float(sum(hel)/(N-2)) ; entr=-R*T*sum(S1)
        sig=[] 
        for i in ((K_Z_vse).atoms(Mul)):
            I=lambdify((w, x, y), i, modules='numpy')
            jj=(I(wwt,xx,yy)) ; p=jj/Z_vse  
            II=diff(i, w) ; II=II*w ; nh=II/i ; fH=nh/(N-2)
            sig.append(p*(fH-abwtbound)**2)
    
        ssigma=float(sum(sig)**0.5)    
        dGbindWT=-R*T*log(Z_vse/ZuWT) ; KbindWT=Z_vse/ZuWT
        
        pol=sum(polHEL) ;  celi_dol=sum(celiHEL_dol) ; celi_kratki=sum(celiHEL_kratki) 

        ENTROPIJE_2.append(entr) ; CELI_dol_2.append(celi_dol) ; CELI_kr_2.append(celi_kratki) ; POL_2.append(pol) ; POPULACIJE1_2.append(p_eni) ; POPULACIJE2_2.append(p_dva) 
        SIGME_2.append(ssigma) ; HELICNOSTI_2.append(abwtbound) ; DG_INT.append(G_i) ; SEPARACIJE_2.append(separation) ; AFINITETE_2.append(KbindWT)
    
#savetxt('rezultat_dGch_32_NOVA2b.txt', column_stack([SEPARACIJE,DGCH,SIGME,HELICNOSTI,ENTROPIJE,AFINITETE,POPULACIJE1,POPULACIJE2,CELI_dol,CELI_kr,POL]))
#savetxt('rezultat_dGint_32_NOVA2.txt', column_stack([SEPARACIJE_2,DG_INT,SIGME_2,HELICNOSTI_2,AFINITETE_2,ENTROPIJE_2,POPULACIJE1_2,POPULACIJE2_2,CELI_dol_2,CELI_kr_2,POL_2]))


dGch_start=DGCH[0] ; dGint_start=DG_INT[0]
print dGch_start, dGint_start

St_celX_DOL=[] ; St_celY_DOL=[]  ; St_celX_KR=[] ; St_celY_KR=[] ; St_polX=[] ; St_polY=[] ; St_eniX=[] ; St_eniY=[]
for i,j,k,l,m,n,o in zip(SEPARACIJE,DGCH,POPULACIJE1,POPULACIJE2,CELI_dol,CELI_kr,POL):
    if k>l: 
        St_eniX.append(i) ; St_eniY.append(j)
    else:
        if (m+n)>o: 
            if m>n:
                St_celX_DOL.append(i) ; St_celY_DOL.append(j)
            else:
                St_celX_KR.append(i) ; St_celY_KR.append(j)
        else: 
            St_polX.append(i) ; St_polY.append(j)



St_celX_DOL2=[] ; St_celY_DOL2=[]  ; St_celX_KR2=[] ; St_celY_KR2=[] ; St_polX2=[] ; St_polY2=[] ; St_eniX2=[] ; St_eniY2=[]
for i,j,k,l,m,n,o in zip(SEPARACIJE_2,DG_INT,POPULACIJE1_2,POPULACIJE2_2,CELI_dol_2,CELI_kr_2,POL_2):
    if k>l:
        St_eniX2.append(i) ; St_eniY2.append(j)
    else:
        if (m+n)>o: 
            if m>n: 
                St_celX_DOL2.append(i) ; St_celY_DOL2.append(j)
            else:
                St_celX_KR2.append(i) ; St_celY_KR2.append(j)
        else: 
            St_polX2.append(i) ; St_polY2.append(j)        

fig=plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(111)
ax1.plot(St_eniX, St_eniY, 'ro',alpha=1)
ax1.plot(St_celX_DOL, St_celY_DOL, 'ko',alpha=1)
ax1.plot(St_celX_KR, St_celY_KR, 'yo',alpha=1)
ax1.plot(St_polX, St_polY, 'bo',alpha=1)
ax1.set_xlabel('hotspot separation ; N')
ax1.set_ylabel('residue helix  \n propensity ; cal/mol')

fig=plt.figure(figsize=(10,8))
ax2 = fig.add_subplot(111)
ax2.plot(St_eniX2, St_eniY2, 'ro',alpha=1)
ax2.plot(St_celX_DOL2, St_celY_DOL2, 'ko',alpha=1)
ax2.plot(St_celX_KR2, St_celY_KR2, 'yo',alpha=1)
ax2.plot(St_polX2, St_polY2, 'bo',alpha=1)
ax2.set_xlabel('hotspot separation ; N')
ax2.set_ylabel('hotspot interaction  \n energy ; cal/mol')

print("--- %s seconds ---" % (time.time() - start_time))   
