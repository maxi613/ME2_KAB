import numpy as np
import sympy as sy
import pandas as pd
from math import *
from matplotlib import pyplot as plt
#Eingangsdaten variante 2

P=18400  #W
n=1540  #1/min
i=3.2   #Übersetzung

def pow2(basis):
    ergebnis=basis**2
    return ergebnis

def pown(basis,n):
    ergebnis=basis**n
    return ergebnis

def cosd(winkel):
    ergebnis=cos(winkel*pi/180)
    return ergebnis

def sind(winkel):
    ergebnis=sin(winkel*pi/180)
    return ergebnis

def involut (x):
  alphan=tan(x)-x
  return alphan

def wurzel_3(x):
    ergebnis=(x)**(1/3)
    return ergebnis

def tand(winkel):
    ergebnis=tan(winkel*pi/180)
    return ergebnis

def atand(winkel):
    ergebnis=atan(winkel)*180/pi
    return ergebnis

def acosd(winkel):
    ergebnis=acos(winkel)*180/pi
    return ergebnis

#---Zahnzahl/Modul/Übersetzung/Achsabstand---#

T=P*60/(2*pi*n)*1000 
print("Torsionsmoment: ",round(T,2)) #Nmm
psi_d=1.1



#Werkstoff: 34CrMo4
sigma_flim=525  #Mpa
sigma_hlim=1650 #Mpa
z1=31
yfs=4
Sf=1.5
Ka=1.25
beta=14

#Berechnung Modul
m_n=1.85*wurzel_3(Ka*T*pow2(cosd(beta))/(pow2(z1)*psi_d*sigma_flim)) #Normalmodul
print("Normalmodul gerechnet: ",round(m_n,2))
m_n=1.5 #gerundet
print("Normalmodul gewählt ~: ",round(m_n,2))
print("")

#Zähnezahl und Übersetzung
z2=z1*i
print("Zähnenzahl Z1 =",z1)
print("Zähnenzahl Z2 =", z2)
z2=99 #gerundet
print("Zähnenzahl Z2 ~", z2)
i=z2/z1
print("Tatsächliches Übersetzungsverhältnis: ",round(i,2))
print("")

#Nullachsabstand und Teilkreisdurchmeser
d1=z1*m_n/cosd(beta) #Teilkreisdurchmesser
d2=z2*m_n/cosd(beta) #Teilkreisdurchmesser
ad=(d1+d2)/2

print("Teilkreisdurchmesser d1 =",round(d1,2),"mm")
print("Teilkreisdurchmesser d2 =",round(d2,2),"mm")
print("Nullachsenabstand ad =",round(ad,2),"mm")
a=100 #gewählt aus Normreihe R40
print ("Achsabstand a:", a,"mm")
n2=n/i
print("Drehzahl n2:", round(n2,2),"1/min")
T_2=P*60/(2*pi*n2)*1000
print("T_2 :",round(T_2,2))
print("")


#---Abmessung Zahnräder---#

mt=m_n/cosd(beta) #Stirnmodul
alpha=20
alpha_t=atand(tand(alpha)/cosd(beta))    #Strineingriffswinkel
alpha_wt=acosd(ad/a*cosd(alpha_t))    #Betriebseingriffswinkel im Stirnschnitt
alpha_w=acosd(ad/a*cosd(alpha))
print("Stirnmodul: ", round(mt,2))
print("Strineingriffswinkel alpha_t: ",round(alpha_t,2))
print("Betriebseingriffswinkel im Stirnschnitt alpha_wt: ", round(alpha_wt,2))
print("Betriebseingriffswinkel im Normalschnitt alpha_w: ", round(alpha_w,2))
print("")

#Verschiebungsfaktoren
x=(involut(alpha_wt*pi/180)-involut(alpha_t*pi/180))/(2*tand(alpha))*(z1+z2)
z1n=z1/pown(cosd(beta),3)
z2n=z2/pown(cosd(beta),3)
x1=x/2+(0.5-x/2)*log(i)/(log((z2n*z1n)/100))
x2=x-x1
print("Ersatzzähnezahl Z1: ", round(z1n,2))
print("Ersatzzähnezahl Z2: ", round(z2n,2))
print("")
#Verschiebung
V1=m_n*x1
V2=m_n*x2
print("Verschiebung1 V1",round(V1,2))
print("Verschiebung2 V2",round(V2,2))
print("")

#Kopfhöhenänderung
k=a-ad-m_n*x
print("Kopfhöhenänderung k:",round(k,2))
print("")
#Kopfkreisdurchmesser
da1=d1+2*(m_n+V1+k)
da2=d2+2*(m_n+V2+k)

#Fußkreisdurchmesser
hf=1.25*m_n
df1=d1-2*hf+2*V1
df2=d2-2*hf+2*V2

#Grundkreisdurchmesser
db1=d1*cosd(alpha_t)
db2=d2*cosd(alpha_t)

#Betriebswälzdurchmesser
dw1=d1*cosd(alpha)/cosd(alpha_w)
dw2=d2*cosd(alpha)/cosd(alpha_w)


print("Gesamtverschiebung x =",round(x,2))
print("Verschiebungsfaktor x1 =",round(x1,2))
print("Verschiebungsfaktor x2 =",round(x2,2))
print("")
print("Verschiebung V1 =", round(V1,2),"mm")
print("Verschiebung V2 =", round(V2,2),"mm")
print("")
print("Kopfkreisdurchmesser da1 =", round(da1,2),"mm")
print("Kopfkreisdurchmesser da2 =", round(da2,2),"mm")
print("")
print("Fußkreisdurchmesser df1 =",round(df1,2),"mm")
print("Fußkreisdurchmesser df2 =",round(df2,2),"mm")
print("")
print("Grundkreisdurchmesser db1 =",round(db1,2),"mm")
print("Grundkreisdurchmesser db2 =",round(db2,2),"mm")
print("")
print("Betriebswälzdurchmesser dw1 =",round(dw1,2),"mm")
print("Betriebswälzdurchmesser dw2 =",round(dw2,2),"mm")
print("")


#---Zahnradbreite/Überdeckung und Tragfähigkeit---#
b1=d1*psi_d
b1=floor(b1)      #abgerundet
b2=b1-2

print("Breite Zahnrad b1 =", b1," mm")
print("Breite Zahnrad b2 =", b2," mm")
print("")

#Überdeckung
epsilon_alpha=(0.5*(sqrt(pow2(da1)-pow2(db1))+sqrt(pow2(da2)-pow2(db2)))-a*sind(alpha_wt))/(pi*mt*cosd(alpha_t))
print("Profilüberdeckung im Stirnschnitt epsilon_alpha (21.57): ", round(epsilon_alpha,4))
epsilon_beta=b2*sind(beta)/(pi*m_n)
print("Profilüberdeckung im Normalschnitt epsilon_beta (21.44): ", round(epsilon_beta,2))
epsilon_gamma= epsilon_alpha+epsilon_beta
print("Gesamtüberdeckung epsilon_gamma (21.46): ", round(epsilon_gamma,2))
print("")

#Kräfte am Zahnrad
F_t1=2*T/dw1
F_t2=2*T_2/dw2
F_t=F_t1=F_t2
print("Umfangskraft im Stirnschnitt F_t1 (21.70): ", round(F_t1,2))
print("Umfangskraft im Stirnschnitt F_t2 (21.70): ", round(F_t2,2))
Fr=F_t1*tand(alpha_wt)
print("Radialkraft Fr (21.71): ",round(Fr,2))
print("")

#Zahnfußtraffähigkeit
print("Zahnfußtragfähigkeit:")
Y_epsilon=0.25+0.75*pow2(cosd(beta))/epsilon_alpha
Y_beta=1-epsilon_beta*beta/120
if Y_beta<0.75:
    Y_beta=0.75


print("Y_epsilon: ", round(Y_epsilon,2))
print("Y_beta: ",round(Y_beta,2))
print("")

Y_sa1 = 1.54  # Tabelle 21-19b
Y_sa2 = 1.89

Y_Fa1=2.94  # Tabelle 21-19a
Y_Fa2=2.05

sigma_F01=F_t1/(m_n*b2)*Y_Fa1*Y_sa1*Y_beta*Y_epsilon
sigma_F02=F_t1/(m_n*b2)*Y_Fa2*Y_sa2*Y_beta*Y_epsilon
print("örtliche Zahnfußspannung simga_F01 (21.82): ",round(sigma_F01,2))
print("örtliche Zahnfußspannung simga_F02 (21.82): ",round(sigma_F02,2))
print("")

KFges=1.25*1.1*1.25*1.25 #Vorlesung
sigma_F1=sigma_F01*KFges
sigma_F2=sigma_F02*KFges
print("örtliche Zahnfußnennspannung sigma_F1 (21.83): ", round(sigma_F1,2))
print("örtliche Zahnfußnennspannung sigma_F2 (21.83): ", round(sigma_F2,2))
print("")

Y_st=2
Y_nt=1
Y_delta=1
Y_relT=1
Y_x=1

sigma_FG=sigma_flim*Y_st*Y_nt*Y_delta*Y_x
print("Zahnfuß Grenzfestigkeit sigma_FG (21.84a): ", round(sigma_FG,2))

if sigma_F1>=sigma_F2:

    S_F=sigma_FG/sigma_F1
    if S_F<1.5:
        print("Sicherheit ist nicht gewährleiset", round(S_F,2))
    else:
        print("Berechnung mit sigma_F1 Sicherheit ist gewährleistet S_F (21.85): ", round(S_F,2))

elif sigma_F2>=sigma_F1:
    
    S_F=sigma_FG/sigma_F2
    if S_F<1.5:
        print("Sicherheit ist nicht gewährleiset", round(S_F,2))
    else:
        print("Berechnung mit simga_F2. Sicherheit ist gewährleistet S_F (21.85): ", round(S_F,2))
print("")

#Grübchentrafgähigkeit

betab=acosd(sind(alpha)/sind(alpha_t))
Z_H=sqrt(2*cosd(betab)/(pow2(cosd(alpha_t))*tand(alpha_wt)))
Z_E=sqrt(0.175*210000)
Z_e=sqrt(1/epsilon_alpha)
Z_b=sqrt(cosd(beta))
print("Winkel beta_b: ", round(betab,2))
print("Zonenfaktor Z_H: ",round(Z_H,2))
print("Elastizitätsfaktor Z_E: ", round(Z_E,2))
print("Überdeckungsfaktor Z_epsilon: ", round(Z_e,2))
print("Schrägenfakor Z_beta: ", round(Z_b,2))
print("")

sigma_HC=sqrt(F_t1/(b2*d1)*(i+1)/i)*Z_H*Z_E
sigma_H0=Z_e*Z_b*sigma_HC
print("Pressung im Wälzpunkt sigma_HC (27.87): ", round(sigma_HC,2))
print("Flankenpressung sigma_H0 (21.88): ",round(sigma_H0,2))

K_Hges=sqrt(1.25*1.1*1.25*1.25)
print("Belastungseinflussfaktoren K_Hges (21.81): ", round(K_Hges,2))
sigma_H=sigma_H0*K_Hges
print("Flankenpressung am Wälzkreis sigma_H (21.89): ", round(sigma_H,2))
sigma_HG=sigma_hlim
print("Flankengrenzfestigkeit sigma_HG: ", sigma_HG)
print("")

S_H=sigma_HG/sigma_H

if S_H > 1.6:
    print("Festigkeit gewährleistet S_H (21.90a): ", round(S_H))
else:
    print("Sicherheit ist nicht gewährleistet S_H: ", round(S_H))
print("")

#---Welle---#

#Werkstoff E295
sigma_bwn=245 
M_v1=Ka*1.17*T #Welle Ritzel
M_v2=Ka*1.17*T_2 #Welle Rad

#Vorauslegung der Durchmesser
d_w1=3.4*wurzel_3(M_v1/(sigma_bwn))
d_w2=3.4*wurzel_3(M_v2/(sigma_bwn))
print("Vorausgelegter Durhmessser Welle Ritzel d_w1 (11.8): ", round(d_w1,2))
print("Vorausgelegter Durhmessser Welle Ritzel d_w2 (11.8): ", round(d_w2,2))
print("")

s_r1=(df1-d_w1)/2 #Kranzdicke
s_r2=(df2-d_w2)/2

if s_r1>3.5*m_n and s_r2>3.5*m_n:
    print("Kranzdicke Welle Ritzel s_r1: ", round(s_r1,2),"mm")
    print("Kranzdicke Welle Ritzel s_r2: ", round(s_r2,2),"mm")
print("")


#Überprüfung
F_a=F_t*tand(beta)
print("Axialkraft F_a (21.72): ", round(F_a,2))

l=90 #Lagerlänge l=90 mm
l_1=l_2=l/2
print("")

#--X-Z Ebene--
""" 
Indize 1 => Welle Ritzel
Indize 2 => Welle Rad
"""
print("X-Z Ebene:")
#Kräfte im Lager B

Fb1=sy.S('Fb1')
sigma_M_B1=sy.Eq(Fb1*l+Fr*l_1+F_a*dw1/2,0) #Summe aller Momente um Lager A an Welle Ritzel = 0
F_b1=str(sy.solve(sigma_M_B1)).replace("[","")
F_b1=round(float(F_b1.replace("]","")),2)
print("Lagerkraft im Lager B Welle Ritzel F_b1: ", F_b1)


Fb2=sy.S('Fb2')
sigma_M_B2=sy.Eq(Fb2*l-Fr*l_2+F_a*dw2/2,0) #Summer aller Momente um Lager A an Welle Rad = 0
F_b2=str(sy.solve(sigma_M_B2)).replace("[","")
F_b2=round(float(F_b2.replace("]","")),2)
print("Lagerkraft im Lager B Welle Rad F_b2: ", F_b2)
print("")


#Kräfter im Lager A
Fa1=sy.S('Fa1')
sigma_Fy1=sy.Eq(Fr+F_b1+Fa1,0) #Summe aller Kräfte Welle Ritzel = 0
F_a1=str(sy.solve(sigma_Fy1)).replace("[","")
F_a1=round(float(F_a1.replace("]","")),2)
print("Lagerkraft im Lager A Welle Ritzel F_a1: ", F_a1)

Fa2=sy.S('Fa2')
sigma_Fy2=sy.Eq(-Fr+F_b2+Fa2,0) #Summe aller Kräfte Welle Rad = 0
F_a2=str(sy.solve(sigma_Fy2)).replace("[","")
F_a2=round(float(F_a2.replace("]","")),2)
print("Lagerkraft im Lager A Welle Rad F_a2: ", F_a2)
print("")

#Momente im Lager A
M_a1=F_a1*l_1
M_a2=F_a2*l_2
print("Moment im Lager A Welle Ritzel M_a1: ", round(M_a1,2))
print("Moment im Lager A Welle Rad M_a2: ", round(M_a2,2))

#Momente im Lager B
M_b1=F_b1*l_1
M_b2=F_b2*l_2
print("Moment im Lager B Welle Ritzel M_b1: ", round(M_b1,2))
print("Moment im Lager B Welle Rad M_b2: ", round(M_b2,2))
print("")

#--Y-Z Ebene--
print("Y-Z Ebene:")
F_y12=F_t*l_1/l
M_y12=F_y12*l_1
print("Momente Welle Rad und Ritzel in beiden Lager M_y12: ", round(M_y12,2))
print("")

#Welle Rad
if abs(M_b1)>abs(M_a1):
    M_eq1=sqrt(pow2(M_b1)+pow2(M_y12))*Ka
    print("größeres Moment M_b1 => M_eq1: ", round(M_eq1,2))
elif abs(M_a1)>abs(M_b1):
    M_eq1=sqrt(pow2(M_a1)+pow2(M_y12))*Ka
    print("größeres Moment M_a1 => M_eq1: ", round(M_eq1,2))

#Welle Ritzel
if M_b2>M_a2:
    M_eq1=sqrt(pow2(M_b2)+pow2(M_y12))*Ka
    print("größeres Moment M_b2 => M_eq2: ", round(M_eq2,2))
elif M_a2>M_b2:
    M_eq2=sqrt(pow2(M_a2)+pow2(M_y12))*Ka
    print("größeres Moment M_a2 => M_eq2: ", round(M_eq2,2))
print("")

T_eq1=Ka*T
T_eq2=Ka*T_2
print("equivalentes Torsionsmoment T_eq1: ", round(T_eq1,2))
print("equivalentes Torsionsmoment T_eq2: ", round(T_eq2,2))
print("")

M_v11=sqrt(pow2(M_eq1)+0.75*pow2(0.7*T_eq1))
M_v22=sqrt(pow2(M_eq2)+0.75*pow2(0.7*T_eq2))
print("Vergleichsmoment Welle Rad M_v11: ", round(M_v11))
print("Vergleichsmoment Welle Rad M_v22: ", round(M_v22))
print("")

d_1=3.4*wurzel_3(M_v11/sigma_bwn)
d_2=3.4*wurzel_3(M_v22/sigma_bwn)
print("tatsächlich benötigter Durchmesser d_1: ", round(d_1,2))
print("tatsächlich benötiger Duchmesser d_2: ", round(d_2,2))
print("")

if d_1<d_w1:
    print("Vorausgewählter Durchmesser an Welle Ritzel kann verwendet werden")
else:
    print("Neuer Durchmesser an Welle Ritzel muss gewählt werden. Gegebenfalls Rechnung überprüfen")

if d_2<d_w2:
    print("Vorausgewählter Durchmesser an Welle Rad kann verwendet werden")
else:
    print("Neuer Durchmesser an Welle Rad muss gewählt werden. Gegebenfalls Rechnung überprüfen")

input("Press enter to exit ;)")