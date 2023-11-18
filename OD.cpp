import numpy as np
from math import cos,sin,asin,acos,atan2,sqrt,pi,radians,degrees
k = 2.9812312e-5 #GD/D
c = 173.144632674240 #AU/day
c2 = -1
obliquity = radians(0) #degrees -> radians
#IMPORTANT FUNCTIONS TO BE USED IN THE PROGRAM
#converts values to julian date
def julian(year,month,day,UTC):
return(367 * year - int(7 * (year + int( (month + 9) / 12 ) ) / 4) + int(275 * month / 9) + day + 1721013.5 + UTC / 24)
  
#takes raw line from input file and converts it into usuable values in usuable units
def cleanObs(rawObs):
rawObs = rawObs.split()
UTCList = rawObs[3].split(":")
UTC = float(UTCList[0]) + float(UTCList[1]) / 60 + float(UTCList[2]) / 3600
julianDate = julian(float(rawObs[0]),float(rawObs[1]),float(rawObs[2]),UTC)
raList = rawObs[4].split(":")
ra = radians(15 * (float(raList[0]) + float(raList[1]) / 60 + float(raList[2]) / 3600))
decList = rawObs[5].split(":")
if float(decList[0]) > 0 and float(decList[1]) > 0 and float(decList[2]) > 0:
dec = radians(float(decList[0]) + float(decList[1]) / 60 + float(decList[2]) / 3600)
else:
dec = radians(float(decList[0]) - float(decList[1]) / 60 - float(decList[2]) / 3600)
sunvec = [float(rawObs[6]),float(rawObs[7]),float(rawObs[8])]
return([ra,dec,julianDate,sunvec])
#takes in cleaned observation data and returns the rhohats
def getRhoHat(obsData):
return([cos(obsData[1]) * cos(obsData[0]),
cos(obsData[1]) * sin(obsData[0]),
sin(obsData[1])])
#takes in cleaned data and rhohats to return Ds
def getD(obsData,rhoHatList):
D0 = np.dot(rhoHatList[0],np.cross(rhoHatList[1],rhoHatList[2]))
D1 = np.dot(np.cross(obsData[3],rhoHatList[1]),rhoHatList[2])
D2 = np.dot(np.cross(rhoHatList[0],obsData[3]),rhoHatList[2])
D3 = np.dot(np.cross(rhoHatList[1],obsData[3]),rhoHatList[0])
return(D0,D1,D2,D3)
#6 orbital element functions
def get_a(radius,speed):
return((2/radius - speed**2)**(-1))
def get_e(a,radius,velocity):
return(sqrt( 1 - (np.linalg.norm(np.cross(radius,velocity)))**2 / a))
def get_i(radius,velocity):
return(acos((np.cross(radius,velocity))[2] /(np.linalg.norm(np.cross(radius,velocity)))))
def get_Omega(h,i):
omega = atan2(h[0] / np.linalg.norm(h) / sin(i), -h[1] / np.linalg.norm(h) / sin(i) )
if omega < 0:
omega += 2*pi
return(omega)
def get_fw(R,i,Omega,e,a,RDOT):
fwsin = R[2] / np.linalg.norm(R) / sin(i)
fwcos = (R[0] / np.linalg.norm(R) + cos(i) * fwsin * sin(Omega) ) / cos(Omega)
fw = atan2(fwsin,fwcos)
fcos = (a * (1-e**2) / np.linalg.norm(R) - 1) / e
fsin = sqrt(a * (1-e**2)) * np.dot(R,RDOT) / e / np.linalg.norm(R)
f = atan2(fsin,fcos)
if fw - f < 0:
  f = f - 2*pi
  return(f,fw-f)
def get_M(e,R,a,f):
  E = acos((1 - np.linalg.norm(R) / a) / e)
  if f < 0:
    f = f + 2 * pi
  if f > pi and E < pi:
    E = 2 * pi - E
  elif f < pi and E > pi:
    E = 2 * pi - E
  return(E - e * sin(E))
#Scalar equation of lagrange - takes in taus, rhohats, Ds, and sun vectors in order to generate the first guess of r2
def SEL(taus,sun2,rhohat2,Ds):
#absurd constants
  A1 = taus[1]/taus[2]
  B1 = A1 * (taus[2]**2 - taus[1]**2) / 6
  A3 = - taus[0] / taus[2]
  B3 = A3 * (taus[2]**2 - taus[0]**2) / 6
  A = (A1 * Ds[1] - Ds[2] + A3 * Ds[3] ) / -Ds[0]
  B = (B1 * Ds[1] + B3 * Ds[3]) / -Ds[0]
  E = -2 * (np.dot(rhohat2,sun2))
  F = np.dot(sun2,sun2)
  a = -(A**2 + A * E + F)
  b = -(2 * A * B + B * E)
  c = -1 * (B**2)
  print(Ds)
  print(a,b,c)
#calculates roots
roots = np.roots([1,0,a,0,0,b,0,0,c])
roots2 = []
print(roots)
#filters complex roots
for item in roots:
  if np.real(item) > 0 and np.imag(item) == 0:
    roots2.append(np.real(item))
  rhos = []
  realRoots = []
  print(Ds)
#filters negative rhos
for item in roots2:
  if A + B / item**3 > 0:
  realRoots.append(item)
  rhos.append(A + B / item**3)
return(realRoots,rhos)
#computes fg values through function or taylor series.
def fg(tau,r2,r2dot,flag):
  r2scal = np.linalg.norm(r2)
  r2dotscal = np.linalg.norm(r2dot)
  a = get_a(r2scal,r2dotscal)
  n = sqrt(1 / a**3)
#series
if flag == 0:
  ecosE0 = 1 - r2scal / a
  esinE0 = np.dot(r2, r2dot) / n / a**2
  Eguess1 = n * tau
  Eguess2 = Eguess1 + 100
while abs(Eguess1 - Eguess2) > 1e-12:
  fE = Eguess1 - ecosE0 * sin(Eguess1) + esinE0 * (1 - cos(Eguess1) ) - n * tau
  fdE = 1 - ecosE0 * cos(Eguess1) + esinE0 * sin(Eguess1)
  Eguess2 = Eguess1 - fE/fdE
  temporary_swap = Eguess1
  Eguess1 = Eguess2
  Eguess2 = temporary_swap
  f = 1 - a * (1 - cos(Eguess1)) / r2scal
  g = tau + (sin(Eguess1) - Eguess1) / n
  return(f,g)
u = 1 / r2scal**3
z = np.dot(r2,r2dot) / r2scal**2
q = r2dotscal**2 / r2scal**2 - u
#3rd order
if flag == 3:
  f = 1 - u * tau**2 / 2 + u * z * tau**3 / 2
  g = tau - u * tau**3 / 6
  return(f,g)
#4th order
if flag == 4:
  f = 1 - u * tau**2 / 2 + u * z * tau**3 / 2 + (3*u*q - 15 * u * z**2 + u**2) * tau**4 / 24
  g = tau - u * tau**3 / 6 + u * z * tau**4 / 4
  return(f,g)
def MOG(obsline1,obsline2,obsline3,filename):
  #dumps textfile into a list
  inputFile = open(filename)
  allObs = []
  for line in inputFile.readlines():
    allObs.append(line)
    inputFile.close()
    #makes each observation line into a neat list, with format ra dec time sunvec
    obs1 = cleanObs(allObs[obsline1])
    obs2 = cleanObs(allObs[obsline2])
    obs3 = cleanObs(allObs[obsline3])
    #calculates initial values of tau
    tau1 = k * (obs1[2] - obs2[2])
    tau3 = k * (obs3[2] - obs2[2])
    tau = tau3 - tau1
    print("taus",tau1, tau3, tau)
#adds rhohats to each observation information
obs1.append(getRhoHat(obs1))
obs2.append(getRhoHat(obs2))
obs3.append(getRhoHat(obs3))
#makes own rhohat list for no good reason lol, but I had code already that needed a list with all rhohats
rhoHatList = [np.array(obs1[4]),np.array(obs2[4]),np.array(obs3[4])]
#makes list of D
Dj1 = getD(obs1,rhoHatList)
Dj2 = getD(obs2,rhoHatList)
Dj3 = getD(obs3,rhoHatList)
#SCALAR EQUATION OF LAGRANGE to find first initial guesses of r2
r2candidates = SEL([tau1,tau3,tau],obs2[3],obs2[4],[Dj1[0],Dj1[2],Dj2[2],Dj3[2]])
print("******************************************")
print(len(r2candidates[0]),"roots found")
print("******************************************")
#loops over each root found
for i in range(len(r2candidates[0])):
  r2guess = r2candidates[0][i]
  rho2 = r2candidates[1][i]
  print()
  print("Running case number",i+1,", where r2 :",r2candidates[0][i],"au")
fgflag = int(input("Enter 3 for 3rd order series, 4 for 4th, and 0 for function "))
print()
print()
#initial r2dot guess + more absurd constants
f1 = 1 - tau1**2 / 2 / r2guess**3
f3 = 1 - tau3**2 / 2 / r2guess**3
g3 = tau3 - tau3**3 / 6 / r2guess**3
g1 = tau1 - tau1**3 / 6 / r2guess**3
c1 = g3 / (f1*g3 - g1 * f3)
c3 = -g1 / (f1 *g3 - g1 * f3)
d1 = -f3 / (f1 * g3 - f3 * g1)
d3 = f1 / (f1 * g3 - f3 * g1)
rho3 = np.array(rhoHatList[2]) * (c1 * Dj1[3] - Dj2[3] + c3 * Dj3[3]) / c3 / Dj1[0]
rho1 = np.array(rhoHatList[0]) * (c1 * Dj1[1] - Dj2[1] + c3 * Dj3[1]) / c1 / Dj1[0]
rho2 = np.array(rhoHatList[1]) * (c1 * Dj1[2] - Dj2[2] + c3 * Dj3[2]) / c2 / Dj2[0]
r3 = rho3 - obs3[3]
r1 = rho1 - obs1[3]
r2dotguess = d1 * r1 + d3 * r3
print("HIII",np.linalg.norm(r2dotguess))
r2prev = r2guess + 100
r2guess = c1 * r1 + c3 * r3
#iteration number tracker
i = 0
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#meat of the MoG, where it iterates r2 until it changes by <1e-20
while abs(np.linalg.norm(r2guess) - np.linalg.norm(r2prev)) > 1e-20:
#light travel correction
tau1 = k * ((obs1[2]- np.linalg.norm(rho1) / c) - (obs2[2]- np.linalg.norm(rho2) / c))
tau3 = k * ((obs3[2]- np.linalg.norm(rho3) / c) - (obs2[2]- np.linalg.norm(rho2) / c))
tau = tau3 - tau1
print("iteration",i,", rho2:",np.linalg.norm(rho2),"AU. Light travel time =",np.linalg.norm(rho2)/c * 24 * 3600,"seconds")
#constants + standard mog stuff
f1 = fg(tau1,r2guess,r2dotguess,fgflag)[0]
f3 = fg(tau3,r2guess,r2dotguess,fgflag)[0]
g1 = fg(tau1,r2guess,r2dotguess,fgflag)[1]
g3 = fg(tau3,r2guess,r2dotguess,fgflag)[1]
c1 = g3 / (f1*g3 - g1 * f3)
c3 = -g1 / (f1 * g3 - g1 * f3)
d1 = -f3 / (f1 * g3 - f3 * g1)
d3 = f1 / (f1 * g3 - f3 * g1)
rho3 = np.array(rhoHatList[2]) * (c1 * Dj1[3] - Dj2[3] + c3 * Dj3[3]) / c3 / Dj1[0]
rho2 = np.array(rhoHatList[1]) * (c1 * Dj1[2] - Dj2[2] + c3 * Dj3[2]) / -1 / Dj1[0]
rho1 = np.array(rhoHatList[0]) * (c1 * Dj1[1] - Dj2[1] + c3 * Dj3[1]) / c1 / Dj1[0]
r2prev = r2guess
r2guess = rho2 - obs2[3]
r3 = rho3 - obs3[3]
r1 = rho1 - obs1[3]
r2dotguess = d1 * r1 + d3 * r3
i += 1
if i > 1000:
break
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print()
print("Converged!!!!!!")
print("r2 =",r2guess,"-------->",np.linalg.norm(r2guess),"AU" )
print("r2 dot =",r2dotguess * k,"-------->",np.linalg.norm(r2dotguess * k),"AU / day" )
print("rho2 =",np.linalg.norm(rho2),"AU")
print()
print()
#orbital elements time!!!!
print("ORBITAL ELEMENTS")
print("____________")
print()
ecliptic_rotation = np.array([[1,0,0],
[0,cos(obliquity),sin(obliquity)],
[0,-sin(obliquity),cos(obliquity)]])
#rotates to equatorial + calculates values
R = ecliptic_rotation @ r2guess
RDOT = ecliptic_rotation @ r2dotguess
h = np.cross(R,RDOT)
a = get_a(np.linalg.norm(R),np.linalg.norm(RDOT))
e = get_e(a,R,RDOT)
i = get_i(R,RDOT)
Omega = get_Omega(h,i)
w = get_fw(R,i,Omega,e,a,RDOT)[1]
f = get_fw(R,i,Omega,e,a,RDOT)[0]
M = get_M(e,R,a,f)
print("a =",a,"AU")
print("e =",e)
print("i =",degrees(i),"degrees")
print("Longitude of ascending node =",degrees(Omega),"degrees")
print("Argument of perihelion =",degrees(w),"degrees")
print("Mean anomaly =",degrees(M),"degrees at JD =",obs2[2])
print("M0 =",(degrees(M) + degrees(k * (DUEDATE - obs2[2]) / a**(3/2))) % 360,"degrees at 2022 July 24 6:00 UTC")
print("____________")
DUEDATE = julian(2022,7,24,6)
#startup. first reads the number of lines, to appropriately run permutations of observations as shown in test cases
inputFile = input("Please enter input file name: ")
countFileLines = open(inputFile)
numLines = len(countFileLines.readlines())
countFileLines.close()
print("==================================================")
#implements what was shown in test case
if numLines > 3:
print(numLines,"lines detected. MoG will do",numLines - 2,"permutations, keeping observations 1 and",numLines,"in each permutation.")
print("==================================================")
print()
print(numLines)
#runs MoG for each permutation as shown in test case
for i in range(numLines - 2):
print("==================================================")
print("Running observations 1,",i+2,"and",numLines)
print("==================================================")
MOG(0,i+1,numLines-1,inputFile)
