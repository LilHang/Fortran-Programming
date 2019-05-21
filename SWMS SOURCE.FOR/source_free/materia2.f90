! Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

Real Function fk(h, par)

Implicit Double Precision (A-H, O-Z)
Double Precision n, m, ks, kr, kk
Real h, par(10)
Integer ppar

bpar = .5D0
ppar = 2
qr = par(1)
qs = par(2)
qa = par(3)
qm = par(4)
alfa = par(5)
n = par(6)
ks = par(7)
kk = par(8)
qk = par(9)
m = 1.D0 - 1.D0/n
hmin = -1.D300**(1.D0/n)/max(alfa, 1.D0)
hh = max(dble(h), hmin)
qees = dmin1((qs-qa)/(qm-qa), .999999999999999D0)
qeek = dmin1((qk-qa)/(qm-qa), qees)
hs = -1.D0/alfa*(qees**(-1.D0/m)-1.D0)**(1.D0/n)
hk = -1.D0/alfa*(qeek**(-1.D0/m)-1.D0)**(1.D0/n)
If (dble(h)<hk) Then
  qee = (1.D0+(-alfa*hh)**n)**(-m)
  qe = (qm-qa)/(qs-qa)*qee
  qek = (qm-qa)/(qs-qa)*qeek
  ffq = 1.D0 - (1.D0-qee**(1.D0/m))**m
  ffqk = 1.D0 - (1.D0-qeek**(1.D0/m))**m
  If (ffq<=0.D0) ffq = m*qee**(1.D0/m)
  kr = (qe/qek)**bpar*(ffq/ffqk)**ppar*kk/ks
  fk = sngl(max(ks*kr,1.D-37))
  Return
End If
If (dble(h)>=hk .And. dble(h)<hs) Then
  kr = (1.D0-kk/ks)/(hs-hk)*(dble(h)-hs) + 1.D0
  fk = sngl(ks*kr)
End If
If (dble(h)>=hs) fk = sngl(ks)
Return
End Function fk

!***********************************************************************

Real Function fc(h, par)

Implicit Double Precision (A-H, O-Z)
Double Precision n, m
Real h, par(9)

qr = par(1)
qs = par(2)
qa = par(3)
qm = par(4)
alfa = par(5)
n = par(6)
m = 1.D0 - 1.D0/n
hmin = -1.D300**(1.D0/n)/max(alfa, 1.D0)
hh = max(dble(h), hmin)
qees = dmin1((qs-qa)/(qm-qa), .999999999999999D0)
hs = -1.D0/alfa*(qees**(-1.D0/m)-1.D0)**(1.D0/n)
If (dble(h)<hs) Then
  c1 = (1.D0+(-alfa*hh)**n)**(-m-1.D0)
  c2 = (qm-qa)*m*n*(alfa**n)*(-hh)**(n-1.D0)*c1
  fc = sngl(max(c2,1.D-37))
  Return
Else
  fc = 0.0
End If
Return
End Function fc

!***********************************************************************

Real Function fq(h, par)

Implicit Double Precision (A-H, O-Z)
Double Precision n, m
Real h, par(9)

qr = par(1)
qs = par(2)
qa = par(3)
qm = par(4)
alfa = par(5)
n = par(6)
m = 1.D0 - 1.D0/n
hmin = -1.D300**(1.D0/n)/max(alfa, 1.D0)
hh = max(dble(h), hmin)
qees = dmin1((qs-qa)/(qm-qa), .999999999999999D0)
hs = -1.D0/alfa*(qees**(-1.D0/m)-1.D0)**(1.D0/n)
If (dble(h)<hs) Then
  qee = (1.D0+(-alfa*hh)**n)**(-m)
  fq = sngl(max(qa+(qm-qa)*qee,1.D-37))
  Return
Else
  fq = sngl(qs)
End If
Return
End Function fq

!***********************************************************************

Real Function fh(qe, par)

Implicit Double Precision (A-H, O-Z)
Double Precision n, m
Real qe, par(9)

qr = par(1)
qs = par(2)
qa = par(3)
qm = par(4)
alfa = par(5)
n = par(6)

m = 1.D0 - 1.D0/n
hmin = -1.D300**(1.D0/n)/max(alfa, 1.D0)
qeem = (1.D0+(-alfa*hmin)**n)**(-m)
qee = dmin1(dmax1(qe*(qs-qa)/(qm-qa),qeem), .999999999999999D0)
fh = sngl(max(-1.D0/alfa*(qee**(-1.D0/m)-1.D0)**(1.D0/n),-1.D37))
Return
End Function fh

! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
