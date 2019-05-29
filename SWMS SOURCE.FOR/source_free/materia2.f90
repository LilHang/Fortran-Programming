! Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

Real Function fk(h, par)

Implicit Double Precision (A-H, O-Z)
Double Precision n, m, ks, kr, kk
Real h, par(10)
Integer ppar

bpar = .5D0 ! 幂
ppar = 2    ! 幂
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
hs = -1.D0/alfa*(qees**(-1.D0/m)-1.D0)**(1.D0/n) ! 根据vg模型计算hs，下同计算hk
hk = -1.D0/alfa*(qeek**(-1.D0/m)-1.D0)**(1.D0/n)
If (dble(h)<hk) Then
  qee = (1.D0+(-alfa*hh)**n)**(-m) ! vg模型右端项
  qe = (qm-qa)/(qs-qa)*qee   ! 计算结果为式2.16
  qek = (qm-qa)/(qs-qa)*qeek ! 计算结果为式2.17
  ffq = 1.D0 - (1.D0-qee**(1.D0/m))**m ! 计算结果为式2.14
  ffqk = 1.D0 - (1.D0-qeek**(1.D0/m))**m ! 同上
  If (ffq<=0.D0) ffq = m*qee**(1.D0/m) ! 若小于0 则 ？？？
  kr = (qe/qek)**bpar*(ffq/ffqk)**ppar*kk/ks ! 计算结果为式2.13，kr（式中f（theta r/a项结果为0））
  fk = sngl(max(ks*kr,1.D-37)) ! 得到fk，限制其最小值为1.d-37
  Return
End If
If (dble(h)>=hk .And. dble(h)<hs) Then
  kr = (1.D0-kk/ks)/(hs-hk)*(dble(h)-hs) + 1.D0 ! 式2.12提取ks变形，形式与上面保持一致
  fk = sngl(ks*kr)
End If
If (dble(h)>=hs) fk = sngl(ks)
Return
End Function fk

!***********************************************************************
! 比水容量c与h的函数，土壤水分特征曲线theta(h)求导而得
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
  c1 = (1.D0+(-alfa*hh)**n)**(-m-1.D0)     ! c1 为求导中间量，c2为最终求导结果
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
hh = max(dble(h), hmin) ! 限制h不超过hmin
qees = dmin1((qs-qa)/(qm-qa), .999999999999999D0) ! 限制有效饱和度小于“1”
hs = -1.D0/alfa*(qees**(-1.D0/m)-1.D0)**(1.D0/n)  ! 计算饱和含水量对应的水头hs
If (dble(h)<hs) Then  ! 小于hs是非饱和，大于是饱和
  qee = (1.D0+(-alfa*hh)**n)**(-m)
  fq = sngl(max(qa+(qm-qa)*qee,1.D-37)) ! 根据vg模型计算土壤含水率
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
hmin = -1.D300**(1.D0/n)/max(alfa, 1.D0) ! 确定一个最小负压即最大土壤吸力
qeem = (1.D0+(-alfa*hmin)**n)**(-m)      ! 根据hmin确定一个最小有效饱和度
qee = dmin1(dmax1(qe*(qs-qa)/(qm-qa),qeem), .999999999999999D0) ! min为了限制有效饱和度小于“1”，max为了限制有效饱和度大于qeem(最小有效饱和度)
fh = sngl(max(-1.D0/alfa*(qee**(-1.D0/m)-1.D0)**(1.D0/n),-1.D37)) ! 根据vg模型计算计算h并转换为单精度，max是为了限制土壤吸力小于1d37
Return
End Function fh

! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
