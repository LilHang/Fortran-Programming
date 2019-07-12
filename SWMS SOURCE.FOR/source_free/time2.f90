! Source file TIME2.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

Subroutine tmcont(dt, dtmaxw, dtopt, dmul, dmul2, dtmin, iter, tprint, tatm, t, tmax, dtmaxc, lminstep, dtinit)
    Logical lminstep

    If (lminstep) Then
      dtmax = amin1(dtmaxw, dtmaxc, dtinit)
      lminstep = .False.
    Else
      dtmax = amin1(dtmaxw, dtmaxc)
    End If
    tfix = amin1(tprint, tatm, tmax)
    If (iter<=3 .And. (tfix-t)>=dmul*dtopt) dtopt = amin1(dtmax, dmul*dtopt)
    If (iter>=7) dtopt = amax1(dtmin, dmul2*dtopt)
    dt = amin1(dtopt, tfix-t)
    dt = amin1((tfix-t)/anint((tfix-t)/dt), dtmax)
    If (tfix-t/=dt .And. dt>(tfix-t)/2.) dt = (tfix-t)/2.
    Return
  End Subroutine tmcont

  !***********************************************************************

  Subroutine setatm(tatm, rtop, rroot, hcrita, width, kxb, numbp, kode, hnew, q, numnp, &
                  gwl0l, qgwlf, freed, cprec, cht, crt, lminstep)

    Logical qgwlf, freed, lminstep
    Dimension width(numbp), kxb(numbp), kode(numnp), hnew(numnp), q(numnp)

    rtopold = rtop
    Read (31, *) tatm, prec, cprec, rsoil, rroot, hcrita, rgwl, gwl, crt, cht
    prec = abs(prec)
    rsoil = abs(rsoil)
    rroot = abs(rroot)
    hcrita = -abs(hcrita)
    hgwl = gwl + gwl0l ! prescribed value of pressure head
    rtop = rsoil - prec ! potential fluid flux
    Do i = 1, numbp
      n = kxb(i)
      k = kode(n)
      If (k==4 .Or. k==-4) Then ! 大气边界
        kode(n) = -4
        q(n) = -width(i)*rtop
        Goto 11
      End If
      If (k==3) Then  ! 变水头边界
        If (abs(hnew(n)-hgwl)>1.E-8) lminstep = .True.
        hnew(n) = hgwl
      End If
      If (k==-3 .And. .Not. qgwlf .And. .Not. freed) Then ! 变流量边界
        If (width(i)>0.) rgwlold = -q(n)/width(i)
        If (abs(rgwlold-rgwl)>1.E-8) lminstep = .True.
        q(n) = -width(i)*rgwl
      End If
    11 End Do
    If ((prec-rsoil)>0.) Then
      cprec = prec/(prec-rsoil)*cprec
    Else
      cprec = 0.
    End If
    If (abs(rtop-rtopold)>1.E-8 .And. abs(rtop)>0.) lminstep = .True.
    Return
  End Subroutine setatm

  !***********************************************************************

  Real Function fqh(gwl, aqh, bqh)
    fqh = -aqh*exp(bqh*abs(gwl))
    Return
  End Function fqh

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
