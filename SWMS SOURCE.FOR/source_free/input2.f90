! Source file INPUT2.FOR |||||||||||||||||||||||||||||||||||||||||||||||

Subroutine basinf(kat, maxit, tolth, tolh, lwat, lchem, atminf, shortf, seepf, checkf, fluxf, freed, drainf)

    Character *72 hed
    Character *5 lunit, tunit, munit
    Logical lwat, lchem, checkf, atminf, shortf, seepf, fluxf, freed, drainf
    Dimension iu(11)
    Data iu/50, 71, 72, 75, 76, 77, 78, 79, 80, 81, 82/

    Read (30, *)
    Read (30, *)
    Read (30, *) hed
    Read (30, *)
    Read (30, *) lunit, tunit, munit
    Read (30, *)
    Read (30, *) kat
    Read (30, *)
    Read (30, *) maxit, tolth, tolh
    Read (30, *)
    Read (30, *) lwat, lchem, checkf, shortf, fluxf, atminf, seepf, freed, drainf
    Do i = 1, 11
      Write (iu(i), *) hed
      Write (iu(i), *)
      Write (iu(i), *) 'Program SWMS_2D'
  !        call getdat(ii,imonth,iday)
  !        call gettim(ihours,mins,isecs,ii)
  !        write(IU(i),100) iday,imonth,ihours,mins,isecs
      If (atminf) Then
        Write (iu(i), *) 'Time dependent boundary conditions'
      Else
        Write (iu(i), *) 'Time independent boundary conditions'
      End If
      If (kat==0) Write (iu(i), 110)
      If (kat==1) Write (iu(i), 120)
      If (kat==2) Write (iu(i), 130)
      Write (iu(i), *) 'Units: L = ', lunit, ', T = ', tunit, ', M = ', munit
    End Do
    Write (*, *) '-----------------------------------------------------'
    Write (*, *) '|                                                   |'
    Write (*, *) '|                     SWMS_2D                       |'
    Write (*, *) '|                                                   |'
    Write (*, *) '|     Code for simulating water flow and solute     |'
    Write (*, *) '|       transport in two-dimensional variably       |'
    Write (*, *) '|             saturated porous media                |'
    Write (*, *) '|                                                   |'
    Write (*, *) '|                  version 1.22                     |'
    Write (*, *) '|          Last modified: January, 1994             |'
    Write (*, *) '|                                                   |'
    Write (*, *) '-----------------------------------------------------'
    Write (*, *)
    Write (*, *) hed
    If (kat==0) Write (*, 110)
    If (kat==1) Write (*, 120)
    If (kat==2) Write (*, 130)
    Write (50, 140) maxit, tolth, tolh


    100 Format (' Date: ', I3, '.', I2, '.', '    Time: ', I3, ':', I2, ':', I2)
    110 Format (' Horizontal plane flow, V = L*L')
    120 Format (' Axisymmetric flow, V = L*L*L')
    130 Format (' Vertical plane flow, V = L*L')
    140 Format (/' Max. number of iterations           ', I4/' Absolute water content tolerance [-]',&
               F8.5/' Absolute pressure head tolerance [L]', F8.5/)
    Return
End Subroutine basinf

  !***********************************************************************

  Subroutine matin(nmatd, nmat, nlay, par, htab1, htabn)

    Real k
    Dimension par(10, nmatd), qe(10)
    Data qe/1., .99, .90, .85, .75, .65, .50, .35, .20, .10/

    Write (*, *) 'reading material information'
    imax = 10
    Read (30, *)
    Read (30, *)
    Read (30, *) nmat, nlay, htab1, htabn, npar
    If (nmat>nmatd) Then
      Write (*, *) 'Dimension in NMatD is exceeded'
      Stop
    End If
    htab1 = -amin1(abs(htab1), abs(htabn))
    htabn = -amax1(abs(htab1), abs(htabn))
    Read (30, *)
    Write (50, 110)
    Do m = 1, nmat
      Read (30, *)(par(i,m), i=1, npar)
      Write (50, 120) m, (par(i,m), i=1, npar)
    End Do
    Write (50, 130)
    Do m = 1, nmat
      Write (50, *)
      Do i = 1, imax
        h = fh(qe(i), par(1,m))
        k = fk(h, par(1,m))
        c = fc(h, par(1,m))
        q = fq(h, par(1,m))
        Write (50, 140) m, qe(i), q, h, c, k
      End Do
    End Do

    110 Format (/' MatNum,        thR    thS    tha    thm      alpha         n          Ks          Kk          thk'/)
    120 Format (I5, 8X, 4F7.3, 16E12.3)
    130 Format (//' MatNum          Qe     Q        h         C         K')
    140 Format (I5, 8X, 2F7.3, F10.3, E10.2, E12.3)
    Return
  End Subroutine matin

  !***********************************************************************

  Subroutine genmat(ntab, ntabd, nmat, thr, hsat, par, htab, contab, captab, consat, thetab, thsat)

    Dimension par(10, nmat), contab(ntabd, nmat), captab(ntabd, nmat), thetab(ntabd, nmat), htab(ntab), consat(nmat), &
                 hsat(nmat), thr(nmat), thsat(nmat)

    Write (*, *) 'generating materials'
    htab1 = htab(1)
    htabn = htab(ntab)
    dlh = (alog10(-htabn)-alog10(-htab1))/(ntab-1)
    Do i = 1, ntab
      alh = alog10(-htab1) + (i-1)*dlh
      htab(i) = -10**alh
    End Do
    Do m = 1, nmat
      hr = fh(0.0, par(1,m))
      hsat(m) = fh(1.0, par(1,m))
      consat(m) = fk(0.0, par(1,m))
      thr(m) = fq(hr, par(1,m))
      thsat(m) = fq(0.0, par(1,m))
      Do i = 1, ntab
        contab(i, m) = fk(htab(i), par(1,m))
        captab(i, m) = fc(htab(i), par(1,m))
        thetab(i, m) = fq(htab(i), par(1,m))
      End Do
    End Do
    Return
  End Subroutine genmat

  !***********************************************************************

  Subroutine tmin(tinit, tmax, tatm, told, dt, dtmax, dmul, dmul2, dtmin, tprint, t, dtopt, dtold, atminf)

    Logical atminf
    Dimension tprint(50)

    Write (*, *) 'reading time information'
    Read (30, *)
    Read (30, *)
    Read (30, *) dt, dtmin, dtmax, dmul, dmul2, mpl
    Read (30, *)
    Read (30, *)(tprint(i), i=1, mpl)
    dtopt = dt
    dtold = dt
    If (.Not. atminf) Then
      tmax = tprint(mpl)
      tatm = tmax
    End If
    tprint(mpl+1) = tmax
    told = tinit
    t = tinit + dt
    Return
  End Subroutine tmin

  !***********************************************************************

  Subroutine seepin(nseepd, numspd, nseep, nsp, np)

    Dimension nsp(nseepd), np(nseepd, numspd)

    Write (*, *) 'reading seepage face information'
    Read (30, *)
    Read (30, *)
    Read (30, *) nseep
    If (nseep>nseepd) Then
      Write (*, *) 'Dimension in NSeepD is exceeded'
      Stop
    End If
    Read (30, *)
    Read (30, *)(nsp(i), i=1, nseep)
    Do i = 1, nseep
      If (nsp(i)>numspd) Then
        Write (*, *) 'Dimension in NumSPD is exceeded'
        Stop
      End If
    End Do
    Read (30, *)
    Do i = 1, nseep
      Read (30, *)(np(i,j), j=1, nsp(i))
    End Do
    Return
  End Subroutine seepin

  !***********************************************************************

  Subroutine drainin(ndr, ndrd, neldrd, numel, nd, ned, keldr, efdim, conaxx, conaxz, conazz)

    Integer e
    Dimension nd(ndrd), ned(ndrd), keldr(ndrd, neldrd), efdim(2, ndrd), conaxx(numel), conazz(numel), conaxz(numel)

    Write (*, *) 'reading drainage information'
    Read (30, *)
    Read (30, *)
    Read (30, *) ndr, drcorr
    If (ndr>ndrd) Then
      Write (*, *) 'Dimension in NDrD is exceeded'
      Stop
    End If
    Read (30, *)
    Read (30, *)(nd(i), i=1, ndr)
    Read (30, *)
    Read (30, *)(ned(i), i=1, ndr)
    Do i = 1, ndr
      If (ned(i)>neldrd) Then
        Write (*, *) 'Dimension in NElDrD is exceeded'
        Stop
      End If
    End Do
    Read (30, *)
    Do i = 1, ndr
      Read (30, *)(efdim(i,j), j=1, 2)
    End Do
    Read (30, *)
    Do i = 1, ndr
      Read (30, *)(keldr(i,j), j=1, ned(i))
    End Do
    Do i = 1, ndr
      rho = efdim(i, 2)/efdim(i, 1)
      a = (1.+0.405*rho**(-4))/(1.-0.405*rho**(-4))
      b = (1.+0.163*rho**(-8))/(1.-0.163*rho**(-8))
      c = (1.+0.067*rho**(-12))/(1.-0.067*rho**(-12))
      red = 376.7/(138.*alog10(rho)+6.48-2.34*a-0.48*b-0.12*c)/drcorr
      Do j = 1, ned(i)
        e = keldr(i, j)
        conaxx(e) = conaxx(e)*red
        conaxz(e) = conaxz(e)*red
        conazz(e) = conazz(e)*red
      End Do
    End Do
    Return
  End Subroutine drainin

  !***********************************************************************

  Subroutine nodinf(numnp, numel, ij, numbp, numnpd, numeld, numbpd, numkd, nobs, nobsd, kode, q, conc, hnew, hold, htemp,&
                   x, y, matnum, beta, axz, bxz, dxz, checkf)

    Logical checkf
    Dimension kode(numnpd), q(numnpd), hold(numnpd), x(numnpd), y(numnpd), hnew(numnpd), htemp(numnpd), matnum(numnpd), &
                  beta(numnpd), axz(numnpd), bxz(numnpd), dxz(numnpd), conc(numnpd)

    Write (*, *) 'reading nodal information'
    Read (32, *)
    Read (32, *)
    Read (32, *) numnp, numel, ij, numbp, nobs
    If (numnp>numnpd) Then
      Write (*, *) 'Dimension in NumNPD is exceeded'
      Stop
    Else If (numel>numeld) Then
      Write (*, *) 'Dimension in NumElD is exceeded'
      Stop
    Else If (numbp>numbpd) Then
      Write (*, *) 'Dimension in NumBPD is exceeded'
      Stop
    Else If (nobs>nobsd) Then
      Write (*, *) 'Dimension in NObsD is exceeded'
      Stop
    End If
    Read (32, *)
    npr = 0
    k = 0

    11 k = k + 1
    Read (32, *) n, kode(n), x(n), y(n), hold(n), conc(n), q(n), matnum(n), beta(n), axz(n), bxz(n), dxz(n)
    If (kode(n)>numkd) Then
      Write (*, *) 'Dimension in NumKD is exceeded'
      Stop
    End If
    If (n-k) 12, 15, 13
    12 Write (*, 130) n
    Stop
    13 deno = n - k + 1
    dx = (x(n)-x(npr))/deno
    dy = (y(n)-y(npr))/deno
    dp = (hold(n)-hold(npr))/deno
    dconc = (conc(n)-conc(npr))/deno
    dbeta = (beta(n)-beta(npr))/deno
    da = (axz(n)-axz(npr))/deno
    db = (bxz(n)-bxz(npr))/deno
    dd = (dxz(n)-dxz(npr))/deno
    14 x(k) = x(k-1) + dx
    y(k) = y(k-1) + dy
    hold(k) = hold(k-1) + dp
    conc(k) = conc(k-1) + dconc
    beta(k) = beta(k-1) + dbeta
    axz(k) = axz(k-1) + da
    bxz(k) = bxz(k-1) + db
    dxz(k) = dxz(k-1) + dd
    matnum(k) = matnum(k-1)
    kode(k) = kode(k-1)
    q(k) = q(k-1)
    k = k + 1
    If (k<n) Goto 14
    15 npr = n
    If (k<numnp) Goto 11

    Do n = 1, numnp
      hnew(n) = hold(n)
      htemp(n) = hold(n)
    End Do
    If (checkf) Then
      Write (50, 110)
      Do n = 1, numnp
        Write (50, 120) n, kode(n), x(n), y(n), hold(n), q(n), conc(n), matnum(n), beta(n), axz(n), bxz(n), dxz(n)
      End Do
    End If

    110 Format (////' NODAL POINT INFORMATION'//'  NODE NO.', 6X, 'KODE', 7X, 'X,R', 12X, 'Y,Z', 11X, &
              '.PSI.', 12X, 'Q', 11X, 'Conc'/)
    120 Format (2I10, 5E15.6, I10, 4F7.3)
    130 Format (' ERROR IN NodInf AT N=', I5)
    Return
  End Subroutine nodinf

  !***********************************************************************

  Subroutine elemin(numel, numeld, numnp, kx, laynum, conaxx, conazz, conaxz, checkf, listne, ij, mband, mbandd, lchem, lort)

    Logical checkf, lchem, lconst, lort
    Integer e
    Dimension kx(numeld, 4), conaxx(numel), conazz(numel), conaxz(numel), laynum(numel), listne(numnp)

    Write (*, *) 'reading element information'
    num = 0
    Read (32, *)
    Read (32, *)
    Do e = 1, numel
      If (num-e) 11, 14, 12
      11 Read (32, *) num, (kx(num,i), i=1, 4), conaxz(num), conaxx(num), conazz(num), laynum(num)
      If (kx(num,4)==0) kx(num, 4) = kx(num, 3)
      If (num==e) Goto 14
      12 Do i = 1, 4
        kx(e, i) = kx(e-1, i) + 1
      End Do
      conaxx(e) = conaxx(e-1)
      conazz(e) = conazz(e-1)
      conaxz(e) = conaxz(e-1)
      laynum(e) = laynum(e-1)
    14 End Do
    aa = 3.141592654/180.
    Do e = 1, numel
      ang = aa*conaxz(e)
      caxx = conaxx(e)
      cazz = conazz(e)
      conaxx(e) = caxx*cos(ang)*cos(ang) + cazz*sin(ang)*sin(ang)
      conazz(e) = caxx*sin(ang)*sin(ang) + cazz*cos(ang)*cos(ang)
      conaxz(e) = (caxx-cazz)*sin(ang)*cos(ang)
    End Do
    If (checkf) Then
      Write (50, 110)
      Do e = 1, numel
        Write (50, 120) e, (kx(e,i), i=1, 4), conaxz(e), conaxx(e), conazz(e), laynum(e)
      End Do
    End If

    Do i = 1, numnp
      listne(i) = 0
    End Do
    Do e = 1, numel
      ncorn = 4
      If (kx(e,3)==kx(e,4)) ncorn = 3
      Do n = 1, ncorn - 2
        i = kx(e, 1)
        j = kx(e, n+1)
        k = kx(e, n+2)
        listne(i) = listne(i) + 1
        listne(j) = listne(j) + 1
        listne(k) = listne(k) + 1
      End Do
    End Do

    lort = .False.
    lconst = .True.
    mband = 1
    Do e = 1, numel
      nus = 4
      If (kx(e,3)==kx(e,4)) nus = 3
      Do kk = 1, nus - 2
        mb = 1
        i = kx(e, 1)
        j = kx(e, kk+1)
        k = kx(e, kk+2)
        If (abs(i-j)>mb) mb = abs(i-j)
        If (abs(i-k)>mb) mb = abs(i-k)
        If (abs(j-k)>mb) mb = abs(j-k)
        If (mb>mband) mband = mb
        If (e==1 .And. kk==1) Then
          mb1 = mb
        Else
          If (mb1/=mb) lconst = .False.
        End If
      End Do
    End Do
    mband = mband + 1
    If (mband>mbandd .Or. (lchem .And. 2*mband-1>mbandd)) lort = .True.
    If (.Not. lconst) ij = numnp
    If (mband>10 .Or. numnp>200) lort = .True.

    110 Format (////' Element Information'//' Element    C O R N E R    N O D E S   ConAxz    ConAxx    ConAzz  LayNum'/)
    120 Format (I6, I9, 3I6, E14.3, 2F8.3, I5)
    Return
  End Subroutine elemin

  !***********************************************************************

  Subroutine geomin(numkd, numnp, numbp, nobs, nobsd, swidth, width, kode, kxb, rlen, node)

    Character *50 text1
    Dimension kxb(numbp), width(numbp), swidth(numkd), kode(numnp), node(nobsd)

    Write (*, *) 'reading geometric information'
    Read (32, *)
    Read (32, *)
    Read (32, *)(kxb(i), i=1, numbp)
    Read (32, *)
    Read (32, *)(width(i), i=1, numbp)
    Read (32, *)
    Read (32, *) rlen
    Do i = 1, numkd
      swidth(i) = 0.
    End Do
    Do i = 1, numbp
      n = kxb(i)
      j = iabs(kode(n))
      If (j==0) Goto 12
      swidth(j) = swidth(j) + width(i)
    12 End Do
    If (nobs>0) Then
      Read (32, *)
      Read (32, *)(node(i), i=1, nobs)
      text1 = '    hNew     theta     conc    '
      Write (92, 110)(node(j), j=1, nobs)
      Write (92, 120)(text1, i=1, nobs)
    End If
    110 Format (///14X, 5(11X,'Node(',I3,')',11X))
    120 Format (/'     time     ', 5(A31))
    Return
  End Subroutine geomin

  !***********************************************************************

  Subroutine atmin(gwl0l, sinkf, qgwlf, tinit, tmax, aqh, bqh, hcrits, maxal)

    Logical sinkf, qgwlf

    Write (*, *) 'reading atmospheric information'
    Read (31, *)
    Read (31, *)
    Read (31, *)
    Read (31, *)
    Read (31, *) sinkf, qgwlf
    Read (31, *)
    Read (31, *) gwl0l, aqh, bqh
    Read (31, *)
    Read (31, *) tinit, maxal
    Read (31, *)
    Read (31, *) hcrits
    Read (31, *)
    Do i = 1, maxal - 1
      Read (31, *)
    End Do
    Read (31, *) tmax
    Rewind 31
    Do i = 1, 12
      Read (31, *)
    End Do
    Return
  End Subroutine atmin

  !***********************************************************************

  Subroutine sinkin(nmat, numel, numnp, numeld, kat, kx, x, y, p0, poptm, p2h, p2l, p3, r2h, r2l, beta)

    Integer e
    Dimension poptm(nmat), beta(numnp), kx(numeld, 4), x(numnp), y(numnp)

    Write (*, *) 'reading sink information'
    Read (30, *)
    Read (30, *)
    Read (30, *) p0, p2h, p2l, p3, r2h, r2l
    Read (30, *)
    Read (30, *)(poptm(i), i=1, nmat)
    p0 = -abs(p0)
    p2l = -abs(p2l)
    p2h = -abs(p2h)
    p3 = -abs(p3)
    xmul = 1.
    sbeta = 0.
    Do e = 1, numel
      nus = 4
      If (kx(e,3)==kx(e,4)) nus = 3
      Do k = 1, nus - 2
        i = kx(e, 1)
        j = kx(e, k+1)
        l = kx(e, k+2)
        cj = x(i) - x(l)
        ck = x(j) - x(i)
        bj = y(l) - y(i)
        bk = y(i) - y(j)
        ae = (ck*bj-cj*bk)/2.
        If (kat==1) xmul = 2.*3.1416*(x(i)+x(j)+x(l))/3.
        betae = (beta(i)+beta(j)+beta(l))/3.
        sbeta = sbeta + xmul*ae*betae
      End Do
    End Do
    Do i = 1, numnp
      beta(i) = beta(i)/sbeta
    End Do
    Return
  End Subroutine sinkin

  !***********************************************************************

  Subroutine chemin(nmat, numbp, cbound, chpar, epsi, tpulse, kodcb, nlevel, lupw, lartd, pecr)

    Logical lupw, lartd
    Dimension chpar(10, nmat), kodcb(numbp), cbound(6)

    Write (*, *) 'reading transport information'
    nlevel = 1
    Write (50, 110)
    Read (30, *)
    Read (30, *)
    Read (30, *) epsi, lupw, lartd, pecr
    pecr = amax1(pecr, 0.01)
    If (epsi<0.999) nlevel = 2
    Read (30, *)
    If (lupw) Then
      Write (50, 120)
    Else
      Write (50, 130)
      If (lartd) Write (50, 140) pecr
    End If
    Write (50, 150)
    Do m = 1, nmat
      Read (30, *)(chpar(j,m), j=1, 9)
      Write (50, 160)(chpar(j,m), j=1, 9)
    End Do
    Read (30, *)
    Read (30, *)(kodcb(i), i=1, numbp)
    Read (30, *)
    Read (30, *)(cbound(i), i=1, 6)
    Write (50, 170)(cbound(i), i=1, 6)
    Read (30, *)
    Read (30, *) tpulse
    Write (50, 180) tpulse

    110 Format (/' Solute transport information'/1X, 28('='))
    120 Format (/' Upstream weighting finite-element method')
    130 Format (/' Galerkin finite-element method')
    140 Format (/' Artificial dispersion is added when Peclet number is', ' higher than', F10.3)
    150 Format (/'  Bulk.d.   Difus.         Disper.        Adsorp.    ', 'SinkL1    SinkS1    SinkL0    SinkS0')
    160 Format (10E10.3)
    170 Format (/'   Conc1     Conc2     Conc3     Conc4     cSink     ', 'cWell'/8E10.3)
    180 Format (/' tPulse =   ', F15.3)
    Return
  End Subroutine chemin

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
