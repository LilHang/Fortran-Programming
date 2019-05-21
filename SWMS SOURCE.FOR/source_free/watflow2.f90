! Source file WATFLOW2.FOR |||||||||||||||||||||||||||||||||||||||||||||

Subroutine watflow(numnp, numel, numeld, ntab, ntabd, mband, mbandd, nmat, nseep, nseepd, &
                numspd, nsp, np, numbp, itcum, maxit, iter, kode, kat, t, dt, dtmin, dtopt,&
                 dtold, told, hcrita, hcrits, tolth, tolh, rlen, width, rtop, vmeanr, hmeanr, &
                 atminf, sinkf, seepf, qgwlf, freed, par, htab, contab, captab, thetab, hnew, hold,&
                  htemp, thr, thsat, thnew, thold, matnum, con, cap, consat, axz, bxz, dxz, hsat, &
                  a, b, q, f, x, y, kx, sink, ds, beta, conaxx, conazz, conaxz, kxb, explic, gwl0l, &
                  aqh, bqh, lwat, tlevel, lort, drainf, nd, ndr, ndrd, rroot, p0, poptm, p2h, p2l, p3,&
                   r2h, r2l, cono, a1, b1, numnpd, iad, iadn, iadd, vrv, res, rqi, rq, qq, qi, rqidot, &
                   ecnvrg, rcnvrg, acnvrg, mnorth, maxito)

    Logical atminf, sinkf, seepf, explic, itcrit, lwat, qgwlf, freed, lort, drainf
    Double Precision a, b, a1, b1, vrv, res, rqi, rq, qq, qi, rqidot, ecnvrg, rcnvrg, acnvrg
    Integer tlevel
    Dimension a(mbandd, numnp), b(numnp), kode(numnp), q(numnp), f(numnp), hnew(numnp), htemp(numnp), hold(numnp),&
             consat(nmat), htab(ntab), contab(ntabd, nmat), captab(ntabd, nmat), con(numnp), cap(numnp), x(numnp), y(numnp),&
             matnum(numnp), kx(numeld, 4), kxb(numbp), par(10, nmat), width(numbp), conaxx(numel), conazz(numel), conaxz(numel),&
             hsat(nmat), np(nseepd, numspd), nsp(nseepd), ds(numnp), beta(numnp), axz(numnp), bxz(numnp), dxz(numnp), sink(numnp),&
             thr(nmat), thsat(nmat), thetab(ntabd, nmat), thnew(numnp), nd(ndrd), thold(numnp), poptm(nmat), cono(numnp), &
             a1(mbandd, numnp), b1(numnp), res(numnp), iad(mbandd, numnp), iadn(numnp), iadd(numnp), vrv(numnp),&
              rqi(numnpd, mnorth),rq(numnp), qq(numnp), rqidot(mnorth), qi(numnpd, mnorth)

    If (lwat .And. tlevel/=1) Then
      Do i = 1, numnp
        hold(i) = hnew(i)
        thold(i) = thnew(i)
        cono(i) = con(i)
        If (kode(i)<1) Then
          htemp(i) = hnew(i) + (hnew(i)-hold(i))*dt/dtold
          hnew(i) = htemp(i)
        Else
          htemp(i) = hnew(i)
        End If
      End Do
    End If

    11 Continue

    iter = 0
    explic = .False.

  ! --- Beginning of iteration loop --------------------------------------

    12 Continue
    If (sinkf) & !.and..not.lWat.and.Iter.ne.0)
      Call setsnk(numnp, nmat, matnum, hnew, rroot, sink, p0, poptm, p2h, p2l, p3, r2h, r2l, beta, rlen)
    Call setmat(numnp, ntab, ntabd, nmat, htab, contab, captab, hnew, hold, matnum, par, con, cap, consat, axz, bxz, dxz, hsat,&
                htemp, explic, thetab, thsat, thr, thnew)
    Call reset(a, b, kode, q, hnew, f, con, cap, x, y, kx, numnp, numel, numeld, mband, mbandd, kat, dt, iter, sinkf, sink, ds,&
               beta, rlen, vmeanr, hmeanr, conaxx, conazz, conaxz, thnew, thold, lwat, p3, lort, iad, iadn, iadd, b1)
    Call shift(numnp, numbp, nseepd, numspd, nseep, nsp, np, hnew, hold, q, kode, rtop, width, kxb, hcrita, hcrits, seepf, atminf, &
               explic, qgwlf, freed, gwl0l, aqh, bqh, con, cono, drainf, nd, ndr, ndrd, lwat, iter)
    Call dirich(a, b, kode, hnew, numnp, mband, mbandd, lort, iadd)
    If (lort) Then
      Call ilu(a, numnp, mbandd, iad, iadn, iadd, a1)
      north = 0
      Call orthomin(a, b1, b, numnp, mbandd, numnpd, iad, iadn, iadd, a1, vrv, res, rqi, rq, qq, qi,&
                   rqidot, ecnvrg, rcnvrg, acnvrg,north, mnorth, maxito)
    Else
      Call solve(a, b, numnp, mband, mbandd)
    End If

    Do i = 1, numnp
      htemp(i) = hnew(i)
      If (lort) b(i) = b1(i)
      hnew(i) = sngl(b(i))
      If (abs(hnew(i))>1.E+10) hnew(i) = sign(1.E+10, hnew(i))
    End Do
    iter = iter + 1
    itcum = itcum + 1
    If (explic) Goto 18

  !     Test for convergence
    itcrit = .True.
    Do i = 1, numnp
      m = matnum(i)
      epsth = 0.
      epsh = 0.
      If (htemp(i)<hsat(m) .And. hnew(i)<hsat(m)) Then
        th = thnew(i) + cap(i)*(hnew(i)-htemp(i))/(thsat(m)-thr(m))/dxz(i)
        epsth = abs(thnew(i)-th)
      Else
        epsh = abs(hnew(i)-htemp(i))
      End If
      If (epsth>tolth .Or. epsh>tolh) Then
        itcrit = .False.
        Goto 15
      End If
    End Do
    15 Continue

    If (.Not. itcrit) Then
      If (iter<maxit .Or. (.Not. lwat .And. iter<5*maxit)) Then
        Goto 12
      Else If (dt<=dtmin) Then
        explic = .True.
        Do i = 1, numnp
          hnew(i) = hold(i)
          htemp(i) = hold(i)
        End Do
        Goto 12
      Else If (.Not. lwat) Then
        Write (*, *) ' No steady state solution found'
        Stop
      Else
        Do i = 1, numnp
          hnew(i) = hold(i)
          htemp(i) = hold(i)
        End Do
        dt = amax1(dt/3., dtmin)
        dtopt = dt
        t = told + dt
        Goto 11
      End If
    End If

  ! --- End of iteration loop --------------------------------------------

    18 Continue

    If (.Not. lwat .And. tlevel==1) Then
      Write (*, 101) iter
      Do i = 1, numnp
        hold(i) = hnew(i)
        thold(i) = thnew(i)
      End Do
    End If
    Return

    101 Format (' Steady state was reached after', I6, ' iterations.')
  End Subroutine watflow

  !***********************************************************************

  Subroutine reset(a, b, kode, q, hnew, f, con, cap, x, y, kx, numnp, numel, numeld, mband, mbandd,&
                  kat, dt, iter, sinkf, sink, ds, beta, rlen, vmeanr, hmeanr, conaxx, conazz, conaxz, &
                  thnew, thold, lwat, p3, lort, iad, iadn, iadd, b1)

    Logical sinkf, lwat, lort
    Double Precision a, b, b1
    Dimension a(mbandd, numnp), b(numnp), q(numnp), hnew(numnp), f(numnp), con(numnp), cap(numnp), x(numnp), &
              y(numnp), kx(numeld, 4), conaxx(numel), conazz(numel), conaxz(numel), kode(numnp), sink(numnp), ds(numnp),&
              beta(numnp), thnew(numnp), thold(numnp), e(3, 3), iloc(3), bi(3), ci(3), b1(numnp), iad(mbandd, numnp),&
              iadn(numnp), iadd(numnp)

  !     Initialisation
    xmul = 1.
    If (iter==0) Then
      vmeanr = 0.
      hmeanr = 0.
      arear = 0.
    End If
    Do i = 1, numnp
      b(i) = 0.D0
      If (lort) b1(i) = hnew(i)
      f(i) = 0.
      If (iter==0) ds(i) = 0.
      Do j = 1, mbandd
        a(j, i) = 0.D0
      End Do
    End Do

  !     Loop on elements
    Do n = 1, numel
      condi = conaxx(n)
      condj = conazz(n)
      condk = conaxz(n)
      nus = 4
      If (kx(n,3)==kx(n,4)) nus = 3

  !       Loop on subelements
      Do k = 1, nus - 2
        i = kx(n, 1)
        j = kx(n, k+1)
        l = kx(n, k+2)
        iloc(1) = 1
        iloc(2) = k + 1
        iloc(3) = k + 2
        ci(1) = x(l) - x(j)
        ci(2) = x(i) - x(l)
        ci(3) = x(j) - x(i)
        bi(1) = y(j) - y(l)
        bi(2) = y(l) - y(i)
        bi(3) = y(i) - y(j)
        ae = (ci(3)*bi(2)-ci(2)*bi(3))/2.
        cone = (con(i)+con(j)+con(l))/3.
        If (kat==1) xmul = 2.*3.1416*(x(i)+x(j)+x(l))/3.
        amul = xmul*cone/4./ae
        bmul = xmul*cone/2.
        fmul = xmul*ae/12.
        betae = (beta(i)+beta(j)+beta(l))/3.
        If (sinkf .And. betae>0. .And. iter==0) Then
          sinke = (sink(i)+sink(j)+sink(l))/3.
          If (hnew(i)>p3) ds(i) = ds(i) + fmul*(3.*sinke+sink(i))
          If (hnew(j)>p3) ds(j) = ds(j) + fmul*(3.*sinke+sink(j))
          If (hnew(l)>p3) ds(l) = ds(l) + fmul*(3.*sinke+sink(l))
          hnewe = (hnew(i)+hnew(j)+hnew(l))/3.
          vmeanr = vmeanr + xmul*ae*sinke/rlen
          hmeanr = hmeanr + xmul*hnewe*ae
          arear = arear + xmul*ae
        End If
        Do i = 1, 3
          ig = kx(n, iloc(i))
          f(ig) = f(ig) + fmul*4.
          If (kat>=1) b(ig) = b(ig) + bmul*(condk*bi(i)+condj*ci(i))
          Do j = 1, 3
            jg = kx(n, iloc(j))
            e(i, j) = condi*bi(i)*bi(j) + condk*(bi(i)*ci(j)+ci(i)*bi(j)) + condj*ci(i)*ci(j)
            If (lort) Then
              Call find(ig, jg, kk, numnp, mbandd, iad, iadn)
              a(kk, ig) = a(kk, ig) + amul*e(i, j)
            Else
              ib = ig - jg + 1
              If (ib>=1) a(ib, jg) = a(ib, jg) + amul*e(i, j)
            End If
          End Do
        End Do
      End Do

    End Do
    If (arear>0. .And. iter==0) hmeanr = hmeanr/arear

  !     Determine boundary fluxes
    Do n = 1, numnp
      If (kode(n)<1) Goto 19
      qn = sngl(b(n)) + ds(n) + f(n)*(thnew(n)-thold(n))/dt
      If (lort) Then
        Do j = 1, iadn(n)
          qn = qn + sngl(a(j,n))*hnew(iad(j,n))
        End Do
      Else
        qn = qn + sngl(a(1,n))*hnew(n)
        Do j = 2, mband
          k = n - j + 1
          If (k>=1) Then
            qn = qn + sngl(a(j,k))*hnew(k)
          End If
          k = n + j - 1
          If (k<=numnp) Then
            qn = qn + sngl(a(j,n))*hnew(k)
          End If
        End Do
      End If
      q(n) = qn
    19 End Do

  !     Complete construction of RHS vector and form effective matrix
    Do i = 1, numnp
      If (.Not. lwat) f(i) = 0.
      j = 1
      If (lort) j = iadd(i)
      a(j, i) = a(j, i) + f(i)*cap(i)/dt
      b(i) = f(i)*cap(i)*hnew(i)/dt - f(i)*(thnew(i)-thold(i))/dt + q(i) - b(i) - ds(i)
    End Do
    Return
  End Subroutine reset

  !***********************************************************************

  Subroutine dirich(a, b, kode, hnew, numnp, mband, mbandd, lort, iadd)

    Logical lort
    Double Precision a, b
    Dimension a(mbandd, numnp), b(numnp), kode(numnp), hnew(numnp), iadd(numnp)

    Do n = 1, numnp
      If (kode(n)<1) Goto 12
      If (lort) Then
        a(iadd(n), n) = 10.D30
        b(n) = 10.D30*hnew(n)
      Else
        Do m = 2, mband
          k = n - m + 1
          If (k>0) Then
            b(k) = b(k) - a(m, k)*hnew(n)
            a(m, k) = 0.D0
          End If
          l = n + m - 1
          If (numnp-l>=0) Then
            b(l) = b(l) - a(m, n)*hnew(n)
          End If
          a(m, n) = 0.D0
        End Do
        a(1, n) = 1.D0
        b(n) = hnew(n)
      End If
    12 End Do
    Return
  End Subroutine dirich

  !***********************************************************************

  Subroutine shift(numnp, numbp, nseepd, numspd, nseep, nsp, np, hnew, hold, q, kode, ei, &
                   width, kxb, hcrita, hcrits, seepf, atminf, explic, qgwlf, freed, gwl0l, aqh, bqh, con, &
                   cono, drainf, nd, ndr, ndrd, lwat, iter)

    Logical seepf, atminf, explic, qgwlf, freed, drainf, lwat
    Dimension kode(numnp), q(numnp), hnew(numnp), hold(numnp), width(numbp),&
             kxb(numbp), np(nseepd, numspd), nsp(nseepd), con(numnp), cono(numnp), nd(ndrd)

  !     Modify conditions on seepage faces
    If (seepf) Then
      Do i = 1, nseep
  !          iCheck=0
        ns = nsp(i)
        Do j = 1, ns
          n = np(i, j)
          If (kode(n)==-2) Then
            If (hnew(n)<0.) Then
  !                iCheck=1
            Else
              kode(n) = 2
              hnew(n) = 0.
            End If
          Else
  !              if(iCheck.gt.0.or.Q(n).ge.0.) then
            If (q(n)>=0.) Then
              kode(n) = -2
              q(n) = 0.
  !                iCheck=1
            End If
          End If
        End Do
      End Do
    End If

  !     Modify conditions in drainage node
    If (drainf) Then
      Do i = 1, ndr
        n = nd(i)
        If (kode(n)==-5) Then
          If (hnew(n)>=0.) Then
            kode(n) = 5
            hnew(n) = 0.
          End If
        Else
          If (q(n)>=0.) Then
            kode(n) = -5
            q(n) = 0.
          End If
        End If
      End Do
    End If

  !     Modify potential surface flux boundaries
    If (atminf) Then
      Do i = 1, numbp
        n = kxb(i)
        k = kode(n)
        If (explic .And. iabs(k)==4) Then
          kode(n) = -iabs(k)
          Goto 14
        End If

  !         Critical surface pressure on ...
        If (k==4) Then
          If (abs(q(n))>abs(-ei*width(i)) .Or. q(n)*(-ei)<=0) Then
            kode(n) = -4
            q(n) = -ei*width(i)
          End If
          Goto 14
        End If

  !         Surface flux on ...
        If (k==-4 .And. iter/=0) Then
          If (hnew(n)<=hcrita) Then
            kode(n) = 4
            hnew(n) = hcrita
            Goto 14
          End If
          If (hnew(n)>=hcrits) Then
            kode(n) = 4
            hnew(n) = hcrits
          End If
        End If

  !         Time variable flux
        If (k==-3) Then
          If (lwat) Then
            If (qgwlf) q(n) = -width(i)*fqh(hold(n)-gwl0l, aqh, bqh)
          Else
            If (qgwlf) q(n) = -width(i)*fqh(hnew(n)-gwl0l, aqh, bqh)
          End If
        End If
      14 End Do
    End If

  !     Free Drainage
    If (freed) Then
      Do i = 1, numbp
        n = kxb(i)
        k = kode(n)
        If (lwat) Then
          If (k==-3) q(n) = -width(i)*cono(n)
        Else
          If (k==-3) q(n) = -width(i)*con(n)
        End If
      End Do
    End If
    Return
  End Subroutine shift

  !***********************************************************************

  Subroutine setmat(numnp, ntab, ntabd, nmat, htab, contab, captab, hnew, hold, matnum, par,&
                    con, cap, consat, axz, bxz, dxz, hsat, htemp, explic, thetab, thsat, thr, theta)

    Logical explic
    Dimension htab(ntab), contab(ntabd, nmat), captab(ntabd, nmat), hnew(numnp), hold(numnp), matnum(numnp), &
              par(10, nmat), con(numnp), cap(numnp), consat(nmat), hsat(nmat), axz(numnp), bxz(numnp), dxz(numnp),&
             htemp(numnp), thetab(ntabd, nmat), thsat(nmat), thr(nmat), theta(numnp)

    alh1 = alog10(-htab(1))
    dlh = (alog10(-htab(ntab))-alh1)/(ntab-1)
    Do i = 1, numnp
      m = matnum(i)
      hi1 = amin1(hsat(m), htemp(i)/axz(i))
      hi2 = amin1(hsat(m), hnew(i)/axz(i))
      If (explic) hi2 = hi1
      him = 0.1*hi1 + 0.9*hi2
      If (hi1>=hsat(m) .And. hi2>=hsat(m)) Then
        ci = consat(m)
      Else If (him>=htab(ntab) .And. him<=htab(1)) Then
        it = int((alog10(-him)-alh1)/dlh) + 1
        s1 = (contab(it+1,m)-contab(it,m))/(htab(it+1)-htab(it))
        ci = contab(it, m) + s1*(him-htab(it))
      Else
        ci = fk(him, par(1,m))
      End If
      con(i) = bxz(i)*ci
      hi1 = hold(i)/axz(i)
      hi2 = hnew(i)/axz(i)
      If (explic) hi2 = hi1
      If (hi2>=hsat(m)) Then
        ci = 0
        ti = thsat(m)
      Else If (hi2>=htab(ntab) .And. hi2<=htab(1)) Then
        it = int((alog10(-hi2)-alh1)/dlh) + 1
        s1 = (captab(it+1,m)-captab(it,m))/(htab(it+1)-htab(it))
        s2 = (thetab(it+1,m)-thetab(it,m))/(htab(it+1)-htab(it))
        ci = captab(it, m) + s1*(hi2-htab(it))
        ti = thetab(it, m) + s2*(hi2-htab(it))
      Else
        ci = fc(hi2, par(1,m))
        ti = fq(hi2, par(1,m))
      End If
      cap(i) = ci*dxz(i)/axz(i)
      theta(i) = thr(m) + (ti-thr(m))*dxz(i)
    End Do
    Return
  End Subroutine setmat

  !***********************************************************************

  Subroutine solve(a, b, numnp, mband, mbandd)

    Double Precision a, b, c
    Dimension a(mbandd, numnp), b(numnp)

  !     Reduction
    Do n = 1, numnp
      Do m = 2, mband
        If (dabs(a(m,n))<1.D-30) Goto 12
        c = a(m, n)/a(1, n)
        i = n + m - 1
        If (i>numnp) Goto 12
        j = 0
        Do k = m, mband
          j = j + 1
          a(j, i) = a(j, i) - c*a(k, n)
        End Do
        a(m, n) = c
        b(i) = b(i) - a(m, n)*b(n)
      12 End Do
      b(n) = b(n)/a(1, n)
    End Do

  !     Back substitution
    n = numnp
    14 Do k = 2, mband
      l = n + k - 1
      If (l>numnp) Goto 16
      b(n) = b(n) - a(k, n)*b(l)
    End Do
    16 n = n - 1
    If (n>0) Goto 14
    Return
  End Subroutine solve

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
