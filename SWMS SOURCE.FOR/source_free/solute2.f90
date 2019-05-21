! Source file SOLUTE2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

Subroutine solute(numnp, numel, numeld, mband, mbandd, nmat, t, kode, a, b,&
                 q, hnew, hold, f, x, y, kx, kat, dt, ds, sink, matnum, con, cono, &
                conaxx, conazz, conaxz, vx, vz, dispxx, dispzz, dispxz, chpar, thnew, &
                thold, thsat, ac, fc, gc, qc, conc, listne, cbound, tpulse, numbp, kodcb,&
                kxb, nlevel, cprec, crt, cht, lwat, lupw, wetab, epsi, cumch0, cumch1, &
                cumchr, dtmaxc, peclet, courant, lartd, pecr, lort, a1, b1, numnpd, iad, &
                iadn, iadd, vrv, res, rqi, rq, qq, qi, rqidot, ecnvrg, rcnvrg, acnvrg, mnorth, maxito)

    Double Precision a, b, a1, b1, vrv, res, rqi, rq, qq, qi, rqidot, ecnvrg, rcnvrg, acnvrg
    Logical lwat, lupw, lort, lartd
    Dimension a(mbandd, numnp), b(numnp), q(numnp), hnew(numnp), f(numnp), kx(numeld, 4),&
              matnum(numnp), sink(numnp), ds(numnp), x(numnp), y(numnp), kode(numnp),&
              vx(numnp), vz(numnp), thnew(numnp), thold(numnp), kxb(numbp), conc(numnp),&
              con(numnp), conaxx(numel), conazz(numel), conaxz(numel), qc(numnp), ac(numnp),&
              fc(numnp), gc(numnp), kodcb(numbp), dispxx(numnp), dispzz(numnp), dispxz(numnp),&
              listne(numnp), chpar(10, nmat), cbound(6), thsat(nmat), s(3, 3), bi(3), ci(3),&
              list(3), wx(3), wz(3), wetab(3, 2*numel), vxe(3), vze(3), cono(numnp), &
              hold(numnp), a1(mbandd, numnp), b1(numnp), res(numnp), iad(mbandd, numnp),&
              iadn(numnp), iadd(numnp), vrv(numnp), rqi(numnpd, mnorth), rq(numnp), qq(numnp),&
              rqidot(mnorth), qi(numnpd, mnorth)

  !     Initialisation
    xmul = 1.
    alf = 1. - epsi
    jjj = mband
    If (t>tpulse) Then
      Do i = 1, 4
        cbound(i) = 0.
      End Do
    End If
    Do i = 1, numnp
      b(i) = 0.D0
      If (lort) b1(i) = conc(i)
      qc(i) = 0.
      If (epsi<0.001) Then
        If (lort) jjj = iadd(i)
        a(jjj, i) = 0.D0
      Else
        Do j = 1, mbandd
          a(j, i) = 0.D0
        End Do
      End If
    End Do

    Do level = 1, nlevel
      If (level==nlevel) Then
        eps = epsi
        If (lwat) Call veloc(kat, numnp, numel, numeld, hnew, x, y, kx, listne, con, conaxx, conazz, conaxz, vx, vz)
        Call disper(numnp, nmat, dispxx, dispzz, dispxz, vx, vz, thnew, thsat, chpar, matnum, lartd, pecr, dt)
        Call pecour(numnp, numel, numeld, nmat, x, y, vx, vz, kx, matnum, dispxx, dispzz, chpar, thnew, dt, dtmaxc,&
                   peclet, courant, lupw, lartd, pecr)
        If (lupw .And. lwat) Call wefact(numnp, numel, numeld, x, y, kx, wetab, vx, vz, dispxx, dispzz, dispxz)
      Else
        eps = 1. - epsi
        Call disper(numnp, nmat, dispxx, dispzz, dispxz, vx, vz, thnew, thsat, chpar, matnum, lartd, pecr, dt)
      End If
      Do i = 1, numnp
        m = matnum(i)
        If (level/=nlevel) Then
          If (.Not. lartd .And. .Not. lupw) Then
            dpom = dt/6./(thold(i)+chpar(1,m)*chpar(5,m))
            dispxx(i) = dispxx(i) + vx(i)*vx(i)*dpom
            dispzz(i) = dispzz(i) + vz(i)*vz(i)*dpom
            dispxz(i) = dispxz(i) + vx(i)*vz(i)*dpom
          End If
        Else
          ac(i) = -(thold(i)*alf+thnew(i)*epsi) - chpar(1, m)*chpar(5, m)
          If (.Not. lartd .And. .Not. lupw) Then
            dpom = dt/6./(thnew(i)+chpar(1,m)*chpar(5,m))
            dispxx(i) = dispxx(i) - vx(i)*vx(i)*dpom
            dispzz(i) = dispzz(i) - vz(i)*vz(i)*dpom
            dispxz(i) = dispxz(i) - vx(i)*vz(i)*dpom
          End If
          cs = cbound(5)
          If (cs>conc(i)) cs = conc(i)
          gc(i) = chpar(8, m)*thnew(i) + chpar(1, m)*chpar(9, m) - sink(i)*cs
          fc(i) = chpar(6, m)*thnew(i) + chpar(1, m)*chpar(7, m)*chpar(5, m) + sink(i)
        End If
      End Do
      Do i = 1, numnp
        f(i) = 0.
        If (level==nlevel) ds(i) = 0.
      End Do

  !       Loop on elements
      numsel = 0
      Do n = 1, numel
        caxx = conaxx(n)
        cazz = conazz(n)
        caxz = conaxz(n)
        nus = 4
        If (kx(n,3)==kx(n,4)) nus = 3

  !         Loop on subelements
        Do k = 1, nus - 2
          numsel = numsel + 1
          i = kx(n, 1)
          j = kx(n, k+1)
          l = kx(n, k+2)
          list(1) = i
          list(2) = j
          list(3) = l
          ci(1) = x(l) - x(j)
          ci(2) = x(i) - x(l)
          ci(3) = x(j) - x(i)
          bi(1) = y(j) - y(l)
          bi(2) = y(l) - y(i)
          bi(3) = y(i) - y(j)
          ae = (ci(3)*bi(2)-ci(2)*bi(3))/2.

  !           Calculate Velocities
          ae1 = 1./ae/2.
          ai = caxx*bi(1) + caxz*ci(1)
          aj = caxx*bi(2) + caxz*ci(2)
          ak = caxx*bi(3) + caxz*ci(3)
          If (level==nlevel) Then
            vxx = ae1*(ai*hnew(i)+aj*hnew(j)+ak*hnew(l))
          Else
            vxx = ae1*(ai*hold(i)+aj*hold(j)+ak*hold(l))
          End If
          If (kat>0) vxx = vxx + caxz
          ai = caxz*bi(1) + cazz*ci(1)
          aj = caxz*bi(2) + cazz*ci(2)
          ak = caxz*bi(3) + cazz*ci(3)
          If (level==nlevel) Then
            vzz = ae1*(ai*hnew(i)+aj*hnew(j)+ak*hnew(l))
          Else
            vzz = ae1*(ai*hold(i)+aj*hold(j)+ak*hold(l))
          End If
          If (kat>0) vzz = vzz + cazz
          If (level/=nlevel) Then
            cone = (cono(i)+cono(j)+cono(l))/3.
            vxe(1) = -cono(i)*vxx
            vxe(2) = -cono(j)*vxx
            vxe(3) = -cono(l)*vxx
            vze(1) = -cono(i)*vzz
            vze(2) = -cono(j)*vzz
            vze(3) = -cono(l)*vzz
          Else
            cone = (con(i)+con(j)+con(l))/3.
            vxe(1) = -con(i)*vxx
            vxe(2) = -con(j)*vxx
            vxe(3) = -con(l)*vxx
            vze(1) = -con(i)*vzz
            vze(2) = -con(j)*vzz
            vze(3) = -con(l)*vzz
          End If
          vxee = -cone*vxx
          vzee = -cone*vzz

          If (kat==1) xmul = 2.*3.1416*(x(i)+x(j)+x(l))/3.
          If (level==nlevel) Then
            cs = cbound(5)
            If (abs(cbound(5))>1.E-30) Then
              If (cbound(5)>conc(i)) cs = cs + (conc(i)-cbound(5))/3.
              If (cbound(5)>conc(j)) cs = cs + (conc(j)-cbound(5))/3.
              If (cbound(5)>conc(l)) cs = cs + (conc(l)-cbound(5))/3.
            End If
            rootch = xmul*ae*dt*cs*(sink(i)+sink(j)+sink(l))/3.
            cumch0 = cumch0 - xmul*ae*dt*(gc(i)+gc(j)+gc(l))/3. + rootch
            cumch1 = cumch1 - xmul*ae*dt*((fc(i)-sink(i))*conc(i)+(fc(j)-sink(j))*conc(j)+(fc(l)-sink(l))*conc(l))/3.
            cumchr = cumchr + rootch
          End If
          fmul = xmul*ae/4.
          gce = (gc(i)+gc(j)+gc(l))/3.
          ec1 = (dispxx(i)+dispxx(j)+dispxx(l))/3.
          ec2 = (dispxz(i)+dispxz(j)+dispxz(l))/3.
          ec3 = (dispzz(i)+dispzz(j)+dispzz(l))/3.
          If (level==nlevel) ace = (ac(i)+ac(j)+ac(l))/3.
          fce = (fc(i)+fc(j)+fc(l))/3.
          smul1 = -1./ae/4.*xmul
          smul2 = ae/20.*xmul
          If (lupw) Then
            ns = numsel
            w1 = wetab(1, ns)
            w2 = wetab(2, ns)
            w3 = wetab(3, ns)
            wx(1) = 2*vxe(1)*(w2-w3) + vxe(2)*(w2-2.*w3) + vxe(3)*(2.*w2-w3)
            wx(2) = vxe(1)*(2.*w3-w1) + 2*vxe(2)*(w3-w1) + vxe(3)*(w3-2.*w1)
            wx(3) = vxe(1)*(w1-2.*w2) + vxe(2)*(2.*w1-w2) + 2*vxe(3)*(w1-w2)
            wz(1) = 2*vze(1)*(w2-w3) + vze(2)*(w2-2.*w3) + vze(3)*(2.*w2-w3)
            wz(2) = vze(1)*(2.*w3-w1) + 2*vze(2)*(w3-w1) + vze(3)*(w3-2.*w1)
            wz(3) = vze(1)*(w1-2.*w2) + vze(2)*(2.*w1-w2) + 2*vze(3)*(w1-w2)
          End If
          Do j1 = 1, 3
            i1 = list(j1)
            f(i1) = f(i1) + fmul*(gce+gc(i1)/3.)
            If (level==nlevel) ds(i1) = ds(i1) + fmul*(ace+ac(i1)/3.)
            ibound = 0
            If (kode(i)/=0) Then
              Do id = 1, numbp
                If (kxb(id)==i1 .And. kodcb(id)>0) ibound = 1
              End Do
            End If
            If (ibound==1) qc(i1) = qc(i1) - eps*fmul*(gce+gc(i1)/3.)
            Do j2 = 1, 3
              i2 = list(j2)
              s(j1, j2) = smul1*(ec1*bi(j1)*bi(j2)+ec3*ci(j1)*ci(j2)+ec2*(bi(j1)*ci(j2)+ci(j1)*bi(j2)))
              s(j1, j2) = s(j1, j2) - (bi(j2)/8.*(vxee+vxe(j1)/3.)+ci(j2)/8.*(vzee+vze(j1)/3.))*xmul
              If (lupw) s(j1, j2) = s(j1, j2) - xmul*(bi(j2)/40.*wx(j1)+ci(j2)/40.*wz(j1))
              ic = 1
              If (i1==i2) ic = 2
              s(j1, j2) = s(j1, j2) + smul2*ic*(fce+(fc(i1)+fc(i2))/3.)
              If (level/=nlevel) Then
                b(i1) = b(i1) - alf*s(j1, j2)*conc(i2)
              Else
                If (lort) Then
                  Call find(i1, i2, kk, numnp, mbandd, iad, iadn)
                  ib = kk
                Else
                  ib = mband + i2 - i1
                End If
                a(ib, i1) = a(ib, i1) + epsi*s(j1, j2)
              End If
              If (ibound==1) qc(i1) = qc(i1) - eps*s(j1, j2)*conc(i2)
            End Do
          End Do
        End Do
      End Do

      Do i = 1, numnp
        m = matnum(i)
        If (level/=nlevel) Then
          b(i) = b(i) - alf*f(i)
        Else
          If (lort) jjj = iadd(i)
          a(jjj, i) = a(jjj, i) + ds(i)/dt
          b(i) = b(i) + ds(i)/dt*conc(i) - epsi*f(i)
        End If
      End Do
    End Do

  !     Boundary condition
    Call c_bound(numnp, mband, mbandd, numbp, a, b, q, qc, conc, kode, &
                kxb, kodcb, cbound, cprec, crt, cht, epsi, dt, ds, lort, iadd)

  !     Solve the global matrix equation for transport
    If (epsi<0.001) Then
      Do i = 1, numnp
        If (lort) jjj = iadd(i)
        b(i) = b(i)/a(jjj, i)
      End Do
    Else
      If (lort) Then
        Call ilu(a, numnp, mbandd, iad, iadn, iadd, a1)
        north = 4
        Call orthomin(a, b1, b, numnp, mbandd, numnpd, iad, iadn, iadd, a1, vrv,&
                      res, rqi, rq, qq, qi, rqidot, ecnvrg, rcnvrg, acnvrg, north, mnorth, maxito)
      Else
        Call solvet(a, b, mband, mbandd, numnp)
      End If
    End If
    Do i = 1, numnp
      If (lort) b(i) = b1(i)
      conc(i) = sngl(b(i))
      If (abs(conc(i))<1.E-38) conc(i) = 0.
    End Do
    Return
  End Subroutine solute

  !***********************************************************************

  Subroutine c_bound(numnp, mband, mbandd, numbp, a, b, q, qc, conc, kode, kxb, kodcb, &
                    cbound, cprec, crt, cht, epsi, dt, ds, lort, iadd)

    Double Precision a, b
    Integer ckod
    Logical lort
    Dimension a(mbandd, numnp), b(numnp), q(numnp), conc(numnp), qc(numnp), kode(numnp),&
              kxb(numbp), kodcb(numbp), cbound(6), ds(numnp), iadd(numnp)

    alf = 1. - epsi
    jjj = mband
    Do i = 1, numnp
      If (kode(i)/=0) Then
        Do j = 1, numbp
          If (kxb(j)==i) Then
            If (kodcb(j)>0) Then
              ckod = 1
              If (abs(kode(i))<=2 .Or. abs(kode(i))>=5) cbnd = cbound(kodcb(j))
              If (abs(kode(i))==3) cbnd = cht
              If (abs(kode(i))==4) cbnd = cprec
            Else
              If (q(i)>0.) Then
                ckod = 3
                If (abs(kode(i))==1 .Or. abs(kode(i))>=5) cbnd = cbound(-kodcb(j))
                If (abs(kode(i))==3) cbnd = crt
                If (abs(kode(i))==4) cbnd = cprec
              Else
                ckod = 2
                If (kode(i)==-4) Then
                  ckod = 3
                  cbnd = 0.
                End If
              End If
            End If
            If (abs(kode(i))==2) ckod = 2
            Goto 12
          End If
        End Do

  !     Point source or sink
        If (q(i)<0.) Then
          ckod = 2
        Else
          cbnd = cbound(6)
          ckod = 3
        End If

        12 Continue

  !     Dirichlet boundary condition
        If (ckod==1) Then
          qc(i) = qc(i) + q(i)*(epsi*cbnd+alf*conc(i)) - ds(i)*(cbnd-conc(i))/dt
          If (lort) Then
            a(iadd(i), i) = 1.D30
            b(i) = 1.D30*cbnd
          Else
            Do j = 1, 2*mband - 1
              a(j, i) = 0.D0
            End Do
            a(mband, i) = 1.D0
            b(i) = cbnd
          End If
        End If

  !     Neumann boundary condition
        If (ckod==2) Then
          qc(i) = q(i)*conc(i)
        End If

  !     Cauchy boundary condition
        If (ckod==3) Then
          b(i) = b(i) - q(i)*(cbnd-alf*conc(i))
          If (lort) jjj = iadd(i)
          a(jjj, i) = a(jjj, i) - epsi*q(i)
          qc(i) = q(i)*cbnd
        End If

      End If
    End Do
    Return
  End Subroutine c_bound

  !***********************************************************************

  !     Initial values for solute transport calculation

  Subroutine chinit(numnp, numel, numeld, nmat, x, y, kx, matnum, nlevel, con, hnew, &
                   sink, csink, vx, vz, conaxx, conazz, conaxz, dispxx, dispzz, dispxz,&
                   chpar, theta, thsat, conc, fc, gc, listne, lupw, wetab, dt, dtmaxc,&
                    peclet, courant, kat, lartd, pecr, cono)

    Logical lupw, lartd
    Dimension hnew(numnp), x(numnp), y(numnp), kx(numeld, 4), theta(numnp), sink(numnp),&
              chpar(10, nmat), vx(numnp), vz(numnp), matnum(numnp), con(numnp), conaxx(numel),&
              conazz(numel), conaxz(numel), dispxx(numnp), dispzz(numnp), dispxz(numnp), &
            listne(numnp), fc(numnp), gc(numnp), wetab(3, 2*numel), thsat(nmat), conc(numnp), cono(numnp)

    Do i = 1, numnp
      m = matnum(i)
      If (nlevel==2) Then
        cs = csink
        If (cs>conc(i)) cs = conc(i)
        gc(i) = chpar(8, m)*theta(i) + chpar(1, m)*chpar(9, m) - sink(i)*cs
        fc(i) = chpar(6, m)*theta(i) + chpar(1, m)*chpar(7, m)*chpar(5, m) + sink(i)
      End If
      cono(i) = con(i)
    End Do
    Call veloc(kat, numnp, numel, numeld, hnew, x, y, kx, listne, con, conaxx, conazz, conaxz, vx, vz)
    Call disper(numnp, nmat, dispxx, dispzz, dispxz, vx, vz, theta, thsat, chpar, matnum, lartd, pecr, dt)
    Call pecour(numnp, numel, numeld, nmat, x, y, vx, vz, kx, matnum, dispxx, dispzz, chpar, theta, dt, &
               dtmaxc, peclet, courant, lupw, lartd, pecr)
    If (lupw) Call wefact(numnp, numel, numeld, x, y, kx, wetab, vx, vz, dispxx, dispzz, dispxz)
    Return
  End Subroutine chinit

  !***********************************************************************

  !     Calculate  velocities

  Subroutine veloc(kat, numnp, numel, numeld, hnew, x, y, kx, listne, con, conaxx, conazz, conaxz, vx, vz)

    Integer e
    Dimension hnew(numnp), x(numnp), y(numnp), listne(numnp), con(numnp), kx(numeld, 4), &
              vx(numnp), vz(numnp), conaxx(numel), conazz(numel), conaxz(numel), list(3)

    Do i = 1, numnp
      vx(i) = 0.
      vz(i) = 0.
    End Do
    Do e = 1, numel
      caxx = conaxx(e)
      cazz = conazz(e)
      caxz = conaxz(e)
      ncorn = 4
      If (kx(e,3)==kx(e,4)) ncorn = 3
      Do n = 1, ncorn - 2
        i = kx(e, 1)
        j = kx(e, n+1)
        k = kx(e, n+2)
        list(1) = i
        list(2) = j
        list(3) = k
        vi = y(j) - y(k)
        vj = y(k) - y(i)
        vk = y(i) - y(j)
        wi = x(k) - x(j)
        wj = x(i) - x(k)
        wk = x(j) - x(i)
        area = .5*(wk*vj-wj*vk)
        a = 1./area/2.
        ai = caxx*vi + caxz*wi
        aj = caxx*vj + caxz*wj
        ak = caxx*vk + caxz*wk
        vxx = a*(ai*hnew(i)+aj*hnew(j)+ak*hnew(k))
        If (kat>0) vxx = vxx + caxz
        ai = caxz*vi + cazz*wi
        aj = caxz*vj + cazz*wj
        ak = caxz*vk + cazz*wk
        vzz = a*(ai*hnew(i)+aj*hnew(j)+ak*hnew(k))
        If (kat>0) vzz = vzz + cazz
        Do m = 1, 3
          l = list(m)
          vx(l) = vx(l) - con(l)*vxx
          vz(l) = vz(l) - con(l)*vzz
        End Do
      End Do
    End Do
    Do i = 1, numnp
      vx(i) = vx(i)/listne(i)
      vz(i) = vz(i)/listne(i)
    End Do
    Return
  End Subroutine veloc

  !***********************************************************************

  !     Calculate the dispersion coefficient

  Subroutine disper(numnp, nmat, dispxx, dispzz, dispxz, vx, vz, theta, thsat, chpar, matnum, lartd, pecr, dt)

    Logical lartd
    Dimension vx(numnp), vz(numnp), theta(numnp), chpar(10, nmat), dispxx(numnp), dispzz(numnp), dispxz(numnp),&
              matnum(numnp), thsat(nmat)

    Do i = 1, numnp
      m = matnum(i)
      tau = theta(i)**(7./3.)/thsat(m)**2
      vabs = sqrt(vx(i)*vx(i)+vz(i)*vz(i))
      dif = theta(i)*chpar(2, m)*tau
      displ = chpar(3, m)
      dispt = chpar(4, m)
      If (lartd .And. vabs>1.E-20) displ = amax1(displ, vabs*dt/(theta(i)+chpar(1,m)*chpar(5,m))/pecr-dif/vabs)
      dispxx(i) = dif
      dispzz(i) = dif
      dispxz(i) = 0.
      If (vabs>1.E-20) Then
        dispxx(i) = displ*vx(i)*vx(i)/vabs + dispt*vz(i)*vz(i)/vabs + dif
        dispzz(i) = displ*vz(i)*vz(i)/vabs + dispt*vx(i)*vx(i)/vabs + dif
        dispxz(i) = (displ-dispt)*vx(i)*vz(i)/vabs
      End If
    End Do
    Return
  End Subroutine disper


  !***********************************************************************

  !     Calculate upstream weighing factors

  Subroutine wefact(numnp, numel, numeld, x, y, kx, wetab, vx, vz, dispxx, dispzz, dispxz)

    Integer e
    Dimension x(numnp), y(numnp), kx(numeld, 4), vx(numnp), vz(numnp), dispxx(numnp), dispzz(numnp),&
              dispxz(numnp), wetab(3, 2*numel), beta(3), list(3)

    tanh(z) = (exp(z)-exp(-z))/(exp(z)+exp(-z))

    numsel = 0
    Do e = 1, numel
      ncorn = 4
      If (kx(e,3)==kx(e,4)) ncorn = 3
      Do n = 1, ncorn - 2
        numsel = numsel + 1
        m1 = kx(e, 1)
        m2 = kx(e, n+1)
        m3 = kx(e, n+2)
        a = y(m2) - y(m1)
        b = x(m2) - x(m1)
        beta(1) = atan2(a, b)
        a = y(m3) - y(m2)
        b = x(m3) - x(m2)
        beta(2) = atan2(a, b)
        a = y(m1) - y(m3)
        b = x(m1) - x(m3)
        beta(3) = atan2(a, b)
        list(1) = m1
        list(2) = m2
        list(3) = m3
        Do j = 1, 3
          k = j - 1
          If (k==0) k = 3
          wetab(k, numsel) = 0.
          m1 = list(j)
          jp1 = j + 1
          If (j==3) jp1 = 1
          m2 = list(jp1)
          vxx = (vx(m1)+vx(m2))/2.
          vzz = (vz(m1)+vz(m2))/2.
          If (abs(vxx)<1.E-30 .And. abs(vzz)<1.E-30) Goto 11
          betav = atan2(vzz, vxx)
          delta = abs(betav-beta(j))
          If (delta>0.314 .And. abs(delta-3.1416)>0.314) Goto 11
          aleng = sqrt((x(m2)-x(m1))**2+(y(m2)-y(m1))**2)
          cbeta = cos(beta(j))
          sbeta = sin(beta(j))
          val = vxx*cbeta + vzz*sbeta
          vv = sqrt(vxx*vxx+vzz*vzz)
          dll = (dispxx(m1)+dispxx(m2))/2.
          dlt = (dispxz(m1)+dispxz(m2))/2.
          dtt = (dispzz(m1)+dispzz(m2))/2.
          dal = abs(dll*cbeta*cbeta+2.0*cbeta*sbeta*dlt+dtt*sbeta*sbeta)
          vel = val*aleng
          disp = 2.0*dal
          aa = 11.
          If (abs(disp)>1.E-30) aa = abs(vel/disp)
          If (abs(disp)<1.E-30 .Or. abs(vel)<0.001*vv .Or. abs(aa)>10.) Then
            If (abs(vel)<0.001*vv) wetab(k, numsel) = 0.0
            If (vel>0.001*vv) wetab(k, numsel) = 1.0
            If (vel<-0.001*vv) wetab(k, numsel) = -1.0
          Else
            wetab(k, numsel) = 1.0/tanh(vel/disp) - disp/vel
          End If
        11 End Do
      End Do
    End Do
    Return
  End Subroutine wefact

  !************************************************************************

  !     Calculate the maximum local Peclet and Courant numbers

  Subroutine pecour(numnp, numel, numeld, nmat, x, y, vx, vz, kx, matnum, dispxx, dispzz, chpar,&
                    theta, dt, dtmaxc, peclet, courant, lupw, lartd, pecr)

    Logical lupw, lartd
    Dimension kx(numeld, 4), x(numnp), y(numnp), vx(numnp), vz(numnp), matnum(numnp), dispxx(numnp),&
              dispzz(numnp), theta(numnp), chpar(10, nmat), bi(3), ci(3)

    peclet = 0.
    courant = 0.
    dtmaxc = 1.E+30
    Do n = 1, numel
      nus = 4
      If (kx(n,3)==kx(n,4)) nus = 3
      Do k = 1, nus - 2
        pecx = 99999.
        pecy = 99999.
        dt1 = 1.E+30
        dt2 = 1.E+30
        i = kx(n, 1)
        j = kx(n, k+1)
        l = kx(n, k+2)
        ci(1) = x(l) - x(j)
        ci(2) = x(i) - x(l)
        ci(3) = x(j) - x(i)
        bi(1) = y(j) - y(l)
        bi(2) = y(l) - y(i)
        bi(3) = y(i) - y(j)
        delx = amax1(abs(ci(1)), abs(ci(2)), abs(ci(3)))
        dely = amax1(abs(bi(1)), abs(bi(2)), abs(bi(3)))
        dxe = (dispxx(i)+dispxx(j)+dispxx(l))/3.
        dze = (dispzz(i)+dispzz(j)+dispzz(l))/3.
        vxe = abs(vx(i)+vx(j)+vx(l))/3.
        vze = abs(vz(i)+vz(j)+vz(l))/3.
        If (dxe>1.E-20) pecx = vxe*delx/dxe
        If (dze>1.E-20) pecy = vze*dely/dze
        If (pecx/=99999.) peclet = amax1(peclet, pecx)
        If (pecy/=99999.) peclet = amax1(peclet, pecy)
        peclet = amin1(peclet, 99999.)

        vxmax = amax1(abs(vx(i))/theta(i), abs(vx(j))/theta(j), abs(vx(l))/theta(l))
        vzmax = amax1(abs(vz(i))/theta(i), abs(vz(j))/theta(j), abs(vz(l))/theta(l))
        r1 = 1. + chpar(1, matnum(i))*chpar(5, matnum(i))/theta(i)
        r2 = 1. + chpar(1, matnum(j))*chpar(5, matnum(j))/theta(j)
        r3 = 1. + chpar(1, matnum(l))*chpar(5, matnum(l))/theta(l)
        rmin = amin1(r1, r2, r3)
        courx = vxmax*dt/delx/rmin
        coury = vzmax*dt/dely/rmin
        courant = amax1(courant, courx, coury)

        cour1 = 1.0
        cour2 = 1.0
        If (.Not. lupw .And. .Not. lartd) Then
          If (pecx/=99999.) cour1 = amin1(1., pecr/amax1(0.5,pecx))
          If (pecy/=99999.) cour2 = amin1(1., pecr/amax1(0.5,pecy))
        End If
        If (vxmax>1.E-20) dt1 = cour1*delx*rmin/vxmax
        If (vzmax>1.E-20) dt2 = cour2*dely*rmin/vzmax
        dtmaxc = amin1(dtmaxc, dt1, dt2)

      End Do
    End Do
    Return
  End Subroutine pecour

  !***********************************************************************

  !     Solve the global matrix equation for transport

  Subroutine solvet(a, b, mband, mbandd, numnp)

    Double Precision a, b, p, c, sum
    Dimension a(mbandd, numnp), b(numnp)

    n1 = numnp - 1
    Do k = 1, n1
      p = 1.D0/a(mband, k)
      kk = k + 1
      kc = mband
      Do i = kk, numnp
        kc = kc - 1
        If (kc<=0) Goto 12
        c = -p*a(kc, i)
        a(kc, i) = c
        ii = kc + 1
        l = kc + mband - 1
        Do j = ii, l
          jj = j + mband - kc
          a(j, i) = a(j, i) + c*a(jj, k)
        End Do
      End Do
    12 End Do
    Do i = 2, numnp
      jj = mband + 1 - i
      ii = 1
      If (jj<=0) Then
        jj = 1
        ii = i - mband + 1
      End If
      sum = 0.
      Do j = jj, mband - 1
        sum = sum + a(j, i)*b(ii)
        ii = ii + 1
      End Do
      b(i) = b(i) + sum
    End Do
    b(numnp) = b(numnp)/a(mband, numnp)
    Do k = 1, n1
      i = numnp - k
      jj = i
      m = min0(2*mband-1, mband+k)
      sum = 0.
      Do j = mband + 1, m
        jj = jj + 1
        sum = sum + a(j, i)*b(jj)
      End Do
      b(i) = (b(i)-sum)/a(mband, i)
    End Do
    Return
  End Subroutine solvet

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
