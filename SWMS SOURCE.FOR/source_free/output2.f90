! Source file OUTPUT2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

Subroutine tlinf(numnp, numbp, kode, q, hnew, cumq, width, swidth, kxb, t, dt, tlevel, &
               shortf, tprint, iter, itcum, rtop, rroot, vmeanr, hmeant, hmeanr, hmeang,&
                atminf, sinkf, cumqrt, cumqrr, cumqvr, numkd, hmean, vmean, lwat, lchem, &
                rlen, peclet, courant, wcumt, wcuma)

    Integer tlevel
    Logical shortf, sinkf, atminf, lwat, lchem
    Dimension q(numnp), kode(numnp), kxb(numbp), swidth(numkd), hnew(numnp), width(numbp), cumq(numkd), hmean(numkd), vmean(numkd)

    If (tlevel==1) Then
      If (.Not. lchem) Then
        Write (70, 110)
      Else
        If (lwat) Then
          Write (70, 120)
        Else
          Write (70, 130)
        End If
      End If
      Write (71, 140)
      Write (77, 150)
      If (lwat) Write (78, 160)
    End If
    If (lwat .Or. tlevel==1) Then
      Do i = 1, numkd
        vmean(i) = 0.
        hmean(i) = 0.
      End Do
      Do i = 1, numbp
        n = kxb(i)
        j = iabs(kode(n))
        If (j==0) Goto 12
        hmean(j) = hmean(j) + hnew(n)*width(i)/swidth(j)
        If (j==4) vmean(j) = vmean(j) - q(n)/swidth(j)
      12 End Do
      hmeant = hmean(4)
      hmeang = hmean(3)
      Do i = 1, numnp
        j = iabs(kode(i))
        wcuma = wcuma + abs(q(i))*dt
        If (j/=0 .And. j/=4) Then
          vmean(j) = vmean(j) - q(i)
        End If
      End Do
      If (.Not. lwat .And. tlevel==1) Then
        Write (71, 170) t, rtop, rroot, vmean(4), vmeanr, vmean(3), vmean(1), vmean(2), (vmean(i), i=5, numkd)
        Write (77, 180) t, hmean(4), hmeanr, hmean(3), hmean(1), hmean(2), (hmean(i), i=5, numkd)
      End If
    End If
    If (lwat) Then
      wcuma = wcuma + abs(vmeanr*dt*rlen)
      wcumt = cumqvr
      Do j = 1, numkd
        If (j==4) Then
          cumq(j) = cumq(j) + vmean(j)*dt*swidth(4)
        Else
          cumq(j) = cumq(j) + vmean(j)*dt
        End If
        wcumt = wcumt + cumq(j)
      End Do
      cumqrt = cumqrt + rtop*dt*swidth(4)
      cumqrr = cumqrr + rroot*dt*rlen
      cumqvr = cumqvr + vmeanr*dt*rlen
    End If
    If (.Not. shortf .Or. abs(tprint-t)<0.001*dt) Then
      If (lwat) Then
        Write (71, 170) t, rtop, rroot, vmean(4), vmeanr, vmean(3), vmean(1), vmean(2), (vmean(i), i=5, numkd)
        Write (77, 180) t, hmean(4), hmeanr, hmean(3), hmean(1), hmean(2), (hmean(i), i=5, numkd)
        Write (78, 190) t, cumqrt, cumqrr, cumq(4), cumqvr, cumq(3), cumq(1), cumq(2), (cumq(i), i=5, numkd)
      End If
      If (lchem) Then
        If (lwat) Then
          Write (70, 200) tlevel, t, dt, iter, itcum, peclet, courant
        Else
          Write (70, 210) tlevel, t, dt, peclet, courant
        End If
      Else
        Write (70, 200) tlevel, t, dt, iter, itcum
      End If
    End If
    If (lwat) Then
      If (sinkf) Then
        Write (*, 220) t, iter, itcum, cumq(4), cumqvr, cumq(3), hmean(4), hmeanr, hmean(3)
      Else
        If (atminf) Then
          Write (*, 220) t, iter, itcum, cumq(4), cumq(1), cumq(3), hmean(4), hmean(1), hmean(3)
        Else
          Write (*, 220) t, iter, itcum, vmean(1), cumq(1), cumq(2), hmean(1), hmean(2)
        End If
      End If
    End If

    110 Format (//' TLevel   Time         dt      Iter  ItCum'/)
    120 Format (//' TLevel   Time         dt      Iter  ItCum  Peclet   ', 'Courant'/)
    130 Format (//' TLevel   Time         dt        Peclet   Courant'/)
    140 Format (' All fluxes (v) are positive out of the region'//&
                '     Time       rAtm       rRoot      vAtm       vRoot      ',&
                'vKode3     vKode1     vSeep      vKode5     vKode6 ...'/&
                '      [T]      [L/T]       [L/T]     [L/T]       [L/T]      ',&
                ' [V/T]      [V/T]     [V/T]       [V/T]      [V/T]'/)
    150 Format (//&
                '     Time         hAtm       hRoot     hKode3     hKode1 ',&
                '     hSeep     hKode5     hKode6 ... '/&
                '      [T]          [L]        [L]        [L]        [L]  ',&
                '      [L]        [L]        [L]'/)
    160 Format (' All cumulative fluxes (CumQ) are positive out of the region'//&
                '     Time     CumQAP      CumQRP     CumQA      CumQR     CumQ3',&
                '       CumQ1      CumQS      CumQ5       CumQ6 ....'/&
                '      [T]       [V]         [V]       [V]        [V]       [V] ',&
                '        [V]        [V]        [V]         [V]'/)
    170 Format (F12.4, 9E11.3)
    180 Format (F12.4, 9F11.1)
    190 Format (F12.4, 9E11.3)
    200 Format (I5, 2E12.3, I5, I6, 2F10.3)
    210 Format (I5, 2E12.3, 2F10.3)
    220 Format (F14.4, I3, I6, 1X, 3E11.3, 3F7.0)
    Return
  End Subroutine tlinf

  !***********************************************************************

  Subroutine alinf(t, cumq, hmeant, hmeanr, hmeang, alevel, cumqrt, cumqrr, cumqvr, numkd)

    Integer alevel
    Dimension cumq(numkd)

    If (alevel==1) Write (72, 110)
    Write (72, 120) t, cumqrt, cumqrr, cumq(4), cumqvr, cumq(3), hmeant, hmeanr, hmeang, alevel
    Return

    110 Format (' All cumulative fluxes (CumQ) are positive out of the region'//&
                '      Time      CumQAP     CumQRP     CumQA      CumQR ',&
                '     CumQ3        hAtm       hRoot     hKode3    A-level'/&
                '      [T]         [V]        [V]       [V]        [V]  ',&
                '      [V]          [L]        [L]        [L] '/)
    120 Format (F12.4, 5E11.3, 3F11.1, I8)
  End Subroutine alinf

  !***********************************************************************

  Subroutine subreg(numel, numeld, numnp, nmat, hnew, tho, thn, x, y, matnum, laynum, kx,&
                     kat, t, dt, nlay, plevel, lwat, lchem, conc, chpar, wcuma, wcumt, &
                     ccuma, ccumt, wvoli, cvoli, watin, solin)

    Logical lwat, lchem
    Integer plevel, e
    Dimension hnew(numnp), x(numnp), y(numnp), matnum(numnp), conc(numnp), kx(numeld, 4), &
                 chpar(10, nmat), laynum(numel), tho(numnp), thn(numnp), watin(numel), solin(numel),&
                 area(10), hmean(10), subvol(10), subcha(10), cmean(10), consub(10)

    xmul = 1.
    atot = 0.
    If (lwat .Or. plevel<=1) Then
      volume = 0.
      change = 0.
      htot = 0.
      deltw = 0.
    End If
    If (lchem) Then
      ctot = 0.
      convol = 0.
      deltc = 0.
    End If
    Do i = 1, nlay
      area(i) = 0.
      If (lwat .Or. plevel<=1) Then
        subvol(i) = 0.
        subcha(i) = 0.
        hmean(i) = 0.
      End If
      If (lchem) Then
        consub(i) = 0.
        cmean(i) = 0.
      End If
    End Do
    Do e = 1, numel
      lay = laynum(e)
      wel = 0.
      cel = 0.
      nus = 4
      If (kx(e,3)==kx(e,4)) nus = 3
      Do k = 1, nus - 2
        i = kx(e, 1)
        j = kx(e, k+1)
        l = kx(e, k+2)
        mi = matnum(i)
        mj = matnum(j)
        mk = matnum(l)
        cj = x(i) - x(l)
        ck = x(j) - x(i)
        bj = y(l) - y(i)
        bk = y(i) - y(j)
        If (kat==1) xmul = 2.*3.1416*(x(i)+x(j)+x(l))/3.
        ae = xmul*(ck*bj-cj*bk)/2.
        area(lay) = area(lay) + ae
        If (lwat .Or. plevel<=1) Then
          he = (hnew(i)+hnew(j)+hnew(l))/3.
          vnewe = ae*(thn(i)+thn(j)+thn(l))/3.
          volde = ae*(tho(i)+tho(j)+tho(l))/3.
          volume = volume + vnewe
          wel = wel + vnewe
          change = change + (vnewe-volde)/dt
          subvol(lay) = subvol(lay) + vnewe
          subcha(lay) = subcha(lay) + (vnewe-volde)/dt
          htot = htot + he*ae
          hmean(lay) = hmean(lay) + he*ae
        End If
        If (lchem) Then
          ce = (conc(i)+conc(j)+conc(l))/3.
          cnewe = ae*((thn(i)+chpar(1,mi)*chpar(5,mi))*conc(i)+(thn(j)+chpar(1,mj)*chpar(5,mj))&
                  *conc(j)+(thn(l)+chpar(1,mk)*chpar(5,mk))*conc(l))/3.
          convol = convol + cnewe
          cel = cel + cnewe
          consub(lay) = consub(lay) + cnewe
          ctot = ctot + ce*ae
          cmean(lay) = cmean(lay) + ce*ae
        End If
        If (k==nus-2) Then
          If (plevel==0) Then
            If (lwat) watin(e) = wel
            If (lchem) solin(e) = cel
          Else
            If (lwat) deltw = deltw + abs(watin(e)-wel)
            If (lchem) deltc = deltc + abs(solin(e)-cel)
          End If
        End If
      End Do
    End Do
    Do lay = 1, nlay
      If (lwat .Or. plevel<=1) hmean(lay) = hmean(lay)/area(lay)
      If (lchem) cmean(lay) = cmean(lay)/area(lay)
      atot = atot + area(lay)
    End Do
    If (lwat .Or. plevel<=1) htot = htot/atot
    If (lchem) ctot = ctot/atot
    If (plevel==0) Write (80, 110)
    Write (80, 120) t, (i, i=1, nlay)
    Write (80, 130) atot, (area(i), i=1, nlay)
    If (lwat .Or. plevel<=1) Then
      Write (80, 140) volume, (subvol(i), i=1, nlay)
      Write (80, 150) change, (subcha(i), i=1, nlay)
      Write (80, 160) htot, (hmean(i), i=1, nlay)
    End If
    If (lchem) Then
      Write (80, 170) convol, (consub(i), i=1, nlay)
      Write (80, 180) ctot, (cmean(i), i=1, nlay)
    End If

  !     Mass balance calculation
    If (plevel==0) Then
      wvoli = volume
      If (lchem) cvoli = convol
    Else
      If (lwat) Then
        wbalt = volume - wvoli + wcumt
        Write (80, 190) wbalt
        ww = amax1(deltw, wcuma)
        If (ww>=1.E-25) Then
          wbalr = abs(wbalt)/ww*100.
          Write (80, 200) wbalr
        End If
      End If
      If (lchem) Then
        cbalt = convol - cvoli + ccumt
        Write (80, 210) cbalt
        cc = amax1(deltc, ccuma)
        If (cc>=1.E-25) Then
          cbalr = abs(cbalt)/cc*100.
          Write (80, 220) cbalr
        End If
      End If
    End If


    110 Format (/' Time [T]             Total     Sub-region number ...')
    120 Format (/F12.4, 16X, 10(I7,4X))
    130 Format (' Area    [V]       ', E11.3, 10E11.3)
    140 Format (' Volume  [V]       ', E11.3, 10E11.3)
    150 Format (' InFlow  [V/T]     ', E11.3, 10E11.3)
    160 Format (' hMean   [L]       ', E11.3, 10F11.1)
    170 Format (' ConcVol [VM/L3]   ', E11.3, 10E11.3)
    180 Format (' cMean   [M/L3]    ', E11.3, 10E11.3)
    190 Format (' WatBalT [V]       ', E11.3)
    200 Format (' WatBalR [%]       ', F11.3)
    210 Format (' CncBalT [VM/L3]   ', E11.3)
    220 Format (' CncBalR [%]       ', F11.3)
    Return
  End Subroutine subreg

  !***********************************************************************

  Subroutine bouout(numnp, numbp, t, hnew, theta, q, width, kxb, kode, x, y, conc)

    Dimension hnew(numnp), q(numnp), width(numbp), theta(numnp), kxb(numbp), kode(numnp), x(numnp), y(numnp), conc(numnp)

    Write (79, 110) t
    ii = 0
    Do i = 1, numnp
      If (kode(i)/=0) Then
        Do j = 1, numbp
          n = kxb(j)
          If (n==i) Then
            ii = ii + 1
            v = -q(i)/width(j)
            Write (79, 120) ii, i, x(i), y(i), kode(i), q(i), v, hnew(i), theta(i), conc(i)
            Goto 12
          End If
        End Do
        ii = ii + 1
        Write (79, 130) ii, i, x(i), y(i), kode(i), q(i), hnew(i), theta(i), conc(i)
      End If
    12 End Do


    110 Format (//' Time:', F12.4//&
            '    i    n    x      z    Code     Q          v       ',&
            '    h       th      Conc'/                              &
            '                                 [V/T]      [L/T]     ',&
            '   [L]     [-]     [M/L3]'/)
    120 Format (2I5, 2F7.1, I5, 2E11.3, F11.1, F7.3, E10.3)
    130 Format (2I5, 2F7.1, I5, E11.3, 11X, F11.1, F7.3, E10.3)
    Return
  End Subroutine bouout

  !***********************************************************************

  Subroutine solinf(numnp, kode, qc, t, dt, tlevel, shortf, tprint, numkd, smean, &
                  chems, cumch0, cumch1, cumchr, ccuma, ccumt, lwat)

    Integer tlevel
    Logical shortf, lwat
    Dimension qc(numnp), kode(numnp), chems(numkd), smean(numkd)

    Do i = 1, numkd
      smean(i) = 0.
    End Do
    Do i = 1, numnp
      j = iabs(kode(i))
      If (j/=0) Then
        smean(j) = smean(j) - qc(i)
      End If
    End Do
    ccuma = abs(cumch0) + abs(cumch1) + abs(cumchr)
    ccumt = cumch0 + cumch1 + cumchr
    Do j = 1, numkd
      chems(j) = chems(j) + smean(j)*dt
      ccumt = ccumt + chems(j)
      ccuma = ccuma + abs(chems(j))
    End Do
    If (tlevel==1) Write (74, 110)
    If (.Not. shortf .Or. abs(tprint-t)<0.001*dt) Write (74, 120) t, cumch0, cumch1, cumchr, &
       (chems(j), j=1, numkd), (smean(i), i=1, numkd)
    If (.Not. lwat) Write (*, 130) t, tlevel, cumch0, cumch1, chems(1)
    Return

    110 Format (' All solute fluxes (SMean) and cumulative solute fluxes',&
     ' (ChemS) are positive out of the region'//&
     '     Time     CumCh0     CumCh1     CumChR   ', 20('-'), &
     '  ChemS(i),i=1,NumKD  ', 22('-'), '  ', 21('-'), '  SMean(j),j=1,NumKD ', 22('-')&
     /'      [T]    [VM/L3]    [VM/L3]    [VM/L3]', 31(' '), '[VM/L3]', 59(' '),&
      '[VM/T/L3]'/)
    120 Format (F10.2, 15E11.3)
    130 Format (F14.4, I6, 1X, 3E11.3, 2X, 2E11.3)
  End Subroutine solinf

  !**********************************************************************

  Subroutine obsnod(t, numnp, nobs, nobsd, node, hnew, thnew, conc)

    Dimension node(nobsd), hnew(numnp), thnew(numnp), conc(numnp)

    Write (92, 110) t, (hnew(node(i)), thnew(node(i)), conc(node(i)), i=1, nobs)
    Return

    110 Format (F11.3, 5(F11.3,F9.4,E11.3))
  End Subroutine obsnod

  !***********************************************************************

  Subroutine hout(hnew, x, y, numnp, t, ij)

    Dimension hnew(numnp), x(numnp), y(numnp)

    Write (75, 110) t
    l1 = (ij-1)/10 + 1
    Do n = 1, numnp, ij
      Do l = 1, l1
        m = n + (l-1)*10
        k = m + 9
        If (l==l1) k = n + ij - 1
        Write (75, 120) m, x(m), y(m), (hnew(j), j=m, k)
      End Do
    End Do
    Return

    110 Format (//' Time  ***', F12.4, ' ***'//'    n    x(n)   z(n)       h(n)      h(n+1) ...'/)
    120 Format (I5, 2F8.1, 10F10.1)
  End Subroutine hout

  !**********************************************************************

  Subroutine qout(q, x, y, numnp, t, ij)

    Dimension q(numnp), x(numnp), y(numnp)

    Write (73, 110) t
    l1 = (ij-1)/10 + 1
    Do n = 1, numnp, ij
      Do l = 1, l1
        m = n + (l-1)*10
        k = m + 9
        If (l==l1) k = n + ij - 1
        Write (73, 120) m, x(m), y(m), (q(j), j=m, k)
      End Do
    End Do
    Return

    110 Format (//' Time  ***', F12.4, ' ***'//'    n    x(n)   z(n)       Q(n)      Q(n+1) ...'/)
    120 Format (I5, 2F8.1, 10E11.3)
  End Subroutine qout

  !***********************************************************************

  Subroutine thout(theta, x, y, numnp, t, ij)

    Dimension theta(numnp), x(numnp), y(numnp)

    Write (76, 110) t
    l1 = (ij-1)/16 + 1
    Do n = 1, numnp, ij
      Do l = 1, l1
        m = n + (l-1)*16
        k = m + 15
        If (l==l1) k = n + ij - 1
        Write (76, 120) m, x(m), y(m), (theta(j), j=m, k)
      End Do
    End Do
    Return

    110 Format (//' Time  ***', F12.4, ' ***'//'    n    x(n)   z(n)      th(n)     th(n+1) ...'/)
    120 Format (I5, 2F8.1, 16F6.3)
  End Subroutine thout

  !***********************************************************************

  Subroutine flxout(vx, vz, x, y, numnp, t, ij)

    Dimension x(numnp), y(numnp), vx(numnp), vz(numnp)

    Write (81, 110) t
    Write (82, 120) t
    l1 = (ij-1)/10 + 1
    Do n = 1, numnp, ij
      Do l = 1, l1
        m = n + (l-1)*10
        k = m + 9
        If (l==l1) k = n + ij - 1
        Write (81, 130) m, x(m), y(m), (vz(j), j=m, k)
        Write (82, 130) m, x(m), y(m), (vx(j), j=m, k)
      End Do
    End Do
    Return

    110 Format (//' Time  ***', F12.4, ' ***'//'    n    x(n)  z(n)     vz(n)     vz(n+1) ...'/)
    120 Format (//' Time  ***', F12.4, ' ***'//'    n    x(n)  z(n)     vx(n)     vx(n+1) ...'/)
    130 Format (I5, 2F8.1, 10E10.2)
  End Subroutine flxout

  !***********************************************************************

  Subroutine cout(numnp, conc, x, y, t, ij)

    Dimension conc(numnp), x(numnp), y(numnp)

    Write (83, 110) t
    l1 = (ij-1)/10 + 1
    Do n = 1, numnp, ij
      Do l = 1, l1
        m = n + (l-1)*10
        k = m + 9
        If (l==l1) k = n + ij - 1
        Write (83, 120) m, x(m), y(m), (conc(j), j=m, k)
      End Do
    End Do
    Return

    110 Format (//' Time  ***', F12.4, ' ***'//'    n    x(n)   z(n)      Conc(n)   Conc(n+1)  ...'/)
    120 Format (I5, 2F8.1, 10E11.3)
  End Subroutine cout

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
