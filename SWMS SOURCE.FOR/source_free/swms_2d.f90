!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
!                                                                      *
!     SWMS_2D  - Numerical model of two-dimensional flow and solute    *
!                transport in a variably saturated porous medium       *
!                Conjugate gradient solver for symmetric matrix        *
!                ORTHOMIN solver for asymmetric matrix                 *
!                version 1.22                                          *
!                                                                      *
!     Updated by J.Simunek (1994)                                      *
!     Based on model SWMS_2D (Simunek et al., 1992)                    *
!                                                                      *
!                                         Last modified: January, 1996 *
!                                                                      *
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

Program swms_2d

    Parameter (numnpd=5000, numeld=6000, numbpd=250, mbandd=20, nseepd=2, numspd=50, ndrd=2, &
    neldrd=8, nmatd=20, ntabd=100, numkd=6,nobsd=4, mnorth=4)

    Double Precision a, b
  !      double precision RTime1,RTime2
    Double Precision a1, b1, vrv, res, rqi, rq, qq, qi, rqidot, ecnvrg, rcnvrg, acnvrg
    Logical lwat, lchem, sinkf, qgwlf, atminf, shortf, seepf, checkf, fluxf, explic, &
             lupw, freed, drainf, lartd, lort, lminstep
    Integer plevel, alevel, tlevel

    Dimension a(mbandd, numnpd), b(numnpd), kode(numnpd), q(numnpd), hnew(numnpd), &
             htemp(numnpd), hold(numnpd), consat(nmatd), f(numnpd), htab(ntabd), &
             contab(ntabd, nmatd), captab(ntabd, nmatd), con(numnpd), cap(numnpd), &
             x(numnpd), y(numnpd), matnum(numnpd), laynum(numeld), kx(numeld, 4), &
             kxb(numbpd), tprint(50), par(10, nmatd), width(numbpd), conaxx(numeld), &
             conazz(numeld), conaxz(numeld), swidth(numkd), np(nseepd, numspd), nsp(nseepd), &
             hsat(nmatd), watin(numeld), axz(numnpd), bxz(numnpd), dxz(numnpd), thr(nmatd), &
             thsat(nmatd), thetab(ntabd, nmatd), thnew(numnpd), thold(numnpd), listne(numnpd), &
             sink(numnpd), poptm(nmatd), beta(numnpd), ds(numnpd), cumq(numkd), vmean(numkd), &
             hmean(numkd), kodcb(numbpd), qc(numnpd), vx(numnpd), vz(numnpd), chpar(10, nmatd), &
             dispxx(numnpd), dispzz(numnpd), dispxz(numnpd), cbound(6), ac(numnpd), fc(numnpd), &
             solin(numeld), conc(numnpd), smean(numkd), chems(numkd), wetab(3, 2*numeld), &
             gc(numnpd), nd(ndrd), ned(ndrd), efdim(2, ndrd), keldr(ndrd, neldrd), cono(numnpd), &
             node(nobsd), b1(numnpd), iad(mbandd, numnpd), iadn(numnpd), iadd(numnpd), a1(mbandd, numnpd), &
             res(numnpd), vrv(numnpd), rqi(numnpd, mnorth), rq(numnpd), qq(numnpd), rqidot(mnorth), qi(numnpd, mnorth)

    ! 打开输入文件，创建输出文件
    Open (30, File='SWMS_2D.IN\Selector.in', Status='old')
    Open (32, File='SWMS_2D.IN\Grid.in', Status='old')
    Open (50, File='SWMS_2D.OUT\Check.out', Status='unknown')
    Open (71, File='SWMS_2D.OUT\v_Mean.out', Status='unknown')
    Open (72, File='SWMS_2D.OUT\A_Level.out', Status='unknown')
    Open (75, File='SWMS_2D.OUT\h.out', Status='unknown')
    Open (76, File='SWMS_2D.OUT\th.out', Status='unknown')
    Open (77, File='SWMS_2D.OUT\h_Mean.out', Status='unknown')
    Open (78, File='SWMS_2D.OUT\Cum_Q.out', Status='unknown')
    Open (79, File='SWMS_2D.OUT\Boundary.out', Status='unknown')
    Open (80, File='SWMS_2D.OUT\Balance.out', Status='unknown')
    Open (81, File='SWMS_2D.OUT\vz.out', Status='unknown')
    Open (82, File='SWMS_2D.OUT\vx.out', Status='unknown')
    Open (92, File='SWMS_2D.OUT\ObsNod.out', Status='unknown')

    Data sinkf, qgwlf, tinit, ntab, itcum, iter, tlevel, alevel, plevel/.False., .False., 0., 100, 0, 0, 1, 1, 1/&
    cumq, sink, cumqrt, cumqrr, cumqvr, chems, rroot, rtop/numkd*0., numnpd*0., 0., 0., 0., numkd*0., 0., 0./&
    cumch0, cumch1, cumchr, dtmaxc, wcuma, ccuma, explic, lminstep/0., 0., 0., 1.E+30, 0., 0., .False., .False./

    Data ecnvrg, acnvrg, rcnvrg, maxito/1.0D-8, 1.0D-8, 1.0D-8, 200/

  ! --- Reading of the input files and initial calculations --------------

    ! 读入basic information赋值给相关变量；创建输出文件并在里写入相关标题；在终端输出标题
    Call basinf(kat, maxit, tolth, tolh, lwat, lchem, atminf, shortf, seepf, checkf, fluxf, freed, drainf)
    ! 读入节点信息 根据checkf写入check.out
    Call nodinf(numnp, numel, ij, numbp, numnpd, numeld, numbpd, numkd, nobs, nobsd, kode, q, conc, hnew, &
               hold, htemp, x, y, matnum, beta, axz, bxz, dxz, checkf)
    ! 读入单元信息， 判断节点周围单元个数，计算带宽
    Call elemin(numel, numeld, numnp, kx, laynum, conaxx, conazz, conaxz, checkf, listne, ij, mband, mbandd, lchem, lort)
    ! 读入几何信息，观测点信息
    Call geomin(numkd, numnp, numbp, nobs, nobsd, swidth, width, kode, kxb, rlen, node)
    ! 得到iad，iadn，iadd数组
    Call iadmake(kx, numnp, numel, numeld, mbandd, iad, iadn, iadd)
    ! 关闭网格文件
    Close (32)
    ! 读入materia information，包括土壤种类数、子区数，土壤水分特征曲线参数等，子程序内调用material2模块土壤水力参数函数，计算每类土壤10个有效饱和度（qe）
    ! 下的土壤水力参数并写入check.out
    Call matin(nmatd, nmat, nlay, par, htab(1), htab(ntab))
    ! 生成相关土壤水力参数数组
    Call genmat(ntab, ntabd, nmat, thr, hsat, par, htab, contab, captab, consat, thetab, thsat)
    ! 计算节点theta k  C，得到con  cap  thold数组
    Call setmat(numnp, ntab, ntabd, nmat, htab, contab, captab, hnew, hold, matnum, par, con, cap, consat, axz, bxz, dxz, hsat, &
               htemp, explic, thetab, thsat, thr, thold)

    If (atminf) Then
      Open (31, File='SWMS_2D.IN\Atmosph.in', Status='old')
      Call atmin(gwl0l, sinkf, qgwlf, tinit, tmax, aqh, bqh, hcrits, maxal)
      ! update time-dependent boundary conditon
      Call setatm(tatm, rtop, rroot, hcrita, width, kxb, numbp, kode, hnew, q, numnp, gwl0l, qgwlf, &
                 freed, cprec, cht, crt, lminstep)
    End If
    Call tmin(tinit, tmax, tatm, told, dt, dtmax, dmul, dmul2, dtmin, tprint, t, dtopt, dtold, atminf)
    dtinit = dt
    If (sinkf) Then
      Call sinkin(nmat, numel, numnp, numeld, kat, kx, x, y, p0, poptm, p2h, p2l, p3, r2h, r2l, beta)
      Call setsnk(numnp, nmat, matnum, hnew, rroot, sink, p0, poptm, p2h, p2l, p3, r2h, r2l, beta, rlen)
    End If
    If (seepf) Call seepin(nseepd, numspd, nseep, nsp, np)
    ! 读入drain信息，计算调整的K
    If (drainf) Call drainin(ndr, ndrd, neldrd, numel, nd, ned, keldr, efdim, conaxx, conaxz, conazz)
    If (lchem) Then
      Call chemin(nmat, numbp, cbound, chpar, epsi, tpulse, kodcb, nlevel, lupw, lartd, pecr)
      If (lwat) Call chinit(numnp, numel, numeld, nmat, x, y, kx, matnum, nlevel, con, hnew, sink, cbound(5),&
                           vx, vz, conaxx, conazz, conaxz, dispxx, dispzz, dispxz, chpar, thold, thsat, conc,&
                           fc, gc, listne, lupw, wetab, dt, dtmaxc, peclet, courant, kat, lartd, pecr, cono)
      Open (83, File='SWMS_2D.OUT\Conc.out', Status='unknown')
      Open (74, File='SWMS_2D.OUT\Solute.out', Status='unknown')
      Call cout(numnp, conc, x, y, tinit, ij)
    End If
    Close (30)
    Close (50)
    Open (73, File='SWMS_2D.OUT\Q.out', Status='unknown')
    Open (70, File='SWMS_2D.OUT\Run_Inf.out', Status='unknown')

    ! 在h.out写入节点初始水头
    Call hout(hnew, x, y, numnp, tinit, ij)
    ! 在th.out写入节点初始含水率
    Call thout(thold, x, y, numnp, tinit, ij)
    Call subreg(numel, numeld, numnp, nmat, hnew, thold, thold, x, y, matnum, laynum, kx, kat, tinit, dt, nlay, 0,&
                lwat, lchem, conc, chpar, wcuma, wcumt, ccuma, ccumt, wvoli, cvoli, watin, solin)
    If (nobs>0) Call obsnod(tinit, numnp, nobs, nobsd, node, hnew, thold, conc)

    Write (*, *) 'beginning of numerical solution'
  !      call getdat(i,i,iday)
  !      call gettim(ihours,mins,isecs,i)
  !      Rtime1=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs

  ! --- Beginning of time loop -------------------------------------------

    11 Continue

  !     Calculate water flow
    If (lwat .Or. tlevel==1) Call watflow(numnp, numel, numeld, ntab, ntabd, mband, mbandd, nmat, nseep, nseepd,&
                                        numspd, nsp, np, numbp, itcum, maxit, iter, kode, kat, t, dt, dtmin, &
                                        dtopt, dtold, told, hcrita, hcrits, tolth, tolh, rlen, width, rtop, &
                                        vmeanr, hmeanr, atminf, sinkf, seepf, qgwlf, freed, par, htab, contab, &
                                        captab, thetab, hnew, hold, htemp, thr, thsat, thnew, thold, matnum, &
                                        con, cap, consat, axz, bxz, dxz, hsat, a, b, q, f, x, y, kx, sink, ds, &
                                        beta, conaxx, conazz, conaxz, kxb, explic, gwl0l, aqh, bqh, lwat, tlevel, &
                                        lort, drainf, nd, ndr, ndrd, rroot, p0, poptm, p2h, p2l, p3, r2h, r2l, &
                                        cono, a1, b1, numnpd, iad, iadn, iadd, vrv, res, rqi, rq, qq, qi, rqidot, &
                                        ecnvrg, rcnvrg, acnvrg, mnorth, maxito)
    If (.Not. lwat .And. tlevel==1) Then ! 稳定流解法
      If (lchem) Then
        Call chinit(numnp, numel, numeld, nmat, x, y, kx, matnum, nlevel, con, hnew, sink, cbound(5), vx, vz, &
                   conaxx, conazz, conaxz, dispxx, dispzz, dispxz, chpar, thnew, thsat, conc, fc, gc, listne,&
                    lupw, wetab, dt, dtmaxc, peclet, courant, kat, lartd, pecr, cono)
      Else
        Call veloc(kat, numnp, numel, numeld, hnew, x, y, kx, listne, con, conaxx, conazz, conaxz, vx, vz)
      End If
      iter = 1
    End If

  !     Calculate solute transport
    If (lchem) Call solute(numnp, numel, numeld, mband, mbandd, nmat, t, kode, a, b, q, hnew, hold, f, x, y, kx, &
                          kat, dt, ds, sink, matnum, con, cono, conaxx, conazz, conaxz, vx, vz, dispxx, dispzz, &
                          dispxz, chpar, thnew, thold, thsat, ac, fc, gc, qc, conc, listne, cbound, tpulse, numbp, &
                          kodcb, kxb, nlevel, cprec, crt, cht, lwat, lupw, wetab, epsi, cumch0, cumch1, cumchr, &
                          dtmaxc, peclet, courant, lartd, pecr, lort, a1, b1, numnpd, iad, iadn, iadd, vrv, res, &
                          rqi, rq, qq, qi, rqidot, ecnvrg, rcnvrg, acnvrg, mnorth, maxito)

  !     T-Level information
    Call tlinf(numnp, numbp, kode, q, hnew, cumq, width, swidth, kxb, t, dt, tlevel, shortf, tprint(plevel), &
              iter, itcum, rtop, rroot, vmeanr, hmeant, hmeanr, hmeang, atminf, sinkf, cumqrt, cumqrr, cumqvr,&
               numkd, hmean, vmean, lwat, lchem, rlen, peclet, courant, wcumt, wcuma)
    If (lchem) Call solinf(numnp, kode, qc, t, dt, tlevel, shortf, tprint(plevel), numkd, smean, chems, cumch0,&
                         cumch1, cumchr, ccuma, ccumt, lwat)
    If (nobs>0) Call obsnod(t, numnp, nobs, nobsd, node, hnew, thnew, conc)

  !     P-Level information
    If (abs(tprint(plevel)-t)<0.001*dt .Or. (.Not. lwat .And. .Not. lchem)) Then
      If (lwat .Or. (.Not. lwat .And. plevel==1)) Then
        Call hout(hnew, x, y, numnp, t, ij)
        Call thout(thnew, x, y, numnp, t, ij)
        If (fluxf) Then
          If (.Not. lchem) Call veloc(kat, numnp, numel, numeld, hnew, x, y, kx, listne, con, conaxx, conazz, conaxz, vx, vz)
          Call flxout(vx, vz, x, y, numnp, t, ij)
          Call qout(q, x, y, numnp, t, ij)
        End If
      End If
      Call subreg(numel, numeld, numnp, nmat, hnew, thold, thnew, x, y, matnum, laynum, kx, kat, t, dt, nlay, plevel, &
                 lwat, lchem, conc, chpar, wcuma, wcumt, ccuma, ccumt, wvoli, cvoli, watin, solin)
      Call bouout(numnp, numbp, t, hnew, thnew, q, width, kxb, kode, x, y, conc)
      If (lchem) Call cout(numnp, conc, x, y, t, ij)
      plevel = plevel + 1
    End If

  !     A-level information
    If (abs(t-tatm)<=0.001*dt .And. atminf) Then
      If (lwat) Call alinf(t, cumq, hmeant, hmeanr, hmeang, alevel, cumqrt, cumqrr, cumqvr, numkd)
      If (alevel<maxal) Then
        Call setatm(tatm, rtop, rroot, hcrita, width, kxb, numbp, kode, hnew, q, numnp, gwl0l, qgwlf,&
                   freed, cprec, cht, crt, lminstep)
        alevel = alevel + 1
      End If
    End If

  !     Root extraction
  !      if(SinkF)
  !     !  call SetSnk(NumNP,NMat,MatNum,hNew,rRoot,Sink,P0,POptm,P2H,P2L,
  !     !              P3,r2H,r2L,Beta,rLen)

  !     Time governing
    If (abs(t-tmax)<=0.001*dt .Or. (.Not. lwat .And. .Not. lchem)) Then
  !        call getdat(i,i,iday)
  !        call gettim(ihours,mins,isecs,i)
  !        Rtime2=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
  !        write(70,*)
  !        write(70,*) 'Real time [sec]',Rtime2-RTime1
  !        write( *,*) 'Real time [sec]',Rtime2-RTime1
      Stop
    End If
    told = t
    dtold = dt
    Call tmcont(dt, dtmax, dtopt, dmul, dmul2, dtmin, iter, tprint(plevel), tatm, t, tmax, dtmaxc, lminstep, dtinit)
    tlevel = tlevel + 1
    t = t + dt
    Goto 11

  ! --- end of time loop -------------------------------------------------

  End Program swms_2d

  !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
