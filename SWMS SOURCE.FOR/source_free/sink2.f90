! Source file SINK2.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

Subroutine setsnk(numnp, nmat, matnum, hnew, tpot, sink, p0, poptm, p2h, p2l, p3, r2h, r2l, beta, length)

    Real length
    Dimension matnum(numnp), hnew(numnp), poptm(nmat), beta(numnp), sink(numnp)

    Do i = 1, numnp
      If (beta(i)>0.) Then
        m = matnum(i)
        alfa = falfa(tpot, hnew(i), p0, poptm(m), p2h, p2l, p3, r2h, r2l)
        sink(i) = alfa*beta(i)*length*tpot
      End If
    End Do
    Return
  End Subroutine setsnk

  !***********************************************************************

  Real Function falfa(tpot, h, p0, p1, p2h, p2l, p3, r2h, r2l)

    If (tpot<r2l) p2 = p2l
    If (tpot>r2h) p2 = p2h
    If ((tpot>=r2l) .And. (tpot<=r2h)) p2 = p2h + (r2h-tpot)/(r2h-r2l)*(p2l-p2h)
    falfa = 0.0
    If ((h>p3) .And. (h<p2)) falfa = (h-p3)/(p2-p3)
    If ((h>=p2) .And. (h<=p1)) falfa = 1.0
    If ((h>p1) .And. (h<p0)) falfa = (h-p0)/(p1-p0)
    Return
  End Function falfa

  ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
