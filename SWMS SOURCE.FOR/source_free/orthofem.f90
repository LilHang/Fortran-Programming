! Source file ORTHOFEM.FOR |||||||||||||||||||||||||||||||||||||||||||||
!
!                            ORTHOFEM
!
!                           VERSION 1.02
!
!                      FORTRAN SUBROUTINES FOR
!                   ORTHOMIN OR CONJUGATE GRADIENT
!               MATRIX SOLUTION ON FINITE-ELEMENT GRIDS
!
!     ***********************************************************
!
!                           CARL A. MENDOZA
!
!                      WITH CONTRIBUTIONS FROM:
!                            RENE THERRIEN
!
!                     BASED ON AN ORIGINAL CODE BY:
!                         FRANK W. LETNIOWSKI
!
!              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
!                       UNIVERSITY OF WATERLOO
!                         WATERLOO, ONTARIO
!                          CANADA, N2L 3G1
!
!                    LATEST UPDATE: JANUARY 1991
!
!
!     ***********************************************************
!
!                   COPYRIGHT (c) 1989, 1990, 1991
!                            E.A. SUDICKY
!                            C.A. MENDOZA
!              WATERLOO CENTRE FOR GROUNDWATER RESEARCH
!
!          DUPLICATION OF THIS PROGRAM, OR ANY PART THEREOF,
!          WITHOUT THE EXPRESS PERMISSION OF THE COPYRIGHT
!          HOLDERS IS STRICTLY FORBIDDEN
!
!     ***********************************************************
!
!               Modified for SWMS_2D code by Jirka Simunek
!                            april 1993
!
!     ***********************************************************
!
!                             DISCLAIMER
!
!     ALTHOUGH GREAT CARE HAS BEEN TAKEN IN PREPARING THIS CODE
!     AND THE ACCOMPANYING DOCUMENTATION, THE AUTHORS CANNOT BE
!     HELD RESPONSIBLE FOR ANY ERRORS OR OMISSIONS.  THE USER IS
!     EXPECTED TO BE FAMILIAR WITH THE FINITE-ELEMENT METHOD,
!     PRECONDITIONED ITERATIVE TECHNIQUES AND FORTRAN PROGRAMING.
!
!     ***********************************************************
!
!             A USER'S GUIDE IS AVAILABLE -> CONSULT IT!
!
!     THESE SUBROUTINES SOLVE A BANDED (OR SPARSE) MATRIX USING:
!       - PRECONDITIONED ORTHOMIN FOR ASYMMETRIC MATRICES, OR
!       - PRECONDITIONED CONJUGATE GRADIENT FOR SYMMETRIC MATRICES
!         (FULL MATRIX STORAGE REQUIRED)
!
!     PRECONDITIONING IS BY INCOMPLETE LOWER-UPPER DECOMPOSITION
!       - ONLY ONE FACTORIZATION (GAUSSIAN ELIMINATION) IS PERFORMED
!       - EQUIVALENT TO DKR FACTORIZATION
!
!     THE SUBROUTINES ARE DESIGNED FOR FINITE-ELEMENT GRIDS
!       - ARBITRARY ELEMENT SHAPES AND NUMBERING MAY BE USED
!         - NUMBERING MAY, HOWEVER, AFFECT EFFICIENCY
!           - TRY TO MINIMIZE THE BANDWIDTH AS MUCH AS POSSIBLE
!         - ALL ELEMENTS MUST HAVE THE SAME NUMBER OF LOCAL NODES
!
!
!     THE FOLLOWING ROUTINES ARE CALLED FROM THE SOURCE PROGRAM:
!       IADMAKE (IN,NN,NE,NLN,MNLN,MAXNB,IAD,IADN,IADD)
!         -> ASSEMBLE ADJACENCY MATRIX
!       FIND (I,J,K,NN,MAXNB,IAD,IADN)
!         -> LOCATE MATRIX POSITION FOR A NODAL PAIR (ASSEMBLY)
!       ILU (R,NN,MAXNB,IAD,IADN,IADD,B)
!         -> DECOMPOSE GLOBAL MATRIX
!       ORTHOMIN (R,C,GT,NNR,MAXNB,MAXNN,IAD,IADN,IADD,B,VRV,
!                 RES,RQI,RQ,Q,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
!                 NORTH,MNORTH,MAXIT)
!         -> SOLVE DECOMPOSED MATRIX
!
!     THESE ROUTINES CALL OTHER ROUTINES (LOCATED DIRECTLY BELOW THE
!     APPROPRIATE PRIMARY ROUTINE IN THE CODE)
!
!     THE FOLLOWING ARRAYS MUST BE DEFINED IN THE SOURCE PROGRAM
!     (THESE ARRAYS ARE PASSED TO THE SOLVER SUBROUTINES):
!
!     IN(MNLN,MAXNE) - INCIDENCE MATRIX (ELEMENTAL NODE DEFINITION)
!
!     GT(MAXNN)      - RIGHT-HAND-SIDE VECTOR
!     C(MAXNN)       - SOLUTION VECTOR
!     R(MAXNB,MAXNN) - GLOBAL MATRIX TO BE SOLVED
!
!     ARRAY DIMENSIONING PARAMETERS
!
!     MAXNN  - MAXIMUM NUMBER OF NODES
!     MAXNE  - MAXIMUM NUMBER OF ELEMENTS
!     MNLN   - MAXIMUM NUMBER OF LOCAL NODES (IN AN ELEMENT)
!     MAXNB  - MAXIMUM NUMBER OF NODES ADJACENT TO A PARTICULAR NODE
!              (INCLUDING ITSELF).
!            - IE. THE MAXIMUM NUMBER OF INDEPENDENT NODES THAT A
!              PARTICULAR NODE SHARES AN ELEMENT WITH.
!            - THIS WILL BE IDENTICALLY EQUIVALENT TO THE MAXIMUM
!              NUMBER OF NONZERO ENTRIES IN A ROW OF THE FULL MATRIX.
!     MNORTH - MAXIMUM NUMBER OF ORTHOGONALIZATIONS PERFORMED
!              (AT LEAST MNORTH = 1 REQUIRED FOR CONJUGATE GRADIENT)
!
!
!     ORTHOMIN ARRAY SPACE/VARIABLES
!
!     NORTH  - NUMBER OF ORTHOGONALIZATIONS TO PERFORM
!            - SET NORTH=0 FOR SYMMETRIC MATRICES (CONJUGATE GRADIENT)
!     ECNVRG - RESIDUAL CONVERGENCE TOLERANCE
!     ACNVRG - ABSOLUTE CONVERGENCE TOLERANCE
!     RCNVRG - RELATIVE CONVERGENCE TOLERANCE
!     MAXIT  - MAXIMUM NUMBER OF ITERATIONS TO PERFORM
!     ITERP  - NUMBER OF ITERATIONS PERFORMED
!
!     B(MAXNB,MAXNN) - ILU DECOMPOSED MATRIX
!     Q(MAXNN)   - SEARCH DIRECTION Q
!     RQ(MAXNN)  - PRODUCT OF R AND Q
!     VRV(MAXNN) - EITHER V OR PRODUCT OF R AND V
!     RES(MAXNN) - RESIDUAL
!
!     QI(MAXNN,MNORTH)  - STORAGE OF Q'S
!     RQI(MAXNN,MNORTH) - STORAGE OF PRODUCTS OF R AND Q
!     RQIDOT(MNORTH)    - STORAGE OF DOT PRODUCTS OF RQ AND RQ
!
!     RESV   - PREVIOUS VALUE OF RES V DOT PRODUCT (CONJUGATE GRADIENT)
!
!     IAD(MAXNB,MAXNN) - ADJACENCY MATRIX (NODAL CONNECTIONS)
!
!     IADN(MAXNN) - NUMBER OF ADJACENT NODES IN IAD (SELF-INCLUSIVE)
!
!     IADD(MAXNN) - POSITION OF DIAGONAL IN ADJACENCY MATRIX
!
!
!     OTHER PARAMETERS PASSED FROM SOURCE PROGRAM
!
!     NN  - NUMBER OF NODES
!     NE  - NUMBER OF ELEMENTS
!     NLN - NUMBER OF LOCAL NODES IN AN ELEMENT
!
!
!     APPROXIMATE REAL STORAGE SPACE FOR ORTHOMIN AND MATRIX EQUATION
!
!       ((6 + 2*MAXNB + 2*MNORTH)*MAXNN)*(8 BYTES)
!
!***********************************************************************

Subroutine iadmake(kx, numnp, numel, numeld, maxnb, iad, iadn, iadd)

    !     Generate the adjacency matrix for nodes from the element
    !     indidence matrix

    !     Requires subroutine Insert

      Dimension kx(numeld, 4), iad(maxnb, numnp), iadn(numnp), iadd(numnp)
      Integer e

    !     Determine independent adjacency within each element
    !     version for SWMS_2D

      Do i = 1, numnp  ! 初始化iad，iadn，iadd数组
        iadn(i) = 0
        iadd(i) = 0
        Do j = 1, maxnb
          iad(j, i) = 0
        End Do
      End Do

      Do e = 1, numel   ! 逐个扫描，建立iad数组
        ncorn = 4
        If (kx(e,3)==kx(e,4)) ncorn = 3
        Do n = 1, ncorn - 2
          i = kx(e, 1)
          j = kx(e, n+1)
          k = kx(e, n+2)

          Call insert(i, j, kk, numnp, maxnb, iad, iadn) ! 将节点j插入iad的第列，下同
          Call insert(j, i, kk, numnp, maxnb, iad, iadn)
          Call insert(i, k, kk, numnp, maxnb, iad, iadn)
          Call insert(k, i, kk, numnp, maxnb, iad, iadn)
          Call insert(j, k, kk, numnp, maxnb, iad, iadn)
          Call insert(k, j, kk, numnp, maxnb, iad, iadn)
        End Do
      End Do  ! 已生成iad，iadn数组

    !     Determine self-adjacency terms

      Do i = 1, numnp
        Call insert(i, i, kk, numnp, maxnb, iad, iadn) ! 在此调用insert，传入i， i，kk，得到i的位置
    !     Store self-adjacency position
        iadd(i) = kk
      End Do ! 已生成iadd数组
      Return
    End Subroutine iadmake

    !***********************************************************************

    Subroutine insert(i, j, k, nn, maxnb, iad, iadn)

    !     ADD J TO THE ADJACENCY LIST FOR I

    !     RETURNS THE POSITION K WHERE IT HAS BEEN ADDED, OR WHERE IT
    !     WAS ALREADY IN THE LIST.

      Dimension iad(maxnb, nn), iadn(nn)

      Logical found
      found = .False.

    !     DETERMINE NUMBER OF NODES ALREADY IN ADJACENCY LIST

      n = iadn(i)
      k = n + 1

    !     DETERMINE WHETHER ALREADY IN LIST

      Do l = 1, n  ! 判断当前iad中，节点i相邻节点数组是否包含j
        inode = iad(l, i)
        If (inode>=j) Then
          k = l
          If (inode==j) found = .True.
          Goto 15
        End If
      End Do

      15 Continue

    !     PLACE IN LIST (NUMERICAL ORDER)

      If (found) Then ! 若j在数组，则跳过
        Continue
      Else            ! 否则，增加iad行数插入节点，若大于最大带宽则输出错误信息并终止程序
        If ((n+1)>maxnb) Then
          Write (*, 601) i, maxnb
          Write (50, 601) i, maxnb
          Stop
        End If

        iadn(i) = n + 1
        Do l = (n+1), (k+1), (-1)
          iad(l, i) = iad(l-1, i)
        End Do
        iad(k, i) = j
      End If

      Return
      601 Format (//5X, 'ERROR IN IADMAKE: NODE ', I5, ' HAS > ', I5, ' ADJACENCIES')
    End Subroutine insert

    !***********************************************************************

    Subroutine find(i, j, k, nn, maxnb, iad, iadn)

    !     FOR NODE I, DETERMINE THE 'BAND' (K) RELATED TO ITS ADJACENCY TO
    !     NODE J.

    !     IF NODE NOT ADJACENT, RETURN 0 AS THE 'BAND'

      Dimension iad(maxnb, nn), iadn(nn)

      k = 0
      n = iadn(i)

      Do l = 1, n
        inode = iad(l, i)

    !     EXIT THE LOOP IF AT OR PAST THE REQUIRED POSITION

        If (inode>=j) Then
          If (inode==j) k = l ! j存在则返回位置退出循环，否则j不存在则返回初始值0
          Goto 20
        End If
      End Do

      20 Continue

      Return
    End Subroutine find

    !***********************************************************************

    Subroutine ilu(r, nn, maxnb, iad, iadn, iadd, b)

    !     INCOMPLETE LOWER-UPPER DECOMPOSITION OF MATRIX R INTO B
    !     ONE STEP OF GAUSSIAN ELIMINATION PERFORMED
    !     DIAGONAL DOMINANCE IS ASSUMED - NO PIVOTING PERFORMED
    !     REQUIRES FUNCTION DU

      Implicit Real *8(A-H, O-Z)
      Dimension r(maxnb, nn), iad(maxnb, nn), iadn(nn), iadd(nn), b(maxnb, nn)

    !     INITIALIZE B

      Do i = 1, nn
        Do j = 1, maxnb
          b(j, i) = 0.0D0
        End Do
      End Do

    !     LOOP OVER NODES

      Do i = 1, nn

    !     DETERMINE NUMBER OF BANDS/POSITION OF DIAGONAL IN THIS ROW

        n = iadn(i)
        k = iadd(i)

    !     LOWER TRIANGULAR MATRIX

        Do j = 1, (k-1)
          sum = r(j, i)
          icur = iad(j, i)
          Do l = 1, (j-1)
            inode = iad(l, i)
            sum = sum - b(l, i)*du(inode, icur, nn, maxnb, iad, iadn, iadd, b)
          End Do
          b(j, i) = sum
        End Do

    !     DIAGONAL

        sum = r(k, i)
        Do l = 1, (k-1)
          inode = iad(l, i)
          sum = sum - b(l, i)*du(inode, i, nn, maxnb, iad, iadn, iadd, b)
        End Do
        d = 1.0D0/sum
        b(k, i) = d

    !     UPPER TRIANGULAR MATRIX
    !       - ACTUALLY D*U TO OBTAIN UNIT DIAGONAL

        Do j = (k+1), n
          sum = r(j, i)
          icur = iad(j, i)
          Do l = 1, (k-1)
            inode = iad(l, i)
            sum = sum - b(l, i)*du(inode, icur, nn, maxnb, iad, iadn, iadd, b)
          End Do
          b(j, i) = d*sum
        End Do

      End Do

      Return
    End Subroutine ilu

    !***********************************************************************

    Function du(i, inode, nn, maxnb, iad, iadn, iadd, b)

    !     SEARCHES THE I'TH ROW OF THE UPPER DIAGONAL MATRIX
    !     FOR AN ADJACENCY TO THE NODE 'INODE'

    !     RETURNS CORRESPONDING VALUE OF B (OR ZERO)

      Implicit Real *8(A-H, O-Z)
      Dimension iad(maxnb, nn), iadn(nn), iadd(nn), b(maxnb, nn)

      temp = 0.0D0
      n = iadn(i)
      k = iadd(i)

      If (i==inode) Then
        temp = 1.0D0
        Goto 20
      End If

      Do j = (k+1), n
        If (inode==iad(j,i)) Then
          temp = b(j, i)
          Goto 20
        End If
      End Do

      20 Continue
      du = temp

      Return
    End Function du

    !***********************************************************************

    Subroutine orthomin(r, c, gt, nn, maxnb, maxnn, iad, iadn, iadd, b, vrv, res, rqi, rq, q, &
                       qi, rqidot, ecnvrg, rcnvrg, acnvrg, north, mnorth, maxit)

    !     ORTHOMIN OR CONJUGATE GRADIENT ACCELERATION/SOLUTION

    !     CONJUGATE GRADIENT (SYMMETRIC MATRIX) IF NORTH=0
    !     (HOWEVER, NOTE THAT MNORTH MUST BE AT LEAST 1)

    !     REQUIRES FUNCTIONS SDOT,SDOTK,SNRM
    !     REQUIRES SUBROUTINES LUSOLV,MATM2,SAXPYK,SCOPY,SCOPYK

      Implicit Real *8(A-H, O-Z)
      Dimension r(maxnb, nn), c(nn), gt(nn), iad(maxnb, nn), iadn(nn), iadd(nn), b(maxnb, nn), &
                vrv(nn), res(nn), rqi(maxnn, mnorth), rq(nn), q(nn), rqidot(mnorth), qi(maxnn, mnorth)

    !     INITIALIZE RESIDUAL VECTOR

      Call matm2(res, r, c, nn, iad, iadn, maxnb)

    !     Solution for homogeneous system of equations - Modified by Simunek
      ihomog = 0
      Do i = 1, nn
        res(i) = gt(i) - res(i)
        If (abs(gt(i))>1.D-300) ihomog = 1
      End Do
      If (ihomog==0) Then
        Do i = 1, nn
          c(i) = gt(i)
        End Do
        Return
      End If

    !     LOOP OVER A MAXIMUM OF MAXIT ITERATIONS

      norcur = 0

      Do iter = 1, maxit

    !     INVERT LOWER/UPPER MATRICES

        Call scopy(nn, res, vrv)

        Call lusolv(nn, maxnb, iad, iadn, iadd, b, vrv)

    !     COPY V INTO Q

        Call scopy(nn, vrv, q)

    !     CALCULATE PRODUCT OF R AND V

        Call matm2(vrv, r, q, nn, iad, iadn, maxnb)

    !     COPY RV INTO RQ

        Call scopy(nn, vrv, rq)

    !     RES V DOT PRODUCT (CONJUGATE GRADIENT)

        If (north==0) Then

          dot = sdot(nn, res, q)
          If (norcur==0) resv = dot

        End If

    !     LOOP OVER PREVIOUS ORTHOGONALIZATIONS

        k = 1

        20 If (k>norcur) Goto 30

    !     DETERMINE WEIGHTING FACTOR (CONJUGATE GRADIENT)

        If (north==0) Then

          alpha = dot/resv
          resv = dot

    !     DETERMINE WEIGHTING FACTOR (ORTHOMIN)

        Else

          dot = sdotk(nn, k, rqi, vrv, maxnn, mnorth)
          alpha = -dot/rqidot(k)

        End If

    !     SUM TO OBTAIN NEW Q AND RQ

        Call saxpyk(nn, alpha, k, qi, q, maxnn, mnorth)
        Call saxpyk(nn, alpha, k, rqi, rq, maxnn, mnorth)

        k = k + 1
        Goto 20
        30 Continue

    !     CALCULATE WEIGHTING FACTOR (CONJUGATE GRADIENT)

        If (north==0) Then

          dot = sdot(nn, q, rq)
          omega = resv/dot

    !     CALCULATE WEIGHTING FACTOR (ORTHOMIN)

        Else

          dot = sdot(nn, res, rq)
          rqnorm = sdot(nn, rq, rq)
          omega = dot/rqnorm

        End If

    !     SAVE VALUES FOR FUTURE ORTHOGONALIZATIONS

        norcur = norcur + 1
        If (norcur>north) norcur = 1

        Call scopyk(nn, norcur, q, qi, maxnn, mnorth)
        Call scopyk(nn, norcur, rq, rqi, maxnn, mnorth)
        rqidot(norcur) = rqnorm

    !     UPDATE SOLUTION/RESIDUAL VECTORS

        Call saxpyk(nn, omega, norcur, qi, c, maxnn, mnorth)
        Call saxpyk(nn, -omega, norcur, rqi, res, maxnn, mnorth)

    !     DETERMINE CONVERGENCE PARAMETERS

        resmax = snrm(nn, res)
        dxnorm = dabs(omega)*snrm(nn, q)
        cnorm = snrm(nn, c)
        xratio = dxnorm/cnorm

    !     ITERATION (DEBUG) OUTPUT

    !     STOP ITERATING IF CONVERGED

        If (resmax<ecnvrg) Goto 200
        If (xratio<rcnvrg) Goto 200
        If (dxnorm<acnvrg) Goto 200

      End Do

    !     TERMINATE IF TOO MANY ITERATIONS

      Write (*, 602) maxit, resmax, xratio, dxnorm
      Write (70, 602) maxit, resmax, xratio, dxnorm

      Stop

      200 Continue

    !     RETURN NUMBER OF ITERATIONS REQUIRED

      iterp = iter

      Return

      602 Format (///5X, 'ORTHOMIN TERMINATES -- TOO MANY ITERATIONS', /8X, 'MAXIT  = ', I5, /8X, 'RESMAX = ', &
                  E12.4, /8X, 'XRATIO = ', E12.4, /8X, 'DXNORM = ', E12.4)
    End Subroutine orthomin

    !***********************************************************************

    Subroutine lusolv(nn, maxnb, iad, iadn, iadd, b, vrv)

    !     LOWER DIAGONAL MATRIX INVERSION BY FORWARD SUBSTITUTION
    !     UPPER DIAGONAL MATRIX INVERSION BY BACKWARD SUBSTITUTION
    !     LOWER/UPPER MATRICES ARE IN B
    !     RIGHT-HAND-SIDE VECTOR IS IN VRV AT START
    !     'SOLUTION' IS RETURNED IN VRV UPON EXIT

      Implicit Real *8(A-H, O-Z)
      Dimension iad(maxnb, nn), iadn(nn), iadd(nn), b(maxnb, nn), vrv(nn)

    !     LOWER INVERSION

      Do i = 1, nn
        sum = vrv(i)
        k = iadd(i)
        Do j = 1, (k-1)
          inode = iad(j, i)
          sum = sum - b(j, i)*vrv(inode)
        End Do

        vrv(i) = b(k, i)*sum

      End Do

    !     UPPER INVERSION

      Do i = nn, 1, -1
        sum = vrv(i)
        n = iadn(i)
        k = iadd(i)
        Do j = (k+1), n
          inode = iad(j, i)
          sum = sum - b(j, i)*vrv(inode)
        End Do

        vrv(i) = sum

      End Do

      Return
    End Subroutine lusolv

    !***********************************************************************

    Subroutine matm2(s1, r, p, nn, iad, iadn, maxnb)

    !     MULTIPLY MATRIX R BY VECTOR P TO OBTAIN S1

      Implicit Real *8(A-H, O-Z)
      Dimension s1(nn), p(nn), r(maxnb, nn), iad(maxnb, nn), iadn(nn)

      Do i = 1, nn
        sum = 0.0D0
        n = iadn(i)
        Do j = 1, n
          inode = iad(j, i)
          sum = sum + r(j, i)*p(inode)
        End Do

        s1(i) = sum
      End Do

      Return
    End Subroutine matm2

    !***********************************************************************

    Function sdot(nn, r, b)

    !     OBTAIN DOT PRODUCT OF R AND B

      Implicit Real *8(A-H, O-Z)
      Dimension r(nn), b(nn)

      sdot = 0.0D0
      Do l = 1, nn
        sdot = sdot + r(l)*b(l)
      End Do

      Return
    End Function sdot

    !***********************************************************************

    Function sdotk(nn, k, r, b, maxnn, mnorth)

    !     OBTAIN DOT PRODUCT OF R AND B

      Implicit Real *8(A-H, O-Z)
      Dimension r(maxnn, mnorth), b(nn)

      sdotk = 0.0D0
      Do l = 1, nn
        sdotk = sdotk + r(l, k)*b(l)
      End Do

      Return
    End Function sdotk

    !***********************************************************************

    Function snrm(nn, r)

    !     COMPUTE MAXIMUM NORM OF R

      Implicit Real *8(A-H, O-Z)
      Dimension r(nn)

      snrm = 0.0D0
      Do l = 1, nn
        temp = dabs(r(l))
        If (temp>snrm) snrm = temp
      End Do

      Return
    End Function snrm

    !***********************************************************************

    Subroutine saxpyk(nn, sa, k, fx, fy, maxnn, mnorth)

    !     MULTIPLY VECTOR FX BY SCALAR SA AND ADD TO VECTOR FY

      Implicit Real *8(A-H, O-Z)
      Dimension fx(maxnn, mnorth), fy(nn)

      If (nn>0) Then
        Do i = 1, nn
          fy(i) = sa*fx(i, k) + fy(i)
        End Do
      End If

      Return
    End Subroutine saxpyk

    !***********************************************************************

    Subroutine scopy(nn, fx, fy)

    !     COPY VECTOR FX INTO VECTOR FY

      Implicit Real *8(A-H, O-Z)
      Dimension fx(nn), fy(nn)

      If (nn>0) Then
        Do i = 1, nn
          fy(i) = fx(i)
        End Do
      End If

      Return
    End Subroutine scopy

    !***********************************************************************

    Subroutine scopyk(n, k, fx, fy, maxnn, mnorth)

    !     COPY VECTOR FX INTO VECTOR FY

      Implicit Real *8(A-H, O-Z)
      Dimension fx(n), fy(maxnn, mnorth)

      If (n>0) Then
        Do i = 1, n
          fy(i, k) = fx(i)
        End Do
      End If

      Return
    End Subroutine scopyk

    ! ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
