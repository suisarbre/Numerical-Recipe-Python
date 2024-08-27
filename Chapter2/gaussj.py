"""
    SUBROUTINE GAUSSJ (A,N,NP,B,M,MP)
        Linier equation solution by Gauss-Jordan elimination, equation (2.1.1) above.
        A is an input matrix of N by N elements, stored in an array of physical dimensions NP by NP.
        B is an input matrix of N by M containing the M right-hand side vectors, stored in an array of physical dimensions NP by MP.
        On output, A is replaced by its matrix inverse, and B is replaced by the corresponding set of solution vectors.
        
    PARAMETER (NMAX = 50)
    DIMENSION A(NP,NP), B(NP,MP), IPIV(NMAX), INDXR(NMAX), INDXC(NMAX)
        The integer arrays IPIV, INDXR, and INDXC are used for bookkeeping on the pivoting.
        NMAX should be as large as the largest anticipated value of N.
    DO 11 J = 1, N
        IPIV(J) = 0
    11 CONTINUE
    DO 22 I = 1, N  #This is the main loop over the columns to be reduced
        BIG = 0
        DO 13 J = 1, N  #This is the outer loop of the search for a pivot element
            IF(IPIV(J) .NE. 1) THEN
                DO 12 K = 1, N
                    IF (IPIV(K) .EQ. 0) THEN
                        IF (ABS(A(J,K)) .GE. BIG) THEN
                            BIG = ABS(A(J,K))
                            IROW = J
                            ICOL = K
                        ENDIF
                    ELSE IF (IPIV(K) .GT. 1) THEN
                        PAUSE 'Singular matrix'
                    ENDIF
                12 CONTINUE
            ENDIF
        13 CONTINUE
        IPIV(ICOL) = IPIV(ICOL) + 1
            We now have the pivot element, so we interchange rows, if needed, to put the pivot element on the diagonal.
            The columns are not physically interchanged, only relabeld: INDX(I), the column of the Ith pivot element, is the Ith colum that is reduced, while INDXR(I) is the row in which that pivot element was originally located.
            If  INDXR(I) != INDXC(I) there is an implied column interchange.
            With this form of bookkeeping, the solution B's will end up in the correct order, and the inverse matrix will be scrambled by columns.
        IF (IROW .NE. ICOL) THEN
            DO 14 L = 1, N
                DUM = A(IROW,L)
                A(IROW,L) = A(ICOL,L)
                A(ICOL,L) = DUM
            14 CONTINUE
            
            DO 15 L = 1, M
                DUM = B(IROW,L)
                B(IROW,L) = B(ICOL,L)
                B(ICOL,L) = DUM
            15 CONTINUE
        ENDIF
        INDXR(I) = IROW   #We are now ready to divide the pivot row by the pivot element, located at IROW and ICOL
        INDXC(I) = ICOL
        IF (A(ICOL,ICOL) .EQ. 0) PAUSE 'Singular matrix'
        PIVINV = 1 ./A(ICOL,ICOL)
        A(ICOL,ICOL) = 1
        
        DO 16 L = 1, N
            A(ICOL,L) = A(ICOL,L)*PIVINV
        16 CONTINUE
        
        DO 17 L = 1,N
            B(ICOL,L) = B(ICOL,L)*PIVINV
        17 CONTINUE
        
        DO 21 LL = 1, N #Next we reduce the rows...except for the pivot one, of course
            IF (LL .NE. ICOL) THEN
                DUM = A(LL,ICOL)
                A(LL,ICOL) = 0
                DO 18 L = 1, N
                    A(LL,L) = A(LL,L) -A(ICOL,L) * DUM
                18 CONTINUE
                
                DO 19 L = 1, M
                    B(LL,L) = B(LL,L) - B(ICOL,L) * DUM
                19 CONTINUE
            ENDIF
        21 CONTINUE
    22 CONTINUE 
        This is the end of the main loop over columns of the reduction.
        It only remains to unscramble the solution in view of the column interchanges.
        We do this by interchanging pairs of columns in the reverse order that the permutation was built up.
    
    DO 24 L = N, 1, -1
        IF (INDXR(L) .NE. INDXC(L)) THEN
            DO 23 K = 1, N
                DUM = A(K,INDXR(L))
                A(K,INDXR(L)) = A(K,INDXC(L))
                A(K,INDXC(L)) = DUM
            23 CONTINUE
        ENDIF
    24 CONTINUE
    RETURN
    END
"""

def gaussj(A,B):
    NMAX = 50
    IPIV = [0]*NMAX
    INDXR = [0]*NMAX
    INDXC = [0]*NMAX
    icol, irow = 0, 0
    
    for i in range(len(A)):
        big = 0
        for j in range(len(A)):
            if IPIV[j] != 1:
                for k in range(len(A)):
                    if IPIV[k] == 0:
                        if abs(A[j][k]) >= big:
                            big = abs(A[j][k])
                            irow = j
                            icol = k
                    elif IPIV[k] > 1:
                        print('Singular matrix')
        IPIV[icol] = IPIV[icol] + 1
        if irow != icol:
            for l in range(len(A)):
                dum = A[irow][l]
                A[irow][l] = A[icol][l]
                A[icol][l] = dum
            for l in range(len(B)):
                dum = B[irow][l]
                B[irow][l] = B[icol][l]
                B[icol][l] = dum
        
        INDXR[i] = irow
        INDXC[i] = icol
        if A[icol][icol] == 0:
            print('Singular matrix')
        PIVINV = 1/A[icol][icol]
        A[icol][icol] = 1
        
        for l in range(len(A)):
            A[icol][l] = A[icol][l]*PIVINV
            
        for l in range(len(A)):
            B[icol][l] = B[icol][l]*PIVINV
        
        for ll in range(len(A)):
            if ll != icol:
                dum = A[ll][icol]
                A[ll][icol] = 0
                for l in range(len(A)):
                    A[ll][l] = A[ll][l] - A[icol][l]*dum
                for l in range(len(B)):
                    B[ll][l] = B[ll][l] - B[icol][l]*dum
    for l in range(len(A)-1, -1, -1):
        if INDXR[l] != INDXC[l]:
            for k in range(len(A)):
                dum = A[k][INDXR[l]]
                A[k][INDXR[l]] = A[k][INDXC[l]]
                A[k][INDXC[l]] = dum
                
