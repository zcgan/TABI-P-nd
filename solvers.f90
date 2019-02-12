!####################################################################################################
SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err)
use bicg
Implicit Double Precision(A-H, O-Z)
INTEGER*4 iter,itmax,itol,n
DOUBLE PRECISION err,tol,b(n),x(n),EPS 
PARAMETER (EPS=1.d-20)

!C USES atimes,asolve,snrm
!Solves A ¡¤ x = b for x(1:n), given b(1:n), by the iterative biconjugate gradient method.
!On input x(1:n) should be set to an initial guess of the solution (or all zeros); itol is
!1,2,3, or 4, specifying which convergence test is applied (see text); itmax is the maximum
!number of allowed iterations; and tol is the desired convergence tolerance. On output,
!x(1:n) is reset to the improved solution, iter is the number of iterations actually taken,
!and err is the estimated error. The matrix A is referenced only through the user-supplied
!routines atimes, which computes the product of either A or its transpose on a vector; and
!asolve, which solves .A ¡¤ x = b or .AT
!¡¤ x = b for some preconditioner matrix .A (possibly
!the trivial diagonal part of A).
real*8,allocatable :: p(:),pp(:),r(:),rr(:),z(:),zz(:)
INTEGER*4 j
DOUBLE PRECISION	ak,akden,bk,bkden,bknum,bnrm,dxnrm, &
					xnrm,zm1nrm,znrm,snrm
allocate(p(n),pp(n),r(n),rr(n),z(n),zz(n))

iter=0								!Calculate initial residual.
call atimes(n,x,r,0)				!Input to atimes is x(1:n), output is r(1:n);
									!the final 0 indicates that the matrix (not
									!its transpose) is to be used.
do j=1,n
	r(j)=b(j)-r(j)
	rr(j)=r(j)
enddo

call atimes(n,r,rr,0)				!Uncomment this line to get the ¡°minimum
									!residual¡± variant of the algorithm. 
if (itol.eq.1) then
	bnrm=snrm(n,b,itol)
	call asolve(n,r,z,0)			!Input to asolve is r(1:n), output is z(1:n);
									!the final 0 indicates that the matrix .A
									!(not its transpose) is to be used.
elseif (itol.eq.2) then
	
	call asolve(n,b,z,0)
	bnrm=snrm(n,z,itol)
	call asolve(n,r,z,0)
elseif (itol.eq.3.or.itol.eq.4) then
	call asolve(n,b,z,0)
	bnrm=snrm(n,z,itol)
	call asolve(n,r,z,0)
	znrm=snrm(n,z,itol)
else
		pause 'illegal itol in linbcg'
endif

100 if (iter.le.itmax) then			!Main loop.
		iter=iter+1
		call asolve(n,rr,zz,1)		!Final 1 indicates use of transpose matrix .AT.
		bknum=0.d0
		do j=1,n					!Calculate coefficient bk and direction vectors p and pp. 
			bknum=bknum+z(j)*rr(j)
		enddo
		if (iter.eq.1) then
		do j=1,n
			p(j)=z(j)
			pp(j)=zz(j)
		enddo
	else
		bk=bknum/bkden
		do j=1,n
			p(j)=bk*p(j)+z(j)
			pp(j)=bk*pp(j)+zz(j)
		enddo
	endif
bkden=bknum							!Calculate coefficient ak, new iterate x, and
									!new residuals r and rr. 
call atimes(n,p,z,0)

akden=0.d0
do j=1,n
	akden=akden+z(j)*pp(j)
enddo
ak=bknum/akden

call atimes(n,pp,zz,1)

do j=1,n
	x(j)=x(j)+ak*p(j)
	r(j)=r(j)-ak*z(j)
	rr(j)=rr(j)-ak*zz(j)
enddo

call asolve(n,r,z,0)				!Solve .A¡¤z = r and check stopping criterion.


if (itol.eq.1) then
	err=snrm(n,r,itol)/bnrm
elseif(itol.eq.2)then
	err=snrm(n,z,itol)/bnrm
elseif(itol.eq.3.or.itol.eq.4)then
	zm1nrm=znrm
	znrm=snrm(n,z,itol)
	if(abs(zm1nrm-znrm).gt.EPS*znrm) then
		dxnrm=abs(ak)*snrm(n,p,itol)
		err=znrm/abs(zm1nrm-znrm)*dxnrm
	else
		err=znrm/bnrm				!Error may not be accurate, so loop again.
		goto 100
	endif
	xnrm=snrm(n,x,itol)
	if(err.le.0.5d0*xnrm) then
		err=err/xnrm
	else
		err=znrm/bnrm				!Error may not be accurate, so loop again.
	goto 100
	endif
endif

!print *,iter,err
if(err.gt.tol) goto 100
!write (*,*) 'iter=', iter, 'err=',err, tol
endif

deallocate(p,pp,r,rr,z,zz)
return
END

!----------------------------------------------------------------------------------------------------

!The routine linbcg uses this short utility for computing vector norms:
FUNCTION snrm(n,sx,itol)
Implicit Double Precision(A-H, O-Z)
INTEGER n,itol,i,isamax
DOUBLE PRECISION sx(n),snrm
!Compute one of two norms for a vector sx(1:n), as signaled by itol. Used by linbcg.

if (itol.le.3)then
	snrm=0.d0
	do i=1,n						!Vector magnitude norm.
		snrm=snrm+sx(i)**2
	enddo
	snrm=dsqrt(snrm)
else
	isamax=1
	do i=1,n						!Largest component norm.
		if(dabs(sx(i)).gt.dabs(sx(isamax))) isamax=i
	enddo
	snrm=dabs(sx(isamax))
endif
return
END


SUBROUTINE atimes(n,x,r,itrnsp)
use bicg
Implicit Double Precision(A-H, O-Z)
INTEGER*4 n,itrnsp
DOUBLE PRECISION x(n),r(n)
!PARAMETER (nmax=10000000)
!COMMON /mat/ sa(NMAX),ijka(NMAX)		!The matrix is stored somewhere.
!C USES dsprsax,dsprstx DOUBLE PRECISION versions of sprsax and sprstx.

if (itrnsp.eq.0) then
	call dsprsax(x,r,n)
else
	call dsprstx(x,r,n)
endif
return
END


SUBROUTINE asolve(n,b,x,itrnsp)
use bicg
Implicit Double Precision(A-H, O-Z)
INTEGER n,itrnsp,i
DOUBLE PRECISION x(n),b(n)
!PARAMETER (nmax=10000000)
!COMMON /mat/ sa(NMAX),ijka(NMAX)		!The matrix is stored somewhere.

do i=1,n
	x(i)=b(i)/sa(i)					!The matrix .A is the diagonal part of A, stored in
									!the first n elements of sa. Since the transpose
									!matrix has the same diagonal, the flag itrnsp is not used.

enddo
!print *,b
return
END


!------------------------------------------------------------------------------------
SUBROUTINE dsprsax(x,b,n)
use bicg
Implicit Double Precision(A-H, O-Z)
INTEGER n
REAL*8 b(n),x(n)
!Multiply a matrix in row-index sparse storage arrays sa and ijka by a vector x(1:n), giving
!a vector b(1:n).
INTEGER*4 i,k

if (ijka(1).ne.n+2) pause 'mismatched vector and matrix in sprsax'
do i=1,n
	b(i)=sa(i)*x(i)					!Start with diagonal term.
	do k=ijka(i),ijka(i+1)-1			!Loop over off-diagonal terms.
		b(i)=b(i)+sa(k)*x(ijka(k))
	enddo
enddo
return
END


SUBROUTINE dsprstx(x,b,n)
use bicg
Implicit Double Precision(A-H, O-Z)
INTEGER*4 n
REAL*8 b(n),x(n)
!Multiply the transpose of a matrix in row-index sparse storage arrays sa and ijka by a
!vector x(1:n), giving a vector b(1:n).
INTEGER*4 i,j,k

if (ijka(1).ne.n+2) pause 'mismatched vector and matrix in sprstx'
	do i=1,n						!Start with diagonal terms.
		b(i)=sa(i)*x(i)
	enddo
	do i=1,n						!Loop over off-diagonal terms.
		do k=ijka(i),ijka(i+1)-1
			j=ijka(k)
			b(j)=b(j)+sa(k)*x(i)
		enddo
	enddo
return
END



!--------------------------------------------------------------------------
SUBROUTINE dsprsin(a,n,np,thresh,nmax,sa,ijka)
Implicit Double Precision(A-H, O-Z)
INTEGER*4 n,nmax,np,ijka(nmax)
REAL*8 thresh,a(np,np),sa(nmax)
!Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
!storage mode. Only elements of a with magnitude ¡Ýthresh are retained. Output is in
!two linear arrays with physical dimension nmax (an input parameter): sa(1:) contains
!array values, indexed by ijka(1:). The logical sizes of sa and ijka on output are both
!ijka(ijka(1)-1)-1 (see text).
INTEGER i,j,k
do j=1,n								!Store diagonal elements.
	sa(j)=a(j,j)
enddo 
ijka(1)=n+2								!Index to 1st row off-diagonal element, if any.
k=n+1
do i=1,n								!Loop over rows.
	do j=1,n							!Loop over columns.
		if(dabs(a(i,j)).ge.thresh)then
			if(i.ne.j)then				!Store off-diagonal elements and their columns.
				k=k+1
				if(k.gt.nmax) pause 'nmax too small in sprsin'
				sa(k)=a(i,j)
				ijka(k)=j
			endif
		endif
	enddo
	ijka(i+1)=k+1						!As each row is completed, store index to next.
enddo
return
END

