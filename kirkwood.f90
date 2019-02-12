!****************************************************************************
!
!  PROGRAM: kirkwood
!
!  PURPOSE:  Calculating the Electrostatic Potential at any point for given 
!  sphere with multiple changes 
!
!****************************************************************************

!program kirkwood
!use Molecule
!implicit double precision(a-h,o-z)
!real*8 pxyz(3)
!
!! Set Parameters
!eps0=1.d0
!eps1=20.d0
!rds=1.d0
!kappa=3.d0
!nchr=1
!nt=60
!
!! Setup Charge positions and amount
!allocate(chrpos(3,nchr),atmchr(nchr))
!chrpos(:,1)=(/0.d0,0.d0,0.9d0/)
!atmchr=1.d0

!a2=1.5d0
!chrpos(:,1)=(/a2,0.d0,0.d0/)
!chrpos(:,2)=(/0.d0,0.d0,a2/)

!chrpos(:,2)=(/0.d0,0.5d0,0.d0/)
!chrpos(:,3)=(/0.d0,0.d0,1.d0/)
!chrpos(:,4)=(/0.5d0,0.5d0,0.5d0/)
!chrpos(:,5)=(/-0.5d0,0.5d0,0.5d0/)
!chrpos(:,6)=(/0.5d0,-0.5d0,0.5d0/)

!chrpos(:,1)=(/0.4d0,0.d0,0.d0/)
!chrpos(:,2)=(/0.d0,0.8d0,0.d0/)
!chrpos(:,3)=(/0.d0,0.d0,1.2d0/)
!chrpos(:,4)=(/0.d0,0.d0,-0.4d0/)
!chrpos(:,5)=(/-0.8d0,0.d0,0.d0/)
!chrpos(:,6)=(/0.d0,-1.2d0,0.d0/)

!chrpos(:,1)=(/0.2d0,0.2d0,0.2d0/)
!chrpos(:,2)=(/0.5d0,0.5d0,0.5d0/)
!chrpos(:,3)=(/0.8d0,0.8d0,0.8d0/)
!chrpos(:,4)=(/-0.2d0,0.2d0,-0.2d0/)
!chrpos(:,5)=(/0.5d0,-0.5d0,0.5d0/)
!chrpos(:,6)=(/-0.8d0,-0.8d0,-0.8d0/)

!call lagendini
!
!pxyz=chrpos(:,1)
!call kirk_potential(-1,pxyz,ptl,pmt,psv)
!print *,ptl/4.0/pi,pmt/4.d0/pi,psv/4.d0/pi
!	
!
!end program kirkwood

subroutine lagendini
use Molecule
implicit double precision(a-h,o-z)
integer ierrpos(3)
real*8 pxyz(3)

allocate(chrpos_sph(3,nchr))
chrpos_sph=0.d0
do i=1,nchr
	call ct2sph(chrpos(:,i),chrpos_sph(:,i))
enddo

!truncated at nt_th term
allocate(chgmnx(-nt:nt,0:nt,1:nchr))
chgmnx=0.d0
do k=1,nchr
	call MLPMN(nt,nt,chgmnx(:,:,k),cos(chrpos_sph(2,k)))
enddo

end

!------------------------------------------------------------------------------
subroutine kirk_potential(iside,pxyz,ptl,pmt,psv)
!Input: The position where the energy is required to be calculated in cartition coordinates

use Molecule

implicit double precision(a-h,o-z)
real*8 prpt(3),pxyz(3), test(3) !r, theta, and phi of the charges
real*8 Enm(2), phirik(2), Bnm(2), Cnm(2), Dnm(2), Enmtp(2), phi2(2), phi1(2)
real*8 posmnx(-nt:nt,0:nt)

call ct2sph(pxyz,prpt)
call MLPMN(nt,nt,posmnx,cos(prpt(2)))

r=prpt(1)
ptl=0.d0
pmt=0.d0
psv=0.d0

do n=0,nt
	
	! Coefficient for Poisson Equation
	Bcf_pe=-((eps1/eps0-1)*(n+1))/((n*eps0+(n+1)*eps1)*rds**(2*n+1))
	Ccf_pe=(2*n+1)/(n*eps0+(n+1)*eps1)

	! Coefficient for Poisson Boltzmann Equation
	if (n==0) then
		Bcf_pbe=1/rds*(1/((1+kappa*rds)*eps1)-1/eps0)
	else
		xx=(n+1)+((kappa*rds)**2*fkn(kappa*rds,n-1))/(2*n-1)/fkn(kappa*rds,n)
		fenzi=(n+1)-eps1/eps0*xx
		fenmu=rds**(2*n+1)*(n*eps0+eps1*xx)
		Bcf_pbe=fenzi/fenmu
	endif
	Ccf_pbe=exp(kappa*rds)*(1/eps0+rds**(2*n+1)*Bcf_pbe)/fkn(kappa*rds,n)
	
	!Bcf=Bcf_pe;		Ccf=Ccf_pe
	Bcf=Bcf_pbe;	Ccf=Ccf_pbe

	ptltp=0.d0
	pmttp=0.d0
	psvtp=0.d0
	do m=-n,n
		Enm=0.d0
		do k=1,nchr
			pnm0=chgmnx(m,n,k)
			phi2(1)=cos(-m*chrpos_sph(3,k))
			phi2(2)=sin(-m*chrpos_sph(3,k))
			Enmtp=atmchr(k)*chrpos_sph(1,k)**n*pnm0*phi2
			Enm=Enm+Enmtp
		enddo
		
		Enm=Enm*dble(fact(n-abs(m)))/dble(fact(n+abs(m)))
		Bnm=Bcf*Enm
		Cnm=Ccf*Enm
		Dnm=Enm/eps0/(rds**(2*n+1))

		pnm1=posmnx(m,n)
		phi1(1)=cos(m*prpt(3))
		phi1(2)=sin(m*prpt(3))
		
		aa=pnm1*phi1(1)
		bb=pnm1*phi1(2)
		
		if (iside<=0) then
			cc=Enm(1)/eps0/r**(n+1)+Bnm(1)*r**n
			dd=Enm(2)/eps0/r**(n+1)+Bnm(2)*r**n
			ee=Bnm(1)*r**n
			ff=Bnm(2)*r**n
			gg=Dnm(1)*r**n
			hh=Dnm(2)*r**n
		else
			cc=Cnm(1)*exp(-kappa*r)*fkn(kappa*r,n)/r**(n+1)
			dd=Cnm(2)*exp(-kappa*r)*fkn(kappa*r,n)/r**(n+1)
		endif
		ptltp1=aa*cc-bb*dd
		pmttp1=aa*ee-bb*ff
		psvtp1=aa*gg-bb*hh
		
		ptltp=ptltp+ptltp1
		pmttp=pmttp+pmttp1
		psvtp=psvtp+psvtp1
	enddo
	ptl=ptl+ptltp	!phi
	pmt=pmt+pmttp	!phi0+phi~
	psv=psv+psvtp	!-phi0
enddo

end
!--------------------------------------------------------------------

! find Kn
function fkn(x,n)
implicit double precision(a-h,o-z)
real*8 x
integer n
fkn=0.d0
do is=0,n
	fkn=fkn+x**is*(2**is)*fact(n)*fact(2*n-is)/(fact(is)*fact(2*n)*fact(n-is))
enddo
end

!-----------------------------------------------------------------------------
subroutine MLPMN(mm,nn,flp,xx)
IMPLICIT DOUBLE PRECISION (a-h,o-z)
DIMENSION PM(0:100,0:100),PD(0:100,0:100), flp(-nn:nn,0:nn)
CALL LPMN(100,mm,nn+1,xx,PM,PD)
flp=0.d0
DO J=0,nn
    do m=0,j
		flp(m,j)=pm(m,j)
		flp(-m,j)=pm(m,j)
		!flp(-m,j)=(-1)**m*dble(fact(j-m))/dble(fact(j+m))*pm(m,j)
		!flp(m,j)=(-1)**m*flp(m,j)
		!flp(-m,j)=(-1)**(-m)*flp(-m,j)
	enddo
enddo
END

!----------------------------------------------------------------------------
    SUBROUTINE LPMN(MM,M,N,X,PM,PD)
!C
!C       =====================================================
!C       Purpose: Compute the associated Legendre functions 
!C                Pmn(x) and their derivatives Pmn'(x)
!C       Input :  x  --- Argument of Pmn(x)
!C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
!C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
!C                mm --- Physical dimension of PM and PD
!C       Output:  PM(m,n) --- Pmn(x)
!C                PD(m,n) --- Pmn'(x)
!C       =====================================================
!C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=1.0D+300
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,M
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-(I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)*(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END

!-----------------------------------------------------------------------------
function fact(n) 
implicit double precision(a-h, o-z)
real*8 fact
fact=1.d0
if (n==0) then
	fact=1.d0
else
	do i=1,n
		fact=fact*dble(i)
	enddo
endif
end

!-------------------------------------------------------------------------
subroutine ct2sph(ctpos,sphpos)
implicit double precision(a-h,o-z)
real*8 ctpos(3), sphpos(3)
x=ctpos(1); y=ctpos(2); z=ctpos(3)
pi=acos(-1.d0)
roh=0.d0
phi=0.d0
theta=0.d0

roh=sqrt(x**2+y**2+z**2)
if (roh<1.d-10) then
	phi=0.d0
	theta=0.d0
else
	!phi
	if(abs(z)<1.d-10) then
		phi=pi/2.d0
	else
		if (z<1.d-10) then
			phi=atan(sqrt(x**2+y**2)/z)+pi
		else
			phi=atan(sqrt(x**2+y**2)/z)
		endif
	endif
	
	!theta
	if (abs(y)<1.d-10) then
		if (x>1.d-10 .or. abs(x)<1.d-10) then
			theta=0.d0
		else
			theta=pi
		endif

	elseif (abs(x)<1.d-10) then
		if (y>1.d-10) then
			theta=pi/2.d0
		elseif (y+1.d-10<0.d0) then
			theta=3.d0*pi/2.d0
		endif	
	elseif (x>1.d-10 .and. y>1.d-10) then				!I		Quadrant
		theta=atan(y/x)
	elseif (x+1.d-10<0.d0 .and. y>1.d-10) then			!II		Quadrant
		theta=atan(y/x)+pi
	elseif (x+1.d-10<0.d0 .and.  y+1.d-10<0.d0) then	!III	Quadrant
		theta=atan(y/x)+pi 
	elseif (x>1.d-10 .and. y+1.d-10<0.d0) then			!IV		Quadrant
		theta=atan(y/x)+2*pi
	endif
endif

sphpos=(/roh, phi, theta/)
End

subroutine sph2ct(sphpos,ctpos)
implicit double precision(a-h,o-z)
real*8 ctpos(3), sphpos(3)
ctpos(1)=sphpos(1)*sin(sphpos(2))*cos(sphpos(3))
ctpos(2)=sphpos(1)*sin(sphpos(2))*sin(sphpos(3))
ctpos(3)=sphpos(1)*cos(sphpos(2))
end

