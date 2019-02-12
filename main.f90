! This program solves poisson-boltzmann equation with the following treatment
! 1. Boundary Intergral Formulation derived by Juff et. al
! 2. MSMS interface extended from one sphere to molecules with multiple atoms
! 3. Curved triagular elements as derived by Zauhar et. al ---- No
! 4. Singularities are removed by transformation derived by Atkinson  ---- NO
! 5. Linear Intepolation of the collocation method derived by Atkinson  ------No
! 6. Simply the one-particle-per-element treecode as what Lu has done with FMM PB
! 7. Extension from spherical cavity to real molecules
! 8. Preconditioning: Diagnal Scalling 
program bimpb_pb_treecode_2
use molecule
use comdata
use bicg
use treecode
use treecode3d_procedures
implicit double precision(a-h,o-z)
real*8 r0(3), Pxyz(3), err_surf(10,6), err_reaction(10,6), err_reaction_rel(10,6)
real*8 pi,one_over_4pi, center(3), kappa2

external MATVEC, MSOLVE

common // pi,one_over_4pi
pi=acos(-1.d0)
one_over_4pi=0.25d0/pi

open(10,file='result.dat')
!do idens=2,7
!do level=3,3
!do order=9,9
!if (order==5 .or. order==7 .or. order==9) goto 1022 
!print *,order,level

isolver=2 ! 1: LU decomposition; 2: Bicg Gradient

eps0=1.d0;  eps1=40.d0; eps=eps1/eps0;
bulk_strength=0.d0                      !ion_strength in M
kappa2=8.430325455*bulk_strength/eps1     !kappa2 in 300K
kappa=sqrt(kappa2)                        !kappa
nt=60; !nt for order of harmonic osillators
para=1.0 !332.0716d0
call cpu_time(cpu1)

!####################################################
! read the protein structure and charges
! 1. From MSMS
!  call readin(idens,ichrpos)
! 2. From Icosahedron
center=(/0.d0,0.d0,0.d0/)
level=1.0d0
call readin_sphere_tri(level,center)

!call lagendini ! initialization for harmonic osillators
print *,'Begin to form linear algebraic matrix...'

call cpu_time(cpu2)

! write the discritized integral equation
call form_matrix

print *,'Begin to initialize treecode...'
call treecode_initialization

call cpu_time(cpu23)
print *,'it takes ',cpu23-cpu2,'seconds to form the matrix'
! To solve by GMRES

! parameter for biconjugate method
thresh=1.d-20
itol=2
itmax=1000000
tol=1.d-6
ndim=nface
print *,'Begin to allocate varibles for the solver...'

allocate(sb(ndim),sx(ndim),STAT=jerr)
if (jerr .ne. 0) then
    write(*,*) 'Error allocating sb and sx'
endif

! setup the parameters of the GMRES solver
MAXL=10		! Maximum dimension of Krylov subspace in which X - X0 is to be found
LRGW=1 + ndim*(MAXL+6) + MAXL*(MAXL+3)		! Length of the double precision workspace, RGWK.
JSCAL=0		! Flag indicating whether the scaling arrays SB and SX are to be used
JPRE=-1		! Flag indicating whether preconditioning is being used
NRMAX=10	! Maximum number of restarts of the Krylov iteration
LIGW=20
MLWK=LIGW	! Required minimum length of RGWK array
NMS=0		! The total number of calls to MSOLVE
ISYM=0		! If the symmetric matrix is stored only in half way
!lenw =  10*ndim; leniw = 10*ndim 
lenw =  10; leniw = 10 

allocate(RGWK(LRGW), IGWK(LIGW), RWORK(lenw), IWORK(leniw), STAT=jerr)
if (jerr .ne. 0) then
	Write(*,*) 'Error allocating RGWK, IGWK, RWORK, IWORK, jerr= ', jerr 
	write(*,*) 'LRGW=',LRGW,'LIGW=',LIGW,'lenw= ',lenw,'leniw=', leniw
	stop
endif
RGWK=0.d0;	IGWK=0; RWORK=0.d0;	IWORK=0
IGWK(1:7)=(/MAXL, MAXL, JSCAL, JPRE, NRMAX, MLWK, NMS/)

print *,'Begin to call the solver...'
call DGMRES(	ndim, bvct, xvct, MATVEC, MSOLVE, ITOL, TOL, ITMAX, & 
				ITER, ERR, IERR, 0, SB, SX, RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)

print *,'err=',err,'ierr=',ierr,'iter=',iter

call cpu_time(cpu3)


print *, 'surpot normal derivative=', xvct
         
! Calculate the error of potential at a specific point
soleng=0.d0; soleng_exa=0.d0
do iatm=1,nchr
    r0=chrpos(:,iatm)
    call potential_molecule(r0, ptl,0)
    soleng=soleng+atmchr(iatm)*0.5d0*para*4*pi*ptl
!    soleng_exa=soleng_exa+atmchr(iatm)*0.5d0*para*4*pi*
enddo


print *,'Solvation energy=:' ,soleng

call cpu_time(cpu4)
print *,'setup cpu=', real(cpu2-cpu1), '  solving cpu=', real(cpu3-cpu2), &
'Total cpu= ', real(cpu4-cpu1)
stop

call output_potential

do i=1,nface
    write(level*10+order,*) i,xvct(i),xvct(nface+i)
enddo


! Calculate the interface error
call exact_solution_analytic
!call exact_solution_kirkwood

call ERROR(xvct(1:nface),F_EXA(1:nface),e1,e2,e3,nface)
print *,'potential  errors= ',real(e1),real(e2),real(e3)
!call ERROR(xvct(nface+1:2*nface),F_EXA(nface+1:2*nface),en1,en2,en3,nface)
!print *,'norm field errors= ',real(en1),real(en2),real(en3)


!write(10,'(2i4,3f16.8)') order,idens,real(cpu3-cpu2),real(cpu4-cpu1),soleng
write(10,'(2i4,4f16.8)') order,level,real(cpu3-cpu2),soleng,e2,en2

! Kirkwood solutions
!pxyz=chrpos(:,1)
!call kirk_potential(-1,pxyz,ptl_kirk,pmt_kirk,psv_kirk)

!print *,'kirk solvation energy: ',0.5d0*pmt_kirk*332.0716d0

!ptl_exa=pmt_kirk/4.d0/pi
!print *,ptl_kirk/4.0/pi,pmt_kirk/4.d0/pi,psv_kirk/4.d0/pi

!call potential_exact_analytic(r0, ptl_exa)

!print *, 'at ', real(r0), 'true value: ', real(ptl_exa)
!write (*,*)  real(ptl),real(abs(ptl-ptl_exa)),real(abs(ptl_semiexa-ptl_exa))

! deallocate memory
deallocate(bvct, xvct, F_exa, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating bvct, xvct, F_exa'
    stop
endif

deallocate(tr_xyz,tr_q,tchg,schg,tr_area,kk,der_cof, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating tr_xyz,tr_q,tchg,schg,tr_area,kk'
    stop
endif


deallocate(SB,SX, RGWK, IGWK, RWORK, IWORK, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating SB,SX, RGWK, IGWK, RWORK, IWORK'
    stop
endif 

deallocate(atmpos,atmrad,atmchr,chrpos, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating atmpos,atmrad,atmchr,chrpos'
    stop
endif 

DEALLOCATE(x,y,z,q,orderind,STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(*,*) 'Error deallocating x, y, z, q, or orderind!'
    STOP
END IF

DEALLOCATE(SPTPOS, SPTNRM, NATMAFF, NSFTYPE, NVERT, MFACE, STAT= ierr)
IF (ierr .NE. 0) THEN
    WRITE(*,*) 'Error deallocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE, NVERT, MFACE !'
    STOP
END IF

DEALLOCATE(cf, cf1, cf2, cf3, a, b,STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(*,*) 'Error allocating Taylor variables cf, cf1, cf2, cf3, a, b ! '
    STOP
END IF

DEALLOCATE(orderarr,STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(6,*) 'Error deallocating copy variables orderarr! '
    STOP
END IF

!deallocate(chrpos_sph,chgmnx)
!err_surf(ichrpos,idens)=e1
!err_reaction(ichrpos,idens)=abs(ptl-ptl_exa)
!err_reaction_rel(ichrpos,idens)=abs((ptl-ptl_exa)/ptl_exa)
1022 continue
!enddo
!enddo
!write(21,'(6ES12.5)') err_surf(ichrpos,:)
!write(22,'(6ES12.5)') err_reaction(ichrpos,:) 
!write(23,'(6ES12.5)') err_reaction_rel(ichrpos,:) 
close(10)
end program bimpb_pb_treecode_2

!#######################################################################################
Subroutine output_potential
use comdata
use molecule
use treecode
use treecode3d_procedures
implicit none
real*8, dimension(:), allocatable :: x_temp
real*8, dimension(:,:), allocatable :: xyz_temp
integer i,j,k,jerr,nface_vert
real*8 tot_length,loc_length,aa(3),pi,para_temp,one_over_4pi,phi_star


common // pi,one_over_4pi

allocate(xtemp(numpars),xyz_temp(3,numpars),STAT=jerr)
if (jerr .ne. 0) then
    write(*,*) 'Error allocating xtemp, xyz_temp'
endif
xtemp=0.d0; xyz_temp=0.d0

do i=1,numpars
  xtemp(orderarr(i))=xvct(i)           !put things back
  xyz_temp(:,orderarr(i))=tr_xyz(:,i)
enddo
xvct=xtemp
tr_xyz=xyz_temp

deallocate(xtemp,xyz_temp, STAT=jerr)
if (jerr .ne. 0) then
    write(*,*) 'Error deallocating xtemp, xyz_temp'
endif

End








!-----------------------------------------------------------------------
subroutine treecode_initialization
use molecule
use bicg
use comdata
use treecode
use treecode3d_procedures
implicit none

real*8 pi, one_over_4pi
common // pi,one_over_4pi

! local variables

INTEGER :: level,ierr,err,i,j,k,mm,nn,idx,ijk(3)

! variables needed for cpu time

REAL*8 :: totaltime,timetree
real*8, dimension(:), allocatable:: temp_a,temp_b
real*8, dimension(:,:), allocatable:: temp_q

! Setup Treecode parameters
order=0		! The order of taylor expansion
iflag=1			! force needs to be evaluated or not 
maxparnode=5000000 ! maximum particles per leaf
theta=0.0d0     ! MAC, rc/R<MAC, the bigger MAC, the more treecode is applied

allocate(kk(3,4), der_cof(0:order,0:order,0:order,4), STAT=ierr)	
if (ierr .ne. 0) then
	Write(*,*) 'Error allocating auxilary Taylor coefficients kk and der_ncf'
	stop
endif

! The adjustment of k for the recurrance relation 
kk(:,1)=(/0,0,0/);        ! Original Kernel

kk(:,2)=(/1,0,0/);        ! 1st Order Derivative:	partial x
kk(:,3)=(/0,1,0/);        !		                    partial y           
kk(:,4)=(/0,0,1/);        !                         partial z

! The adjustment of der_cof for the recurrance relation
der_cof=1.d0

DO idx=1,4
    DO k=0,order
        DO j=0,order-k
            DO i=0,order-k-j
                ijk=(/i,j,k/)
                DO mm=1,3
                    IF (kk(mm,idx) .ne. 0) THEN
                        DO nn=1,kk(mm,idx)
                            der_cof(i,j,k,idx)=der_cof(i,j,k,idx)*(ijk(mm)+nn)
                        ENDDO
                    ENDIF
                ENDDO
            ENDDO
         ENDDO
     ENDDO
ENDDO

der_cof=der_cof*one_over_4pi
numpars=nface

ALLOCATE(x(numpars),y(numpars),z(numpars),q(numpars),orderind(numpars),STAT=err)
IF (err .NE. 0) THEN
    WRITE(6,*) 'Error allocating x, y, z, q, or orderind!'
    STOP
END IF


allocate(temp_a(numpars),temp_b(numpars),temp_q(3,numpars), STAT=err)
IF (err .NE. 0) THEN
    WRITE(6,*) 'Error allocating temp_a, temp_b, temp_q!'
    STOP
END IF

      
x=tr_xyz(1,:)
y=tr_xyz(2,:)
z=tr_xyz(3,:)
q=1.d0


! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables. Also, copy variables into global copy arrays. 
CALL SETUP(x,y,z,q,numpars,order,iflag,xyzminmax)

! nullify pointer to root of tree (TROOT) and create tree
NULLIFY(troot)  

! creating tree

level=0
minlevel=50000
maxlevel=0

WRITE(6,*) ' '
WRITE(6,*) 'Creating tree for ',numpars,' particles with max ', maxparnode, ' per node...'

CALL CPU_TIME(timebeg)
CALL CREATE_TREE(troot,1,numpars,x,y,z,q,maxparnode,xyzminmax,level,numpars)

temp_a=tr_area
temp_b=bvct
temp_q=tr_q

do i=1,numpars
  tr_area(i)=temp_a(orderarr(i))
  tr_q(:,i)=temp_q(:,orderarr(i))
  bvct(i)=temp_b(orderarr(i))
  !bvct(i+numpars)=temp_b(orderarr(i)+numpars)
  tr_xyz(:,i)=(/x(i),y(i),z(i)/)
enddo

CALL CPU_TIME(timeend)
totaltime=timeend-timebeg
WRITE(6,*) 'Time to create tree (secs):',totaltime      
      
deallocate(temp_a,temp_b,temp_q, STAT=err)
IF (err .NE. 0) THEN
    WRITE(6,*) 'Error deallocating temp_a, temp_b, temp_q!'
    STOP
END IF
End subroutine



!----------------------------------
subroutine MATVEC(N, XX, bb)
use bicg
use molecule
use comdata
use treecode3d_procedures
implicit double precision(a-h,o-z)
integer N
real*8 xx(N),bb(N),timebeg,timeend

call cpu_time(timebeg)
if (sum(abs(xx))<1.d-10) goto 1022
CALL TREE_COMPP_PB(troot,kappa,eps,xx)
!stop
1022 bb=xx
call cpu_time(timeend)
print *,'time to compute AX=: ',timeend-timebeg 
CALL REMOVE_MMT(troot)

return
end subroutine

!-------------------------------------
subroutine MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
use molecule
implicit double precision(a-h,o-z)
real*8 R(N),Z(N),A(N*N),RWORK(*)
integer IA(N*N), JA(N*N), IWORK
scale1=0.5d0*(1.d0+1.d0/eps)
!scale2=0.5d0*(1.d0+1.d0/eps)
Z(1:N)=R(1:N)/scale1
!Z((N/2+1):N)=R((N/2+1):N)/scale2
end subroutine

!-----------------------------------
! input the data into the matrix
subroutine form_matrix
use molecule
use comdata
use treecode
implicit double precision(a-h,o-z)
integer idx(3), istag, NGR
real*8 r0(3), v0(3),v(3,3), r(3,3), r1(3), v1(3), uv(2,10), x10(3,10),v10(3,10)

! tr_xyz: The position of the particles on surface
! tr_q:	  The normail direction at the particle location
! bvct:	  The right hand side of the pb equation in BIM form
! xvct:	  The vector composed of potential and its normal derivative
! tchg:	  Charges on Target particles
! schg:   Charges on Source particles
! tr_area: the triangular area of each element
allocate(tr_xyz(3,nface),tr_q(3,nface), bvct(nface), xvct(nface))
allocate(tchg(nface,4),schg(nface,4),tr_area(nface))
tr_xyz=0.d0;	tr_q=0.d0;	bvct=0.d0;	xvct=0.d0
tchg=0.d0;		schg=0.d0;	tr_area=0.d0

 
do i=1,nface    ! for phi on each element
    idx=nvert(1:3,i) ! vertices index of the specific triangle
    r0=0.d0;    v0=0.d0
    do k=1,3 
        r0=r0+1.d0/3.d0*sptpos(1:3,idx(k))	!centriod
        v0=v0+1.d0/3.d0*sptnrm(1:3,idx(k))	
	    r(:,k)=sptpos(1:3,idx(k))
	    v(:,k)=sptnrm(1:3,idx(k))
    enddo
    
    ! Way 1: normlize the midpoint v0
    v0=v0/sqrt(dot_product(v0,v0))
	
    !#######################################################################################
    !modification if it is not a sphere
    !r0=r0/sqrt(dot_product(r0,r0))*rds     ! project to the sphere surface
    !v0=r0/rds                              !for sphere only, need to change if for molecule
    !######################################################################################
    tr_xyz(:,i)=r0			! Get the position of particles
    tr_q(:,i)=v0			! Get the normal of the paricles, acting as charge in treecode
    
    aa=sqrt(dot_product(r(:,1)-r(:,2),r(:,1)-r(:,2)))
    bb=sqrt(dot_product(r(:,1)-r(:,3),r(:,1)-r(:,3)))
    cc=sqrt(dot_product(r(:,2)-r(:,3),r(:,2)-r(:,3)))
    tr_area(i)=triangle_area(aa,bb,cc)
    							
    ! setup the right hand side of the system of equations
    do j=1,nchr ! for each atom
    !bvct(i)=bvct(i)+atmchr(j)*G0(chrpos(1:3,j),r0)/eps0
 	bvct(i)=bvct(i)+atmchr(j)*G1(v0,chrpos(1:3,j),r0)/eps0 !gan: right-handside for surpot normal derivative
    enddo
enddo
end subroutine

!--------------------------------------------------------------------------------------------------------------
SUBROUTINE ERROR(F,F_EXA, e1,e2,e3,n)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION F(1:n),F_EXA(1:n)
! e1: L infinity error
! e2: Relative L infinity error
! e3: L2 error
e1=0.0d0
e2=0.0d0
e3=0.0d0
DO I=1,n
!	if (abs(f_exa(i)) < 1.d-15) print *,'relative error failed',f_exa(i)
	err0=ABS(F(I)-F_EXA(I))
	err1=ABS((F(I)-F_EXA(I))/f_exa(i))
	e3=e3+err0**2
    IF (e1<err0) e1=err0
	IF (e2<err1) e2=err1
END DO

e3=sqrt(e3/n)
RETURN

END

!------------------------------------------------------------------------------------------
subroutine exact_solution_analytic
use molecule
use comdata
implicit double precision(a-h,o-z)

common // pi,one_over_4pi

real*8 r(3),r0(3),v0(3),v(3)
integer idx(3)

allocate(F_exa(2*nface))
F_exa=0.d0

do i=1,nface
	idx=nvert(1:3,i) ! vertices index of the specific triangle
    r0=0.d0;    v0=0.d0
    do k=1,3 
        r0=r0+1.d0/3.d0*sptpos(1:3,idx(k))	!centriod
		!r(:,k)=sptpos(1:3,idx(k))
		!v(:,k)=sptnrm(1:3,idx(k))
    enddo
    r0=r0/sqrt(dot_product(r0,r0))*rds          ! project to the sphere surface
    
	rs=sqrt(dot_product(r0,r0))
	if (rs-rds>1.d-10) then	!outside
		F_exa(i)=(1/4.d0/pi)/eps/(1.d0+kappa*rds)/rs*exp(-kappa*(rs-rds))
		F_exa(i+nface)=(1.d0/4.d0/pi)/eps/(1.d0+kappa*rds)*(exp(-kappa*(rs-rds))*(-kappa*rs-1.d0))/rs**2
	    F_exa(i+nface)=F_exa(i+nface)*eps
!	    print *,'outside',i,-1.d0/4.d0/pi/rs**2,F_exa(i+nface),pi,rds,rs,exp(-kappa*(rs-rds)),(kappa*rs+1.d0),(1.d0+kappa*rds),eps
	else								!inside
		F_exa(i)=(1/4.d0/pi)*(1.d0/eps/(1.d0+kappa*rds)/rds-1.d0/rds+1.d0/rs)
		F_exa(i+nface)=-1.d0/4.d0/pi/rs**2
	!	print *,'in side',i,F_exa(i),F_exa(i+nface)
	endif
enddo

! Careful, charge could be different from one
F_exa=F_exa*atmchr(1)
close(13)


End


!---------------------------------------------------------------------------------------------------------------
subroutine exact_solution_kirkwood
use molecule
use comdata
implicit double precision(a-h,o-z)

real*8 r(3),r0(3)
integer iside

allocate(F_exa(2*nspt))
F_exa=0.d0

do i=1,nspt
	r=sptpos(1:3,i)
	rs=sqrt(dot_product(r,r))
	if (rs-rds>1.d-10) then	!outside
	    iside=1
	    call kirk_potential(iside,r,ptl_kirk,pmt_kirk,psv_kirk)
		F_exa(i)=ptl_kirk/4.d0/pi
		
		F_exa(i+nspt)=(1.d0/4.d0/pi)/eps/(1.d0+kappa*rds)*(exp(-kappa*(rs-rds))*(-kappa*rs-1.d0))/rs**2
	    F_exa(i+nspt)=F_exa(i+nspt)*eps
	    !print *,'outside',i,-1.d0/4.d0/pi/rs**2,F_exa(i+nspt),pi,rds,rs,exp(-kappa*(rs-rds)),(kappa*rs+1.d0),(1.d0+kappa*rds),eps
	else								!inside
		iside=-1
		call kirk_potential(iside,r,ptl_kirk,pmt_kirk,psv_kirk)
		F_exa(i)=ptl_kirk/4.d0/pi
		
		F_exa(i+nspt)=-1.d0/4.d0/pi/rs**2
		!print *,'in side',i,F_exa(i),F_exa(i+nspt)
	endif
enddo
close(13)


End

!---------------------------------------------------------------------------------------------------------------
! This subroutine calculate the potential inside or outside given those on surface
subroutine potential_molecule(r0,ptl,iexa)
use molecule
use comdata
use treecode
implicit double precision(a-h,o-z)
real*8 ptl,r0(3),r(3),v(3),s(3)
common // pi,one_over_4pi
ptl=0.d0

do j=1,nface ! for each triangle
    r=tr_xyz(1:3,j)
    v=tr_q(:,j)
          
    s=r0
    rs=sqrt(dot_product(r-s,r-s))
    G0=one_over_4pi/rs
    cos_theta=dot_product(v,r-s)/rs
    tp1=G0/rs
    G1=cos_theta*tp1
    area=tr_area(j)
    ptl=ptl-(eps-1.d0)*G1*xvct(j)*tr_area(j)
  
enddo

End
!------------------------------------------------------------------------------------------------------------------

subroutine potential_exact_analytic(r0,ptl)
use molecule
use comdata
implicit double precision(a-h,o-z)

real*8 r0(3)


common // pi,one_over_4pi

rs=sqrt(dot_product(r0,r0))

! Getting the reaction potential
if (rs-rds>1.d-10) then	!outside
	ptl=(1/4.d0/pi)/eps/(1.d0+kappa*rds)/rs*exp(-kappa*(rs-rds))
else								!inside
	ptl=(1/4.d0/pi)*(1.d0/eps/(1.d0+kappa*rds)/rds-1.d0/rds) !+1.d0/rs)
endif 

End


!----------------------------------------------------------------------
! This subroutine computer the source charge and target charge for the treecode
! Total Number of Kernel = 2*(1+3*2+3*3)=32
! Refer to the table in the paper for detail

subroutine pb_kernel(phi)
use treecode
use treecode3d_procedures
implicit double precision(a-h,o-z)
integer ikp,ixyz,jxyz,indx !ixyz: source; jxyz target;
real*8 phi(numpars) 


	indx=0
	indx=indx+1
	tchg(:,indx)=1.d0
	schg(:,indx)=tr_area*phi
	do ixyz=1,3
		indx=indx+1
		tchg(:,indx)=1.d0
		schg(:,indx)=tr_q(ixyz,:)*tr_area*phi
	enddo
end
