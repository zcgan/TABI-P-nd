!Read protein structure from msms
subroutine readin(idens,ichrpos)
!use IFPORT 
!use mesh_procedures
use molecule
use comdata
implicit double precision(a-h,o-z)
real*8 pos(3),vector(3)
integer nind(5)
CHARACTER(100) :: FHEAD, dens(10)

!dens=(/'1','2','5','10','20','40','80','160','320','640'/)
!dens=(/'1','2','4','8','16','32','64','128','256','512'/)

!Obtain path
	!pathname='e:\My_Programs\BIMPB\bimpb_poisson_treecode_1\'
	!pathname='~/Desktop/My_Programs/bimpb_pb_treecode_1/'
	pathname='~/bimpb/bimpb_poisson_treecode_1/'
	!pathname='/home/geng/bimpb/bimpb_pb_treecode_21/'
        lenpath = len(pathname)
	do while (pathname(lenpath:lenpath) .eq. ' ')
		lenpath = lenpath - 1
	enddo  


!Obtain protein name
	!write(*,*) 'Please input the protein name:'
	!read(*,*) fname
	fname='1a63'
	!fname='1ajj'
	
	write(*,*) 'Please input MSMS triangulation density in # per astrong^2:'
	read(*,*) den
!        den='1'
	!den=dens(idens)
		
	lenfname = len(fname)
	do while (fname(lenfname:lenfname) .eq. ' ')
		lenfname = lenfname - 1
	enddo 


!Read atom coordinates and partial charges
	OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")            
	NATM = 0
		DO 
			READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
			IF(MEOF .LT. 0) EXIT
			NATM = NATM + 1
		END DO         
	CLOSE(1)
	
	OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".pqr")            
	nchr = 0
		DO 
			READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
			IF(MEOF .LT. 0) EXIT
			nchr = nchr + 1
		END DO         
	CLOSE(1)
	
	if (nchr .ne. natm) print *,'atoms completely inside atoms are deleted'
	
    !nchr=natm
	allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
	IF (ierr .NE. 0) THEN
    		WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
    		STOP
	END IF


	open(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")
	do i=1,natm
		read(1,*) atmpos(:,i), atmrad(i)
	enddo
	close(1)
    
    ! For sphere only
    rds=atmrad(1)
	open(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".pqr")
	do i=1,nchr
		read(1,*) chrpos(:,i),atmchr(i)
	enddo
	close(1)
	
	!chrpos(1,1)=(ichrpos-1)*0.1d0                               ! For multiple runs
	

    rslt=system('msms -if '//fname(1:lenfname)//".xyzr"//' -prob 1.4 -de '//den(1:5)//' -of '//fname(1:lenfname))    
! read the surface points
      OPEN(2,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")

      READ(2,*) FHEAD
      READ(2,*) FHEAD
      READ(2,*) NSPT, ppp, qqq, rrr

      ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT), STAT= ierr)
      IF (ierr .NE. 0) THEN
          WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE!'
      STOP
      END IF

	SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	
      
      DO I=1,NSPT
         READ(2,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
		 
	 !radial projection to improve msms accuracy, ONLY FOR SPHERE!!!!!!!!
	 !pos=pos/sqrt(dot_product(pos,pos))*rds;	vector=pos/rds;
         
		 SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
         NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
      END DO
	 
      CLOSE(2)

! read the surface triangulization

      OPEN(3,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")

      READ(3,*) FHEAD
      READ(3,*) FHEAD
      READ(3,*) NFACE, PPP, QQQ, RRR

      ALLOCATE(NVERT(3,NFACE), MFACE(NFACE), STAT=ierr)
      IF (ierr .NE. 0) THEN
          WRITE(6,*) 'Error allocating NVERT, MFACE'
      STOP
      END IF

      NVERT=0; MFACE=0
      
      DO I=1,NFACE 
         READ(3,*) NIND(1:5) 
         NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
      END DO
      CLOSE(3)
	  call surface_area(s_area)
	  !s2=4.d0*pi*rds**2
	  print *,'surface area=', real(s_area)


	  

End
!------------------------------------------------------------------------
subroutine surface_area(s_area)
use molecule
implicit double precision(a-h,o-z)
integer iface(3),jface(3),nfacenew,ialert
real*8 face(3,3),s_area,face_old(3,3),xx(3),yy(3),cpu1,cpu2
real*8,allocatable:: nvert_copy(:,:)
        print *,'# of surfaces=',nface,' # of surface points=',nspt
        call cpu_time(cpu1)
        s_area=0.d0
        nfacenew=nface
        allocate(nvert_copy(3,nface))
        nvert_copy=nvert
        do i=1,nface
                iface=nvert(:,i)
                xx=0.d0
                ialert=0;
                do j=1,3
                    face(:,j)=sptpos(:,iface(j))
                    xx=xx+1/3.d0*(face(:,j))
                enddo
                aa=sqrt(dot_product(face(:,1)-face(:,2),face(:,1)-face(:,2)))
                bb=sqrt(dot_product(face(:,1)-face(:,3),face(:,1)-face(:,3)))
                cc=sqrt(dot_product(face(:,2)-face(:,3),face(:,2)-face(:,3)))
                area_local=triangle_area(aa,bb,cc)
                !goto 1022
                if (2==2) then 
                do ii=max(1,i-10),i-1
                    jface=nvert(:,ii)
                    yy=0.d0
                    do j=1,3
                        face_old(:,j)=sptpos(:,jface(j))
                        yy=yy+1/3.d0*(face_old(:,j))
                    enddo
                    dist_local=dot_product(xx-yy,xx-yy)
                    if (dist_local<1.d-6) then
                       ialert=1
                       print *,i,ii,'particles are too close',dist_local
                    endif
                enddo
                endif
                if (area_local < 1.d-5 .or. ialert==1) then
                    print *,i,j,'small area=', area_local
                    ichanged=nface-nfacenew
                    nvert_copy(:,i-ichanged:nface-1)=nvert_copy(:,i-ichanged+1:nface)
                    nfacenew=nfacenew-1
                endif
1022 continue
                s_area=s_area+area_local
        enddo
        print *,nface-nfacenew,' ugly faces are deleted'
        nface=nfacenew
        deallocate(nvert)
        allocate(nvert(3,nface))
        nvert=nvert_copy
        deallocate(nvert_copy)
        call cpu_time(cpu2)
        print *, 'total MSMS post-processing time =',cpu2-cpu1
        !print *,i,j,real(triangle_area(aa,bb,cc)),surface_area
        !pause

end



!------------------------------------------------------------------------
function triangle_area(aa,bb,cc)
implicit double precision(a-h,o-z)
	s=0.5d0*(aa+bb+cc)
	triangle_area=sqrt(s*(s-aa)*(s-bb)*(s-cc))
end 


!#############################################################################################
!This subroutine gives the triangulation of a sphere starting from icosahedron
Subroutine readin_sphere_tri(level,center)
use molecule
use comdata
!program sphere
implicit none
real*8 rad, center(3), xxx,yyy,zzz,rrr
integer init_type,level,iwind,nedge,err,i,meof,ierr
!integer nface,nspt
!real*8, allocatable:: sptpos(:,:),sptnrm(:,:)
!integer, allocatable:: nvert(:,:)
integer, allocatable:: nefv(:,:)
real*8 surface_area

	!Obtain path
	!pathname='e:\My_Programs\BIMPB\bimpb_poisson_treecode_1\'
	!pathname='~/Desktop/My_Programs/bimpb_pb_treecode_1/'
	pathname='/home/gzc/bie_only_poisson_for_normal_derivative/'
	!pathname='/home/gengweih/bimpb/bimpb_poisson_treecode_1/'


	lenpath = len(pathname)
	do while (pathname(lenpath:lenpath) .eq. ' ')
		lenpath = lenpath - 1
	enddo  

	!Obtain protein name
	!write(*,*) 'Please input the protein name:'
	!read(*,*) fname

	fname='oneb'
	lenfname = len(fname)
	do while (fname(lenfname:lenfname) .eq. ' ')
		lenfname = lenfname - 1
	enddo 
	
	!Read atom coordinates and partial charges
	OPEN(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")            
	NATM = 0
		DO 
			READ(1,*,IOSTAT = MEOF) xxx, yyy, zzz, rrr
			IF(MEOF .LT. 0) EXIT
			NATM = NATM + 1
		END DO         
	CLOSE(1)

    nchr=natm
	
	allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
	IF (ierr .NE. 0) THEN
    		WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
    		STOP
	END IF


	open(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".xyzr")
	do i=1,natm
		read(1,*) atmpos(:,i), atmrad(i)
	enddo
	close(1)
        rad=atmrad(1); rds=rad
	open(1,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".pqr")
	do i=1,nchr
		read(1,*) chrpos(:,i),atmchr(i)
	enddo
	close(1)
! These two are for future usage
init_type=1;	iwind=1

!level=3
!center=(/0.d0, 0.d0, 0.d0/)
!rad=1.d0

! For the icosahedron, initially 20 faces, 12 vertex and 32 edges
! With the relation V+F=E+2 and 12+20=30+2
nface=20;	nspt=12;	nedge=nface+nspt-2

! nefv contains the number of edges,faces and vertexes on each level
allocate(nefv(3,0:level),STAT=err)
IF (err .NE. 0) THEN
	WRITE(6,*) 'Error allocating number of edges,vertex,and face'
	STOP
END IF
nefv=0;	nefv(:,0)=(/nedge,nface,nspt/)

do i=1,level
	nedge=nedge*2+3*nface
	nface=nface*4
	nspt=nedge+2-nface
	nefv(:,i)=(/nedge,nface,nspt/)
enddo

! allocate varibles for triangulation
allocate(sptpos(3,nspt),sptnrm(3,nspt),nvert(3,nface), STAT=err)
IF (err .NE. 0) THEN
	WRITE(6,*) 'Error allocating vertex and faces'
	STOP
END IF
sptpos=0.d0;	nvert=0

ALLOCATE(NATMAFF(NSPT), NSFTYPE(NSPT), MFACE(NFACE), STAT= ierr)
      IF (ierr .NE. 0) THEN
          WRITE(6,*) 'Error allocating NATMAFF, NSFTYPE, MFACE!'
      STOP
      END IF

call sphere_tri(init_type,level,rad,center,iwind,nspt,nface,sptpos,nvert,nefv)
sptnrm=sptpos

!pi=acos(-1.d0)
!print *,surface_area(nface,nspt,nvert,sptpos)-4.d0*pi*rad**2

End


!--------------------------------------------------------------------------------------
Subroutine sphere_tri(init_type,level,rad,center,iwind,nspt,nface,sptpos,nvert,nefv)
! init_type: 1: tetrahedron; 2:
implicit none
integer init_type,level,iwind,nvert(3,nface),nspt,nface,ico_f0(3,20),i,nefv(3,0:level)
real*8 rad,center(3),sptpos(3,nspt),t,tau,one,ico_v0(3,12)

!Twelve vertices of icosahedron on unit sphere

!tau = 0.8506508084; one = 0.5257311121
t=(1.d0+sqrt(5.d0))/2.d0
tau=t/sqrt(1.d0+t**2)
one=1/sqrt(1.d0+t**2)

! 12 vertex
ico_v0(:,1) = (/  tau,  one,	0.d0 /); ! ZA
ico_v0(:,2) = (/ -tau,  one,	0.d0 /); ! ZB
ico_v0(:,3) = (/ -tau, -one,	0.d0 /); ! ZC
ico_v0(:,4) = (/  tau, -one,	0.d0 /); ! ZD
ico_v0(:,5) = (/  one, 0.d0,	tau  /); ! YA
ico_v0(:,6) = (/  one, 0.d0,	-tau /); ! YB
ico_v0(:,7) = (/ -one, 0.d0,	-tau /); ! YC
ico_v0(:,8) = (/ -one, 0.d0,	tau  /); ! YD
ico_v0(:,9) = (/ 0.d0,  tau,	one  /); ! XA
ico_v0(:,10)= (/ 0.d0, -tau,	one  /); ! XB
ico_v0(:,11)= (/ 0.d0, -tau,	-one /); ! XC
ico_v0(:,12)= (/ 0.d0,  tau,	-one /); ! XD

   
! Structure for unit icosahedron
! 20 faces
ico_f0(:,1) =	(/5,  8,  9/)
ico_f0(:,2) =	(/5, 10,  8/)
ico_f0(:,3) =	(/6, 12,  7/)
ico_f0(:,4) =	(/6,  7, 11/)
ico_f0(:,5) =	(/1,  4,  5/)
ico_f0(:,6) =	(/1,  6,  4/)
ico_f0(:,7) =	(/3,  2,  8/)
ico_f0(:,8) =	(/3,  7,  2/)
ico_f0(:,9) =	(/9, 12,  1/)
ico_f0(:,10)=	(/9,  2, 12/)
ico_f0(:,11)=	(/10, 4, 11/)
ico_f0(:,12)=	(/10, 11, 3/)
ico_f0(:,13)=	(/9,  1,  5/)
ico_f0(:,14)=	(/12, 6,  1/)
ico_f0(:,15)=	(/5,  4, 10/)
ico_f0(:,16)=	(/6, 11,  4/)
ico_f0(:,17)=	(/8,  2,  9/)
ico_f0(:,18)=	(/7, 12,  2/)
ico_f0(:,19)=	(/8, 10,  3/)
ico_f0(:,20)=	(/7,  3, 11/)

sptpos(:,1:12)=ico_v0
nvert(:,1:20)=ico_f0


do i=1,level
	call mesh_refine_tri4(i,nspt,nface,sptpos,nvert)
	call sphere_project(nspt,sptpos,rad,nefv(3,i),center);
enddo

End


!---------------------------------------------------------------------------------------
Subroutine sphere_project(nspt,sptpos,rad,nspt_level,center)
implicit none
integer nspt,nspt_level,isptpos	
real*8 sptpos(3,nspt),rad,center(3),x(3),xnorm

do isptpos=1,nspt_level
	x=sptpos(:,isptpos)-center
	xnorm=sqrt(dot_product(x,x))
	x=x/xnorm
	x=x*rad
	sptpos(:,isptpos)=x+center
enddo


End

!---------------------------------------------------------------------------------------
Subroutine mesh_refine_tri4(ilevel,nspt,nface,sptpos,nvert)
implicit none
real*8, allocatable :: sptpos0(:,:)
integer, allocatable:: nvert0(:,:)

integer ilevel,i,nspt,nface,nface0,nspt0,nedge0,nvert(3,nface),err,NABC(3)
integer indx_sptpos,iface,NA,NB,NC,ifind
real*8 sptpos(3,nspt),ABC(3,3),A(3),B(3),C(3)

nface0=20
nspt0=12
nedge0=nface0+nspt0-2
do i=1,ilevel-1
	nedge0=nedge0*2+3*nface0
	nface0=nface0*4
	nspt0=nedge0+2-nface0
enddo

allocate(sptpos0(3,nspt0),nvert0(3,nface0), STAT=err)
IF (err .NE. 0) THEN
	WRITE(6,*) 'Error allocating local vertex and faces'
	STOP
END IF
sptpos0=0.d0
nvert0=0

sptpos0(:,1:nspt0)=sptpos(:,1:nspt0)
nvert0(:,1:nface0)=nvert(:,1:nface0)

indx_sptpos=nspt0

do iface = 1,nface0
!Get the triangle vertex indices
    NA = nvert0(1,iface)
    NB = nvert0(2,iface)
    NC = nvert0(3,iface)

!Get the triangle vertex coordinates
    A = sptpos0(:,NA)
    B = sptpos0(:,NB)
    C = sptpos0(:,NC)
    
!Now find the midpoints between vertices
    ABC(:,1) = (A + B) / 2;
    ABC(:,2) = (B + C) / 2;
    ABC(:,3) = (C + A) / 2;
    
!Store the midpoint vertices, while
!checking if midpoint vertex already exists

    do i=1,3
		call mesh_find_vertex(indx_sptpos,ifind,sptpos,nspt,ABC(:,i));
		if (ifind==0) then
			indx_sptpos=indx_sptpos+1
			sptpos(:,indx_sptpos)=ABC(:,i)
			NABC(i)=indx_sptpos
		else
			NABC(i)=ifind
		endif
	enddo
	    
!Create new faces with orig vertices plus midpoints

    nvert(:,iface*4-3) = (/ NA, NABC(1), NABC(3) /)
    nvert(:,iface*4-2) = (/ NABC(1), NB, NABC(2) /)
    nvert(:,iface*4-1) = (/ NABC(3), NABC(2), NC /)
    nvert(:,iface*4-0) = (/ NABC(1), NABC(2), NABC(3)/)
    
enddo

deallocate(sptpos0,nvert0, STAT=err)
IF (err .NE. 0) THEN
	WRITE(6,*) 'Error deallocating local vertex and faces'
	STOP
END IF

End

!---------------------------------------------------------------------------------
! For a new vertex, compared with all the stored vertex
! ifind .ne. 0:	if the vertex is already stored
! ifind=0:		if the vertex is completely new
Subroutine mesh_find_vertex(indx_sptpos,ifind,sptpos,nspt,sptpos_new)
implicit none
real*8 sptpos(3,nspt),sptpos_new(3),diff(3)
integer indx_sptpos,ifind,nspt,isptpos

ifind=0

do isptpos=1,indx_sptpos
	diff=sptpos_new-sptpos(:,isptpos)
	if (sqrt(dot_product(diff,diff))<1.d-10) then
		ifind = isptpos
		return
	endif
enddo

End

!#############################################################################################


!The face file contains three header lines followed by one triangle per line. 
!The first header line provides a comment and the file name of the sphere set. 
!The second header line holds comments about the content of the third line. 
!The third header line provides the number of triangles, the number of spheres in the set, 
!the triangulation density and the probe sphere radius. 

!The first three numbers are (1 based) vertex indices. 

!The next field can be: 
!1 for a triangle in a toric reen trant face, 
!2 for a triangle in a spheric reentrant face and 
!3 for a triangle in a contact face. 

!The last # on the line is the (1 based) face number in the analytical description of the solvent excluded surface. 
!These values are written in the following format ``%6d %6d %6d %2d %6d''.

!The vertex file contains three header lines (similar to the header in the .face file) 
!followed by one vertex per line and provides the coordinates (x,y,z) and the normals (nx,ny,nz) 
!followed by the number of the face (in the analytical description of the solvent excluded surface) 
!to which the vertex belongs. The vertices of the analytical surface have a value 0 in that field 
!and the vertices lying on edges of this surface have negative values. 
!The next field holds the (1 based) index of the closest sphere. 
!The next field is 
!1 for vertices which belong to toric reentrant faces (including ver tices of the analytical surface), 
!2 for vertices inside reentrant faces and 
!3 for vertices inside contact faces. 
!These values are written in the following format ``%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d''.






!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! 
!            POSTPROCESSSING OF MSMS
!            PEIJUN LI, PURDUE UNIVERSITY
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! postprocessing of msms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

MODULE pmsms_procedures

  IMPLICIT NONE

  INTEGER, PARAMETER :: nchar = 1

  INTEGER nnvert
  INTEGER mmvert
  INTEGER jjvert
  INTEGER nnface
  INTEGER mmface
      
  INTEGER, ALLOCATABLE, DIMENSION (:) :: activeface
  INTEGER, ALLOCATABLE, DIMENSION (:) :: inactver
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: label

!  REAL ( KIND = 8 ), ALLOCATABLE, DIMENSION (:) :: q
  REAL ( KIND = 8 ), ALLOCATABLE, DIMENSION (:,:) :: pvert
  REAL ( KIND = 8 ), ALLOCATABLE, DIMENSION (:,:) :: pnorm
!  REAL ( KIND = 8 ), AlLOCATABLE, DIMENSION (:,:) :: pchar
  
END MODULE pmsms_procedures

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

!PROGRAM pmsms

!  USE pmsms_procedures
!  IMPLICIT NONE

!  CALL pmsms_setup
!  CALL check
!  CALL move
!  CALL output
!  CALL pmsms_cleanup

!  STOP

!END PROGRAM pmsms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE pmsms_setup
  USE COMDATA
  USE MOLECULE
  USE pmsms_procedures
  IMPLICIT NONE

  INTEGER i

  !OPEN (unit = 81, file = 'sphere.vert' )
  !READ (81,*) 
  !READ (81,*) 
  !READ (81,*) nvert

  nnvert=nspt
  nnface=nface
  ALLOCATE ( pvert(3,nnvert) ) 
  ALLOCATE ( pnorm(3,nnvert) ) 
 
  pvert=SPTPOS
  pnorm=SPTNRM

  !DO i = 1, nvert
  !   READ (81,*) pvert(1,i), pvert(2,i), pvert(3,i), pnorm(1,i), pnorm(2,i), pnorm(3,i)
  !END DO
  !CLOSE (81)

  !OPEN (unit = 82, file = 'sphere.face' )
  !READ (82,*) 
  !READ (82,*)
  !READ (82,*) nface

  ALLOCATE ( activeface(nnface) )
  ALLOCATE ( inactver(nnvert) )
  ALLOCATE ( label(3,nnface) )
  
  label=NVERT
  print *,pvert(1,1),pnorm(1,1),label(1,1)  
  !DO i = 1, nface
  !   READ (82,*) label(1,i), label(2,i), label(3,i)
  !END DO
  !CLOSE (82)

!  OPEN(unit = 83, file = 'sphere.char' )

!  ALLOCATE ( pchar(3,nchar) )
!  ALLOCATE ( q(nchar) ) 

!  DO i = 1, nchar
!     READ (83,*) pchar(1,i), pchar(2,i), pchar(3,i), q(i)
!  END DO
!  CLOSE (83) 

  RETURN

END SUBROUTINE pmsms_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE check

  USE pmsms_procedures
  IMPLICIT NONE

  INTEGER imin
  INTEGER imax
  INTEGER i
  INTEGER i1
  INTEGER i2
  INTEGER i3
  INTEGER ip
  INTEGER j
  INTEGER j1
  INTEGER j2
  INTEGER j3
  INTEGER jp

  REAL ( KIND = 8 ), PARAMETER :: eps = 1.0E-06

  REAL ( KIND = 8 ), DIMENSION (3) :: x
  REAL ( KIND = 8 ), DIMENSION (3) :: y
  REAL ( KIND = 8 ), DIMENSION (3) :: z
  REAL ( KIND = 8 ), DIMENSION (3) :: d
  REAL ( KIND = 8 ), DIMENSION (3) :: p
  REAL ( KIND = 8 ), DIMENSION (3) :: xx
  REAL ( KIND = 8 ), DIMENSION (3) :: yy
  REAL ( KIND = 8 ), DIMENSION (3) :: zz

  ip = 0
  jjvert = 0

  DO i = 1, nnface
  
     jp = 0

     x = pvert(1,label(:,i))
     y = pvert(2,label(:,i))
     z = pvert(3,label(:,i))
    !!!!!!!!!!!!!!!!!!!! 
     i1=label(1,i)
     i2=label(2,i)
     i3=label(3,i)
    !!!!!!!!!!!!!!!!!!!!
     xx(1) = x(2) - x(1)
     xx(2) = x(3) - x(1)
     xx(3) = x(3) - x(2)

     yy(1) = y(2) - y(1)
     yy(2) = y(3) - y(1)
     yy(3) = y(3) - y(2)

     zz(1) = z(2) - z(1)
     zz(2) = z(3) - z(1)
     zz(3) = z(3) - z(2)

     d = SQRT ( xx * xx + yy * yy + zz * zz ) ! The distance          
     p = SQRT ( x * x + y * y + z * z )
     d = d / p
!     print *,i,real(d)
!     print *,x
!     print *,y
!     print *,z
!     print *,xx
!     print *,yy
!     print *,zz
     IF ( d(3) .LE. eps ) THEN
        jp = 1
        IF ( i2 .LT. i3 ) THEN
           imin = i2
           imax = i3
           CALL insert (i3)
        ELSEIF ( i3 .LT. i2 ) THEN
           imin = i3
           imax = i2
           CALL insert (i2)
        ELSE
           imin = 0
           imax = 0
        END IF
        CALL update ( imin, imax )
     END IF

     IF ( d(2) .LE. eps ) THEN
        jp = 1
        IF ( i1 .LT. i3 ) THEN
           imin = i1
           imax = i3
           CALL insert (i3)
        ELSEIF ( i3 .LT. i1 ) THEN
           imin=i3
           imax=i1
           CALL insert (i1)
        ELSE
           imin=0
           imax=0
        END IF
        CALL update ( imin, imax )
     END IF

     IF ( d(1) .LE. eps ) THEN
        jp = 1
        IF ( i1 .LT. i2 ) THEN
           imin = i1
           imax = i2
           CALL insert (i2)
        ELSEIF ( i2 .LT. i1 ) THEN
           imin = i2
           imax = i1
           CALL insert (i1)
        ELSE
           imin=0
           imax=0
        END IF
        CALL update ( imin, imax )
     END IF

     IF ( jp .EQ. 0 ) THEN
        ip = ip + 1
        activeface(ip) = i     
     END IF

  END DO

  mmface = ip
  mmvert = nnvert - jjvert

  IF ( mmface .EQ. nnface ) THEN
     WRITE(*,*) ''
     WRITE(*,*) 'correct triangulation!'
     WRITE(*,*) ''
     !CALL writeup
     return
  END IF

  CALL sort

  RETURN

END SUBROUTINE check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE insert (imax)

  USE pmsms_procedures
  IMPLICIT NONE

  INTEGER imax

  INTEGER :: i

  DO i = 1, jjvert
     IF ( inactver(i) .EQ. imax ) THEN 
         print *,'inactver(i) equals imax' 
         STOP
     ENDIF
  END DO

  jjvert = jjvert + 1
  inactver(jjvert) = imax

  RETURN
  
END SUBROUTINE insert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE sort

  USE pmsms_procedures
  IMPLICIT NONE
   
  INTEGER i
  INTEGER j
  INTEGER k

  DO j = 2, jjvert
     k = inactver(j)
     DO i = j - 1, 1, -1
        IF ( inactver(i) .GT. k ) GO TO 30
        inactver(i+1) = inactver(i)
     END DO
     i = 0
30   inactver(i+1) = k
   END DO
   
   RETURN

END SUBROUTINE sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE update(imin,imax)

  USE pmsms_procedures
  IMPLICIT NONE

  INTEGER imin
  INTEGER imax

  INTEGER i
  INTEGER j

  IF ( imin .EQ. imax ) THEN
       print *,'imin=imax inside update'
  !     STOP
  ENDIF

  DO i = 1, nnface
     DO j = 1, 3
        IF ( label(j,i) .EQ. imax ) label(j,i) = imin
     END DO
  END DO
   
 RETURN
      
END SUBROUTINE update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE move

  USE pmsms_procedures
  IMPLICIT NONE

  INTEGER i
  INTEGER j
  INTEGER k
  INTEGER ip
  INTEGER jp
  
  DO i = 1, jjvert
     ip = inactver(i)
     DO j = 1, mmface
        jp = activeface(j)
        DO k = 1, 3
           IF ( label(k,jp) .GT. ip ) label(k,jp) = label(k,jp) - 1
        END DO
     END DO
  END DO

  DO i = 1, jjvert
     ip = inactver(i)
     DO j = ip, nnvert - i
        pvert(:,j) = pvert(:,j + 1)
        pnorm(:,j) = pnorm(:,j + 1)
     END DO
  END DO

  RETURN

END SUBROUTINE move

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE output

  USE pmsms_procedures
  IMPLICIT NONE

  WRITE(*,*) ''
  WRITE(*,*) 'defect triangulation!'
  WRITE(*,*) ''
  WRITE(*,*) 'new number of vertex:', mmvert
  WRITE(*,*) ''
  WRITE(*,*) 'new number of face:', mmface
  
!  CALL writeup

  RETURN

END SUBROUTINE output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

!SUBROUTINE writeup

!  USE pmsms_procedures
!  IMPLICIT NONE

!  INTEGER i
!  INTEGER j

!  OPEN ( unit = 1, file = 'sphere.cons' )
!  WRITE(1,*) mvert, mface, nchar
!  CLOSE (1)

!  OPEN ( unit = 2, file = 'sphere.vert' )
!  DO i = 1, mvert
!     WRITE(2,1) pvert(1,i), pvert(2,i), pvert(3,i)
!  END DO
!  CLOSE (2)

!  OPEN ( unit = 3, file = 'sphere.norm' )
!  DO i = 1, mvert
!     WRITE(3,1) pnorm(1,i), pnorm(2,i), pnorm(3,i) 
!  END DO
!  CLOSE (3)

!  OPEN ( unit = 4, file = 'sphere.face' )
!  DO i = 1, mface
!     j = activeface(i)
!     WRITE(4,*) label(1,j), label(2,j), label(3,j)
!  END DO
!  CLOSE(4)

!  OPEN ( unit = 5, file = 'sphere.char' )
!  DO i = 1, nchar
!     WRITE(5,2) pchar(1,i), pchar(2,i), pchar(3,i), q(i)
!  END DO
!  CLOSE(5) 
     
!1 FORMAT(3(F12.6,1x))
!2 FORMAT(4(F12.6,1x))

!  RETURN

!END SUBROUTINE writeup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

SUBROUTINE pmsms_cleanup

  USE pmsms_procedures
  IMPLICIT NONE

  DEALLOCATE (activeface)
  DEALLOCATE(inactver)
  DEALLOCATE(label)
  DEALLOCATE(pvert)
  DEALLOCATE(pnorm)
!  DEALLOCATE(pchar)
!  DEALLOCATE(q)


  RETURN

END SUBROUTINE pmsms_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80







     


