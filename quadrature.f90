!-------------------------------------------------------------------------------
!This subroutine computes the integral of the function on every element
!input uv(2,10) contains the parameters u and v for the 10 quadrature points
subroutine comp_int(uv)
use molecule
use comdata
implicit double precision(a-h,o-z)
real*8 v(3,3),r(3,3),uv(2,10),x10(3,10),v10(3,10),rs(2),xi(3),xj(3),vi(3),vj(3),u,xij(3),vij(3)
real*8 w9(9),rs9(2,9),Xrs(3),dXdr(3),dXds(3),ajacob(3)
real*8 w1(1),rs1(2,1),w3(3),rs3(2,3),w4(4),rs4(2,4),w7(7),rs7(2,7),xx(3),vcent(3)
integer idx(3)

real*8, allocatable :: wx(:),rsx(:,:),SS(:),DSR(:),DSS(:)  


NGR=4	! # of Pts for Gaussian Quadrature (5 places to change: 1:NGR, 2:w-, 3:rs-, 4:function definition and 5:reference)
S_total=0.d0
S_area=0.d0		! Calculate the surface integral via (centroid + area of triangles)
allocate(wx(NGR),rsx(2,NGR),SS(10),DSR(10),DSS(10)) 
rsx=0.d0;	wx=0.d0;	SS=0.d0;	DSR=0.d0;	DSS=0.d0;

call Gauss_Radau_Quad(wx,rsx)													! Quad change

!print *,real(w4),real(rs4)
do j=1,nface ! for each triangle
	idx=nvert(1:3,j) ! vertices index of the specific triangle
	vcent=0.d0
	do k=1,3
		v(1:3,k)=sptnrm(1:3,idx(k))	! normal direction
		r(1:3,k)=sptpos(1:3,idx(k))	! position
		vcent=vcent+1.d0/3.d0*sptpos(1:3,idx(k))	!centriod
	enddo
	
	aa=sqrt(dot_product(r(:,1)-r(:,2),r(:,1)-r(:,2)))
	bb=sqrt(dot_product(r(:,1)-r(:,3),r(:,1)-r(:,3)))
	cc=sqrt(dot_product(r(:,2)-r(:,3),r(:,2)-r(:,3)))
	area=triangle_area(aa,bb,cc)	
	
	! Find the coordinates and normal of the 10 picked point on the jth curved elememt
	
	! Three known points indexed as 1,4 and 9
	x10(:,1)=r(:,1);	v10(:,1)=v(:,1)
	x10(:,4)=r(:,2);	v10(:,4)=v(:,2)
	x10(:,9)=r(:,3);	v10(:,9)=v(:,3)
	
	! 1st side: ponts (1), 2 , 3 and (4)
	xi=r(1:3,1);		vi=v(1:3,1)
	xj=r(1:3,2);		vj=v(1:3,2)
	
	u=uv(1,2)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,2)=xij;		v10(:,2)=vij

	u=uv(1,3)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,3)=xij;		v10(:,3)=vij
	
	
	! 2nd side: points (1), 5 ,7 and (9)
	xi=r(1:3,1);		vi=v(1:3,1)
	xj=r(1:3,3);		vj=v(1:3,3)
	
	u=uv(1,5)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,5)=xij;		v10(:,5)=vij
	
	u=uv(1,7)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,7)=xij;		v10(:,7)=vij
	
	! 3rd side: points (4), 6, 8 and (9)
	xi=r(1:3,2); vi=v(1:3,2)
	xj=r(1:3,3); vj=v(1:3,3)
	
	u=uv(2,6)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,6)=xij;		v10(:,6)=vij

	u=uv(2,8)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,8)=xij;		v10(:,8)=vij

	! center: point 10
	xi=x10(:,7);		vi=v10(:,7)
	xj=x10(:,3); 		vj=v10(:,3)
	
	u=uv(2,10)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,10)=xij;		v10(:,10)=vij

	
	!do iprt=1,10
	!	write(10,'(6f12.6)') x10(:,iprt), v10(:,iprt) 
	!enddo
	
	!now integration
	S_ele=0.d0
	do l=1,NGR
		rs=rsx(:,l)														! Quad change

		Xrs=0.d0
		dXdr=0.d0
		dXds=0.d0
		call tenpts_intp(rs(1),rs(2),SS,DSR,DSS)

		do i=1,3
			Xrs(i)=dot_product(x10(i,:),SS)
			dXdr(i)=dot_product(x10(i,:),DSR)
			dXds(i)=dot_product(x10(i,:),DSS)
		enddo

		call cross_product(dXdr,dXds, ajacob)

		S_quad=wx(l)*func(Xrs)*sqrt(dot_product(ajacob,ajacob))			! Quad change

		S_ele=S_ele+S_quad

	enddo


	deallocate(wx,rsx,SS,DSR,DSS) 
	!print *,j,real(area),real(S_ele),real(abs(area-S_ele))
	S_total=S_total+S_ele

	S_area=S_area+area*func(Vcent)
	!print *,j,real(S_area),real(func(vcent)),real(area)
enddo
xx=(/1.d0,2.d0,3.d0/)
Print *,'S_total= ', S_total, 'error= ',abs(S_total-true_soln(xx))
!Print *,'S_area= ', S_area, 'error= ',abs(S_area-true_soln(xx))
End
!----------------------------------------------------------------------------------------------------

! alphi and beta are functions needed to calculate the trajectory
function alphi(xi,xj,vi,vj)
! xi,xj are poistions and vi,vj are normal vectors
implicit double precision(a-h,o-z)
real*8 xi(3),xj(3),vi(3),vj(3)

aa=dot_product((xj-xi),(vi*dot_product(vj,vj)+0.5d0*vj*(dot_product(vi,vj)))) 
bb=(2.d0/3.d0)*dot_product(vi,vi)*dot_product(vj,vj)-(1.d0/6.d0)*dot_product(vi,vj)
alphi=aa/bb

End 

!----------------------------------------------------------------------------------------------------
function beta(xi,xj,vi,vj)
! xi,xj are poistions and vi,vj are normal vectors
implicit double precision(a-h,o-z)
real*8 xi(3),xj(3),vi(3),vj(3)

aa=-dot_product((xj-xi),vj)-(1.d0/3.d0)*alphi(xi,xj,vi,vj)*dot_product(vi,vj)
bb=(2.d0/3.d0)*dot_product(vi,vi)
beta=aa/bb
End 

!----------------------------------------------------------------------------------------------
subroutine coor_trans(uv)
implicit double precision(a-h,o-z) 
real*8 rs(2,10),uv(2,10)

onethird=1.d0/3.d0
twothird=2.d0/3.d0

! r and s
rs(:,1)=(/0.d0,		0.d0/)
rs(:,2)=(/onethird, 0.d0/)
rs(:,3)=(/twothird, 0.d0/)
rs(:,4)=(/1.d0,		0.d0/)
rs(:,5)=(/0.d0,		onethird/)
rs(:,6)=(/twothird, onethird/)
rs(:,7)=(/0.d0,		twothird/)
rs(:,8)=(/onethird, twothird/)
rs(:,9)=(/0.d0,		1.d0/)
rs(:,10)=(/onethird,onethird/)

! calculate u and v by 1: u=r+s and v=s/(r+s) if r+s <>0; 2: u=0 and v =0 if r=s=0 
do i=1,10
	if (sum(abs(rs(:,i)))<1.d-10) then
		uv(:,i)=0.d0	
	else
		uv(1,i)=sum(rs(:,i))
		uv(2,i)=rs(2,i)/uv(1,i)
	endif
	!write(*,*) i,uv(:,i)
enddo
End

!------------------------------------------------------------------------------------------------
subroutine comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
implicit double precision(a-h,o-z)
real*8 xi(3),xj(3),vi(3),vj(3),u,xij(3),vij(3), c0(3), c1(3), c2(3), c3(3), dt(3), dtt(3)

af=alphi(xi,xj,vi,vj)
bt=beta(xi,xj,vi,vj)

c0=xi
c2=af*vi
c3=1/3.d0*(bt*vj-c2)
c1=xj-c0-c2-c3

! output the position xij corresponding to parameter u
xij=c0+c1*u+c2*u**2+c3*u**3


dt=c1+2.d0*c2*u+3.d0*c3*u**2
dtt=2.d0*c2+6.d0*c3*u

t1=sqrt(dot_product(dt,dt))
t2=dot_product(dt,dtt)

! output the normal conditions vij corresponding to parameter u
vij=dtt-t2/t1**2*dt
vij=vij/t1
vij=-vij/sqrt(dot_product(vij,vij))        
!print *,u,real(vij)
End

!--------------------------------------------------------------------------------------
subroutine tenpts_intp(r,s,SS,DSR,DSS)
implicit double precision(a-h,o-z)
real*8 SS(1:10), T(9), DTR(9), DTS(9), DSR(10), DSS(10)

T(1)=1.d0-r-s;			!L1		
T(2)=r;					!L2
T(3)=s;					!L3	

T(4)=3.d0*T(1)-1.d0		!3L1-1
T(5)=3.d0*T(2)-1.d0		!3L2-1
T(6)=3.d0*T(3)-1.d0		!3L3-1

T(7)=3.d0*T(1)-2.d0		!3L1-2
T(8)=3.d0*T(2)-2.d0		!3L2-2
T(9)=3.d0*T(3)-2.d0		!3L3-2


DTR(1)=-1.d0			!dL1/dr
DTR(2)=1.d0				!dL2/dr
DTR(3)=0.d0				!dL3/dr

DTR(4:6)=DTR(1:3)*3.d0	!dT4/dr-dT6/dr
DTR(7:9)=DTR(4:6)		!dT7/dr-dT9/dr


DTS(1)=-1.d0			!dT1/ds
DTS(2)=0.d0				!dT2/ds
DTS(3)=1.d0				!dT3/ds

DTS(4:6)=DTS(1:3)*3.d0	!dT4/ds-dT6/ds
DTS(7:9)=DTS(4:6)		!dT7/ds-dT9/ds


SS(1)=0.5d0*T(1)*T(4)*T(7)	!N1-N10		
SS(2)=4.5d0*T(1)*T(2)*T(4)
SS(3)=4.5d0*T(1)*T(2)*T(5)

SS(4)=0.5d0*T(2)*T(5)*T(8)
SS(5)=4.5d0*T(1)*T(3)*T(4)
SS(6)=4.5d0*T(2)*T(3)*T(5)

SS(7)=4.5d0*T(1)*T(3)*T(6)
SS(8)=4.5d0*T(2)*T(3)*T(6)
SS(9)=0.5d0*T(3)*T(6)*T(9)
SS(10)=27.d0*T(1)*T(2)*T(3)


DSR(1)=0.5d0*(DTR(1)*T(4)*T(7)+T(1)*DTR(4)*T(7)+T(1)*T(4)*DTR(7)) !dN1/dr-dN10/dr
DSR(2)=4.5d0*(DTR(1)*T(2)*T(4)+T(1)*DTR(2)*T(4)+T(1)*T(2)*DTR(4))
DSR(3)=4.5D0*(DTR(1)*T(2)*T(5)+T(1)*DTR(2)*T(5)+T(1)*T(2)*DTR(5))

DSR(4)=0.5D0*(DTR(2)*T(5)*T(8)+T(2)*DTR(5)*T(8)+T(2)*T(5)*DTR(8))
DSR(5)=4.5D0*(DTR(1)*T(3)*T(4)+T(1)*DTR(3)*T(4)+T(1)*T(3)*DTR(4))
DSR(6)=4.5D0*(DTR(2)*T(3)*T(5)+T(2)*DTR(3)*T(5)+T(2)*T(3)*DTR(5))

DSR(7)=4.5D0*(DTR(1)*T(3)*T(6)+T(1)*DTR(3)*T(6)+T(1)*T(3)*DTR(6))
DSR(8)=4.5D0*(DTR(2)*T(3)*T(6)+T(2)*DTR(3)*T(6)+T(2)*T(3)*DTR(6))
DSR(9)=0.5D0*(DTR(3)*T(6)*T(9)+T(3)*DTR(6)*T(9)+T(3)*T(6)*DTR(9))
DSR(10)=27.D0*(DTR(1)*T(2)*T(3)+T(1)*DTR(2)*T(3)+T(1)*T(2)*DTR(3)) 
	

DSS(1)=0.5d0*(DTS(1)*T(4)*T(7)+T(1)*DTS(4)*T(7)+T(1)*T(4)*DTS(7)) !dN1/ds-dN10/ds
DSS(2)=4.5d0*(DTS(1)*T(2)*T(4)+T(1)*DTS(2)*T(4)+T(1)*T(2)*DTS(4))
DSS(3)=4.5D0*(DTS(1)*T(2)*T(5)+T(1)*DTS(2)*T(5)+T(1)*T(2)*DTS(5))

DSS(4)=0.5D0*(DTS(2)*T(5)*T(8)+T(2)*DTS(5)*T(8)+T(2)*T(5)*DTS(8))
DSS(5)=4.5D0*(DTS(1)*T(3)*T(4)+T(1)*DTS(3)*T(4)+T(1)*T(3)*DTS(4))
DSS(6)=4.5D0*(DTS(2)*T(3)*T(5)+T(2)*DTS(3)*T(5)+T(2)*T(3)*DTS(5))

DSS(7)=4.5D0*(DTS(1)*T(3)*T(6)+T(1)*DTS(3)*T(6)+T(1)*T(3)*DTS(6))
DSS(8)=4.5D0*(DTS(2)*T(3)*T(6)+T(2)*DTS(3)*T(6)+T(2)*T(3)*DTS(6))
DSS(9)=0.5D0*(DTS(3)*T(6)*T(9)+T(3)*DTS(6)*T(9)+T(3)*T(6)*DTS(9))
DSS(10)=27.D0*(DTS(1)*T(2)*T(3)+T(1)*DTS(2)*T(3)+T(1)*T(2)*DTS(3)) 

End

!---------------------------------------------------------------------------------

subroutine Gauss_Legendre_Quad(w4,x4)
implicit double precision(a-h,o-z)
real*8 w1,x1,w2(2),x2(2),w3(3),x3(3),w4(4),x4(4),w5(5),x5(5),w8(8),x8(8),x6(6),w6(6)

x1=0.d0
w1=2.d0	

x2=(/-sqrt(1.d0/3.d0),sqrt(1.d0/3.d0)/)
w2=1.d0

x3=(/-sqrt(3.d0/5.d0),0.d0,sqrt(3.d0/5.d0)/)
w3=(/5.d0/9.d0, 8.d0/9.d0, 5.d0/9.d0/)

x4(1:2)=(/-sqrt((3.d0+2.d0*sqrt(6.d0/5.d0))/7.d0),-sqrt((3.d0-2.d0*sqrt(6.d0/5.d0))/7.d0)/)
x4(4)=-x4(1); x4(3)=-x4(2);
w4(1:4:3)=(18.d0-sqrt(30.d0))/36.d0
w4(2:3)=(18.d0+sqrt(30.d0))/36.d0

x5(3)=0.d0
w5(3)=128.d0/225.d0
x5(1)=-1.d0/3.d0*sqrt(5.d0+2.d0*sqrt(10.d0/7.d0))
x5(2)=-1.d0/3.d0*sqrt(5.d0-2.d0*sqrt(10.d0/7.d0))
x5(4)=-x5(2)
x5(5)=-x5(1)

w5(1:5:4)=(322.d0-13.d0*sqrt(70.d0))/900.d0
w5(2:4:2)=(322.d0+13.d0*sqrt(70.d0))/900.d0


x6=(/   -0.932469514203153, -0.661209386466264, 0.932469514203154, &
        0.661209386466263,  -0.238619186083197, 0.238619186083197/)
w6=(/   0.171324492379168,  0.360761573048139,  0.171324492379166, &
        0.360761573048139,  0.467913934572691,  0.467913934572691/)

x8=(/   -0.960289856497540, -0.796666477413625, 0.960289856497543,  0.796666477413622,   &
        -0.525532409916329, 0.525532409916331,  -0.183434642495650, 0.183434642495650 /)
w8=(/   0.101228536290368,  0.222381034453375,  0.101228536290360,  0.222381034453380,   &
        0.313706645877887,  0.313706645877887,  0.362683783378362,  0.362683783378362/)        
End
!------------------------------------------------------------------------------------------
subroutine Gauss_Radau_Quad(w4,rs4)
implicit double precision(a-h,o-z)
real*8 a4(4),w2(2),a9(9),w6(6),w8(8)
real*8 w4(4),rs4(2,4),w9(9),rs9(2,9),w16(16),rs16(2,16)
real*8 w3(3),rs3(2,3),w7(7),rs7(2,7),w1(1),rs1(2,1)

!n=1
w1=0.5d0
rs1=1.d0/3.d0

!n=3
w3=1.d0/6.d0
rs3(1,:)=(/0.5d0,0.5d0,0.d0/)
rs3(2,:)=(/0.d0,0.5d0,0.5d0/)

!n=4
a4=(/0.0750311102, 0.1785587283, 0.2800199155, 0.6663902460/)
w2=(/0.0909793091, 0.1590206909/)
w4=(/w2(1),w2(2),w2(1),w2(2)/)
rs4(1,:)=(/a4(3),a4(4),a4(1),a4(2)/)
rs4(2,:)=a4

!n=7
w7=(/	1.d0/40.d0, 1.d0/15.d0, 1.d0/40.d0, 1.d0/15.d0, 1.d0/40.d0, 1.d0/15.d0, 9.d0/40.d0	/)
rs7(1,:)=(/	0.d0, 0.5d0, 1.d0, 0.5d0, 0.d0, 0.d0,  1.d0/3.d0 /)
rs7(2,:)=(/	0.d0, 0.d0,  0.d0, 0.5d0, 1.d0, 0.5d0, 1.d0/3.d0 /)

!n=9
a9=(/	0.02393113229,	0.06655406786,	0.10271765483, &
		0.10617026910,	0.188409405913, 0.29526656780, &
		0.45570602025,	0.52397906774,	0.80869438567		/)

w6=(/	0.019396383304, 0.031034213285, 0.055814420490, &
		0.063678085097, 0.089303072783, 0.101884936154		/)

w9=(/w6(1),w6(4),w6(3),w6(2),w6(6),w6(5),w6(1),w6(4),w6(3)/)
rs9(1,:)=(/a9(5),a9(8),a9(9),a9(4),a9(6),a9(7),a9(1),a9(2),a9(3)/)
rs9(2,:)=a9

!n=16
w8=(/	0.005423225910,	0.010167259561,	0.022584049287, 0.023568368199, &
		0.035388067900,	0.042339724518,	0.044185088522, 0.066344216039/)

w16= (/	w8(1),w8(3),w8(5),w8(4),w8(2),w8(6),w8(8),w8(7), &
		w8(2),w8(6),w8(8),w8(7),w8(1),w8(3),w8(5),w8(4) /)


End

subroutine find10pts_pos_norm(v,r,x10,v10,uv) 
implicit none
real*8 v(3,3),r(3,3),x10(3,10),v10(3,10),uv(2,10)
real*8 xi(3),vi(3),xj(3),vj(3),xij(3),vij(3),u


! Find the coordinates and normal of the 10 picked point on the jth curved elememt
	
	! Three known points indexed as 1,4 and 9
	x10(:,1)=r(:,1);	v10(:,1)=v(:,1)
	x10(:,4)=r(:,2);	v10(:,4)=v(:,2)
	x10(:,9)=r(:,3);	v10(:,9)=v(:,3)
	
	! 1st side: ponts (1), 2 , 3 and (4)
	xi=r(1:3,1);		vi=v(1:3,1)
	xj=r(1:3,2);		vj=v(1:3,2)
	
	u=uv(1,2)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,2)=xij;		v10(:,2)=vij

	u=uv(1,3)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,3)=xij;		v10(:,3)=vij
	
	
	! 2nd side: points (1), 5 ,7 and (9)
	xi=r(1:3,1);		vi=v(1:3,1)
	xj=r(1:3,3);		vj=v(1:3,3)
	
	u=uv(1,5)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,5)=xij;		v10(:,5)=vij
	
	u=uv(1,7)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,7)=xij;		v10(:,7)=vij
	
	! 3rd side: points (4), 6, 8 and (9)
	xi=r(1:3,2); vi=v(1:3,2)
	xj=r(1:3,3); vj=v(1:3,3)
	
	u=uv(2,6)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,6)=xij;		v10(:,6)=vij

	u=uv(2,8)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,8)=xij;		v10(:,8)=vij

	! center: point 10
	xi=x10(:,7);		vi=v10(:,7)
	xj=x10(:,3); 		vj=v10(:,3)
	
	u=uv(2,10)
	call comp_pos_nrm(xi,xj,vi,vj,u,xij,vij)
	x10(:,10)=xij;		v10(:,10)=vij

end

!----------------------------------------------
subroutine cross_product(a,b,c)
implicit double precision(a-h,o-z)
real*8 a(3),b(3),c(3)

c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)

End

!-------------------------------------------------------------------
function func(x)
implicit double precision(a-h,o-z)
real*8 x(3), x0(3)
!Calculate the area only
!func=1.d0

!f(r)=1/sqrt(r-r0)
x0=(/1.d0,0.d0,0.d0/)
r=sqrt(dot_product(x-x0,x-x0))
if (abs(r)>1.d-10) then
	func=1.d0/r
else
	func=0.d0
endif
!f(r)=e^z  
!func=exp(x(3))

!f(r)=x^2*z+y^2*z
!func=(x(1)**2+x(2)**2)*x(3)
End


function true_soln(x)
implicit double precision(a-h,o-z)
real*8 x(3)
!area or f(r)=1/sqrt(r-r0)
true_soln=4.d0*acos(-1.d0)

!f(r)=e^z
!true_soln=2.d0*acos(-1.d0)*(exp(1.d0)-exp(-1.d0))

!f(r)=x^2*z+y^2*z
!true_soln=0.d0
End



