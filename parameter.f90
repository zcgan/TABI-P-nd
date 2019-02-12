! The fundermantal solution of poisson equation 
function G0(r,s)
implicit double precision(a-h,o-z)
real*8 r(3),s(3), G0
pi=acos(-1.d0)
rs=sqrt(dot_product(r-s,r-s))
G0=1.d0/(4.d0*pi*rs)

end function 


! The fundermantal solution of linearized PB equation
function Gk(r,s,kappa)
implicit double precision(a-h,o-z)
real*8 r(3),s(3),rs,kappa, Gk
pi=acos(-1.d0)
rs=sqrt(dot_product(r-s,r-s))
Gk=exp(-kappa*rs)/(4.d0*pi*rs)

end function

! find cos(theta) and cos(theta0)
! v is the unit vector
function cos_theta(v,r,s)
implicit double precision(a-h,o-z)
real*8 v(3),r(3),s(3), cos_theta
rs=sqrt(dot_product(r-s,r-s))
cos_theta=dot_product(v,r-s)/rs
end function


! G1, G2, G3 and G4 are the derivatives of the fundermental solution 
function G1(v0,r,s)
implicit double precision(a-h,o-z)
real*8 v0(3),r(3),s(3), G1
pi=acos(-1.d0)
rs=sqrt(dot_product(r-s,r-s))
G1=cos_theta(v0,r,s)/rs**2/4.d0/pi
end function

function G2(v0,r,s,kappa)
implicit double precision(a-h,o-z)
real*8 v0(3),r(3),s(3),kappa, G2
rs=sqrt(dot_product(r-s,r-s))
G2=(1.d0+kappa*rs)*exp(-kappa*rs)*G1(v0,r,s)
end function

function G3(v0,v,r,s)
implicit double precision(a-h,o-z)
real*8 v0(3),v(3),r(3),s(3), G3
pi=acos(-1.d0)
rs=sqrt(dot_product(r-s,r-s))
G3=(dot_product(v0,v)-3.d0*cos_theta(v0,r,s)*cos_theta(v,r,s))/rs**3/4.d0/pi
end function

function G4(v0,v,r,s,kappa)
implicit double precision(a-h,o-z)
real*8 v0(3),v(3),r(3),s(3),kappa, G4
pi=acos(-1.d0)
rs=sqrt(dot_product(r-s,r-s))
G4=(1.d0+kappa*rs)*exp(-kappa*rs)*G3(v0,v,r,s)-kappa**2*exp(-kappa*rs)*cos_theta(v0,r,s)*cos_theta(v,r,s)/rs/4.d0/pi
end function


! H1-H4 are L1-L4 in Juffer's Paper as coefficient in the integral equation
function H1(v,r,s,eps,kappa)
implicit double precision(a-h,o-z)
real*8 v(3),r(3),s(3),eps,kappa, H1
H1=eps*G2(v,r,s,kappa)-G1(v,r,s)
H1=-H1
!################
!H1=0.d0
!################
end function

function H2(r,s,kappa)
implicit double precision(a-h,o-z)
real*8 r(3),s(3),kappa, H2
H2=G0(r,s)-Gk(r,s,kappa)
!###############
!H2=0.d0
!H2=-Gk(r,s,kappa)
!###############
!print *,kappa,G0(r,s),Gk(r,s,kappa)

end function

function H3(v0,v,r,s,kappa)
implicit double precision(a-h,o-z)
real*8 v0(3),v(3),r(3),s(3),kappa, H3
H3=G4(v0,v,r,s,kappa)-G3(v0,v,r,s)
!################
!H3=0.d0
!################
end function

function H4(v0,r,s,kappa,eps)
implicit double precision(a-h,o-z)
real*8 v0(3),r(3),s(3),kappa,eps, H4
H4=G1(v0,r,s)-G2(v0,r,s,kappa)/eps
!###############
!H4=0.d0
!###############
end function
