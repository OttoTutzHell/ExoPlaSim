! ************************************************************************
! * module legmod - direct and indirect Legendre transformation routines *
! * E. Kirk 20-Feb-2009 tested for T21 - T341 resolutions 32 & 64 bit    *
! ************************************************************************

module legmod
! ************************
! * Legendre Polynomials *
! ************************
use pumamod, only:NTRU,NTP1,NCSP,NESP,NLON,NLPP,NLAT,NHOR,NLEV,gwd,sid,plavor,nfilter
use pumamod, only:ngptfilter, nspvfilter,landhoskn0,filterkappa,nfilterexp,nud,mypid,NROOT

real :: qi(NCSP,NLPP) ! P(m,n) = Associated Legendre Polynomials
real :: qj(NCSP,NLPP) ! Q(m,n) = Used for d/d(mu)
real :: qc(NCSP,NLPP) ! P(m,n) * gwd              used in fc2sp
real :: qe(NCSP,NLPP) ! Q(mn,) * gwd / cos2       used in mktend
real :: qm(NCSP,NLPP) ! P(m,n) * gwd / cos2 * m   used in mktend
real :: qq(NCSP,NLPP) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
real :: qu(NCSP,NLPP) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
real :: qv(NCSP,NLPP) ! Q(m,n) / (n*(n+1))        used in dv2uv
real :: skgpsp(NTP1) ! Physics filter for GP -> SP
real :: skspgp(NTP1) ! Physics filter for SP -> GP

end module legmod

! =================
! SUBROUTINE LEGINI
! =================

subroutine legini
use legmod
implicit none

integer :: jlat ! Latitude
integer :: lm
integer :: m
integer :: n

real (kind=8) :: amsq
real (kind=8) :: z1
real (kind=8) :: z2
real (kind=8) :: z3
real (kind=8) :: f1m
real (kind=8) :: f2m
real (kind=8) :: znn1
real (kind=8) :: zsin    ! sin
real (kind=8) :: zcsq    ! cos2
real (kind=8) :: zgwd    ! gw
real (kind=8) :: zgwdcsq ! gw / cos2

real (kind=8) :: zpli(NCSP)
real (kind=8) :: zpld(NCSP)

do jlat = 1 , NLPP

! set p(0,0) and p(0,1)

   zgwd    = gwd(jlat)            ! gaussian weight - from inigau
   zsin    = sid(jlat)            ! sin(phi) - from inigau
   zcsq    = 1.0_8 - zsin * zsin  ! cos(phi) squared
   zgwdcsq = zgwd / zcsq          ! weight / cos squared
   f1m     = sqrt(1.5_8)
   zpli(1) = sqrt(0.5_8)
   zpli(2) = f1m * zsin
   zpld(1) = 0.0
   lm      = 2

! loop over wavenumbers

   do m = 0 , NTRU
      if (m > 0) then
         lm  = lm + 1
         f2m = -f1m * sqrt(zcsq / (m+m))
         f1m =  f2m * sqrt(m+m + 3.0_8)
         zpli(lm) = f2m
         if (lm < NCSP) then
            lm = lm + 1
            zpli(lm  ) =       f1m * zsin
            zpld(lm-1) =  -m * f2m * zsin
         endif ! (lm < NCSP)
      endif ! (m > 0)

      amsq = m * m

      do n = m+2 , NTRU
         lm = lm + 1
         z1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
         z2 = zsin * zpli(lm-1) - z1 * zpli(lm-2)
         zpli(lm  ) = z2 * sqrt((4*n*n-1) / (n*n-amsq))
         zpld(lm-1) = (1-n) * z2 + n * z1 * zpli(lm-2)
      enddo ! n

      if (lm < NCSP) then ! mode (m,NTRU)
         z3 = sqrt((NTRU*NTRU-amsq) / (4*NTRU*NTRU-1))
         zpld(lm)=-NTRU*zsin*zpli(lm) + (NTRU+NTRU+1)*zpli(lm-1)*z3
      else                ! mode (NTRU,NTRU)
         zpld(lm)=-NTRU*zsin*zpli(lm)
      endif
   enddo ! m

   lm = 0
   do m = 0 , NTRU
      do n = m , NTRU
           lm = lm + 1
           znn1 = 0.0
           if (n > 0) znn1 = 1.0_8 / (n*(n+1))
           qi(lm,jlat) = zpli(lm)
           qj(lm,jlat) = zpld(lm)
           qc(lm,jlat) = zpli(lm) * zgwd
           qu(lm,jlat) = zpli(lm) * znn1 * m
           qv(lm,jlat) = zpld(lm) * znn1
           qe(lm,jlat) = zpld(lm) * zgwdcsq
           qq(lm,jlat) = zpli(lm) * zgwdcsq * n * (n+1) * 0.5_8
           qm(lm,jlat) = zpli(lm) * zgwdcsq * m
      enddo ! n
   enddo ! m
   

   
   
enddo ! jlat

skgpsp(:) = 1.0
skspgp(:) = 1.0

   
do n=1,NTP1
   if (nfilter .eq. 0) then
     skgpsp(n) = 1.0    ! Filter for conversion from gp->sp
     skspgp(n) = 1.0    ! Filter for conversion from sp->gp
     
   else if (nfilter .eq. 1) then !Cesaro filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) +  ngptfilter*(1.0 - real(n)/(NTP1))
     skspgp(n) = (1-nspvfilter)*skspgp(n) +  nspvfilter*(1.0 - real(n)/(NTP1))
   
   else if (nfilter .eq. 2) then !Exponential filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) + &
&                 ngptfilter*(exp(-filterkappa*(real(n)/NTRU)**nfilterexp))
     skspgp(n) = (1-nspvfilter)*skspgp(n) + &
&                 nspvfilter*(exp(-filterkappa*(real(n)/NTRU)**nfilterexp))
   
   else if (nfilter .eq. 3) then !Lander-Hoskins physics filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) + ngptfilter*(exp(-(real(n)*real(n+1)/(landhoskn0*(landhoskn0+1)))**2))
     skspgp(n) = (1-nspvfilter)*skspgp(n) + nspvfilter*(exp(-(real(n)*real(n+1)/(landhoskn0*(landhoskn0+1)))**2))
   
   else if (nfilter .eq. 4) then !Riesz-2 filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) + ngptfilter*((1.0 - real(n)/(NTP1))**2)
     skspgp(n) = (1-nspvfilter)*skspgp(n) + nspvfilter*((1.0 - real(n)/(NTP1))**2)
   
   else
     skgpsp(n) = 1.0
     skspgp(n) = 1.0
   endif
   
enddo

! Fold physics filters into the Legendre matrices once so they are not
! re-multiplied on every inner-loop iteration of every timestep.
do jlat = 1 , NLPP
   lm = 0
   do m = 1 , NTP1
      do n = m , NTP1
         lm = lm + 1
         qi(lm,jlat) = qi(lm,jlat) * skspgp(n)
         qj(lm,jlat) = qj(lm,jlat) * skspgp(n)
         qu(lm,jlat) = qu(lm,jlat) * skspgp(n)
         qv(lm,jlat) = qv(lm,jlat) * skspgp(n)
         qc(lm,jlat) = qc(lm,jlat) * skgpsp(n)
         qe(lm,jlat) = qe(lm,jlat) * skgpsp(n)
         qm(lm,jlat) = qm(lm,jlat) * skgpsp(n)
         qq(lm,jlat) = qq(lm,jlat) * skgpsp(n)
      enddo ! n
   enddo ! m
enddo ! jlat

if (mypid==NROOT) then
   write(nud,*)"*************************************************"
   write(nud,*)"* GP->SP Physics Filter Coefficients:           *"
   do n=1,NTP1
      write(nud,*)"* n=",n,":  ",skgpsp(n)," *"
   enddo
   write(nud,*)"* SP->GP Physics Filter Coefficients:           *"
   do n=1,NTP1
      write(nud,*)"* n=",n,":  ",skspgp(n)," *"
   enddo
   write(nud,*)"*************************************************"
endif

return
end


! ================
! SUBROUTINE FC2SP
! ================

subroutine fc2sp(fc,sp)
use legmod
implicit none
real, intent(in ) :: fc(2,NLON/2,NLPP)
real, intent(out) :: sp(2,NESP/2)

integer :: k ! Index for southern latitude
integer :: l ! Index for latitude
integer :: m ! Index for zonal wavenumber
integer :: w, w0, w1 ! Index for spherical harmonic
real :: f1, f2, f1s, f2s, f1a, f2a

sp(:,:) = 0.0

if (NLPP < NLAT) then  ! Universal (parallel executable) version
!----------------------------------------------------------------------
  do l = 1 , NLPP
    w0 = 1
    do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      f1 = fc(1,m,l)
      f2 = fc(2,m,l)
      do w = w0 , w1
        sp(1,w) = sp(1,w) + qc(w,l) * f1
        sp(2,w) = sp(2,w) + qc(w,l) * f2
      enddo ! w
      w0 = w1 + 1
    enddo ! m
  enddo ! l
else                   ! Single CPU version (symmetry conserving)
!----------------------------------------------------------------------
  do l = 1 , NLAT/2
    k = NLAT + 1 - l
    w0 = 1
    do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      f1s = fc(1,m,l) + fc(1,m,k)
      f2s = fc(2,m,l) + fc(2,m,k)
      f1a = fc(1,m,l) - fc(1,m,k)
      f2a = fc(2,m,l) - fc(2,m,k)
      do w = w0 , w1 , 2
        sp(1,w) = sp(1,w) + qc(w,l) * f1s
        sp(2,w) = sp(2,w) + qc(w,l) * f2s
      enddo
      do w = w0 + 1 , w1 , 2
        sp(1,w) = sp(1,w) + qc(w,l) * f1a
        sp(2,w) = sp(2,w) + qc(w,l) * f2a
      enddo
      w0 = w1 + 1
    enddo ! m
  enddo ! l
!----------------------------------------------------------------------
endif ! parallel ?
return
end


! ================
! SUBROUTINE SP2FC
! ================

subroutine sp2fc(sp,fc) ! Spectral to Fourier
use legmod
implicit none

real :: sp(2,NCSP)        ! Coefficients of spherical harmonics
real :: fc(2,NLON/2,NLPP) ! Fourier coefficients

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w, w0, w1 ! Loop index for spectral mode
real :: f1, f2

fc(:,:,:) = 0.0

do l = 1 , NLPP
   w0 = 1
   do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      f1 = 0.0
      f2 = 0.0
      do w = w0 , w1
         f1 = f1 + qi(w,l) * sp(1,w)
         f2 = f2 + qi(w,l) * sp(2,w)
      enddo ! w
      fc(1,m,l) = f1
      fc(2,m,l) = f2
      w0 = w1 + 1
   enddo ! m
enddo ! l
return
end


! ===================
! SUBROUTINE SP2FCDMU
! ===================

subroutine sp2fcdmu(sp,fc) ! Spectral to Fourier d/dmu
use legmod
implicit none

real :: sp(2,NCSP)        ! Coefficients of spherical harmonics
real :: fc(2,NLON/2,NLPP) ! Fourier coefficients

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w, w0, w1 ! Loop index for spectral mode
real :: f1, f2

fc(:,:,:) = 0.0

do l = 1 , NLPP
   w0 = 1
   do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      f1 = 0.0
      f2 = 0.0
      do w = w0 , w1
         f1 = f1 + qj(w,l) * sp(1,w)
         f2 = f2 + qj(w,l) * sp(2,w)
      enddo ! w
      fc(1,m,l) = f1
      fc(2,m,l) = f2
      w0 = w1 + 1
   enddo ! m
enddo ! l
return
end


! ================
! SUBROUTINE SP3FC
! ================

subroutine sp3fc
use pumamod, only:NLEV,sd,st,sz,gd,gt,gz
implicit none
integer :: v ! Loop index for level

do v = 1 , NLEV
   call sp2fc(sd(1,v),gd(1,v))
   call sp2fc(st(1,v),gt(1,v))
   call sp2fc(sz(1,v),gz(1,v))
enddo
return
end


! ================
! SUBROUTINE DV2UV        !SP->GP
! ================

subroutine dv2uv(pd,pz,pu,pv)
use legmod
implicit none

real :: pd(2,NESP/2,NLEV)
real :: pz(2,NESP/2,NLEV)
real :: pu(2,NLON/2,NLPP,NLEV)
real :: pv(2,NLON/2,NLPP,NLEV)
real :: zsave, u1, u2, v1, v2

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: v ! Loop index for level
integer :: w, w0, w1 ! Loop index for spectral mode

pu(:,:,:,:) = 0.0
pv(:,:,:,:) = 0.0

do v = 1 , NLEV
  zsave = pz(1,2,v)
  pz(1,2,v) = zsave - plavor
  do l = 1 , NLPP
    w0 = 1
    do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      u1 = 0.0
      u2 = 0.0
      v1 = 0.0
      v2 = 0.0
      do w = w0 , w1
        u1 = u1 + qv(w,l)*pz(1,w,v) + qu(w,l)*pd(2,w,v)
        u2 = u2 + qv(w,l)*pz(2,w,v) - qu(w,l)*pd(1,w,v)
        v1 = v1 + qu(w,l)*pz(2,w,v) - qv(w,l)*pd(1,w,v)
        v2 = v2 - qu(w,l)*pz(1,w,v) - qv(w,l)*pd(2,w,v)
      enddo ! w
      pu(1,m,l,v) = u1
      pu(2,m,l,v) = u2
      pv(1,m,l,v) = v1
      pv(2,m,l,v) = v2
      w0 = w1 + 1
    enddo ! m
  enddo ! l
  pz(1,2,v) = zsave
enddo ! jv
return
end


! ================
! SUBROUTINE UV2DV        !GP->SP
! ================

subroutine uv2dv(pu,pv,pd,pz)
use legmod
implicit none

real :: pd(2,NESP/2,NLEV)
real :: pz(2,NESP/2,NLEV)
real :: pu(2,NLON/2,NLPP,NLEV)
real :: pv(2,NLON/2,NLPP,NLEV)

integer :: k ! Loop index for southern latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: v ! Loop index for level
integer :: w, w0, w1 ! Loop index for spectral mode
real :: u1, u2, v1, v2
real :: u1p, u2p, u1m, u2m, v1p, v2p, v1m, v2m

pd(:,:,:) = 0.0
pz(:,:,:) = 0.0

if (NLPP < NLAT) then  ! Universal (parallel executable) version
!----------------------------------------------------------------------
do v = 1 , NLEV
  do l = 1 , NLPP
    w0 = 1
    do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      u1 = pu(1,m,l,v)
      u2 = pu(2,m,l,v)
      v1 = pv(1,m,l,v)
      v2 = pv(2,m,l,v)
      do w = w0 , w1
        pz(1,w,v) = pz(1,w,v) + qe(w,l)*u1 - qm(w,l)*v2
        pz(2,w,v) = pz(2,w,v) + qe(w,l)*u2 + qm(w,l)*v1
        pd(1,w,v) = pd(1,w,v) - qe(w,l)*v1 - qm(w,l)*u2
        pd(2,w,v) = pd(2,w,v) - qe(w,l)*v2 + qm(w,l)*u1
      enddo ! w
      w0 = w1 + 1
    enddo ! m
  enddo ! l
enddo ! v
else                   ! Single CPU version (symmetry conserving)
!----------------------------------------------------------------------
do v = 1 , NLEV
  do l = 1 , NLAT/2
    k = NLAT+1-l
    w0 = 1
    do m = 1 , NTP1
      w1 = w0 + NTP1 - m
      u1m = pu(1,m,l,v) - pu(1,m,k,v)
      u2m = pu(2,m,l,v) - pu(2,m,k,v)
      u1p = pu(1,m,l,v) + pu(1,m,k,v)
      u2p = pu(2,m,l,v) + pu(2,m,k,v)
      v1m = pv(1,m,l,v) - pv(1,m,k,v)
      v2m = pv(2,m,l,v) - pv(2,m,k,v)
      v1p = pv(1,m,l,v) + pv(1,m,k,v)
      v2p = pv(2,m,l,v) + pv(2,m,k,v)
      do w = w0 , w1 , 2
        pz(1,w,v) = pz(1,w,v) + qe(w,l)*u1m - qm(w,l)*v2p
        pz(2,w,v) = pz(2,w,v) + qe(w,l)*u2m + qm(w,l)*v1p
        pd(1,w,v) = pd(1,w,v) - qe(w,l)*v1m - qm(w,l)*u2p
        pd(2,w,v) = pd(2,w,v) - qe(w,l)*v2m + qm(w,l)*u1p
      enddo
      do w = w0 + 1 , w1 , 2
        pz(1,w,v) = pz(1,w,v) + qe(w,l)*u1p - qm(w,l)*v2m
        pz(2,w,v) = pz(2,w,v) + qe(w,l)*u2p + qm(w,l)*v1m
        pd(1,w,v) = pd(1,w,v) - qe(w,l)*v1p - qm(w,l)*u2m
        pd(2,w,v) = pd(2,w,v) - qe(w,l)*v2p + qm(w,l)*u1m
      enddo
      w0 = w1 + 1
    enddo ! m
  enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end 


! ================
! SUBROUTINE QTEND
! ================

subroutine qtend(q,qn,uq,vq)
use legmod
implicit none

real, intent(in) :: qn(2,NLON/2,NLPP,NLEV)
real, intent(in) :: uq(2,NLON/2,NLPP,NLEV)
real, intent(in) :: vq(2,NLON/2,NLPP,NLEV)

real, intent(out) :: q(2,NESP/2,NLEV)

integer :: k ! Loop index for southern latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: v ! Loop index for level
integer :: w, w0, w1 ! Loop index for spectral mode
real :: qn1, qn2, uq1, uq2, vq1, vq2
real :: qn1p, qn2p, qn1m, qn2m, uq1p, uq2p, uq1m, uq2m, vq1p, vq2p, vq1m, vq2m

q(:,:,:) = 0.0

if (NLPP < NLAT) then  ! Universal (parallel executable) version
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLPP
  w0 = 1
  do m = 1 , NTP1
   w1 = w0 + NTP1 - m
   qn1 = qn(1,m,l,v)
   qn2 = qn(2,m,l,v)
   uq1 = uq(1,m,l,v)
   uq2 = uq(2,m,l,v)
   vq1 = vq(1,m,l,v)
   vq2 = vq(2,m,l,v)
   do w = w0 , w1
    q(1,w,v) = q(1,w,v) + qe(w,l)*vq1 + qc(w,l)*qn1 + qm(w,l)*uq2
    q(2,w,v) = q(2,w,v) + qe(w,l)*vq2 + qc(w,l)*qn2 - qm(w,l)*uq1
   enddo ! w
   w0 = w1 + 1
  enddo ! m
 enddo ! l
enddo ! v
else                   ! Single CPU version (symmetry conserving)
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLAT/2
  k = NLAT+1-l
  w0 = 1
  do m = 1 , NTP1
   w1 = w0 + NTP1 - m
   vq1m = vq(1,m,l,v) - vq(1,m,k,v)
   vq2m = vq(2,m,l,v) - vq(2,m,k,v)
   vq1p = vq(1,m,l,v) + vq(1,m,k,v)
   vq2p = vq(2,m,l,v) + vq(2,m,k,v)
   qn1p = qn(1,m,l,v) + qn(1,m,k,v)
   qn2p = qn(2,m,l,v) + qn(2,m,k,v)
   qn1m = qn(1,m,l,v) - qn(1,m,k,v)
   qn2m = qn(2,m,l,v) - qn(2,m,k,v)
   uq1p = uq(1,m,l,v) + uq(1,m,k,v)
   uq2p = uq(2,m,l,v) + uq(2,m,k,v)
   uq1m = uq(1,m,l,v) - uq(1,m,k,v)
   uq2m = uq(2,m,l,v) - uq(2,m,k,v)
   do w = w0 , w1 , 2
      q(1,w,v) = q(1,w,v) + qe(w,l)*vq1m + qc(w,l)*qn1p + qm(w,l)*uq2p
      q(2,w,v) = q(2,w,v) + qe(w,l)*vq2m + qc(w,l)*qn2p - qm(w,l)*uq1p
   enddo
   do w = w0 + 1 , w1 , 2
      q(1,w,v) = q(1,w,v) + qe(w,l)*vq1p + qc(w,l)*qn1m + qm(w,l)*uq2m
      q(2,w,v) = q(2,w,v) + qe(w,l)*vq2p + qc(w,l)*qn2m - qm(w,l)*uq1m
   enddo
   w0 = w1 + 1
  enddo ! m
 enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end

! =================
! SUBROUTINE MKTEND
! =================

subroutine mktend(d,t,z,tn,fu,fv,ke,ut,vt)
use legmod
implicit none

real, intent(in) :: tn(2,NLON/2,NLPP,NLEV)
real, intent(in) :: fu(2,NLON/2,NLPP,NLEV)
real, intent(in) :: fv(2,NLON/2,NLPP,NLEV)
real, intent(in) :: ke(2,NLON/2,NLPP,NLEV)
real, intent(in) :: ut(2,NLON/2,NLPP,NLEV)
real, intent(in) :: vt(2,NLON/2,NLPP,NLEV)

real, intent(out) :: d(2,NESP/2,NLEV)
real, intent(out) :: t(2,NESP/2,NLEV)
real, intent(out) :: z(2,NESP/2,NLEV)

integer :: k ! Loop index for southern latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: v ! Loop index for level
integer :: w, w0, w1 ! Loop index for spectral mode
real :: tn1, tn2, fu1, fu2, fv1, fv2, ke1, ke2, ut1, ut2, vt1, vt2
real :: tn1p, tn2p, tn1m, tn2m, fu1p, fu2p, fu1m, fu2m
real :: fv1p, fv2p, fv1m, fv2m, ke1p, ke2p, ke1m, ke2m
real :: ut1p, ut2p, ut1m, ut2m, vt1p, vt2p, vt1m, vt2m

d(:,:,:) = 0.0
t(:,:,:) = 0.0
z(:,:,:) = 0.0

if (NLPP < NLAT) then  ! Universal (parallel executable) version
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLPP
  w0 = 1
  do m = 1 , NTP1
   w1 = w0 + NTP1 - m
   ke1 = ke(1,m,l,v); ke2 = ke(2,m,l,v)
   fu1 = fu(1,m,l,v); fu2 = fu(2,m,l,v)
   fv1 = fv(1,m,l,v); fv2 = fv(2,m,l,v)
   tn1 = tn(1,m,l,v); tn2 = tn(2,m,l,v)
   ut1 = ut(1,m,l,v); ut2 = ut(2,m,l,v)
   vt1 = vt(1,m,l,v); vt2 = vt(2,m,l,v)
   do w = w0 , w1
    d(1,w,v) = d(1,w,v) + qq(w,l)*ke1 - qe(w,l)*fv1 - qm(w,l)*fu2
    d(2,w,v) = d(2,w,v) + qq(w,l)*ke2 - qe(w,l)*fv2 + qm(w,l)*fu1
    t(1,w,v) = t(1,w,v) + qe(w,l)*vt1 + qc(w,l)*tn1 + qm(w,l)*ut2
    t(2,w,v) = t(2,w,v) + qe(w,l)*vt2 + qc(w,l)*tn2 - qm(w,l)*ut1
    z(1,w,v) = z(1,w,v) + qe(w,l)*fu1 - qm(w,l)*fv2
    z(2,w,v) = z(2,w,v) + qe(w,l)*fu2 + qm(w,l)*fv1
   enddo ! w
   w0 = w1 + 1
  enddo ! m
 enddo ! l
enddo ! v
else                   ! Single CPU version (symmetry conserving)
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLAT/2
  k = NLAT+1-l
  w0 = 1
  do m = 1 , NTP1
   w1 = w0 + NTP1 - m
   ke1p = ke(1,m,l,v) + ke(1,m,k,v); ke1m = ke(1,m,l,v) - ke(1,m,k,v)
   ke2p = ke(2,m,l,v) + ke(2,m,k,v); ke2m = ke(2,m,l,v) - ke(2,m,k,v)
   fu1p = fu(1,m,l,v) + fu(1,m,k,v); fu1m = fu(1,m,l,v) - fu(1,m,k,v)
   fu2p = fu(2,m,l,v) + fu(2,m,k,v); fu2m = fu(2,m,l,v) - fu(2,m,k,v)
   fv1p = fv(1,m,l,v) + fv(1,m,k,v); fv1m = fv(1,m,l,v) - fv(1,m,k,v)
   fv2p = fv(2,m,l,v) + fv(2,m,k,v); fv2m = fv(2,m,l,v) - fv(2,m,k,v)
   tn1p = tn(1,m,l,v) + tn(1,m,k,v); tn1m = tn(1,m,l,v) - tn(1,m,k,v)
   tn2p = tn(2,m,l,v) + tn(2,m,k,v); tn2m = tn(2,m,l,v) - tn(2,m,k,v)
   ut1p = ut(1,m,l,v) + ut(1,m,k,v); ut1m = ut(1,m,l,v) - ut(1,m,k,v)
   ut2p = ut(2,m,l,v) + ut(2,m,k,v); ut2m = ut(2,m,l,v) - ut(2,m,k,v)
   vt1p = vt(1,m,l,v) + vt(1,m,k,v); vt1m = vt(1,m,l,v) - vt(1,m,k,v)
   vt2p = vt(2,m,l,v) + vt(2,m,k,v); vt2m = vt(2,m,l,v) - vt(2,m,k,v)
   do w = w0 , w1 , 2
      d(1,w,v) = d(1,w,v) + qq(w,l)*ke1p - qe(w,l)*fv1m - qm(w,l)*fu2p
      d(2,w,v) = d(2,w,v) + qq(w,l)*ke2p - qe(w,l)*fv2m + qm(w,l)*fu1p
      t(1,w,v) = t(1,w,v) + qe(w,l)*vt1m + qc(w,l)*tn1p + qm(w,l)*ut2p
      t(2,w,v) = t(2,w,v) + qe(w,l)*vt2m + qc(w,l)*tn2p - qm(w,l)*ut1p
      z(1,w,v) = z(1,w,v) + qe(w,l)*fu1m - qm(w,l)*fv2p
      z(2,w,v) = z(2,w,v) + qe(w,l)*fu2m + qm(w,l)*fv1p
   enddo
   do w = w0 + 1 , w1 , 2
      d(1,w,v) = d(1,w,v) + qq(w,l)*ke1m - qe(w,l)*fv1p - qm(w,l)*fu2m
      d(2,w,v) = d(2,w,v) + qq(w,l)*ke2m - qe(w,l)*fv2p + qm(w,l)*fu1m
      t(1,w,v) = t(1,w,v) + qe(w,l)*vt1p + qc(w,l)*tn1m + qm(w,l)*ut2m
      t(2,w,v) = t(2,w,v) + qe(w,l)*vt2p + qc(w,l)*tn2m - qm(w,l)*ut1m
      z(1,w,v) = z(1,w,v) + qe(w,l)*fu1p - qm(w,l)*fv2m
      z(2,w,v) = z(2,w,v) + qe(w,l)*fu2p + qm(w,l)*fv1m
   enddo
   w0 = w1 + 1
  enddo ! m
 enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end


! ================
! SUBROUTINE SP2FL
! ================

subroutine sp2fl(psp,pfc,klev)
use legmod
implicit none

integer, intent(in ) :: klev
real,    intent(in ) :: psp(NESP,klev)
real,    intent(out) :: pfc(NHOR,klev)

integer :: jlev

do jlev = 1,klev
   call sp2fc(psp(1,jlev),pfc(1,jlev))
enddo

return
end 


! ==================
! SUBROUTINE INVLEGA
! ==================

subroutine invlega
use pumamod
implicit none

integer :: jlev

call dv2uv(sd,sz,gu,gv)

do jlev = 1,NLEV
  call sp2fc(sd(1,jlev),gd(1,jlev))
  call sp2fc(st(1,jlev),gt(1,jlev))
  call sp2fc(sz(1,jlev),gz(1,jlev))
  if (nqspec == 1) call sp2fc(sq(1,jlev),gq(1,jlev))
enddo

call sp2fc(sp,gp)
call sp2fcdmu(sp,gpj)

return
end


! ==================
! SUBROUTINE INVLEGD
! ==================

subroutine invlegd
use pumamod
implicit none

integer :: jlev

call dv2uv(sd,sz,gu,gv)

do jlev = 1,NLEV
  if (nqspec == 1) call sp2fc(sq(1,jlev),gq(1,jlev))
  call sp2fc(st(1,jlev),gt(1,jlev))
enddo

call sp2fc(sp,gp)

return
end

