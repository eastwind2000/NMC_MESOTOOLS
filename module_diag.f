      module diag
      
      real:: psadiprs(150),psadithte(150),psaditmk(150,150)

      real,parameter :: rgas=287.04  !J/K/kg
      real,parameter :: rgasmd=.608   ! rgas_moist=rgas*(1.+rgasmd*qvp)
      real,parameter :: cp=1004.     ! J/K/kg  Note: not using Bolton's value of 1005.7
      real,parameter :: cpmd=.887   ! cp_moist=cp*(1.+cpmd*qvp)
      real,parameter :: gamma=rgas/cp
      real,parameter :: gammamd=rgasmd-cpmd  ! gamma_moist=gamma*(1.+gammamd*qvp)
      real,parameter :: grav=9.81           ! m/s**2
      real,parameter :: sclht=rgas*256./grav   ! 256 K is avg. trop. temp. from USSA.
      real,parameter :: eps=0.622
      real,parameter :: ezero=6.112  ! hPa
      real,parameter :: xlhc0=3.1484e6   ! J/kg
      real,parameter :: xlhctd=2370.  !
      real,parameter :: xlhf=3.34e5
      real,parameter :: rktpmps=1.94
      real,parameter :: celkel=273.15
      real,parameter :: eslcon1=17.67
      real,parameter :: eslcon2=29.65
      real,parameter :: esicon1=22.514
      real,parameter :: esicon2=6.15e3
      real,parameter :: thtecon1=3376. ! K
      real,parameter :: thtecon2=2.54
      real,parameter :: thtecon3=.81
      real,parameter :: tlclc1=2840.
      real,parameter :: tlclc2=3.5
      real,parameter :: tlclc3=4.805
      real,parameter :: tlclc4=55.
      real,parameter :: rhoice=917.
      real,parameter :: rhowat=1000.
      real,parameter :: pi=4.*atan(1.)
      real,parameter :: rpd=pi/180.
      real,parameter :: abscoef=.145      ! cloud water absorption coefficient in m^2/g
      real,parameter :: abscoefi=.272     ! cloud ice absorption coefficient in m^2/g
      real,parameter :: ussalr=.0065      ! deg C per m
      real,parameter :: rmsg=9.0e+9       ! indicates missing data or specification

     
      contains
!===============================================================      
     
      subroutine wetbulbcalc(prs,tmk,qvp,twb,miy,mjx,mkzh)
c   This routine calculates the wet bulb temperature (in Celsius)
c
      dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   twb(miy,mjx,mkzh)
c
c
      do 80 k=1,mkzh
      do 80 j=1,mjx-1
      do 80 i=1,miy-1
            q=max(qvp(i,j,k),1.e-15)
            t=tmk(i,j,k)
            p=prs(i,j,k)
            e=q*p/(eps+q)
            tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
            eth=t*(1000./p)**(gamma*(1.+gammamd*q))*
     &         exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
            twb(i,j,k)=tonpsadiabat(eth,p)-celkel
  80  continue
c
      return
      endsubroutine
      
      
!===============================================================      

c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine thecalc(prs,tmk,qvp,the,miy,mjx,mkzh)
c
         
        dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   the(miy,mjx,mkzh)
c

c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
         gammam=gamma*(1.+gammamd*qvp(i,j,k))
         the(i,j,k)=tmk(i,j,k)*(1000./prs(i,j,k))**gammam
 1000 continue
      return
      end subroutine

c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine sateqthecalc(tmk,prs,sateth,miy,mjx,mkzh)

      dimension tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),sateth(miy,mjx,mkzh)
c

c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
            t=tmk(i,j,k)
            p=prs(i,j,k)
            esat=ezero*exp(eslcon1*(t-celkel)/
     &         (t-eslcon2))
            qsat=eps*esat/(p-esat)
            sateth(i,j,k)=t*(1000./p)**(gamma*(1.+gammamd*qsat))*
     &         exp((thtecon1/t-thtecon2)*qsat*(1.+thtecon3*qsat))
 1000 continue
      return
      end subroutine

c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine tdpcalc(qvp,prs,tdp,miy,mjx,mkzh)

c   This routine calculates dewpoint temperature (in Celsius)
c
      dimension qvp(miy,mjx,mkzh),prs(miy,mjx,mkzh),tdp(miy,mjx,mkzh)
c

c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
         q=max(qvp(i,j,k),1.e-15)
         p=prs(i,j,k)
         e=q*p/(eps+q)
         elog=alog(e/ezero)
         tdp(i,j,k)=(eslcon2*elog-eslcon1*celkel)/(elog-eslcon1)-
     &      celkel
 1000 continue
      return
      end subroutine
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine qvpcalc(rhu,tmk,prs,icerel,qvp,miy,mjx,mkzh)

      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),rhu(miy,mjx,mkzh)
c
      do 200 k=1,mkzh
      do 200 j=1,mjx-1
      do 200 i=1,miy-1
	    p = prs(i,j,k) 
         rh = rhu(i,j,k)/100.0
          t = tmk(i,j,k)

c          e = q*prs(i,j,k)/(eps+q)
c
c   es should be wrt ice if user has specified icerel and temperature
c   is below freezing.
c
         if (icerel.eq.1.and.t.lt.celkel) then
            es = ezero * exp( esicon1-esicon2/t )
         else
            es = ezero * exp( eslcon1*(t-celkel)/(t-eslcon2) )
         endif

         e = rh*es*p / ( p + rh*es - es   ) 
         
	   qvp(i,j,k)= e*eps/( p-e ) 

  200 continue
      return
      end subroutine
!========================================================
      
      
     
      subroutine pwatcalc( psfc, qs, q, plev, pwatc, nz) 
      integer           :: nz
      real,dimension(nz):: plev, q
      real              :: psfc, qs, pwatc
      real              :: dp0, dp

      pwatc=0

      if( psfc >= plev(1)) then

        dp0=psfc-plev(1)
        q0= ( qs+q(1) )/2.0
        pwatc=q0*dp0

        do k= 1, nz-1
           dp=plev(k)-plev(k+1)
           pwatc=pwatc+ 0.5*( q(k) + q(k+1) )*dp
        enddo
       
       elseif( psfc<plev(1) )then

           do  k=nz,1,-1
              if(psfc>plev(k))then
                k0=k
              endif
          enddo

c          print*, k0, plev(k0), psfc

          dp0=psfc-plev(k0)
          q0=(qs+q(k0))/2.0
          pwatc=q0*dp0      
          
c          print*, k0+2, nz-1     
          do k=k0, nz-1
             dp=plev(k)-plev(k+1)
             pwatc=pwatc+ 0.5*( q(k) + q(k+1) )*dp
          enddo
       endif

       pwatc=pwatc/grav

      end subroutine



    
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine eqthecalc(qvp,tmk,prs,eth,miy,mjx,mkzh)

      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),eth(miy,mjx,mkzh)
c

c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
            q=max(qvp(i,j,k),1.e-15)
            t=tmk(i,j,k)
            p=prs(i,j,k)
            e=q*p/(eps+q)
            tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
            eth(i,j,k)=t*(1000./p)**(gamma*(1.+gammamd*q))*
     &         exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
 1000 continue
      return
      end subroutine

c
c*********************************************************************c
c                                                                     c
      function tonpsadiabat(thte,prs)
c
c   This function gives the temperature (in K) on a moist adiabat
c   (specified by thte in K) given pressure in hPa.  It uses a
c   lookup table, with data that was generated by the Bolton (1980)
c   formula for theta_e.
c
c
c     First check if pressure is less than min pressure in lookup table.
c     If it is, assume parcel is so dry that the given theta-e value can
c     be interpretted as theta, and get temperature from the simple dry
c     theta formula.
c
      if (prs.le.psadiprs(150)) then
         tonpsadiabat=thte*(prs/1000.)**gamma
         return
      endif
c
c   Otherwise, look for the given thte/prs point in the lookup table.
c
      do jtch=1,150-1
         if (thte.ge.psadithte(jtch).and.thte.lt.psadithte(jtch+1)) then
            jt=jtch
            goto 213
         endif
      enddo
      jt=-1
 213  continue
      do ipch=1,150-1
         if (prs.le.psadiprs(ipch).and.prs.gt.psadiprs(ipch+1)) then
            ip=ipch
            goto 215
         endif
      enddo
      ip=-1
 215  continue
      if (jt.eq.-1.or.ip.eq.-1) then
         write(*,*)
     &      'Outside of lookup table bounds. prs,thte=',prs,thte
         stop
      endif
      fracjt=(thte-psadithte(jt))/(psadithte(jt+1)-psadithte(jt))
      fracjt2=1.-fracjt
      fracip=(psadiprs(ip)-prs)/(psadiprs(ip)-psadiprs(ip+1))
      fracip2=1.-fracip
      if (psaditmk(ip,jt).gt.1e9.or.psaditmk(ip+1,jt).gt.1e9.or.
     &    psaditmk(ip,jt+1).gt.1e9.or.psaditmk(ip+1,jt+1).gt.1e9) then
         write(*,*)
     &      'Tried to access missing tmperature in lookup table.'
         write(*,*)
     &      'Prs and Thte probably unreasonable. prs,thte=',prs,thte
         stop
      endif
      tonpsadiabat=fracip2*fracjt2*psaditmk(ip  ,jt  )+
     &       fracip *fracjt2*psaditmk(ip+1,jt  )+
     &       fracip2*fracjt *psaditmk(ip  ,jt+1)+
     &       fracip *fracjt *psaditmk(ip+1,jt+1)
c
      return
      endfunction


      
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine smooth(pslab,work,numpas,mabpl,njx,niy)
c
c   This is a smoothing routine, with several choices:
c
c   If numpas is between 1 and 99, then a 9-point weighted smoother is
c   applied numpas times.  The smoother follows equation 11-107 in
c   Haltiner and Williams. One pass completely removes 2-delta-x waves
c   on the interior.  On the outer row and column, and near missing
c   data points, smoothing is carried out in a manner that preserves
c   the domain average value of the field.
c
c   If numpas is between 101 and 199, then a smoother-desmoother is
c   applied (numpas-100) times.  One pass removes a large fraction
c   of the 2-delta-x component, but is not as harsh on longer
c   wavelengths as the 9-point smoother
c
c   If numpas is between 201 and 299, then the smoother-desmoother is
c   applied (numpas-200) times, and, after each pass, the data field
c   is forced to be non-negative.
c
c   If numpas is between 301 and 399, then a weighted
c   smoother is applied, in which the smoothed value
c   is given by a weighted average of values at
c   surrounding grid points.  The weighting function
c   is the Cressman weighting function:
c
c               w = ( D**2 - d**2 ) / ( D**2 + d**2 )
c
c   In the above, d is the distance (in grid increments)
c   of the neighboring point to the smoothing point, and
c   D is the radius of influence [in grid increments,
c   given by (numpas-300)].
c
c   If numpas is between 401 and 499, then the smoothing
c   is similar for numpas=301-399, except the weighting
c   function is the circular apperture diffraction function
c   (following a suggestion of Barnes et al. 1996):
c
c               w = bessel(3.8317*d/D)/(3.8317*d/D)
c
c   If numpas is between 501 and 599, then the smoothing
c   is similar for numpas=301-399, except the weighting
c   function is the product of the rectangular
c   apperture diffraction function in the x and y directions
c   (the function used in Barnes et al. 1996):
c
c               w = [sin(pi*x/D)/(pi*x/D)]*[sin(pi*y/D)/(pi*y/D)]
c
c   Note, the first index of pslab varies along the abcissa
c   (or x), and the second index varies along the ordinate (or y).
c
      parameter(beszero=3.8317)
c
      dimension pslab(mabpl,niy),work(mabpl,niy),xnu(2),fprint(150,150)
c
!      include 'comconst'
c
      if (mod(numpas,100).eq.0) return
c
      if (numpas.le.99) then   ! 9-point smoother
c
      do ipas=1,numpas
c
      do i=1,niy
      do j=1,njx
         work(j,i)=0.
      enddo
      enddo
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).eq.rmsg) then
            work(j,i)=rmsg
         else
            totgive=0.
            if (i.gt.1) then
               if (pslab(j,i-1).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j,i-1)=work(j,i-1)+give
                  totgive=totgive+give
               endif
               if (j.gt.1) then
                  if (pslab(j-1,i-1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j-1,i-1)=work(j-1,i-1)+give
                     totgive=totgive+give
                  endif
               endif
               if (j.lt.njx) then
                  if (pslab(j+1,i-1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j+1,i-1)=work(j+1,i-1)+give
                     totgive=totgive+give
                  endif
               endif
            endif
            if (i.lt.niy) then
               if (pslab(j,i+1).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j,i+1)=work(j,i+1)+give
                  totgive=totgive+give
               endif
               if (j.gt.1) then
                  if (pslab(j-1,i+1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j-1,i+1)=work(j-1,i+1)+give
                     totgive=totgive+give
                  endif
               endif
               if (j.lt.njx) then
                  if (pslab(j+1,i+1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j+1,i+1)=work(j+1,i+1)+give
                     totgive=totgive+give
                  endif
               endif
            endif
            if (j.gt.1) then
               if (pslab(j-1,i).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j-1,i)=work(j-1,i)+give
                  totgive=totgive+give
               endif
            endif
            if (j.lt.njx) then
               if (pslab(j+1,i).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j+1,i)=work(j+1,i)+give
                  totgive=totgive+give
               endif
            endif
            work(j,i)=work(j,i)+pslab(j,i)-totgive
         endif
      enddo
      enddo
      do i=1,niy
      do j=1,njx
         pslab(j,i)=work(j,i)
      enddo
      enddo
c
      enddo
c
      elseif (numpas.le.299) then   ! smoother-desmoother
c
      if (numpas.ge.200) then
         nump=numpas-200
         inn=1
      else
         nump=numpas-100
         inn=0
      endif
c
      if (nump.lt.1) return
c
c   Check if any data is missing.
c
      imsg=0
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).eq.rmsg) then
            imsg=1
            goto 15
         endif
      enddo
      enddo
 15   continue
c
      if (imsg.eq.1) then
c
c   Get average value of pslab.
c
      nval=0
      total=0.
      do 10 i=1,niy
      do 10 j=1,njx
         if (pslab(j,i).ne.rmsg) then
            total=total+pslab(j,i)
            nval=nval+1
         endif
   10 continue
      if (nval.eq.0) then
         write(*,*)'  All elements of this pslab are rmsg.'
         return
      endif
      avgval=total/nval
c
c   Set each element that is currently rmsg to avgval, and
c   keep track of them with work.
c
      do 20 i=1,niy
      do 20 j=1,njx
         if (pslab(j,i).eq.rmsg) then
            pslab(j,i)=avgval
            work(j,i)=1.
         else
            work(j,i)=0.
         endif
   20 continue
c
      endif
c
c     *** Do calculation and put into pslab array.
c
      xnu(1) = 0.50
      xnu(2) = -0.52
      je = njx - 1
      ie = niy - 1
      do 100 ipass = 1,nump*2
         kp=2-mod(ipass,2)          
c     
c        *** First, smooth in the njx direction.
c     
         do 60 j = 2,je
            asv = pslab(j,1)
            do 50 i = 2,ie
               aplus = pslab(j,i+1)
               cell = pslab(j,i)
               pslab(j,i)= pslab(j,i) + xnu(kp)*
     +            ((asv + aplus)/2.0 - pslab(j,i))
               asv = cell
   50       continue
   60    continue
c     
c        *** Now, smooth in the niy direction.
c     
         do 80 i = 2,ie
            asv = pslab(1,i)
            do 70 j = 2,je
               aplus = pslab(j+1,i)
               cell = pslab(j,i)
               pslab(j,i) = pslab(j,i) + xnu(kp)*
     +            ((asv + aplus)/2.0 - pslab(j,i))
               asv = cell
   70       continue
   80    continue
c
      if (inn.eq.1) then
c
c      Make non-negative.
c
         do i=1,niy
         do j=1,njx
            pslab(j,i)=max(0.,pslab(j,i))
         enddo
         enddo
      endif
c
  100 continue
c
      if (imsg.eq.1) then
c
c      Set rmsg elements back to rmsg
c
         do 200 i=1,niy
         do 200 j=1,njx
            pslab(j,i)=work(j,i)*rmsg + (1.-work(j,i))*pslab(j,i)
  200    continue
      endif
c
      elseif (numpas.le.599) then   ! weighted smoother
c
      idist=mod(numpas,100)
      if (idist.eq.0) return
      nfp=1+2*idist
      npsq=idist*idist
      if (numpas.le.399) then  ! Cressman function
         do i=1,nfp
         do j=1,nfp
            distsq=(i-idist-1.)**2+(j-idist-1.)**2
            fprint(j,i)=max((npsq-distsq)/(npsq+distsq),0.0)
         enddo
         enddo
      elseif (numpas.le.499) then   ! Circular diffraction function
         do i=1,nfp
         do j=1,nfp
            dist=beszero/idist*sqrt((i-idist-1.)**2+(j-idist-1.)**2)
            if (i.eq.idist+1.and.j.eq.idist+1) then
               fprint(j,i)=.5
            else
               fprint(j,i)=max(0.,bes(dist)/dist)
            endif
         enddo
         enddo
      elseif (numpas.le.599) then   ! Rect. diffraction function
         do i=1,nfp
         do j=1,nfp
            if (j.eq.idist+1) then
               xfac=1.
            else
               xdist=pi/idist*(j-idist-1.)
               xfac=sin(xdist)/xdist
            endif
            if (i.eq.idist+1) then
               yfac=1.
            else
               ydist=pi/idist*(i-idist-1.)
               yfac=sin(ydist)/ydist
            endif
            fprint(j,i)=xfac*yfac
         enddo
         enddo
      endif
c
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).ne.rmsg) then
            tot=0.
            totwt=0.
            is=max(1,i-idist)
            ie=min(niy,i+idist)
            js=max(1,j-idist)
            je=min(njx,j+idist)
            do ireg=is,ie
               ifp=ireg-i+idist+1
            do jreg=js,je
               jfp=jreg-j+idist+1
               if (pslab(jreg,ireg).ne.rmsg) then
                  totwt=totwt+fprint(jfp,ifp)
                  tot=tot+fprint(jfp,ifp)*pslab(jreg,ireg)
               endif
            enddo
            enddo
            work(j,i)=tot/totwt
         else
            work(j,i)=rmsg
         endif
      enddo
      enddo
c
      do i=1,niy
      do j=1,njx
         pslab(j,i)=work(j,i)
      enddo
      enddo
c
      endif
c
      return
      endsubroutine

      

      function bes(x)
      rint=0.
      do i=1,1000
         u=i*.001-.0005
         rint=rint+sqrt(1.-u*u)*cos(x*u)*.001
      enddo
      bes=2.*x*rint/(4.*atan(1.))
      return
      endfunction

      
      
      
      
      
      
      end module
      
      
