* SOLID EARTH TIDE
* Dennis Milbert, 2015
* http://home.comcast.net/~dmilbert/softs/solid.htm
* modified by Jose Gama July 2015
* It is no longer a standalone program but a library
* It can create an output file or return an array with the output data

      subroutine initialize()
      implicit double precision(a-h,o-z)
      common/stuff/rad,pi,pi2
      common/comgrs/a,e2
      save /stuff/, /comgrs/
*** constants

      pi=4.d0*datan(1.d0)
      pi2=pi+pi
      rad=180.d0/pi
*** grs80

      a=6378137.d0
      e2=6.69438002290341574957d-03

      end subroutine initialize
*-----------------------------------------------------------------------

      subroutine calctide(iyr, imo, idy, glad, glod, iret, 
     &     retr, ierr)
* iyr	year    [1980-2016]
* imo	month number [1-12]
* idy	day          [1-31]
* glad	Lat. (pos N.) [- 90, +90]
* glod	Lon. (pos E.) [-360,+360]

      implicit double precision(a-h,o-z)
      dimension rsun(3),rmoon(3),etide(3),xsta(3)

      double precision retr(1441, 4)

      common/stuff/rad,pi,pi2
      common/comgrs/a,e2
      call initialize()
      lout=1
*      if(iret .eq.  0.d0) then
*        open(lout,file=fname,form='formatted',status='unknown')
*      end if

*** position of observing point (positive East)

      if(glod.lt.  0.d0) glod=glod+360.d0
      if(glod.ge.360.d0) glod=glod-360.d0

      gla0=glad/rad
      glo0=glod/rad
      eht0=0.d0
      call geoxyz(gla0,glo0,eht0,x0,y0,z0)
      xsta(1)=x0
      xsta(2)=y0
      xsta(3)=z0

*** header
      if(iret .eq.  0.d0) then
*        write(lout,'(a,i5,2i3)') 'year,month,day= ',iyr,imo,idy
*        write(lout,'(a,2f15.9)') 'lat, East lon.= ',glad,glod
      end if
*** here comes the sun  (and the moon)  (go, tide!)

      ihr=   0
      imn=   0
      sec=0.d0                                       !*** GPS time system
      call civmjd(iyr,imo,idy,ihr,imn,sec,mjd,fmjd)
      call mjdciv(mjd,fmjd,iyr,imo,idy,ihr,imn,sec)  !*** normalize civil time
      call setjd0(iyr,imo,idy)

      tdel2=1.d0/60.d0/24.d0                         !*** 1 minute in GPS days
      do iloop=0,60*24
        call sunxyz (mjd,fmjd,rsun)
        call moonxyz(mjd,fmjd,rmoon)
        call detide (xsta,mjd,fmjd,rsun,rmoon,etide)
        xt = etide(1)
        yt = etide(2)
        zt = etide(3)

*** determine local geodetic horizon components (topocentric)

        call rge(gla0,glo0,ut,vt,wt,xt,   yt,   zt)       !*** tide vector

        call mjdciv(mjd,fmjd               +0.001d0/86400.d0,
     *              iyr,imo,idy,ihr,imn,sec-0.001d0)

        tsec=ihr*3600.d0+imn*60.d0+sec
      if(iret .eq.  0.d0) then
*        write(lout,'(f8.1,3f10.6)') tsec,ut,vt,wt
      else
        retr(iloop+1, 1) = tsec
        retr(iloop+1, 2) = ut
        retr(iloop+1, 3) = vt
        retr(iloop+1, 4) = wt
      end if
        fmjd=fmjd+tdel2
        fmjd=(idnint(fmjd*86400.d0))/86400.d0      !*** force 1 sec. granularity
      enddo

      if(iret .eq.  0.d0) then
        close(lout)
      end if
      end








*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine detide(xsta,mjd,fmjd,xsun,xmon,dxtide)

*** computation of tidal corrections of station displacements caused
***    by lunar and solar gravitational attraction

*** step 1 (here general degree 2 and 3 corrections +
***         call st1idiu + call st1isem + call st1l1)
***   + step 2 (call step2diu + call step2lon + call step2idiu)
*** it has been decided that the step 3 un-correction for permanent tide
*** would *not* be applied in order to avoid jump in the reference frame
*** (this step 3 must added in order to get the mean tide station position
*** and to be conformed with the iag resolution.)

*** inputs
***   xsta(i),i=1,2,3   -- geocentric position of the station (ITRF/ECEF)
***   xsun(i),i=1,2,3   -- geoc. position of the sun (ECEF)
***   xmon(i),i=1,2,3   -- geoc. position of the moon (ECEF)
***   mjd,fmjd          -- modified julian day (and fraction) (in GPS time)

****old calling sequence*****************************************************
***   dmjd               -- time in mean julian date (including day fraction)
***   fhr=hr+zmin/60.+sec/3600.   -- hr in the day

*** outputs
***   dxtide(i),i=1,2,3           -- displacement vector (ITRF)

*** author iers 1996 :  v. dehant, s. mathews and j. gipson
***    (test between two subroutines)
*** author iers 2000 :  v. dehant, c. bruyninx and s. mathews
***    (test in the bernese program by c. bruyninx)

*** created:  96/03/23 (see above)
*** modified from dehanttideinelMJD.f by Dennis Milbert 2006sep10
*** bug fix regarding fhr (changed calling sequence, too)
*** modified to reflect table 7.5a and b IERS Conventions 2003
*** modified to use TT time system to call step 2 functions
*** sign correction by V.Dehant to match eq.16b, p.81, Conventions
*** applied by Dennis Milbert 2007may05

      implicit double precision(a-h,o-z)
      double precision xsta(3),xsun(3),xmon(3),dxtide(3),xcorsta(3)
      double precision h20,l20,h3,l3,h2,l2
      double precision mass_ratio_sun,mass_ratio_moon

*** nominal second degree and third degree love numbers and shida numbers

      data h20/0.6078d0/,l20/0.0847d0/,h3/0.292d0/,l3/0.015d0/

*** internal support for new calling sequence
*** also convert GPS time into TT time

      tsecgps=fmjd*86400.d0                       !*** GPS time (sec of day)
      tsectt =gps2tt(tsecgps)                     !*** TT  time (sec of day)
      fmjdtt =tsectt/86400.d0                     !*** TT  time (fract. day)

      dmjdtt=mjd+fmjdtt                           !*** float MJD in TT
*** commented line was live code in dehanttideinelMJD.f
*** changed on the suggestion of Dr. Don Kim, UNB -- 09mar21
*** Julian date for 2000 January 1 00:00:00.0 UT is  JD 2451544.5
*** MJD         for 2000 January 1 00:00:00.0 UT is MJD   51544.0
***** t=(dmjdtt-51545.d0)/36525.d0                !*** days to centuries, TT
      t=(dmjdtt-51544.d0)/36525.d0                !*** days to centuries, TT
      fhr=(dmjdtt-int(dmjdtt))*24.d0              !*** hours in the day, TT

*** scalar product of station vector with sun/moon vector

      call sprod(xsta,xsun,scs,rsta,rsun)
      call sprod(xsta,xmon,scm,rsta,rmon)
      scsun=scs/rsta/rsun
      scmon=scm/rsta/rmon

*** computation of new h2 and l2

      cosphi=dsqrt(xsta(1)*xsta(1) + xsta(2)*xsta(2))/rsta
      h2=h20-0.0006d0*(1.d0-3.d0/2.d0*cosphi*cosphi)
      l2=l20+0.0002d0*(1.d0-3.d0/2.d0*cosphi*cosphi)

*** p2-term

      p2sun=3.d0*(h2/2.d0-l2)*scsun*scsun-h2/2.d0
      p2mon=3.d0*(h2/2.d0-l2)*scmon*scmon-h2/2.d0

*** p3-term

      p3sun=5.d0/2.d0*(h3-3.d0*l3)*scsun**3+3.d0/2.d0*(l3-h3)*scsun
      p3mon=5.d0/2.d0*(h3-3.d0*l3)*scmon**3+3.d0/2.d0*(l3-h3)*scmon

*** term in direction of sun/moon vector

      x2sun=3.d0*l2*scsun
      x2mon=3.d0*l2*scmon
      x3sun=3.d0*l3/2.d0*(5.d0*scsun*scsun-1.d0)
      x3mon=3.d0*l3/2.d0*(5.d0*scmon*scmon-1.d0)

*** factors for sun/moon

      mass_ratio_sun=332945.943062d0
      mass_ratio_moon=0.012300034d0
      re =6378136.55d0
      fac2sun=mass_ratio_sun*re*(re/rsun)**3
      fac2mon=mass_ratio_moon*re*(re/rmon)**3
      fac3sun=fac2sun*(re/rsun)
      fac3mon=fac2mon*(re/rmon)

*** total displacement

      do i=1,3
        dxtide(i)=fac2sun*( x2sun*xsun(i)/rsun + p2sun*xsta(i)/rsta ) +
     *            fac2mon*( x2mon*xmon(i)/rmon + p2mon*xsta(i)/rsta ) +
     *            fac3sun*( x3sun*xsun(i)/rsun + p3sun*xsta(i)/rsta ) +
     *            fac3mon*( x3mon*xmon(i)/rmon + p3mon*xsta(i)/rsta )
      enddo
      call zero_vec8(xcorsta)

*** corrections for the out-of-phase part of love numbers
***     (part h_2^(0)i and l_2^(0)i )

*** first, for the diurnal band

      call st1idiu(xsta,xsun,xmon,fac2sun,fac2mon,xcorsta)
      dxtide(1)=dxtide(1)+xcorsta(1)
      dxtide(2)=dxtide(2)+xcorsta(2)
      dxtide(3)=dxtide(3)+xcorsta(3)

*** second, for the semi-diurnal band

      call st1isem(xsta,xsun,xmon,fac2sun,fac2mon,xcorsta)
      dxtide(1)=dxtide(1)+xcorsta(1)
      dxtide(2)=dxtide(2)+xcorsta(2)
      dxtide(3)=dxtide(3)+xcorsta(3)

*** corrections for the latitude dependence of love numbers (part l^(1) )

      call st1l1(xsta,xsun,xmon,fac2sun,fac2mon,xcorsta)
      dxtide(1)=dxtide(1)+xcorsta(1)
      dxtide(2)=dxtide(2)+xcorsta(2)
      dxtide(3)=dxtide(3)+xcorsta(3)

*** consider corrections for step 2
*** corrections for the diurnal band:

***  first, we need to know the date converted in julian centuries

***  this is now handled at top of code   (also convert to TT time system)
***** t=(dmjd-51545.)/36525.
***** fhr=dmjd-int(dmjd)             !*** this is/was a buggy line (day vs. hr)

***  second, the diurnal band corrections,
***   (in-phase and out-of-phase frequency dependence):

      call step2diu(xsta,fhr,t,xcorsta)
      dxtide(1)=dxtide(1)+xcorsta(1)
      dxtide(2)=dxtide(2)+xcorsta(2)
      dxtide(3)=dxtide(3)+xcorsta(3)

***  corrections for the long-period band,
***   (in-phase and out-of-phase frequency dependence):

      call step2lon(xsta,fhr,t,xcorsta)
      dxtide(1)=dxtide(1)+xcorsta(1)
      dxtide(2)=dxtide(2)+xcorsta(2)
      dxtide(3)=dxtide(3)+xcorsta(3)

*** consider corrections for step 3
*-----------------------------------------------------------------------
* The code below is commented to prevent restoring deformation
* due to permanent tide.  All the code above removes
* total tidal deformation with conventional Love numbers.
* The code above realizes a conventional tide free crust (i.e. ITRF). 
* This does NOT conform to Resolution 16 of the 18th General Assembly
* of the IAG (1983).  This resolution has not been implemented by 
* the space geodesy community in general (c.f. IERS Conventions 2003).
*-----------------------------------------------------------------------

*** uncorrect for the permanent tide  (only if you want mean tide system)

***   pi=3.141592654
***   sinphi=xsta(3)/rsta
***   cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
***   cosla=xsta(1)/cosphi/rsta
***   sinla=xsta(2)/cosphi/rsta
***   dr=-dsqrt(5./4./pi)*h2*0.31460*(3./2.*sinphi**2-0.5)
***   dn=-dsqrt(5./4./pi)*l2*0.31460*3.*cosphi*sinphi
***   dxtide(1)=dxtide(1)-dr*cosla*cosphi+dn*cosla*sinphi
***   dxtide(2)=dxtide(2)-dr*sinla*cosphi+dn*sinla*sinphi
***   dxtide(3)=dxtide(3)-dr*sinphi      -dn*cosphi

      return
      end
*-----------------------------------------------------------------------
      subroutine st1l1(xsta,xsun,xmon,fac2sun,fac2mon,xcorsta)

*** this subroutine gives the corrections induced by the latitude dependence
*** given by l^(1) in mahtews et al (1991)

***  input: xsta,xsun,xmon,fac3sun,fac3mon
*** output: xcorsta

      implicit double precision (a-h,o-z)
      dimension xsta(3),xsun(3),xmon(3),xcorsta(3)
      double precision l1,l1d,l1sd
      data l1d/0.0012d0/,l1sd/0.0024d0/

      rsta=enorm8(xsta)
      sinphi=xsta(3)/rsta
      cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
      sinla=xsta(2)/cosphi/rsta
      cosla=xsta(1)/cosphi/rsta
      rmon=enorm8(xmon)
      rsun=enorm8(xsun)

*** for the diurnal band

      l1=l1d
      dnsun=-l1*sinphi**2*fac2sun*xsun(3)*(xsun(1)*cosla+xsun(2)*sinla)
     *            /rsun**2
      dnmon=-l1*sinphi**2*fac2mon*xmon(3)*(xmon(1)*cosla+xmon(2)*sinla)
     *            /rmon**2
      desun=l1*sinphi*(cosphi**2-sinphi**2)*fac2sun*xsun(3)*
     * (xsun(1)*sinla-xsun(2)*cosla)/rsun**2
      demon=l1*sinphi*(cosphi**2-sinphi**2)*fac2mon*xmon(3)*
     * (xmon(1)*sinla-xmon(2)*cosla)/rmon**2
      de=3.d0*(desun+demon)
      dn=3.d0*(dnsun+dnmon)
      xcorsta(1)=-de*sinla-dn*sinphi*cosla
      xcorsta(2)= de*cosla-dn*sinphi*sinla
      xcorsta(3)=          dn*cosphi

*** for the semi-diurnal band

      l1=l1sd
      costwola=cosla**2-sinla**2
      sintwola=2.d0*cosla*sinla
      dnsun=-l1/2.d0*sinphi*cosphi*fac2sun*((xsun(1)**2-xsun(2)**2)*
     * costwola+2.d0*xsun(1)*xsun(2)*sintwola)/rsun**2
      dnmon=-l1/2.d0*sinphi*cosphi*fac2mon*((xmon(1)**2-xmon(2)**2)*
     * costwola+2.d0*xmon(1)*xmon(2)*sintwola)/rmon**2
      desun=-l1/2.d0*sinphi**2*cosphi*fac2sun*((xsun(1)**2-xsun(2)**2)*
     * sintwola-2.d0*xsun(1)*xsun(2)*costwola)/rsun**2
      demon=-l1/2.d0*sinphi**2*cosphi*fac2mon*((xmon(1)**2-xmon(2)**2)*
     * sintwola-2.d0*xmon(1)*xmon(2)*costwola)/rmon**2
      de=3.d0*(desun+demon)
      dn=3.d0*(dnsun+dnmon)
      xcorsta(1)=xcorsta(1)-de*sinla-dn*sinphi*cosla
      xcorsta(2)=xcorsta(2)+de*cosla-dn*sinphi*sinla
      xcorsta(3)=xcorsta(3)         +dn*cosphi

      return
      end
*-----------------------------------------------------------------------
      subroutine step2diu(xsta,fhr,t,xcorsta)

*** last change:  vd   17 may 00   1:20 pm
*** these are the subroutines for the step2 of the tidal corrections.
*** they are called to account for the frequency dependence
*** of the love numbers.

      implicit double precision (a-h,o-z)
      double precision xsta(3),xcorsta(3),datdi(9,31)
      double precision deg2rad
      data deg2rad/0.017453292519943295769d0/

*** note, following table is derived from dehanttideinelMJD.f (2000oct30 16:10)
*** has minor differences from that of dehanttideinel.f (2000apr17 14:10)
*** D.M. edited to strictly follow published table 7.5a (2006aug08 13:46)

*** cf. table 7.5a of IERS conventions 2003 (TN.32, pg.82)
*** columns are s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op)
*** units of mm

      data ((datdi(i,j),i=1,9),j=1,31)/
     * -3., 0., 2., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
     * -3., 2., 0., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
     * -2., 0., 1.,-1., 0.,-0.02,-0.01, 0.0 , 0.0,
*****-----------------------------------------------------------------------
****** -2., 0., 1., 0., 0.,-0.08,-0.05, 0.01,-0.02,      !*** original entry
     * -2., 0., 1., 0., 0.,-0.08, 0.00, 0.01, 0.01,      !*** table 7.5a
*****-----------------------------------------------------------------------
     * -2., 2.,-1., 0., 0.,-0.02,-0.01, 0.0 , 0.0,
*****-----------------------------------------------------------------------
****** -1., 0., 0.,-1., 0.,-0.10,-0.05, 0.0 ,-0.02,      !*** original entry
     * -1., 0., 0.,-1., 0.,-0.10, 0.00, 0.00, 0.00,      !*** table 7.5a
*****-----------------------------------------------------------------------
****** -1., 0., 0., 0., 0.,-0.51,-0.26,-0.02,-0.12,      !*** original entry
     * -1., 0., 0., 0., 0.,-0.51, 0.00,-0.02, 0.03,      !*** table 7.5a
*****-----------------------------------------------------------------------
     * -1., 2., 0., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0.,-2., 1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 0.,-1., 0., 0., 0.02, 0.01, 0.0 , 0.0,
*****-----------------------------------------------------------------------
******  0., 0., 1., 0., 0., 0.06, 0.02, 0.0 , 0.01,      !*** original entry
     *  0., 0., 1., 0., 0., 0.06, 0.00, 0.00, 0.00,      !*** table 7.5a
*****-----------------------------------------------------------------------
     *  0., 0., 1., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 2.,-1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  1.,-3., 0., 0., 1.,-0.06, 0.00, 0.00, 0.00,      !*** table 7.5a
     *  1.,-2., 0., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
*****-----------------------------------------------------------------------
******  1.,-2., 0., 0., 0.,-1.23,-0.05, 0.06,-0.06,      !*** original entry
     *  1.,-2., 0., 0., 0.,-1.23,-0.07, 0.06, 0.01,      !*** table 7.5a
*****-----------------------------------------------------------------------
     *  1.,-1., 0., 0.,-1., 0.02, 0.0 , 0.0 , 0.0,
     *  1.,-1., 0., 0., 1., 0.04, 0.0 , 0.0 , 0.0,
     *  1., 0., 0.,-1., 0.,-0.22, 0.01, 0.01, 0.00,      !*** table 7.5a
*****-----------------------------------------------------------------------
******  1., 0., 0., 0., 0.,12.02,-0.45,-0.66, 0.17,      !*** original entry
     *  1., 0., 0., 0., 0.,12.00,-0.78,-0.67,-0.03,      !*** table 7.5a
*****-----------------------------------------------------------------------
******  1., 0., 0., 1., 0., 1.73,-0.07,-0.10, 0.02,      !*** original entry
     *  1., 0., 0., 1., 0., 1.73,-0.12,-0.10, 0.00,      !*** table 7.5a
*****-----------------------------------------------------------------------
     *  1., 0., 0., 2., 0.,-0.04, 0.0 , 0.0 , 0.0,
*****-----------------------------------------------------------------------
******  1., 1., 0., 0.,-1.,-0.50, 0.0 , 0.03, 0.0,       !*** original entry
     *  1., 1., 0., 0.,-1.,-0.50,-0.01, 0.03, 0.00,      !*** table 7.5a
*****-----------------------------------------------------------------------
     *  1., 1., 0., 0., 1., 0.01, 0.0 , 0.0 , 0.0,
*****-----------------------------------------------------------------------
******  0., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,       !*** original entry
     *  1., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,       !*** v.dehant 2007
*****-----------------------------------------------------------------------
     *  1., 2.,-2., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
*****-----------------------------------------------------------------------
******  1., 2., 0., 0., 0.,-0.12, 0.01, 0.01, 0.0,       !*** original entry
     *  1., 2., 0., 0., 0.,-0.11, 0.01, 0.01, 0.00,      !*** table 7.5a
*****-----------------------------------------------------------------------
     *  2.,-2., 1., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
     *  2., 0.,-1., 0., 0.,-0.02, 0.02, 0.0 , 0.01,
     *  3., 0., 0., 0., 0., 0.0 , 0.01, 0.0 , 0.01,
     *  3., 0., 0., 1., 0., 0.0 , 0.01, 0.0 , 0.0/

      s=218.31664563d0+481267.88194d0*t-0.0014663889d0*t*t
     * +0.00000185139d0*t**3
      tau=fhr*15.d0+280.4606184d0+36000.7700536d0*t+0.00038793d0*t*t
     * -0.0000000258d0*t**3-s
      pr=1.396971278*t+0.000308889*t*t+0.000000021*t**3
     * +0.000000007*t**4
      s=s+pr
      h=280.46645d0+36000.7697489d0*t+0.00030322222d0*t*t
     * +0.000000020*t**3-0.00000000654*t**4
      p=83.35324312d0+4069.01363525d0*t-0.01032172222d0*t*t
     * -0.0000124991d0*t**3+0.00000005263d0*t**4
      zns=234.95544499d0 +1934.13626197d0*t-0.00207561111d0*t*t
     * -0.00000213944d0*t**3+0.00000001650d0*t**4
      ps=282.93734098d0+1.71945766667d0*t+0.00045688889d0*t*t
     * -0.00000001778d0*t**3-0.00000000334d0*t**4

*** reduce angles to between 0 and 360

      s=  dmod(  s,360.d0)
      tau=dmod(tau,360.d0)
      h=  dmod(  h,360.d0)
      p=  dmod(  p,360.d0)
      zns=dmod(zns,360.d0)
      ps= dmod( ps,360.d0)

      rsta=dsqrt(xsta(1)**2+xsta(2)**2+xsta(3)**2)
      sinphi=xsta(3)/rsta
      cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta

      cosla=xsta(1)/cosphi/rsta
      sinla=xsta(2)/cosphi/rsta
      zla = datan2(xsta(2),xsta(1))
      do i=1,3
        xcorsta(i)=0.d0
      enddo
      do j=1,31
        thetaf=(tau+datdi(1,j)*s+datdi(2,j)*h+datdi(3,j)*p+
     *   datdi(4,j)*zns+datdi(5,j)*ps)*deg2rad
        dr=datdi(6,j)*2.d0*sinphi*cosphi*sin(thetaf+zla)+
     *     datdi(7,j)*2.d0*sinphi*cosphi*cos(thetaf+zla)
        dn=datdi(8,j)*(cosphi**2-sinphi**2)*sin(thetaf+zla)+
     *     datdi(9,j)*(cosphi**2-sinphi**2)*cos(thetaf+zla)
***** following correction by V.Dehant to match eq.16b, p.81, 2003 Conventions
*****   de=datdi(8,j)*sinphi*cos(thetaf+zla)+
        de=datdi(8,j)*sinphi*cos(thetaf+zla)-
     *     datdi(9,j)*sinphi*sin(thetaf+zla)
        xcorsta(1)=xcorsta(1)+dr*cosla*cosphi-de*sinla
     *   -dn*sinphi*cosla
        xcorsta(2)=xcorsta(2)+dr*sinla*cosphi+de*cosla
     *   -dn*sinphi*sinla
        xcorsta(3)=xcorsta(3)+dr*sinphi+dn*cosphi
      enddo

      do i=1,3
         xcorsta(i)=xcorsta(i)/1000.d0
      enddo

      return
      end
*-----------------------------------------------------------------------
      subroutine step2lon(xsta,fhr,t,xcorsta)

      implicit double precision (a-h,o-z)
      double precision deg2rad
      double precision xsta(3),xcorsta(3),datdi(9,5)
      data deg2rad/0.017453292519943295769d0/

*** cf. table 7.5b of IERS conventions 2003 (TN.32, pg.82)
*** columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)
*** IERS cols.= s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op)
*** units of mm

      data ((datdi(i,j),i=1,9),j=1,5)/
     *   0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07,
     *   0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05,
     *   1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04,
     *   2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07,
     *   2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03/

      s=218.31664563d0+481267.88194d0*t-0.0014663889d0*t*t
     * +0.00000185139d0*t**3
      pr=1.396971278*t+0.000308889*t*t+0.000000021*t**3
     * +0.000000007*t**4
      s=s+pr
      h=280.46645d0+36000.7697489d0*t+0.00030322222d0*t*t
     * +0.000000020*t**3-0.00000000654*t**4
      p=83.35324312d0+4069.01363525d0*t-0.01032172222d0*t*t
     * -0.0000124991d0*t**3+0.00000005263d0*t**4
      zns=234.95544499d0 +1934.13626197d0*t-0.00207561111d0*t*t
     * -0.00000213944d0*t**3+0.00000001650d0*t**4
      ps=282.93734098d0+1.71945766667d0*t+0.00045688889d0*t*t
     * -0.00000001778d0*t**3-0.00000000334d0*t**4
      rsta=dsqrt(xsta(1)**2+xsta(2)**2+xsta(3)**2)
      sinphi=xsta(3)/rsta
      cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
      cosla=xsta(1)/cosphi/rsta
      sinla=xsta(2)/cosphi/rsta

*** reduce angles to between 0 and 360

      s=  dmod(  s,360.d0)
***** tau=dmod(tau,360.d0)       !*** tau not used here--09jul28
      h=  dmod(  h,360.d0)
      p=  dmod(  p,360.d0)
      zns=dmod(zns,360.d0)
      ps= dmod( ps,360.d0)

      dr_tot=0.d0
      dn_tot=0.d0
      do i=1,3
        xcorsta(i)=0.d0
      enddo

***             1 2 3 4   5   6      7      8      9
*** columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)

      do j=1,5
        thetaf=(datdi(1,j)*s+datdi(2,j)*h+datdi(3,j)*p+
     *   datdi(4,j)*zns+datdi(5,j)*ps)*deg2rad
        dr=datdi(6,j)*(3.d0*sinphi**2-1.d0)/2.*cos(thetaf)+
     *       datdi(8,j)*(3.d0*sinphi**2-1.d0)/2.*sin(thetaf)
        dn=datdi(7,j)*(cosphi*sinphi*2.d0)*cos(thetaf)+
     *       datdi(9,j)*(cosphi*sinphi*2.d0)*sin(thetaf)
        de=0.d0
        dr_tot=dr_tot+dr
        dn_tot=dn_tot+dn
        xcorsta(1)=xcorsta(1)+dr*cosla*cosphi-de*sinla
     *   -dn*sinphi*cosla
        xcorsta(2)=xcorsta(2)+dr*sinla*cosphi+de*cosla
     *   -dn*sinphi*sinla
        xcorsta(3)=xcorsta(3)+dr*sinphi+dn*cosphi
      enddo

      do i=1,3
        xcorsta(i)=xcorsta(i)/1000.d0
      enddo

      return
      end
*-----------------------------------------------------------------------
      subroutine st1idiu(xsta,xsun,xmon,fac2sun,fac2mon,xcorsta)

*** this subroutine gives the out-of-phase corrections induced by
*** mantle inelasticity in the diurnal band

***  input: xsta,xsun,xmon,fac2sun,fac2mon
*** output: xcorsta

      implicit double precision (a-h,o-z)
      dimension xsta(3),xsun(3),xmon(3),xcorsta(3)
      data dhi/-0.0025d0/,dli/-0.0007d0/

      rsta=enorm8(xsta)
      sinphi=xsta(3)/rsta
      cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
      cos2phi=cosphi**2-sinphi**2
      sinla=xsta(2)/cosphi/rsta
      cosla=xsta(1)/cosphi/rsta
      rmon=enorm8(xmon)
      rsun=enorm8(xsun)
      drsun=-3.d0*dhi*sinphi*cosphi*fac2sun*xsun(3)*(xsun(1)*
     *            sinla-xsun(2)*cosla)/rsun**2
      drmon=-3.d0*dhi*sinphi*cosphi*fac2mon*xmon(3)*(xmon(1)*
     *            sinla-xmon(2)*cosla)/rmon**2
      dnsun=-3.d0*dli*cos2phi*fac2sun*xsun(3)*(xsun(1)*sinla-
     *            xsun(2)*cosla)/rsun**2
      dnmon=-3.d0*dli*cos2phi*fac2mon*xmon(3)*(xmon(1)*sinla-
     *            xmon(2)*cosla)/rmon**2
      desun=-3.d0*dli*sinphi*fac2sun*xsun(3)*
     * (xsun(1)*cosla+xsun(2)*sinla)/rsun**2
      demon=-3.d0*dli*sinphi*fac2mon*xmon(3)*
     * (xmon(1)*cosla+xmon(2)*sinla)/rmon**2
      dr=drsun+drmon
      dn=dnsun+dnmon
      de=desun+demon
      xcorsta(1)=dr*cosla*cosphi-de*sinla-dn*sinphi*cosla
      xcorsta(2)=dr*sinla*cosphi+de*cosla-dn*sinphi*sinla
      xcorsta(3)=dr*sinphi               +dn*cosphi

      return
      end
*-----------------------------------------------------------------------
      subroutine st1isem(xsta,xsun,xmon,fac2sun,fac2mon,xcorsta)

*** this subroutine gives the out-of-phase corrections induced by
*** mantle inelasticity in the diurnal band

***  input: xsta,xsun,xmon,fac2sun,fac2mon
*** output: xcorsta

      implicit double precision (a-h,o-z)
      dimension xsta(3),xsun(3),xmon(3),xcorsta(3)
      data dhi/-0.0022d0/,dli/-0.0007d0/

      rsta=enorm8(xsta)
      sinphi=xsta(3)/rsta
      cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
      sinla=xsta(2)/cosphi/rsta
      cosla=xsta(1)/cosphi/rsta
      costwola=cosla**2-sinla**2
      sintwola=2.d0*cosla*sinla
      rmon=enorm8(xmon)
      rsun=enorm8(xsun)
      drsun=-3.d0/4.d0*dhi*cosphi**2*fac2sun*((xsun(1)**2-xsun(2)**2)*
     * sintwola-2.*xsun(1)*xsun(2)*costwola)/rsun**2
      drmon=-3.d0/4.d0*dhi*cosphi**2*fac2mon*((xmon(1)**2-xmon(2)**2)*
     * sintwola-2.*xmon(1)*xmon(2)*costwola)/rmon**2
      dnsun=1.5d0*dli*sinphi*cosphi*fac2sun*((xsun(1)**2-xsun(2)**2)*
     * sintwola-2.d0*xsun(1)*xsun(2)*costwola)/rsun**2
      dnmon=1.5d0*dli*sinphi*cosphi*fac2mon*((xmon(1)**2-xmon(2)**2)*
     * sintwola-2.d0*xmon(1)*xmon(2)*costwola)/rmon**2
      desun=-3.d0/2.d0*dli*cosphi*fac2sun*((xsun(1)**2-xsun(2)**2)*
     * costwola+2.*xsun(1)*xsun(2)*sintwola)/rsun**2
      demon=-3.d0/2.d0*dli*cosphi*fac2mon*((xmon(1)**2-xmon(2)**2)*
     * costwola+2.d0*xmon(1)*xmon(2)*sintwola)/rmon**2
      dr=drsun+drmon
      dn=dnsun+dnmon
      de=desun+demon
      xcorsta(1)=dr*cosla*cosphi-de*sinla-dn*sinphi*cosla
      xcorsta(2)=dr*sinla*cosphi+de*cosla-dn*sinphi*sinla
      xcorsta(3)=dr*sinphi+dn*cosphi

      return
      end
*-----------------------------------------------------------------------
      subroutine sprod(x,y,scal,r1,r2)

***  computation of the scalar-product of two vectors and their norms

***  input:   x(i),i=1,2,3  -- components of vector x
***           y(i),i=1,2,3  -- components of vector y
***  output:  scal          -- scalar product of x and y
***           r1,r2         -- lengths of the two vectors x and y

      implicit double precision (a-h,o-z)
      double precision x(3),y(3)

      r1=dsqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
      r2=dsqrt(y(1)*y(1) + y(2)*y(2) + y(3)*y(3))
      scal=    x(1)*y(1) + x(2)*y(2) + x(3)*y(3)

      return
      end
*-----------------------------------------------------------------------
      double precision function enorm8(a)

*** compute euclidian norm of a vector (of length 3)

      double precision a(3)

      enorm8=dsqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

      return
      end
*-----------------------------------------------------------------------
      subroutine zero_vec8(v)

*** initialize a vector (of length 3) to zero

      double precision v(3)

      v(1)=0.d0
      v(2)=0.d0
      v(3)=0.d0

      return
      end
*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine moonxyz(mjd,fmjd,rm)

*** get low-precision, geocentric coordinates for moon (ECEF)

*** input:  mjd/fmjd, is Modified Julian Date (and fractional) in GPS time
*** output: rm, is geocentric lunar position vector [m] in ECEF
*** 1."satellite orbits: models, methods, applications" montenbruck & gill(2000)
*** section 3.3.2, pg. 72-73
*** 2."astronomy on the personal computer, 4th ed." montenbruck & pfleger (2005)
*** section 3.2, pg. 38-39  routine MiniMoon

      implicit double precision(a-h,o-z)
      dimension rm(3)

      common/stuff/rad,pi,pi2
      call initialize()

*** use TT for lunar ephemerides

      tsecgps=fmjd*86400.d0                       !*** GPS time (sec of day)
      tsectt =gps2tt(tsecgps)                     !*** TT  time (sec of day)
      fmjdtt =tsectt/86400.d0                     !*** TT  time (fract. day)

*** julian centuries since 1.5 january 2000 (J2000)
***   (note: also low precision use of mjd --> tjd)

      tjdtt = mjd+fmjdtt+2400000.5d0              !*** Julian Date, TT
      t     = (tjdtt - 2451545.d0)/36525.d0       !*** julian centuries, TT

*** el0 -- mean longitude of Moon (deg)
*** el  -- mean anomaly of Moon (deg)
*** elp -- mean anomaly of Sun  (deg)
*** f   -- mean angular distance of Moon from ascending node (deg)
*** d   -- difference between mean longitudes of Sun and Moon (deg)

*** equations 3.47, p.72

      el0=218.31617d0 + 481267.88088d0*t -1.3972*t
      el =134.96292d0 + 477198.86753d0*t
      elp=357.52543d0 +  35999.04944d0*t
      f  = 93.27283d0 + 483202.01873d0*t
      d  =297.85027d0 + 445267.11135d0*t

*** longitude w.r.t. equinox and ecliptic of year 2000

      selond=el0                                      !*** eq 3.48, p.72
     * +22640.d0/3600.d0*dsin((el        )/rad)
     * +  769.d0/3600.d0*dsin((el+el     )/rad)
     * - 4586.d0/3600.d0*dsin((el-d-d    )/rad)
     * + 2370.d0/3600.d0*dsin((d+d       )/rad)
     * -  668.d0/3600.d0*dsin((elp       )/rad)
     * -  412.d0/3600.d0*dsin((f+f       )/rad)
     * -  212.d0/3600.d0*dsin((el+el-d-d )/rad)
     * -  206.d0/3600.d0*dsin((el+elp-d-d)/rad)
     * +  192.d0/3600.d0*dsin((el+d+d    )/rad)
     * -  165.d0/3600.d0*dsin((elp-d-d   )/rad)
     * +  148.d0/3600.d0*dsin((el-elp    )/rad)
     * -  125.d0/3600.d0*dsin((d         )/rad)
     * -  110.d0/3600.d0*dsin((el+elp    )/rad)
     * -   55.d0/3600.d0*dsin((f+f-d-d   )/rad)

*** latitude w.r.t. equinox and ecliptic of year 2000

      q = 412.d0/3600.d0*dsin((f+f)/rad)              !*** temporary term
     *   +541.d0/3600.d0*dsin((elp)/rad)

      selatd=                                         !*** eq 3.49, p.72
     * +18520.d0/3600.d0*dsin((f+selond-el0+q)/rad)
     * -  526.d0/3600.d0*dsin((f-d-d     )/rad)
     * +   44.d0/3600.d0*dsin((el+f-d-d  )/rad)
     * -   31.d0/3600.d0*dsin((-el+f-d-d )/rad)
     * -   25.d0/3600.d0*dsin((-el-el+f  )/rad)
     * -   23.d0/3600.d0*dsin((elp+f-d-d )/rad)
     * +   21.d0/3600.d0*dsin((-el+f     )/rad)
     * +   11.d0/3600.d0*dsin((-elp+f-d-d)/rad)

*** distance from Earth center to Moon (m)

      rse= 385000.d0*1000.d0                          !*** eq 3.50, p.72
     *   -  20905.d0*1000.d0*dcos((el        )/rad)
     *   -   3699.d0*1000.d0*dcos((d+d-el    )/rad)
     *   -   2956.d0*1000.d0*dcos((d+d       )/rad)
     *   -    570.d0*1000.d0*dcos((el+el     )/rad)
     *   +    246.d0*1000.d0*dcos((el+el-d-d )/rad)
     *   -    205.d0*1000.d0*dcos((elp-d-d   )/rad)
     *   -    171.d0*1000.d0*dcos((el+d+d    )/rad)
     *   -    152.d0*1000.d0*dcos((el+elp-d-d)/rad)

*** convert spherical ecliptic coordinates to equatorial cartesian

*** precession of equinox wrt. J2000   (p.71)

      selond=selond + 1.3972d0*t                         !*** degrees

*** position vector of moon (mean equinox & ecliptic of J2000) (EME2000, ICRF)
***                         (plus long. advance due to precession -- eq. above)

      oblir=23.43929111d0/rad        !*** obliquity of the J2000 ecliptic

      sselat=dsin(selatd/rad)
      cselat=dcos(selatd/rad)
      sselon=dsin(selond/rad)
      cselon=dcos(selond/rad)

      t1 = rse*cselon*cselat        !*** meters          !*** eq. 3.51, p.72
      t2 = rse*sselon*cselat        !*** meters          !*** eq. 3.51, p.72
      t3 = rse*       sselat        !*** meters          !*** eq. 3.51, p.72

      call rot1(-oblir,t1,t2,t3,rm1,rm2,rm3)             !*** eq. 3.51, p.72

*** convert position vector of moon to ECEF  (ignore polar motion/LOD)

      call getghar(mjd,fmjd,ghar)                        !*** sec 2.3.1,p.33
      call rot3(ghar,rm1,rm2,rm3,rm(1),rm(2),rm(3))      !*** eq. 2.89, p.37

      return
      end
********************************************************************************
      subroutine getghar(mjd,fmjd,ghar)

*** convert mjd/fmjd in GPS time to Greenwich hour angle (in radians)

*** "satellite orbits: models, methods, applications" montenbruck & gill(2000)
*** section 2.3.1, pg. 33

      implicit double precision(a-h,o-z)

      common/stuff/rad,pi,pi2
      call initialize()

*** need UT to get sidereal time  ("astronomy on the personal computer", 4th ed)
***                               (pg.43, montenbruck & pfleger, springer, 2005)

      tsecgps=fmjd*86400.d0                        !*** GPS time (sec of day)
      tsecutc=gps2utc(tsecgps)                     !*** UTC time (sec of day)
      fmjdutc=tsecutc/86400.d0                     !*** UTC time (fract. day)

***** d = MJD - 51544.5d0                               !*** footnote
      d =(mjd-51544) + (fmjdutc-0.5d0)                  !*** days since J2000

*** greenwich hour angle for J2000  (12:00:00 on 1 Jan 2000)

***** ghad = 100.46061837504d0 + 360.9856473662862d0*d  !*** eq. 2.85 (+digits)
      ghad = 280.46061837504d0 + 360.9856473662862d0*d  !*** corrn.   (+digits)

**** normalize to 0-360 and convert to radians

      i    = ghad/360.d0
      ghar =(ghad-i*360.d0)/rad
    1 if(ghar.gt.pi2) then
        ghar=ghar-pi2
        go to 1
      endif
    2 if(ghar.lt.0.d0) then
        ghar=ghar+pi2
        go to 2
      endif

      return
      end
********************************************************************************
      subroutine sunxyz(mjd,fmjd,rs)

*** get low-precision, geocentric coordinates for sun (ECEF)

*** input, mjd/fmjd, is Modified Julian Date (and fractional) in GPS time
*** output, rs, is geocentric solar position vector [m] in ECEF
*** 1."satellite orbits: models, methods, applications" montenbruck & gill(2000)
*** section 3.3.2, pg. 70-71
*** 2."astronomy on the personal computer, 4th ed." montenbruck & pfleger (2005)
*** section 3.2, pg. 39  routine MiniSun

      implicit double precision(a-h,o-z)
      dimension rs(3)

      common/stuff/rad,pi,pi2
      call initialize()

*** mean elements for year 2000, sun ecliptic orbit wrt. Earth

      obe =23.43929111d0/rad        !*** obliquity of the J2000 ecliptic
      sobe=dsin(obe)
      cobe=dcos(obe)
      opod=282.9400d0               !*** RAAN + arg.peri.  (deg.)

*** use TT for solar ephemerides

      tsecgps=fmjd*86400.d0                       !*** GPS time (sec of day)
      tsectt =gps2tt(tsecgps)                     !*** TT  time (sec of day)
      fmjdtt =tsectt/86400.d0                     !*** TT  time (fract. day)

*** julian centuries since 1.5 january 2000 (J2000)
***   (note: also low precision use of mjd --> tjd)

      tjdtt = mjd+fmjdtt+2400000.5d0              !*** Julian Date, TT
      t     = (tjdtt - 2451545.d0)/36525.d0       !*** julian centuries, TT
      emdeg = 357.5256d0 + 35999.049d0*t          !*** degrees
      em    = emdeg/rad                           !*** radians
      em2   = em+em                               !*** radians

*** series expansions in mean anomaly, em   (eq. 3.43, p.71)

      r=(149.619d0-2.499d0*dcos(em)-0.021d0*dcos(em2))*1.d9      !*** m.
      slond=opod + emdeg + (6892.d0*dsin(em)+72.d0*dsin(em2))/3600.d0

*** precession of equinox wrt. J2000   (p.71)

      slond=slond + 1.3972d0*t                              !*** degrees

*** position vector of sun (mean equinox & ecliptic of J2000) (EME2000, ICRF)
***                        (plus long. advance due to precession -- eq. above)

      slon =slond/rad                                       !*** radians
      sslon=dsin(slon)
      cslon=dcos(slon)

      rs1 = r*cslon              !*** meters             !*** eq. 3.46, p.71
      rs2 = r*sslon*cobe         !*** meters             !*** eq. 3.46, p.71
      rs3 = r*sslon*sobe         !*** meters             !*** eq. 3.46, p.71

*** convert position vector of sun to ECEF  (ignore polar motion/LOD)

      call getghar(mjd,fmjd,ghar)                        !*** sec 2.3.1,p.33
      call rot3(ghar,rs1,rs2,rs3,rs(1),rs(2),rs(3))      !*** eq. 2.89, p.37

      return
      end
********************************************************************************
      subroutine lhsaaz(u,v,w,ra,az,va)

*** determine range,azimuth,vertical angle from local horizon coord.

      implicit double precision(a-h,o-z)
      
      s2=u*u+v*v
      r2=s2+w*w
      
      s =dsqrt(s2)
      ra=dsqrt(r2)
      
      az=datan2(v,u)
      va=datan2(w,s)
      
      return
      end
*-----------------------------------------------------------------------
      subroutine geoxyz(gla,glo,eht,x,y,z)

*** convert geodetic lat, long, ellip ht. to x,y,z

      implicit double precision(a-h,o-z)

      common/comgrs/a,e2
      call initialize()

      sla=dsin(gla)
      cla=dcos(gla)
      w2=1.d0-e2*sla*sla
      w=dsqrt(w2)
      en=a/w

      x=(en+eht)*cla*dcos(glo)
      y=(en+eht)*cla*dsin(glo)
      z=(en*(1.d0-e2)+eht)*sla

      return
      end
*-----------------------------------------------------------------------
      subroutine rge(gla,glo,u,v,w,x,y,z)

*** given a rectangular cartesian system (x,y,z)
*** compute a geodetic h cartesian sys   (u,v,w)

      implicit double precision(a-h,o-z)

      sb=dsin(gla)
      cb=dcos(gla)
      sl=dsin(glo)
      cl=dcos(glo)

      u=-sb*cl*x-sb*sl*y+cb*z
      v=-   sl*x+   cl*y
      w= cb*cl*x+cb*sl*y+sb*z

      return
      end
*-----------------------------------------------------------------------
      subroutine rot1(theta,x,y,z,u,v,w)
      
*** rotate coordinate axes about 1 axis by angle of theta radians
*** x,y,z transformed into u,v,w

      implicit double precision(a-h,o-z)
      
      s=dsin(theta)
      c=dcos(theta)
            
      u=x
      v=c*y+s*z
      w=c*z-s*y
            
      return
      end
*-----------------------------------------------------------------------
      subroutine rot3(theta,x,y,z,u,v,w)
      
*** rotate coordinate axes about 3 axis by angle of theta radians
*** x,y,z transformed into u,v,w

      implicit double precision(a-h,o-z)
      
      s=dsin(theta)
      c=dcos(theta)
            
      u=c*x+s*y
      v=c*y-s*x
      w=z
            
      return
      end
************************************************************************
*** time conversion ****************************************************
************************************************************************
      subroutine setjd0(iyr,imo,idy)

*** set the integer part of a modified julian date as epoch, mjd0
*** the modified julian day is derived from civil time as in civmjd()
*** allows single number expression of time in seconds w.r.t. mjd0

      implicit double precision(a-h,o-z)
      integer y
      save /mjdoff/
      common/mjdoff/mjd0

      if(iyr.lt.1900) then
*stop 34587
        ierr = 34587
        goto 400
      end if

      if(imo.le.2) then
        y=iyr-1
        m=imo+12
      else
        y=iyr
        m=imo
      endif

      it1=365.25d0*y
      it2=30.6001d0*(m+1)
      mjd=it1+it2+idy-679019

*** now set the epoch for future time computations

      mjd0=mjd

      return

*** error - return
  400 return

      end
      subroutine civjts(iyr,imo,idy,ihr,imn,sec,tsec)

*** convert civil date to time in seconds past mjd epoch, mjd0
*** requires initialization of mjd0 by setjd0()

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100     (leap year protocols)
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** adapted from civmjd()

      implicit double precision(a-h,o-z)
      integer y
      save /mjdoff/
      common/mjdoff/mjd0

      if(iyr.lt.1900) then
*stop 34589
        ierr = 34589
        goto 400
      end if

      if(imo.le.2) then
        y=iyr-1
        m=imo+12
      else
        y=iyr
        m=imo
      endif

      it1=365.25d0*y
      it2=30.6001d0*(m+1)
      mjd=it1+it2+idy-679019

      tsec=(mjd-mjd0)*86400.d0+3600*ihr+60*imn+sec

      return

*** error - return
  400 return

      end
      subroutine jtsciv(tsec,iyr,imo,idy,ihr,imn,sec)

*** convert time in seconds past mjd0 epoch into civil date
*** requires initialization of mjd0 by setjd0()

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** adapted from mjdciv()

      implicit double precision(a-h,o-z)
      save /mjdoff/
      common/mjdoff/mjd0

      mjd=mjd0+tsec/86400.d0
*** the following equation preserves significant digits
      fmjd=dmod(tsec,86400.d0)/86400.d0

      rjd=mjd+fmjd+2400000.5d0
      ia=(rjd+0.5d0)
      ib=ia+1537
      ic=(ib-122.1d0)/365.25d0
      id=365.25d0*ic
      ie=(ib-id)/30.6001d0

*** the fractional part of a julian day is (fractional mjd + 0.5)
*** therefore, fractional part of julian day + 0.5 is (fractional mjd)

      it1=ie*30.6001d0
      idy=ib-id-it1+fmjd
      it2=ie/14.d0
      imo=ie-1-12*it2
      it3=(7+imo)/10.d0
      iyr=ic-4715-it3

      tmp=fmjd*24.d0
      ihr=tmp
      tmp=(tmp-ihr)*60.d0
      imn=tmp
      sec=(tmp-imn)*60.d0

      return
      end
************************************************************************
      subroutine civmjd(iyr,imo,idy,ihr,imn,sec,mjd,fmjd)

*** convert civil date to modified julian date

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100     (leap year protocols)
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** operation confirmed against table 3.3 values on pg.34

      implicit double precision(a-h,o-z)
      integer y

      if(iyr.lt.1900)  then
*stop 34588
        ierr = 34588
        goto 400
      end if

      if(imo.le.2) then
        y=iyr-1
        m=imo+12
      else
        y=iyr
        m=imo
      endif

      it1=365.25d0*y
      it2=30.6001d0*(m+1)
      mjd=it1+it2+idy-679019

      fmjd=(3600*ihr+60*imn+sec)/86400.d0

      return

*** error - return
  400 return

      end
      subroutine mjdciv(mjd,fmjd,iyr,imo,idy,ihr,imn,sec)

*** convert modified julian date to civil date

*** imo in range 1-12, idy in range 1-31
*** only valid in range mar-1900 thru feb-2100
*** ref: hofmann-wellenhof, 2nd ed., pg 34-35
*** operation confirmed for leap years (incl. year 2000)

      implicit double precision(a-h,o-z)

      rjd=mjd+fmjd+2400000.5d0
      ia=(rjd+0.5d0)
      ib=ia+1537
      ic=(ib-122.1d0)/365.25d0
      id=365.25d0*ic
      ie=(ib-id)/30.6001d0

*** the fractional part of a julian day is fractional mjd + 0.5
*** therefore, fractional part of julian day + 0.5 is fractional mjd

      it1=ie*30.6001d0
      idy=ib-id-it1+fmjd
      it2=ie/14.d0
      imo=ie-1-12*it2
      it3=(7+imo)/10.d0
      iyr=ic-4715-it3

      tmp=fmjd*24.d0
      ihr=tmp
      tmp=(tmp-ihr)*60.d0
      imn=tmp
      sec=(tmp-imn)*60.d0

      return
      end
************************************************************************
*** supplemental time functions ****************************************
************************************************************************
      double precision function gps2tt(tsec)

*** convert tsec in GPS to tsec in TT

      implicit double precision(a-h,o-z)

      gps2tt=tsec+51.184d0                          !*** fixed offset

      return
      end
      double precision function gps2utc(tsec)

*** convert tsec in GPS to tsec in UTC

      implicit double precision(a-h,o-z)

*** GPS is ahead of UTC  (c.f. USNO)
*** UTC is behind GPS
*** gpsleap() is (so far) positive (and increasing)
*** so, must subtract gpsleap from GPS to get UTC

      gps2utc=tsec-gpsleap(tsec)

      return
      end
      double precision function gpsleap(tsec)

*** return total leap seconds since GPS epoch 1980jan06

*** note: does **NOT** return the full TAI-UTC delta
*** input time is GPS seconds -- initialized by setjd0()
*** Y2K -- only functional between 1980jan06-00:00:00  (GPS time start)
***                            and hard-coded date

      implicit double precision(a-h,o-z)

***** "Julian Date Converter"
***** http://aa.usno.navy.mil/data/docs/JulianDate.php
***** "Bulletin C"
***** http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat
***** parameter(mjdhard=55196)            !*** cut-off date 2009dec31
***** parameter(mjdhard=55377)            !*** cut-off date 2010jun30
***** parameter(mjdhard=55561)            !*** cut-off date 2010dec31
***** parameter(mjdhard=55742)            !*** cut-off date 2011jun30
***** parameter(mjdhard=55926)            !*** cut-off date 2011dec31
***** parameter(mjdhard=56108)            !*** cut-off date 2012jun30
***** parameter(mjdhard=56292)            !*** cut-off date 2012dec31
***** parameter(mjdhard=56473)            !*** cut-off date 2013jun30
***** parameter(mjdhard=56657)            !*** cut-off date 2013dec31
***** parameter(mjdhard=56838)            !*** cut-off date 2014jun30
***** parameter(mjdhard=57022)            !*** cut-off date 2014dec31
***** parameter(mjdhard=57203)            !*** cut-off date 2015jun30
***** parameter(mjdhard=57387)            !*** cut-off date 2015dec31
      parameter(mjdhard=57569)            !*** cut-off date 2016jun30

      save  /mjdoff/
      common/mjdoff/mjd0

*** clone for tests (and do any rollover)

      ttsec=tsec
      mjd0t=mjd0

    1 if(ttsec.ge.86400.d0) then
        ttsec=ttsec-86400.d0
        mjd0t=mjd0t+1
        go to 1
      endif

    2 if(ttsec.lt.0.d0) then
        ttsec=ttsec+86400.d0
        mjd0t=mjd0t-1
        go to 2
      endif

*** test date limits

      if(mjd0t.gt.mjdhard) then
*        write(*,*) 'FATAL ERROR --'
*        write(*,*) 'exceeded cut-off date in gpsleap()'
*        stop 66766
        ierr = 66766
        goto 400
      endif

      if(mjd0t.lt.44244) then             !*** 1980jan06
*        write(*,*) 'FATAL ERROR --'
*        write(*,*) 'cut-off date underflow in gpsleap()'
*        stop 66767
        ierr = 66767
        goto 400
      endif

*** http://maia.usno.navy.mil/ser7/tai-utc.dat
*** 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0s
*** 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0s
*** 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0s
*** 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0s
*** 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0s
*** 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0s
*** 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0s
*** 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0s
*** 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0s
*** 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0s
*** 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0s
*** 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0s
*** 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0s
*** 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0s
*** 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0s
*** 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0s
*** 2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0s
*** 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0s

*** test against newest leaps first

      if    (mjd0t.ge.57204) then       !*** 2015 JUL 1 = 57204
        tai_utc = 36.d0
      elseif(mjd0t.ge.56109) then       !*** 2012 JUL 1 = 56109
        tai_utc = 35.d0
      elseif(mjd0t.ge.54832) then       !*** 2009 JAN 1 = 54832
        tai_utc = 34.d0
      elseif(mjd0t.ge.53736) then       !*** 2006 JAN 1 = 53736
        tai_utc = 33.d0
      elseif(mjd0t.ge.51179) then       !*** 1999 JAN 1 = 51179
        tai_utc = 32.d0
      elseif(mjd0t.ge.50630) then       !*** 1997 JUL 1 = 50630
        tai_utc = 31.d0
      elseif(mjd0t.ge.50083) then       !*** 1996 JAN 1 = 50083
        tai_utc = 30.d0
      elseif(mjd0t.ge.49534) then       !*** 1994 JUL 1 = 49534
        tai_utc = 29.d0
      elseif(mjd0t.ge.49169) then       !*** 1993 JUL 1 = 49169
        tai_utc = 28.d0
      elseif(mjd0t.ge.48804) then       !*** 1992 JUL 1 = 48804
        tai_utc = 27.d0
      elseif(mjd0t.ge.48257) then       !*** 1991 JAN 1 = 48257
        tai_utc = 26.d0
      elseif(mjd0t.ge.47892) then       !*** 1990 JAN 1 = 47892
        tai_utc = 25.d0
      elseif(mjd0t.ge.47161) then       !*** 1988 JAN 1 = 47161
        tai_utc = 24.d0
      elseif(mjd0t.ge.46247) then       !*** 1985 JUL 1 = 46247
        tai_utc = 23.d0
      elseif(mjd0t.ge.45516) then       !*** 1983 JUL 1 = 45516
        tai_utc = 22.d0
      elseif(mjd0t.ge.45151) then       !*** 1982 JUL 1 = 45151
        tai_utc = 21.d0
      elseif(mjd0t.ge.44786) then       !*** 1981 JUL 1 = 44786
        tai_utc = 20.d0
      elseif(mjd0t.ge.44239) then       !*** 1980 JAN 1 = 44239
        tai_utc = 19.d0

*** should never get here

      else
*        write(*,*) 'FATAL ERROR --'
*        write(*,*) 'fell thru tests in gpsleap()'
*        stop 66768
      ierr = 66768
      goto 400
      endif

*** convert TAI-UTC into GPS leap seconds

      gpsleap = tai_utc - 19.d0

*      return

*** error - return
  400 return
      end

