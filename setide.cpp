/*------------------------------------------------------------------------*
 NAME:     setide.cpp

 PURPOSE:  Calculates 3-dimensional solid earth tide displacement vector
           given MJD and station coordinates.

 AUTHOR:   John Gary Sonntag

 DATE:     24 April 2020
 *------------------------------------------------------------------------*/

#include "/home/sonntag/Include/mission.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define C 299792458.0
#define WE 7.292115147e-5
#define MU 3.986005005e14
#define PI (4.0*atan((double)(1.0)))
#define DEG2RAD (PI/180.0)
#define AE 6378137.0
#define FLAT (1.0/298.257223563)
#define RAD2NM (180.0*60.0/PI)
#define RAD2KM (RAD2NM*6076.1*12.0*2.54/100.0/1000.0)


main(int argc, char *argv[])
{
  int mjd,i;
  double mjs,lat,lon,ht;
  double rsta[3],perm_tide[3],tv_tide[3],total_tide[3];
  void setide(int,double,int,double *,double *,double *,double *);
  void geod2cart(double *,double,double,double,double,double);

  //  Check input
  if ( argc != 1 )
  {
    printf("Usage:  setide\n");
    exit(-1);
  }

  //rsta[0] =   538117.2;
  //rsta[1] = -1389036.5;
  //rsta[2] =  6180994.1;

  mjd = 57803;
  mjs = 64800;
  lat = 0.0;
  lon = 0.0;
  ht  = 41.3666;
  for (i=-90;i<=90;i++)
  {
    lat = double(i);
    geod2cart(rsta,lat*DEG2RAD,lon*DEG2RAD,ht,AE,FLAT);
    setide(mjd,mjs,2,rsta,perm_tide,tv_tide,total_tide);
    printf("%lf %lf %lf %lf\n",lat,perm_tide[0],perm_tide[1],perm_tide[2]);
  }
  //printf("x: %7.4lf %7.4lf %7.4lf\n",tv_tide[0],perm_tide[0],total_tide[0]);
  //printf("y: %7.4lf %7.4lf %7.4lf\n",tv_tide[1],perm_tide[1],total_tide[1]);
  //printf("z: %7.4lf %7.4lf %7.4lf\n",tv_tide[2],perm_tide[2],total_tide[2]);

}


void geod2cart(double *rsta,double lat,double lon,double ht,double ae,double flat)
{
  double denom,slat,clat;

  slat = sin(lat);
  clat = cos(lat);
  denom = sqrt(1.0 - flat*(2.0-flat)*slat*slat);
  rsta[0] = (ht + ae/denom)*clat*cos(lon);
  rsta[1] = (ht + ae/denom)*clat*sin(lon);
  rsta[2] = (ht + ae*(1.0-flat)*(1.0-flat)/denom)*slat;

}


/*------------------------------------------------------------------------*
 NAME:      setide(int mjd, double mjs, int tid, double *rsta, 
            double *perm_tide, double *tv_tide,double *total_tide)

 PURPOSE:   Calculate 3-dimensional solid earth tide deformation dtide[3] 
            at time mjd/mjs and at station coordinates sta[3].  User 
            specifies whether input time is GPS or UTC with
            tid = 1: specified time is UTC
            tid = 2: specified time is GPS

 INPUTS:    mjd: Modified Julian Day (interpretation subject to tid)
            mjs: seconds of MJD (interpretation subject to tid)
            tid: timescale flag (see above)
            rsta: ITRF coordinates of crustal location of interest (m)

 OUTPUTS:   all outputs are in ECEF/ITRF XYZ (*not* local north-east-up)
            note that these are algebraically-related:
            total_tide = tv_tide + perm_tide
            perm_tide[0-2]: permanent component of the tidal deformation (m)
            tv_tide[0-2]: time-varying component of the tidal deformation (m)
            total_tide[0-2]: total tidal deformation (m)

 AUTHOR:    John Gary Sonntag

 DATE:      24 April 2020
 *------------------------------------------------------------------------*/
void setide(int mjdin,double mjsin,int tid,double *rsta,
            double *perm_tide,double *tv_tide,double *total_tide)
{
  int mjd,i;
  double leapsecs,mjs,theta,shdg,mhdg;
  double rs[3],rm[3];
  void sunxyz(int,double,double,double *),moonxyz(int,double,double,double *);
  double getleapsecs(int,double,int);
  void detide(double *,int,double,double *,double *,double *,double *);

  // Determine GPS-UTC leapseconds, and MJD/MJS in GPS time scale
  leapsecs = getleapsecs(mjdin,mjsin,tid);
  if (tid!=2)  // input time is UTC
  {
    mjd = mjdin;
    mjs = mjsin+leapsecs;
    while (mjs>=86400.0)
    {
      mjd+=1;
      mjs-=86400.0;
    }
  }
  else  // input time is GPS
  {
    mjd = mjdin;
    mjs = mjsin;
  }

  // Calculate sun and moon position vectors (ECEF)
  sunxyz(mjd,mjs,leapsecs,rs);
  moonxyz(mjd,mjs,leapsecs,rm);
  theta = atan2(rs[1],rs[0]);
  shdg  = ((PI/2.0) - theta)/DEG2RAD;
  theta = atan2(rm[1],rm[0]);
  mhdg  = ((PI/2.0) - theta)/DEG2RAD;
  //printf("sun hdg = %lf deg  moon hdg = %lf deg\n",shdg,mhdg);

  // Calculate the instantaneous 3-d earth tide
  // routine gives total and permanent components
  // so we calculate the time-varying component here
  detide(rsta,mjd,mjs,rs,rm,total_tide,perm_tide);
  for (i=0;i<3;i++)
    tv_tide[i] = total_tide[i] - perm_tide[i];

}



/* This routine adapted from the ICESat-2 production code in Fortran-90
! computation of tidal corrections of station displacements caused
! by lunar and solar gravitational attraction
!
!### References
! * author iers 1996 :  v. dehant, s. mathews and j. gipson
!    (test between two subroutines)
!
! * author iers 2000 :  v. dehant, c. bruyninx and s. mathews
!    (test in the bernese program by c. bruyninx)
!
! * created:  96/03/23 (see above)
!   modified from dehanttideinelMJD.f by Dennis Milbert 2006sep10
!   bug fix regarding fhr (changed calling sequence, too)
!   modified to reflect table 7.5a and b IERS Conventions 2003
!   modified to use TT time system to call step 2 functions
!   sign correction by V.Dehant to match eq.16b, p.81, Conventions
!   applied by Dennis Milbert 2007may05
!
!### Comments
! * step 1 (here general degree 2 and 3 corrections +
!         call st1idiu + call st1isem + call st1l1)
!   + step 2 (call step2diu + call step2lon + call step2idiu)
!
! * it has been decided that the step 3 un-correction for permanent tide
!   would *not* be applied in order to avoid jump in the reference frame
!   (this step 3 must added in order to get the mean tide station position
!   and to be conformed with the iag resolution.)
!
! * inputs
!    rsta[i],i=0,1,2   -- geocentric position of the station (ITRF/ECEF)
!    rsun[i],i=0,1,2   -- geoc. position of the sun (ECEF)
!    rmon[i],i=0,1,2   -- geoc. position of the moon (ECEF)
!    mjd,mjs           -- modified julian day and seconds of day (in GPS time)
!
! * outputs
!    dtide[i],i=0,1,2  -- total displacement vector (ITRF xyz NOT local frame)
!    ptide[i],i=0,1,2  -- permanent component of displacement vector 
!                         (ITRF xyz NOT local frame)
*/

void detide(double *rsta,int mjd,double mjs,double *rsun,double *rmon,double *dtide,double *ptide)
{
  int mjdtt,i;
  double mjstt,dmjdtt,t,fhr,scs,scm,lsta,lsun,lmon;
  double scsun,scmon,sinphi,cosphi,h2,l2,h20,l20,p2sun,p2mon;
  double p3sun,p3mon,h3,l3,x2sun,x2mon,x3sun,x3mon;
  double mass_ratio_sun,mass_ratio_moon,re,dr,dn;
  double fac2sun,fac2mon,fac3sun,fac3mon,cosla,sinla;
  double dcorsta[3];
  void gps2tt(int,double,int *,double *);
  void sprod(double *,double *,double *,double *,double *);
  void st1idiu(double *,double *,double *,double,double,double *dcorsta);
  void st1isem(double *,double *,double *,double,double,double *dcorsta);
  void st1l1(double *,double *,double *,double,double,double *dcorsta);
  void step2diu(double *,double,double,double *);
  void step2lon(double *,double,double *);

  // Constants
  h20 = 0.6078;
  l20 = 0.0847;  
  h3  = 0.292;
  l3  = 0.015;
  mass_ratio_sun  = 332945.943062;
  mass_ratio_moon = 0.012300034;
  re = 6378136.55;

  // Convert GPS time to TT
  gps2tt(mjd,mjs,&mjdtt,&mjstt);     // TT  time (sec of day)
  dmjdtt = mjdtt+mjstt/86400.0;      // double precision MJD in TT

  // commented line was live code in dehanttideinelMJD.f
  // changed on the suggestion of Dr. Don Kim, UNB -- 09mar21
  // Julian date for 2000 January 1 00:00:00.0 UT is  JD 2451544.5
  // MJD         for 2000 January 1 00:00:00.0 UT is MJD   51544.0
  // t=(dmjdtt-51545.d0)/36525.d0                !! days to centuries, TT
  t   = (dmjdtt-51544.0)/36525.0;     // days to centuries, TT
  fhr = (dmjdtt-int(dmjdtt))*24.0;    // hours in the day, TT

  // scalar product of station vector with sun and moon vectors
  sprod(rsta,rsun,&scs,&lsta,&lsun);
  sprod(rsta,rmon,&scm,&lsta,&lmon);
  scsun = scs/lsta/lsun;
  scmon = scm/lsta/lmon;

  // computation of new h2 and l2
  cosphi = sqrt(rsta[0]*rsta[0] + rsta[1]*rsta[1])/lsta;
  h2     = h20-0.0006*(1.0-3.0/2.0*cosphi*cosphi);
  l2     = l20+0.0002*(1.0-3.0/2.0*cosphi*cosphi);

  // p2-term
  p2sun = 3.0*(h2/2.0-l2)*scsun*scsun-h2/2.0;
  p2mon = 3.0*(h2/2.0-l2)*scmon*scmon-h2/2.0;

  // p3-term
  p3sun = 5.0/2.0*(h3-3.0*l3)*pow(scsun,3)+3.0/2.0*(l3-h3)*scsun;
  p3mon = 5.0/2.0*(h3-3.0*l3)*pow(scmon,3)+3.0/2.0*(l3-h3)*scmon;

  // term in direction of sun/moon vector
  x2sun = 3.0*l2*scsun;
  x2mon = 3.0*l2*scmon;
  x3sun = 3.0*l3/2.0*(5.0*scsun*scsun-1.0);
  x3mon = 3.0*l3/2.0*(5.0*scmon*scmon-1.0);

  // factors for sun/moon
  fac2sun = mass_ratio_sun*re*pow((re/lsun),3);
  fac2mon = mass_ratio_moon*re*pow((re/lmon),3);
  fac3sun = fac2sun*(re/lsun);
  fac3mon = fac2mon*(re/lmon);

  // total displacement
  for (i=0;i<3;i++)
  {
    dtide[i]=fac2sun*(x2sun*rsun[i]/lsun + p2sun*rsta[i]/lsta ) + 
       fac2mon*(x2mon*rmon[i]/lmon + p2mon*rsta[i]/lsta ) + 
       fac3sun*(x3sun*rsun[i]/lsun + p3sun*rsta[i]/lsta ) + 
       fac3mon*(x3mon*rmon[i]/lmon + p3mon*rsta[i]/lsta );
  }

  // Apply corrections to total displacement in several steps
  for (i=0;i<3;i++) dcorsta[i] = 0.0;

  // corrections for the out-of-phase part of love numbers
  //     (part h_2^(0)i and l_2^(0)i )

  // first, for the diurnal band
  st1idiu(rsta,rsun,rmon,fac2sun,fac2mon,dcorsta);
  for (i=0;i<3;i++) dtide[i] += dcorsta[i];

  // second, for the semi-diurnal band
  st1isem(rsta,rsun,rmon,fac2sun,fac2mon,dcorsta);
  for (i=0;i<3;i++) dtide[i] += dcorsta[i];

  // corrections for the latitude dependence of love numbers (part l^(1) )
  st1l1(rsta,rsun,rmon,fac2sun,fac2mon,dcorsta);
  for (i=0;i<3;i++) dtide[i] += dcorsta[i];

  // consider corrections for step 2
  // corrections for the diurnal band:
  //  first, we need to know the date converted in julian centuries
  //  this is now handled at top of code   (also convert to TT time system)
  //  second, the diurnal band corrections,
  //   (in-phase and out-of-phase frequency dependence):
  step2diu(rsta,fhr,t,dcorsta);
  for (i=0;i<3;i++) dtide[i] += dcorsta[i];

  //  corrections for the long-period band,
  //  (in-phase and out-of-phase frequency dependence):
  step2lon(rsta,t,dcorsta);
  for (i=0;i<3;i++) dtide[i] += dcorsta[i];
  // dtide[1:3] now represents the total tidal deformation

  // All the code above removes total tidal deformation with conventional 
  // Love numbers.  This realizes a conventional tide free crust (i.e. ITRF).
  // This does NOT conform to Resolution 16 of the 18th General Assembly
  // of the IAG (1983).  This resolution has not been implemented by
  // the space geodesy community in general (c.f. IERS Conventions 2003).

  // Now we separately calculate the "permanent" portion of the total deformation
  sinphi = rsta[2]/lsta;
  cosphi = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1])/lsta;
  cosla  = rsta[0]/cosphi/lsta;
  sinla  = rsta[1]/cosphi/lsta;
  dr     = -sqrt(5.0/4.0/PI)*h2*0.31460*(3.0/2.0*sinphi*sinphi-0.5);
  dn     = -sqrt(5.0/4.0/PI)*l2*0.31460*3.0*cosphi*sinphi;
  //printf("dr,dn: %lf %lf\n",dr,dn);
  ptide[0] = dr*cosla*cosphi+dn*cosla*sinphi;
  ptide[1] = dr*sinla*cosphi+dn*sinla*sinphi;
  ptide[2] = dr*sinphi      -dn*cosphi;

}




// These are the subroutines for the step2 of the tidal corrections.
// they are called to account for the frequency dependence
// of the love numbers.
void step2lon(double *rsta,double t,double *dcorsta)
{
  int i;
  double s,pr,h,p,zns,ps,lsta,sinphi,cosphi,cosla,sinla;
  double dr_tot,dn_tot,thetaf,dr,dn,de;
  double datdi1[5],datdi2[5],datdi3[5],datdi4[5],datdi5[5];
  double datdi6[5],datdi7[5],datdi8[5],datdi9[5];

  // cf. table 7.5b of IERS conventions 2003 (TN.32, pg.82)
  // columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op)
  // IERS cols.= s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op)
  // units of mm

  datdi1[0] =  0.00;  datdi2[0] =  0.00;  datdi3[0] =  0.00;
  datdi4[0] =  1.00;  datdi5[0] =  0.00;  datdi6[0] =  0.47;
  datdi7[0] =  0.23;  datdi8[0] =  0.16;  datdi9[0] =  0.07;

  datdi1[1] =  0.00;  datdi2[1] =  2.00;  datdi3[1] =  0.00;
  datdi4[1] =  0.00;  datdi5[1] =  0.00;  datdi6[1] = -0.20;
  datdi7[1] = -0.12;  datdi8[1] = -0.11;  datdi9[1] = -0.05;

  datdi1[2] =  1.00;  datdi2[2] =  0.00;  datdi3[2] = -1.00;
  datdi4[2] =  0.00;  datdi5[2] =  0.00;  datdi6[2] = -0.11;
  datdi7[2] = -0.08;  datdi8[2] = -0.09;  datdi9[2] = -0.04;

  datdi1[3] =  2.00;  datdi2[3] =  0.00;  datdi3[3] =  0.00;
  datdi4[3] =  0.00;  datdi5[3] =  0.00;  datdi6[3] = -0.13;
  datdi7[3] = -0.11;  datdi8[3] = -0.15;  datdi9[3] = -0.07;

  datdi1[4] =  2.00;  datdi2[4] =  0.00;  datdi3[4] =  0.00;
  datdi4[4] =  1.00;  datdi5[4] =  0.00;  datdi6[4] = -0.05;
  datdi7[4] = -0.05;  datdi8[4] = -0.06;  datdi9[4] = -0.03;

  s   = 218.316645630+481267.881940*t-0.00146638890*t*t
        +0.000001851390*t*t*t;
  pr  = 1.396971278*t+0.000308889*t*t+0.000000021*t*t*t
        +0.000000007*t*t*t*t;
  s   = s+pr;
  h   = 280.466450+36000.76974890*t+0.000303222220*t*t
        +0.000000020*t*t*t-0.00000000654*t*t*t*t;
  p   = 83.353243120+4069.013635250*t-0.010321722220*t*t
        -0.00001249910*t*t*t+0.000000052630*t*t*t*t;
  zns = 234.955444990 +1934.136261970*t-0.002075611110*t*t
        -0.000002139440*t*t*t+0.000000016500*t*t*t*t;
  ps  = 282.937340980+1.719457666670*t+0.000456888890*t*t
        -0.000000017780*t*t*t-0.000000003340*t*t*t*t;

  lsta   = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1]+rsta[2]*rsta[2]);
  sinphi = rsta[2]/lsta;
  cosphi = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1])/lsta;
  cosla  = rsta[0]/cosphi/lsta;
  sinla  = rsta[1]/cosphi/lsta;

  s   = fmod(  s,360.0);
  h   = fmod(  h,360.0);
  p   = fmod(  p,360.0);
  zns = fmod(zns,360.0);
  ps  = fmod( ps,360.0);

  dr_tot = 0.0;
  dn_tot = 0.0;
  for (i=0;i<3;i++) dcorsta[i] = 0.0;
  for (i=0;i<5;i++)
  {
    thetaf     = (datdi1[i]*s+datdi2[i]*h+datdi3[i]*p+ 
                 datdi4[i]*zns+datdi5[i]*ps)*DEG2RAD;
    dr         = datdi6[i]*(3.0*sinphi*sinphi-1.0)/2.0*cos(thetaf)+ 
                 datdi8[i]*(3.0*sinphi*sinphi-1.0)/2.0*sin(thetaf);
    dn         = datdi7[i]*(cosphi*sinphi*2.0)*cos(thetaf)+ 
                 datdi9[i]*(cosphi*sinphi*2.0)*sin(thetaf);
    de         = 0.0;
    dr_tot     = dr_tot+dr;
    dn_tot     = dn_tot+dn;
    dcorsta[0] = dcorsta[0]+dr*cosla*cosphi-de*sinla
                 -dn*sinphi*cosla;
    dcorsta[1] = dcorsta[1]+dr*sinla*cosphi+de*cosla
                 -dn*sinphi*sinla;
    dcorsta[2] = dcorsta[2]+dr*sinphi+dn*cosphi;
  }
  for (i=0;i<3;i++) dcorsta[i] = dcorsta[i]/1000.0;

}



// These are the subroutines for the step2 of the tidal corrections.
// they are called to account for the frequency dependence
// of the love numbers.
void step2diu(double *rsta,double fhr,double t,double *dcorsta)
{
  int i;
  double s,tau,pr,h,p,zns,ps,lsta,sinphi,cosphi,cosla,sinla,zla;
  double thetaf,dr,dn,de;
  double datdi1[31],datdi2[31],datdi3[31],datdi4[31],datdi5[31];
  double datdi6[31],datdi7[31],datdi8[31],datdi9[31];

  datdi1[ 0] = -3.00;  datdi2[ 0] =  0.00;  datdi3[ 0] =  2.00; 
  datdi4[ 0] =  0.00;  datdi5[ 0] =  0.00;  datdi6[ 0] = -0.01;
  datdi7[ 0] = -0.01;  datdi8[ 0] =  0.00;  datdi9[ 0] =  0.00;

  datdi1[ 1] = -3.00;  datdi2[ 1] =  2.00;  datdi3[ 1] =  0.00;
  datdi4[ 1] =  0.00;  datdi5[ 1] =  0.00;  datdi6[ 1] = -0.01;
  datdi7[ 1] = -0.01;  datdi8[ 1] =  0.00;  datdi9[ 1] =  0.00;

  datdi1[ 2] = -2.00;  datdi2[ 2] =  0.00;  datdi3[ 2] =  1.00;
  datdi4[ 2] = -1.00;  datdi5[ 2] =  0.00;  datdi6[ 2] = -0.02;
  datdi7[ 2] = -0.01;  datdi8[ 2] =  0.00;  datdi9[ 2] =  0.00;

  datdi1[ 3] = -2.00;  datdi2[ 3] =  0.00;  datdi3[ 3] =  1.00;
  datdi4[ 3] =  0.00;  datdi5[ 3] =  0.00;  datdi6[ 3] = -0.08;
  datdi7[ 3] =  0.00;  datdi8[ 3] =  0.01;  datdi9[ 3] =  0.01;

  datdi1[ 4] = -2.00;  datdi2[ 4] =  2.00;  datdi3[ 4] = -1.00;
  datdi4[ 4] =  0.00;  datdi5[ 4] =  0.00;  datdi6[ 4] = -0.02;
  datdi7[ 4] = -0.01;  datdi8[ 4] =  0.00;  datdi9[ 4] =  0.00;

  datdi1[ 5] = -1.00;  datdi2[ 5] =  0.00;  datdi3[ 5] =  0.00;
  datdi4[ 5] = -1.00;  datdi5[ 5] =  0.00;  datdi6[ 5] = -0.10;
  datdi7[ 5] =  0.00;  datdi8[ 5] =  0.00;  datdi9[ 5] =  0.00;

  datdi1[ 6] = -1.00;  datdi2[ 6] =  0.00;  datdi3[ 6] =  0.00;
  datdi4[ 6] =  0.00;  datdi5[ 6] =  0.00;  datdi6[ 6] = -0.51;
  datdi7[ 6] =  0.00;  datdi8[ 6] = -0.02;  datdi9[ 6] =  0.03;

  datdi1[ 7] = -1.00;  datdi2[ 7] =  2.00;  datdi3[ 7] =  0.00;
  datdi4[ 7] =  0.00;  datdi5[ 7] =  0.00;  datdi6[ 7] =  0.01;
  datdi7[ 7] =  0.00;  datdi8[ 7] =  0.00;  datdi9[ 7] =  0.00;

  datdi1[ 8] =  0.00;  datdi2[ 8] = -2.00;  datdi3[ 8] =  1.00;
  datdi4[ 8] =  0.00;  datdi5[ 8] =  0.00;  datdi6[ 8] =  0.01;
  datdi7[ 8] =  0.00;  datdi8[ 8] =  0.00;  datdi9[ 8] =  0.00;

  datdi1[ 9] =  0.00;  datdi2[ 9] =  0.00;  datdi3[ 9] = -1.00;
  datdi4[ 9] =  0.00;  datdi5[ 9] =  0.00;  datdi6[ 9] =  0.02;
  datdi7[ 9] =  0.01;  datdi8[ 9] =  0.00;  datdi9[ 9] =  0.00;

  datdi1[10] =  0.00;  datdi2[10] =  0.00;  datdi3[10] =  1.00;
  datdi4[10] =  0.00;  datdi5[10] =  0.00;  datdi6[10] =  0.06;
  datdi7[10] =  0.00;  datdi8[10] =  0.00;  datdi9[10] =  0.00;

  datdi1[11] =  0.00;  datdi2[11] =  0.00;  datdi3[11] =  1.00;
  datdi4[11] =  1.00;  datdi5[11] =  0.00;  datdi6[11] =  0.01;
  datdi7[11] =  0.00;  datdi8[11] =  0.00;  datdi9[11] =  0.00;

  datdi1[12] =  0.00;  datdi2[12] =  2.00;  datdi3[12] = -1.00;
  datdi4[12] =  0.00;  datdi5[12] =  0.00;  datdi6[12] =  0.01;
  datdi7[12] =  0.00;  datdi8[12] =  0.00;  datdi9[12] =  0.00;

  datdi1[13] =  1.00;  datdi2[13] = -3.00;  datdi3[13] =  0.00;
  datdi4[13] =  0.00;  datdi5[13] =  1.00;  datdi6[13] = -0.06;
  datdi7[13] =  0.00;  datdi8[13] =  0.00;  datdi9[13] =  0.00;

  datdi1[14] =  1.00;  datdi2[14] = -2.00;  datdi3[14] =  0.00;
  datdi4[14] =  1.00;  datdi5[14] =  0.00;  datdi6[14] =  0.01;
  datdi7[14] =  0.00;  datdi8[14] =  0.00;  datdi9[14] =  0.00;

  datdi1[15] =  1.00;  datdi2[15] = -2.00;  datdi3[15] =  0.00;
  datdi4[15] =  0.00;  datdi5[15] =  0.00;  datdi6[15] = -1.23;
  datdi7[15] = -0.07;  datdi8[15] =  0.06;  datdi9[15] =  0.01;

  datdi1[16] =  1.00;  datdi2[16] = -1.00;  datdi3[16] =  0.00;
  datdi4[16] =  0.00;  datdi5[16] = -1.00;  datdi6[16] =  0.02;
  datdi7[16] =  0.00;  datdi8[16] =  0.00;  datdi9[16] =  0.00;

  datdi1[17] =  1.00;  datdi2[17] = -1.00;  datdi3[17] =  0.00;
  datdi4[17] =  0.00;  datdi5[17] =  1.00;  datdi6[17] =  0.04;
  datdi7[17] =  0.00;  datdi8[17] =  0.00;  datdi9[17] =  0.00;

  datdi1[18] =  1.00;  datdi2[18] =  0.00;  datdi3[18] =  0.00;
  datdi4[18] = -1.00;  datdi5[18] =  0.00;  datdi6[18] = -0.22;
  datdi7[18] =  0.01;  datdi8[18] =  0.01;  datdi9[18] =  0.00;

  datdi1[19] =  1.00;  datdi2[19] =  0.00;  datdi3[19] =  0.00;
  datdi4[19] =  0.00;  datdi5[19] =  0.00;  datdi6[19] = 12.00;
  datdi7[19] = -0.78;  datdi8[19] = -0.67;  datdi9[19] = -0.03;

  datdi1[20] =  1.00;  datdi2[20] =  0.00;  datdi3[20] =  0.00;
  datdi4[20] =  1.00;  datdi5[20] =  0.00;  datdi6[20] =  1.73;
  datdi7[20] = -0.12;  datdi8[20] = -0.10;  datdi9[20] =  0.00;

  datdi1[21] =  1.00;  datdi2[21] =  0.00;  datdi3[21] =  0.00;
  datdi4[21] =  2.00;  datdi5[21] =  0.00;  datdi6[21] = -0.04;
  datdi7[21] =  0.00;  datdi8[21] =  0.00;  datdi9[21] =  0.00;

  datdi1[22] =  1.00;  datdi2[22] =  1.00;  datdi3[22] =  0.00;
  datdi4[22] =  0.00;  datdi5[22] = -1.00;  datdi6[22] = -0.50;
  datdi7[22] = -0.01;  datdi8[22] =  0.03;  datdi9[22] =  0.00;

  datdi1[23] =  1.00;  datdi2[23] =  1.00;  datdi3[23] =  0.00;
  datdi4[23] =  0.00;  datdi5[23] =  1.00;  datdi6[23] =  0.01;
  datdi7[23] =  0.00;  datdi8[23] =  0.00;  datdi9[23] =  0.00;

  datdi1[24] =  1.00;  datdi2[24] =  1.00;  datdi3[24] =  0.00;
  datdi4[24] =  1.00;  datdi5[24] = -1.00;  datdi6[24] = -0.01;
  datdi7[24] =  0.00;  datdi8[24] =  0.00;  datdi9[24] =  0.00;       // v.dehant 2007

  datdi1[25] =  1.00;  datdi2[25] =  2.00;  datdi3[25] = -2.00;
  datdi4[25] =  0.00;  datdi5[25] =  0.00;  datdi6[25] = -0.01;
  datdi7[25] =  0.00;  datdi8[25] =  0.00;  datdi9[25] =  0.00;

  datdi1[26] =  1.00;  datdi2[26] =  2.00;  datdi3[26] =  0.00;
  datdi4[26] =  0.00;  datdi5[26] =  0.00;  datdi6[26] = -0.11;
  datdi7[26] =  0.01;  datdi8[26] =  0.01;  datdi9[26] =  0.00;

  datdi1[27] =  2.00;  datdi2[27] = -2.00;  datdi3[27] =  1.00;
  datdi4[27] =  0.00;  datdi5[27] =  0.00;  datdi6[27] = -0.01;
  datdi7[27] =  0.00;  datdi8[27] =  0.00;  datdi9[27] =  0.00;

  datdi1[28] =  2.00;  datdi2[28] =  0.00;  datdi3[28] = -1.00;
  datdi4[28] =  0.00;  datdi5[28] =  0.00;  datdi6[28] = -0.02;
  datdi7[28] =  0.02;  datdi8[28] =  0.00;  datdi9[28] =  0.01;

  datdi1[29] =  3.00;  datdi2[29] =  0.00;  datdi3[29] =  0.00;
  datdi4[29] =  0.00;  datdi5[29] =  0.00;  datdi6[29] =  0.00;
  datdi7[29] =  0.01;  datdi8[29] =  0.00;  datdi9[29] =  0.01;

  datdi1[30] =  3.00;  datdi2[30] =  0.00;  datdi3[30] =  0.00;
  datdi4[30] =  1.00;  datdi5[30] =  0.00;  datdi6[30] =  0.00;
  datdi7[30] =  0.01;  datdi8[30] =  0.00;  datdi9[30] =  0.00;

  s   = 218.31664563+481267.88194*t-0.0014663889*t*t
        +0.00000185139*t*t*t;
  tau = fhr*15.0+280.4606184+36000.7700536*t+0.00038793*t*t
        -0.0000000258*t*t*t-s;
  pr  = 1.396971278*t+0.000308889*t*t+0.000000021*t*t*t
        +0.000000007*t*t*t*t;
  s   = s+pr;
  h   = 280.46645+36000.7697489*t+0.00030322222*t*t
        +0.000000020*t*t*t-0.00000000654*t*t*t*t;
  p   = 83.35324312+4069.01363525*t-0.01032172222*t*t
        -0.0000124991*t*t*t+0.00000005263*t*t*t*t;
  zns = 234.95544499 +1934.13626197*t-0.00207561111*t*t
        -0.00000213944*t*t*t+0.00000001650*t*t*t*t;
  ps  = 282.93734098+1.71945766667*t+0.00045688889*t*t
        -0.00000001778*t*t*t-0.00000000334*t*t*t*t;

  // reduce angles to between 0 and 360
  s   = fmod(  s,360.0);
  tau = fmod(tau,360.0);
  h   = fmod(  h,360.0);
  p   = fmod(  p,360.0);
  zns = fmod(zns,360.0);
  ps  = fmod( ps,360.0);

  lsta   = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1]+rsta[2]*rsta[2]);
  sinphi = rsta[2]/lsta;
  cosphi = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1])/lsta;
  cosla  = rsta[0]/cosphi/lsta;
  sinla  = rsta[1]/cosphi/lsta;
  zla    = atan2(rsta[1],rsta[0]);

  for (i=0;i<3;i++) dcorsta[i] = 0.0;
  for (i=0; i<31; i++)
  {
    thetaf = (tau+datdi1[i]*s+datdi2[i]*h+datdi3[i]*p+ 
             datdi4[i]*zns+datdi5[i]*ps)*DEG2RAD;
    dr     = datdi6[i]*2.0*sinphi*cosphi*sin(thetaf+zla)+ 
             datdi7[i]*2.0*sinphi*cosphi*cos(thetaf+zla);
    dn     = datdi8[i]*(cosphi*cosphi-sinphi*sinphi)*sin(thetaf+zla)+ 
             datdi9[i]*(cosphi*cosphi-sinphi*sinphi)*cos(thetaf+zla);

    // following correction by V.Dehant to match eq.16b, p.81, 2003 Conventions
    de          = datdi8[i]*sinphi*cos(thetaf+zla)-
                  datdi9[i]*sinphi*sin(thetaf+zla);
    dcorsta[0]  = dcorsta[0]+dr*cosla*cosphi-de*sinla
                  -dn*sinphi*cosla;
    dcorsta[1]  = dcorsta[1]+dr*sinla*cosphi+de*cosla
                  -dn*sinphi*sinla;
    dcorsta[2]  = dcorsta[2]+dr*sinphi+dn*cosphi;
  }

  for (i=0;i<3;i++) dcorsta[i] = dcorsta[i]/1000.0;

}



// This subroutine gives the corrections induced by the latitude
// dependence given by part l^(1) in Matthews et al (1991)
void st1l1(double *rsta,double *rsun,double *rmon,double fac2sun,double fac2mon,double *dcorsta)
{
  double l1d,l1sd,lsta,sinphi,cosphi,sinla,cosla,lmon,lsun,dnsun,dnmon,l1;
  double desun,demon,de,dn,costwola,sintwola;
  double enorm8(double *);

  //  Constants
  l1d  = 0.0012;
  l1sd = 0.0024;
  lsta   = enorm8(rsta);
  sinphi = rsta[2]/lsta;
  cosphi = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1])/lsta;
  sinla  = rsta[1]/cosphi/lsta;
  cosla  = rsta[0]/cosphi/lsta;
  lmon=enorm8(rmon);
  lsun=enorm8(rsun);

  // for the diurnal band
  l1    = l1d;
  dnsun = -l1*sinphi*sinphi*fac2sun*rsun[2]*(rsun[0]*cosla+rsun[1]*sinla)
          /(lsun*lsun);
  dnmon = -l1*sinphi*sinphi*fac2mon*rmon[2]*(rmon[0]*cosla+rmon[1]*sinla)
          /(lmon*lmon);
  desun = l1*sinphi*(cosphi*cosphi-sinphi*sinphi)*fac2sun*rsun[2]*
          (rsun[0]*sinla-rsun[1]*cosla)/(lsun*lsun);
  demon = l1*sinphi*(cosphi*cosphi-sinphi*sinphi)*fac2mon*rmon[2]*
          (rmon[0]*sinla-rmon[1]*cosla)/(lmon*lmon);
  de    = 3.0*(desun+demon);
  dn    = 3.0*(dnsun+dnmon);
  dcorsta[0] = -de*sinla-dn*sinphi*cosla;
  dcorsta[1] =  de*cosla-dn*sinphi*sinla;
  dcorsta[2] =           dn*cosphi;

  // for the semi-diurnal band
  l1 = l1sd;
  costwola = cosla*cosla-sinla*sinla;
  sintwola = 2.0*cosla*sinla;
  dnsun    = -l1/2.0*sinphi*cosphi*fac2sun*((rsun[0]*rsun[0]-rsun[1]*rsun[1])*
             costwola+2.0*rsun[0]*rsun[1]*sintwola)/(lsun*lsun);
  dnmon    = -l1/2.0*sinphi*cosphi*fac2mon*((rmon[0]*rmon[0]-rmon[1]*rmon[1])*
             costwola+2.0*rmon[0]*rmon[1]*sintwola)/(lmon*lmon);
  desun    = -l1/2.0*sinphi*sinphi*cosphi*fac2sun*((rsun[0]*rsun[0]-rsun[1]*rsun[1])*
             sintwola-2.0*rsun[0]*rsun[1]*costwola)/(lsun*lsun);
  demon    = -l1/2.0*sinphi*sinphi*cosphi*fac2mon*((rmon[0]*rmon[0]-rmon[1]*rmon[1])*
             sintwola-2.0*rmon[0]*rmon[1]*costwola)/(lmon*lmon);
  de       = 3.0*(desun+demon);
  dn       = 3.0*(dnsun+dnmon);

  dcorsta[0] = dcorsta[0]-de*sinla-dn*sinphi*cosla;
  dcorsta[1] = dcorsta[1]+de*cosla-dn*sinphi*sinla;
  dcorsta[2] = dcorsta[2]         +dn*cosphi;

}



// This subroutine gives the out-of-phase corrections induced by
// mantle inelasticity in the semi-diurnal band
void st1isem(double *rsta,double *rsun,double *rmon,double fac2sun,double fac2mon,double *dcorsta)
{
  double dhi,dli,sinphi,cosphi,lsta,costwola,sintwola,cosla,sinla;
  double lmon,lsun,drmon,drsun,dnsun,dnmon,demon,desun,dr,dn,de;
  double enorm8(double *);

  // Constants
  dhi = -0.0022;
  dli = -0.0007;

  lsta     = enorm8(rsta);
  sinphi   = rsta[2]/lsta;
  cosphi   = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1])/lsta;
  sinla    = rsta[1]/cosphi/lsta;
  cosla    = rsta[0]/cosphi/lsta;
  costwola = cosla*cosla-sinla*sinla;
  sintwola = 2.0*cosla*sinla;
  lmon     = enorm8(rmon);
  lsun     = enorm8(rsun);
  drsun    = -3.0/4.0*dhi*cosphi*cosphi*fac2sun*((rsun[0]*rsun[0]-rsun[1]*rsun[1])* 
             sintwola-2.0*rsun[0]*rsun[1]*costwola)/(lsun*lsun);
  drmon    = -3.0/4.0*dhi*cosphi*cosphi*fac2mon*((rmon[0]*rmon[0]-rmon[1]*rmon[1])*
             sintwola-2.0*rmon[0]*rmon[1]*costwola)/(lmon*lmon);
  dnsun    = 1.50*dli*sinphi*cosphi*fac2sun*((rsun[0]*rsun[0]-rsun[1]*rsun[1])*
             sintwola-2.0*rsun[0]*rsun[1]*costwola)/(lsun*lsun);
  dnmon    = 1.50*dli*sinphi*cosphi*fac2mon*((rmon[0]*rmon[0]-rmon[1]*rmon[1])*
             sintwola-2.0*rmon[0]*rmon[1]*costwola)/(lmon*lmon);
  desun    = -3.0/2.0*dli*cosphi*fac2sun*((rsun[0]*rsun[0]-rsun[1]*rsun[1])*
             costwola+2.0*rsun[0]*rsun[1]*sintwola)/(lsun*lsun);
  demon    = -3.0/2.0*dli*cosphi*fac2mon*((rmon[0]*rmon[0]-rmon[1]*rmon[1])*
             costwola+2.0*rmon[0]*rmon[1]*sintwola)/(lmon*lmon);
  dr       = drsun+drmon;
  dn       = dnsun+dnmon;
  de       = desun+demon;

  dcorsta[0] = dr*cosla*cosphi-de*sinla-dn*sinphi*cosla;
  dcorsta[1] = dr*sinla*cosphi+de*cosla-dn*sinphi*sinla;
  dcorsta[2] = dr*sinphi+dn*cosphi;

}



// This subroutine gives the out-of-phase corrections induced by
// mantle inelasticity in the diurnal band
void st1idiu(double *rsta,double *rsun,double *rmon,double fac2sun,double fac2mon,double *dcorsta)
{
  double dhi,dli,lsta,sinphi,cosphi,cos2phi,sinla,cosla,lmon,lsun;
  double drsun,drmon,dnsun,dnmon,desun,demon,dr,dn,de;
  double enorm8(double *);

  // Constants
  dhi = -0.0025;
  dli = -0.0007;

  lsta=enorm8(rsta);
  sinphi  = rsta[2]/lsta;
  cosphi  = sqrt(rsta[0]*rsta[0]+rsta[1]*rsta[1])/lsta;
  cos2phi = cosphi*cosphi-sinphi*sinphi;
  sinla = rsta[1]/cosphi/lsta;
  cosla = rsta[0]/cosphi/lsta;
  lmon = enorm8(rmon);
  lsun = enorm8(rsun);
  drsun = -3.0*dhi*sinphi*cosphi*fac2sun*rsun[2]*(rsun[0]* 
          sinla-rsun[1]*cosla)/(lsun*lsun);
  drmon = -3.0*dhi*sinphi*cosphi*fac2mon*rmon[2]*(rmon[0]*
          sinla-rmon[1]*cosla)/(lmon*lmon);
  dnsun = -3.0*dli*cos2phi*fac2sun*rsun[2]*(rsun[0]*sinla-
          rsun[1]*cosla)/(lsun*lsun);
  dnmon = -3.0*dli*cos2phi*fac2mon*rmon[2]*(rmon[0]*sinla-
          rmon[1]*cosla)/(lmon*lmon);
  desun = -3.0*dli*sinphi*fac2sun*rsun[2]*
          (rsun[0]*cosla+rsun[1]*sinla)/(lsun*lsun);
  demon = -3.0*dli*sinphi*fac2mon*rmon[2]* 
          (rmon[0]*cosla+rmon[1]*sinla)/(lmon*lmon);
  dr = drsun+drmon;
  dn = dnsun+dnmon;
  de = desun+demon;
  dcorsta[0] = dr*cosla*cosphi-de*sinla-dn*sinphi*cosla;
  dcorsta[1] = dr*sinla*cosphi+de*cosla-dn*sinphi*sinla;
  dcorsta[2] = dr*sinphi+dn*cosphi;

}


// Compute euclidian norm of a vector (of length 3)
double enorm8(double *a)
{
  return( sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) );
}



// Computation of the scalar-product of two vectors and their norms
void sprod(double *x,double *y,double *scal,double *r1,double *r2)
{
  *r1   = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  *r2   = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
  *scal = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}


double getleapsecs(int mjd,double mjs,int tid)  //GPS-UTC ***not**** TAI-UTC
{
  char line[85];
  int i,k,mjdutc;
  int mjdls[100];
  double leapsecs,mjsutc;
  double mjsls[100],ls[100];
  double tdiff(int,double,int,double);
  FILE *fptr;

  // Open the leapsecs file
  if ((fptr=fopen("/usr/local/share/leapsecs/leapsecs_mjd.txt","r")) == NULL)
  {
    printf("Leap seconds file /usr/local/share/leapsecs/leapsecs_mjd.txt not found - exiting\n");
    exit(-1);
  }

  // Read the leapsecs file
  k = 0;
  while(fgets(line,85,fptr)!=NULL)
  {
    sscanf(line,"%d %lf %lf",&mjdls[k],&mjsls[k],&ls[k]);
    ++k;
  }

  // Find current leap seconds if input time is UTC
  if (tid!=2)
  {
    //printf("input time is UTC\n");
    i = 0;
    while((tdiff(mjd,mjs,mjdls[i],mjsls[i])>0.0)&&(i<k))
    {
      ++i;
    }
    if (i==0) leapsecs = 0.0;
    else leapsecs = ls[i-1];
  }

  // Find current leap seconds if input time is GPS
  else
  {
    //printf("input time is GPS\n");
    leapsecs = 0.0;
    for (i=0;i<k;i++)
    {
      mjdutc = mjd;
      mjsutc = mjs-leapsecs;
      while (mjsutc<0.0)
      {
        mjsutc += 86400.0;
        mjdutc -= 1;
      }
      if (tdiff(mjdutc,mjsutc,mjdls[i],mjsls[i])>0.0) leapsecs = ls[i];
    }
  }

  // Close the leapsecs file
  fclose(fptr);

  // Return the result
  return(leapsecs);

}


double tdiff(int mjd1,double mjs1,int mjd2,double mjs2)
{
  return((mjd1-mjd2)*86400.0+(mjs1-mjs2));
}



// Get low-precision, geocentric coordinates for sun (ECEF)
//
//### Files Accessed
// * None
//
//### References
// * "satellite orbits: models, methods, applications" montenbruck & gill(2000)
//   section 3.3.2, pg. 70-71
//
// * "astronomy on the personal computer, 4th ed." montenbruck & pfleger (2005)
//   section 3.2, pg. 39  routine MiniSun
//
//### Comments
//  * input, mjd/mjs, is Modified Julian Day (and seconds of day) in GPS time
//   output, rs, is geocentric solar position vector [m] in ECEF
//
void sunxyz(int mjd,double mjs,double leapsecs,double *rs)
{
  int mjdtt;
  double obe,sobe,cobe,opod,tsectt,fmjdtt,rs1,rs2,rs3,ghar;
  double mjstt,tjdtt,t,emdeg,em,em2,r,slond,slon,sslon,cslon;
  void gps2tt(int,double,int *,double *);
  double getghar(int,double,double);

  // mean elements for year 2000, sun ecliptic orbit wrt. Earth
  obe = 23.43929111*DEG2RAD;      // obliquity of the J2000 ecliptic
  sobe = sin(obe);
  cobe = cos(obe);
  opod = 282.9400;                // RAAN + arg.peri.  (deg)

  // use TT for solar ephemerides
  gps2tt(mjd,mjs,&mjdtt,&mjstt);  // TT  time (sec of day)

  // julian centuries since 1.5 january 2000 (J2000)
  // (note: also low precision use of mjd --> tjd)
  tjdtt = mjdtt+(mjstt/86400.0)+2400000.5; // Julian Date, TT
  t     = (tjdtt - 2451545.0)/36525.0;     // julian centuries, TT
  emdeg = 357.52560 + 35999.0490*t;        // degrees
  em    = emdeg*DEG2RAD;                   // radians
  em2   = em+em;                           // radians

  // series expansions in mean anomaly, em   (eq. 3.43, p.71)
  r = (149.619-2.499*cos(em)-0.021*cos(em2))*1.e9;      // m.
  slond = opod + emdeg + (6892.0*sin(em)+72.0*sin(em2))/3600.0;

  // precession of equinox wrt. J2000   (p.71)
  slond = slond + 1.39720*t;                     // degrees

  // position vector of sun (mean equinox & ecliptic of J2000) (EME2000, ICRF)
  //                        (plus long. advance due to precession -- eq. above)
  slon =slond*DEG2RAD;                              // radians
  sslon=sin(slon);
  cslon=cos(slon);
  rs1 = r*cslon;                                // eq. 3.46, p.71 meters
  rs2 = r*sslon*cobe;                           // eq. 3.46, p.71 meters
  rs3 = r*sslon*sobe;                           // eq. 3.46, p.71 meters

  // convert position vector of sun to ECEF  (ignore polar motion/LOD)
  ghar = getghar(mjd,mjs,leapsecs);           // sec 2.3.1,p.33
  rs[0] = rs1*cos(ghar)+rs2*sin(ghar);
  rs[1] = rs2*cos(ghar)-rs1*sin(ghar);
  rs[2] = rs3;

}



// Convert mjd/mjs in GPS time to Greenwich hour angle (in radians)
//
//### Files Accessed
// * None
//
//### References
// * "satellite orbits: models, methods, applications" montenbruck & gill(2000)
//   section 2.3.1, pg. 33
double getghar(int mjd,double mjs,double leapsecs)
{
  int mjdutc;
  double mjsutc,d,ghad,revs;
 
  // Convert to UTC time
  mjdutc = mjd;
  mjsutc = mjs-leapsecs;
  while(mjsutc<0.0)
  {
    mjsutc += 86400.0;
    mjdutc--;
  }
  d =(mjdutc-51544) + (mjsutc/86400.0-0.5);       // days since J2000

  // greenwich hour angle for J2000  (12:00:00 on 1 Jan 2000)
  ghad = 280.460618375040 + 360.98564736628620*d;  // corrn.   (+digits)

  // normalize to 0-360 and convert to radians
  revs = int(ghad/360.0);
  ghad = ghad-double(revs)*360.0;
  while (ghad>=360.0) ghad-=360.0;
  while (ghad<0.0) ghad+=360.0;

  // Return the result
  return(ghad*DEG2RAD);

}



void gps2tt(int mjd,double mjs,int *mjdtt,double *mjstt)
{
  double dtemp;

  *mjdtt = mjd;
  dtemp = mjs+51.184;
  while (dtemp>86400.0)
  {
    *mjdtt++;
    dtemp -= 86400.0;
  }
  *mjstt = dtemp;

}



// get low-precision, geocentric coordinates for moon (ECEF)
//
//
//!### References
// * "satellite orbits: models, methods, applications" montenbruck & gill(2000)
//     section 3.3.2, pg. 72-73
//
// * "astronomy on the personal computer, 4th ed." montenbruck & pfleger (2005)
//     section 3.2, pg. 38-39  routine MiniMoon
//
//### Comments
// * input:  mjd/mjs, is Modified Julian Day (and seconds of day) in GPS time
//
// * output: rm, is geocentric lunar position vector [m] in ECEF
//
void moonxyz(int mjd,double mjs,double leapsecs,double *rm)
{
  int mjdtt;
  double mjstt,tjdtt,t,el0,el,elp,f,d,selond,selatd,q,rse,ghar;
  double oblir,sselat,sselon,cselat,cselon,t1,t2,t3,rm1,rm2,rm3;
  double getghar(int,double,double);

  // use TT for solar ephemerides
  gps2tt(mjd,mjs,&mjdtt,&mjstt);  // TT  time (sec of day)

  // julian centuries since 1.5 january 2000 (J2000)
  // (note: also low precision use of mjd --> tjd)
  tjdtt = mjdtt+(mjstt/86400.0)+2400000.5; // Julian Date, TT
  t     = (tjdtt - 2451545.0)/36525.0;     // julian centuries, TT

  // el0 -- mean longitude of Moon (deg)
  // el  -- mean anomaly of Moon (deg)
  // elp -- mean anomaly of Sun  (deg)
  // f   -- mean angular distance of Moon from ascending node (deg)
  // d   -- difference between mean longitudes of Sun and Moon (deg)
  // equations 3.47, p.72
  el0 = 218.31617 + 481267.88088*t - 1.3972*t;
  el  = 134.96292 + 477198.86753*t;
  elp = 357.52543 +  35999.04944*t;
  f   =  93.27283 + 483202.01873*t;
  d   = 297.85027 + 445267.11135*t;

  // longitude w.r.t. equinox and ecliptic of year 2000  eq 3.48, p.72
  selond = el0                                      
         + 22640.0/3600.0*sin((el        )*DEG2RAD)
         +  769.0/3600.0*sin((el+el     )*DEG2RAD)
         - 4586.0/3600.0*sin((el-d-d    )*DEG2RAD)
         + 2370.0/3600.0*sin((d+d       )*DEG2RAD)
         -  668.0/3600.0*sin((elp       )*DEG2RAD)
         -  412.0/3600.0*sin((f+f       )*DEG2RAD)
         -  212.0/3600.0*sin((el+el-d-d )*DEG2RAD)
         -  206.0/3600.0*sin((el+elp-d-d)*DEG2RAD)
         +  192.0/3600.0*sin((el+d+d    )*DEG2RAD)
         -  165.0/3600.0*sin((elp-d-d   )*DEG2RAD)
         +  148.0/3600.0*sin((el-elp    )*DEG2RAD)
         -  125.0/3600.0*sin((d         )*DEG2RAD)
         -  110.0/3600.0*sin((el+elp    )*DEG2RAD)
         -   55.0/3600.0*sin((f+f-d-d   )*DEG2RAD);

  // latitude w.r.t. equinox and ecliptic of year 2000  eq 3.49, p.72
  q = 412.0/3600.0*sin((f+f)*DEG2RAD)+541.0/3600.0*sin((elp)*DEG2RAD);
  selatd =                                        
         +18520.0/3600.0*sin((f+selond-el0+q)*DEG2RAD)
         -  526.0/3600.0*sin((f-d-d     )*DEG2RAD)
         +   44.0/3600.0*sin((el+f-d-d  )*DEG2RAD)
         -   31.0/3600.0*sin((-el+f-d-d )*DEG2RAD)
         -   25.0/3600.0*sin((-el-el+f  )*DEG2RAD)
         -   23.0/3600.0*sin((elp+f-d-d )*DEG2RAD)
         +   21.0/3600.0*sin((-el+f     )*DEG2RAD)
         +   11.0/3600.0*sin((-elp+f-d-d)*DEG2RAD);

  // distance from Earth center to Moon (m)   eq 3.50, p.72
  rse    = 385000.0*1000.0                 
         -  20905.0*1000.0*cos((el        )*DEG2RAD)
         -   3699.0*1000.0*cos((d+d-el    )*DEG2RAD)
         -   2956.0*1000.0*cos((d+d       )*DEG2RAD)
         -    570.0*1000.0*cos((el+el     )*DEG2RAD)
         +    246.0*1000.0*cos((el+el-d-d )*DEG2RAD)
         -    205.0*1000.0*cos((elp-d-d   )*DEG2RAD)
         -    171.0*1000.0*cos((el+d+d    )*DEG2RAD)
         -    152.0*1000.0*cos((el+elp-d-d)*DEG2RAD);

  // precession of equinox wrt. J2000   (p.71)
  selond = selond + 1.39720*t;                      // degrees

  // position vector of moon (mean equinox & ecliptic of J2000) (EME2000, ICRF)
  //                         (plus long. advance due to precession -- eq. above)
  oblir  = 23.43929111*DEG2RAD;       // obliquity of the J2000 ecliptic
  sselat = sin(selatd*DEG2RAD);
  cselat = cos(selatd*DEG2RAD);
  sselon = sin(selond*DEG2RAD);
  cselon = cos(selond*DEG2RAD);
  t1 = rse*cselon*cselat;        // meters              eq. 3.51, p.72
  t2 = rse*sselon*cselat;        // meters              eq. 3.51, p.72
  t3 = rse*       sselat;        // meters              eq. 3.51, p.72
  rm1 = t1;
  rm2 = t2*cos(-oblir)+t3*sin(-oblir);
  rm3 = t3*cos(-oblir)-t2*sin(-oblir);

  // convert position vector of moon to ECEF  (ignore polar motion/LOD)
  //call getghar(mjd, fmjd, leap_secs, ghar)            !! sec 2.3.1,p.33
  //call rot3(ghar, rm1, rm2, rm3, rm(1), rm(2), rm(3)) !! eq. 2.89, p.37
  ghar = getghar(mjd,mjs,leapsecs);           // sec 2.3.1,p.33
  rm[0] = rm1*cos(ghar)+rm2*sin(ghar);
  rm[1] = rm2*cos(ghar)-rm1*sin(ghar);
  rm[2] = rm3;

}

