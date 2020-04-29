/*------------------------------------------------------------------------*
 NAME:     setide.cpp

 PURPOSE:  Calculates 3-dimensional solid earth tide displacement vector
           given MJD (UTC scale) and station coordinates.

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
  int mjd;
  double mjs;
  double sta[3],dtide[3];
  void setide(int,double,int,double *,double *);

  //  Check input
  if ( argc != 1 )
  {
    printf("Usage:  setide\n");
    exit(-1);
  }

  sta[0] =   538117.2;
  sta[1] = -1389036.5;
  sta[2] =  6180994.1;
  mjd = 57800;
  mjs = 0;
  setide(mjd,mjs,2,sta,dtide);

}



/*------------------------------------------------------------------------*
 NAME:     setide(int mjd, double mjs, int tid, double *sta, double *dtide)

 PURPOSE:  Calculate 3-dimensional solid earth tide deformation dtide[3] 
           at time mjd/mjs and at station coordinates sta[3].  User 
           specifies whether input time is GPS or UTC with
           tid = 1: specified time is UTC
           tid = 2: specified time is GPS

 AUTHOR:   John Gary Sonntag

 DATE:     24 April 2020
 *------------------------------------------------------------------------*/
void setide(int mjdin,double mjsin,int tid,double *sta,double *dtide)
{
  int mjd;
  double leapsecs,mjs;
  double rs[3],rm[3];
  void sunxyz(int,double,double,double *),moonxyz(int,double,double,double *);
  double getleapsecs(int,double,int);

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
  printf("%5d %7.1lf %5d %7.1lf %2.0lf\n",mjdin,mjsin,mjd,mjs,leapsecs);

  // Calculate sun position vector (ECEF)
  sunxyz(mjd,mjs,leapsecs,rs);

  // Calculate moon position vector (ECEF)
  moonxyz(mjd,mjs,leapsecs,rm);
  printf("sun:  %lf %lf %lf\n",rs[0]/1000.0,rs[1]/1000.0,rs[2]/1000.0);
  printf("moon: %lf %lf %lf\n",rm[0]/1000.0,rm[1]/1000.0,rm[2]/1000.0);
 
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

