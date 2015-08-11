


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:    This file contains subroutines required for creating a 
  	disk atmosphere

  Description:	

  Arguments:  


  Returns:

  Notes:

 

  History:
1508	ksl	Begain coding as part of effort to incorporate disk reflection
		into python

 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************

       Space Telescope Science Institute

Synopis:
	From the position calculate the velocity of material
	in a disk atmosphere

Arguments

        double x[]	the position where velocity is 
			to be calculated.       

Returns:
	double v[]     the calculated velocity
	The amplitude of the velocity is returned

Description:

Notes:
	This is just a convenience function.   The routine
	just calls corona_velocity

	v is set to the Keplerian velocity at this radius

History:

1508	ksl	Began coding

************************************************************/


double
disk_atmos_velocity(x, v)
	     double x[], v[];
{
	double speed;
	speed=corona_velocity (x, v);
	return speed;
}


/***********************************************************
	Space Telescope Science Institute

Synopsis:

double disk_atmos_rho(x) calculates the density of a disk atmophere at a position x

Arguments:              
	double x[] 	the position where one desires the density
Returns:

The density at x is returned in gram/cm**3.  

Description:

At present,  the function is just a function of the cylindrical distance of the atmosphere
since we are modelling a uniform slab.  

Notes:

History:
1508	ksl	Coded as part of effort to incoporate a disk atmosphere

**************************************************************/

double disk_atmos_rho(x)
	double x[];
{
	double r;
	double tref, t;
	double v[3],vrot,vsound;
	double h;
	double alpha,nu;
	double sigma,q,rho;


	r=sqrt(x[0]*x[0]+x[1]*x[1]);  // The distance from origin in the xy plane

	tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
	t = teff (tref, r / geo.rstar);  //The disk temperature at r

	vrot=disk_atmos_velocity(x,v);  //The rotatational velocity at r
	vsound=kn_vzero (r);  //Get the sound speed from our existing code for KWD

	h=vsound/vrot*r;  //This is the thickness of the disk

	alpha=1;
	nu=alpha*vsound*h;  //This is the viscosity

	q=1./(3*PI)*geo.disk_mdot*(1.-sqrt(geo.rstar/r));

	sigma=q/nu; // sigma is the surface density of the disk

	rho=sigma/h;

	return rho;



}
