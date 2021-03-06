
/***********************************************************/
/** @file  bilinear.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  As the name indicates the routines solve for 
 * the fractional position of a 2d vector in a 2d cell with
 * 4 corners.
 *
 * There is also a quadratic equation solver
 ***********************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "log.h"


/**********************************************************/
/** 
 * @brief      calculates the fractional position in the coordinate
 *  	grid of a 2d vector x.
 *
 * @param [in] double  x[]   The 2d position within the cell
 * @param [in] double  x00[] The 2d position of the lower-left edge of cell  
 * @param [in] double  x01[] The 2d position of the upper-left edge of the cell
 * @param [in] double  x10[] The 2d position of the lower-right edge of the cell
 * @param [in] double  x11[] The 2d position of the upper-right edge of the cell
 * @param [out] double *  f  The fractional position in the "x" direction
 * @param [out] double *  g  The fractional postion in the "z" direction
 * @return     
 *     The routine itself returns 0 if the postion is inside the cell, and 
 *     -1 otherwise.
 *
 * @details
 * The input 2d input vectors are assumed to have the relationship
 *   shown in real space, that is x00 and x01 are more or less along
 *   the z axis, while x00 and x10 are more or less along the +x acis
 * 
 * * x01 ....... x11
 * * .            .
 * * .            .
 * * x00 ....... x10
 *
 * ### Notes ###
 *
 * Basically
 *         
 * X= (1-g) [ (1-f) X00 + f X10 ] * g [ (1-f)*X01+ f* X11 ]
 * 
 * where the capitalized values of X,etc denote vectors, and
 * f and g are the fractional weightings.  
 * 
 * X = (1-g) A + g B 
 * 
 * It is easy to solve for g in terms of f
 * 
 * g = (X-A)/(B-A)
 * 
 * where A = [ (1-f) X00 + f X10 ] and B =  [ (1-f)*X01+ f* X11 ]
 * 
 * But this gives you two expressions for g depending on whether you
 * consider the 0 or 1 component of the vector, and therefore you
 * can solve for f.  More specifically, 
 * 
 * g= [(X-XOO) + f (XOO-X10)] / [(X01-X10)+f (X11+X00-X01-X10)]
 * 
 * or 
 * 
 * g = (Q + f R) / (S+f T) 
 * 
 * where Q, R, S, and T are defined from the equation above.
 * 
 * So expanding in terms of the component veciots what we have is
 * 
 * (Q[0]+f R[0]) (S[2]+f T[2]) = (Q[2]+f R[2]) (S[0]+f T[0])
 * 
 * which can now be put in the form of a xquadratic
 * 
 * a f**2 + b f +c = 0
 * 
 * where 
 * 
 * * a = r[0]*t[2] - r[2]*t[0]
 * * b=  r[0]*s[2] + q[0]*t[2] - ( r[2]*s[0] + q[2]*t[0] )
 * * c = q[0]*s[2] - q[2]*s[0]
 * 
 * 
 * A problem with this approach is that does not work in the case 
 * where the coordinates are the x axis positions of X01 and X01
 * are the same and the X10 and X11 are the same.  This is because
 * the denominator becomes 0 in the equation for g.  This 
 * particular case is important because it is the case we are
 * trying to solve for.
 * 
 * In this case we need to solve for f by itself. Similarly there 
 * is a degnereate case when
 *               (x00[2] == x10[2] && x01[2] == x11[2])
 *
 **********************************************************/

int
bilin (x, x00, x01, x10, x11, f, g)
     double x[], x00[], x01[], x10[], x11[];
     double *f, *g;

{
  double z;
  double root[2];
  double a, b, c, d;
  double q[3], r[3], s[3], t[3];
  double zz[3];
  int i, xquadratic ();
  void Exit (int error_code);


  z = 0;                        /* Initialize z to prevent warning on compilation */

  if (x00[0] == x01[0] && x10[0] == x11[0])
  {
    *f = z = (x[0] - x00[0]) / (x10[0] - x00[0]);
    a = (1. - z) * x00[2] + z * x10[2];
    b = (1. - z) * x01[2] + z * x11[2];
    *g = (x[2] - a) / (b - a);

  }
  else if (x00[2] == x10[2] && x01[2] == x11[2])
  {
    *g = z = (x[2] - x00[2]) / (x01[2] - x00[2]);
    a = (1. - z) * x00[0] + z * x01[0];
    b = (1. - z) * x10[0] + z * x11[0];
    *f = (x[0] - a) / (b - a);
  }
  else
  {
    q[0] = x[0] - x00[0];
    r[0] = x00[0] - x10[0];
    s[0] = x01[0] - x00[0];
    t[0] = x11[0] - x01[0] - x10[0] + x00[0];

    q[2] = x[2] - x00[2];
    r[2] = x00[2] - x10[2];
    s[2] = x01[2] - x00[2];
    t[2] = x11[2] - x01[2] - x10[2] + x00[2];


    a = r[0] * t[2] - r[2] * t[0];
    b = r[0] * s[2] + q[0] * t[2] - (r[2] * s[0] + q[2] * t[0]);
    c = q[0] * s[2] - q[2] * s[0];


    if (c == 0)
    {
      *f = z = 0;
    }
    else
    {
      i = xquadratic (a, b, c, root);   /* root[i] is the smallest positive root unless i is
                                           negative in which case either both roots were negative or both roots 
                                           were imaginary */
      if (i != (-1))
      {                         // Find the root that closest to 0.5
        zz[0] = fabs (root[0] - 0.5);
        zz[1] = fabs (root[1] - 0.5);
        if (zz[0] < zz[1])
          z = root[0];
        else
          z = root[1];
        *f = z;

      }
      else
      {
        Error ("Bilinear -- imaginary roots from xquadratic\n");
      }

    }


    if ((d = (s[2] + z * t[2])) != 0)
    {
      *g = (q[2] + r[2] * z) / (d);
    }
    else if ((d = (s[0] + z * t[0])) != 0)
    {
      *g = (q[0] + r[0] * z) / (d);
    }
    else
    {
      Error ("bilin: Denominator zero\n");
      Exit (0);
    }
  }

  if (*f < 0. || *f > 1. || *g < 0 || *g > 1.)
  {
    return (-1);
  }

  return (0);

}




/**********************************************************/
/** 
 * @brief      solves a quadratic of the form 0=ax**2+bx+c
 *
 * @param [in, out] double  a   
 * @param [in, out] double  b   
 * @param [in, out] double  c   
 * @param [in, out] double  r[] roots of the quadratic equation   
 * @return    
 * * -1 -> both roots are imaginary
 * * -2 -> both roots are negative or 0
 * *  0 -> the first root is the smallest positive  root 
 * *  1 -> the second root is the smallest positive root
 *
 * @details
 * This solves a simple xquadratic (or if a is zero linear equation).  The return is set up
 * to make it easy to identify the smallest positive root if one exists.  The routine returns
 * a negative number if both roots are negative or imaginary. 
 *
 * ### Notes ###
 *
 * @bug xquadratic looks line for line identical to another routine quadratic which can be
 * found in phot_util.c.  Both versions of the code seem to be called.  One of them
 * should be removed.
 *
 **********************************************************/

int
xquadratic (a, b, c, r)
     double a, b, c, r[];
{
  double q, z;
  double qq;

  if (a == 0.0)
  {                             /* Then it's not really a xquadratic but we can solve it
                                   anyway */
    if (b == 0.0)
    {
      r[0] = r[1] = -99.;
      return (-1);              /* The roots are extremely imaginary, since both a a b were 0 */
    }

    r[0] = r[1] = (-c / b);

    if (r[0] < 0.0)
      return (-2);              /* Generally speaking we are not interested in
                                   negative distances */
    else
      return (0);
  }

  qq = 4. * a * c / (b * b);

  if ((q = 1.0 - qq) < 0.0)
  {
    r[0] = r[1] = -99.;
    return (-1);                /* both roots are imaginary */
  }
  else if (fabs (qq) < 1.e-8)
  {
    r[0] = (-c / b);
    r[1] = (-b / a);
  }
  else
  {

    q = sqrt (q);
    z = 0.5 * b / a;

    r[0] = (-1.0 - q) * z;
    r[1] = (-1.0 + q) * z;
  }



  if (r[0] > 0.0 && (r[0] < r[1] || r[1] <= 0.0))
    return (0);                 /* r[0] is smallest positive root */
  if (r[1] > 0.0 && (r[1] < r[0] || r[0] <= 0.0))
    return (1);                 /* r[1] is smallest positive root */
  return (-2);                  /* both roots are negative */

  /* x1 should be the smallest positive root for most applications */
}
