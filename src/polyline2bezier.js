'use strict';
/**
 * fitCurves.js - Piecewise cubic fitting code
 *
 * original: fitCurves.c
 * https://github.com/erich666/GraphicsGems/blob/master/gems/FitCurves.c
 *
 * ported by ynakajima (https://github.com/ynakajima).
 *
 * THIS SOURCE CODE IS PUBLIC DOMAIN, and
 * is freely available to the entire computer graphics community
 * for study, use, and modification.  We do request that the
 * comment at the top of each file, identifying the original
 * author and its original publication in the book Graphics
 * Gems, be retained in all programs that use these files.
 *
 */
/**
  An Algorithm for Automatically Fitting Digitized Curves
  by Philip J. Schneider
  from "Graphics Gems", Academic Press, 1990
*/
(function() {
  var bezierSegments = [];

  if (typeof module !== 'undefined') {
    module.exports = polyline2bezier;
  } else {
    window.polyline2bezier = polyline2bezier;
  }

  function polyline2bezier(polyline, error) {
    var d = [];
    error = error || 4.0;    // Squared error
    bezierSegments = [];
    polyline.forEach(function(point) {
      d.push(new Point2(point[0], point[1]));
    });
    fitCurve(d, d.length, error); // Fit the Bezier curves
    return bezierSegments;
  }

  /**
   * drawBezierCurve
   * @param {Number} n
   * @param {Array.<Point2>} curve
   */
  function drawBezierCurve(n, curve) {
    bezierSegments.push(curve);
  } 

  /**
   * fitCurve :
   * Fit a Bezier curve to a set of digitized points 
   * @param {Array.<Point2>} d Array of digitized points.
   * @param {Number} nPts Number of digitized points
   * @param {Number} error  User-defined error squared
   */
  function fitCurve(d, nPts, error) {
    /*  Unit tangent vectors at endpoints */
    var tHat1 = new Vector2(), tHat2 = new Vector2();  

    tHat1 = computeLeftTangent(d, 0);
    tHat2 = computeRightTangent(d, nPts - 1);
    fitCubic(d, 0, nPts - 1, tHat1, tHat2, error);
  }

  /**
   * fitCubic :
   * Fit a Bezier curve to a (sub)set of digitized points
   * @param {Array.<Point2>} d Array of digitized points
   * @param {Number} first Indices of first pts in region
   * @param {Number} last Indices of last pts in region
   * @param {Point2} tHat1 Unit tangent vectors at endpoints
   * @param {Point2} tHat2 Unit tangent vectors at endpoints
   * @param {Number} error User-defined error squared
   */
  function fitCubic(d, first, last, tHat1, tHat2, error) {
    var bezCurve,           /*Control points of fitted Bezier curve*/
        u = [],             /*  Parameter values for point  */
        uPrime = [],        /*  Improved parameter values */
        maxError,           /*  Maximum fitting error   */
        splitPoint,         /*  Point to split point set at   */
        nPts,               /*  Number of points in subset  */
        iterationError,     /*Error below which you try iterating  */
        maxIterations = 4,  /*  Max times to try iterating  */
        tHatCenter = new Vector2(),/* Unit tangent vector at splitPoint */
        i;

    iterationError = error * error;
    nPts = last - first + 1;

    /*  Use heuristic if region only has two points in it */
    if (nPts === 2) {
      var dist = v2DistanceBetween2Points(d[last], d[first]) / 3.0;

      bezCurve = [];
      bezCurve[0] = d[first];
      bezCurve[3] = d[last];
      tHat1 = v2Scale(tHat1, dist);
      tHat2 = v2Scale(tHat2, dist);
      bezCurve[1] = v2Add(bezCurve[0], tHat1);
      bezCurve[2] = v2Add(bezCurve[3], tHat2);
      drawBezierCurve(3, bezCurve);
      return;
    }

    /*  Parameterize points, and attempt to fit curve */
    u = chordLengthParameterize(d, first, last);
    bezCurve = generateBezier(d, first, last, u, tHat1, tHat2);

    /*  Find max deviation of points to fitted curve */
    var resultMaxError =
      computeMaxError(d, first, last, bezCurve, u, splitPoint);
    maxError = resultMaxError.maxError;
    splitPoint = resultMaxError.splitPoint;

    if (maxError < error) {
      drawBezierCurve(3, bezCurve);
      return;
    }

    /*  If error not too large, try some reparameterization  */
    /*  and iteration */
    if (maxError < iterationError) {
      for (i = 0; i < maxIterations; i++) {
        uPrime = reparameterize(d, first, last, u, bezCurve);
        bezCurve = generateBezier(d, first, last, uPrime, tHat1, tHat2);
        resultMaxError = computeMaxError(d, first, last,
               bezCurve, uPrime, splitPoint);
        maxError = resultMaxError.maxError;
        splitPoint = resultMaxError.splitPoint;
        if (maxError < error) {
          drawBezierCurve(3, bezCurve);
          return;
        }
        u = uPrime;
      }
    }

    /* Fitting failed -- split at max error point and fit recursively */
    tHatCenter = computeCenterTangent(d, splitPoint);
    fitCubic(d, first, splitPoint, tHat1, tHatCenter, error);
    tHatCenter = v2Negate(tHatCenter);
    fitCubic(d, splitPoint, last, tHatCenter, tHat2, error);
  }


  /**
   * generateBezier :
   * Use least-squares method to find Bezier control points for region.
   * @param {Array.<Point2>} d  Array of digitized points
   * @param {Number} first Indices defining region
   * @param {Number} last Indices defining region
   * @param {Array.<Number>} uPrime Parameter values for region
   * @param {Vector2} tHat1 Unit tangents at endpoints
   * @param {Vector2} tHat2 Unit tangents at endpoints
   * @return {Array.<Point2> BezierCurve
   */
  function generateBezier(d, first, last, uPrime, tHat1, tHat2)
  {
    var i,
        A = [],              /* Precomputed rhs for eqn  */
        nPts,                /* Number of pts in sub-curve */
        C = [[], []],        /* Matrix C    */
        X = [],              /* Matrix X      */
        det_C0_C1,           /* Determinants of matrices  */
        det_C0_X,
        det_X_C1,
        alpha_l,             /* Alpha values, left and right  */
        alpha_r,
        tmp = new Vector2(), /* Utility variable    */
        bezCurve;            /* RETURN bezier curve ctl pts  */

    bezCurve = [];
    nPts = last - first + 1;


    /* Compute the A's  */
    for (i = 0; i < nPts; i++) {
      var v1 = new Vector2(tHat1.x, tHat1.y),
          v2 = new Vector2(tHat2.x, tHat2.y);
      v1 = v2Scale(v1, B1(uPrime[i]));
      v2 = v2Scale(v2, B2(uPrime[i]));
      A[i] = [];
      A[i][0] = v1;
      A[i][1] = v2;
    }

    /* Create the C and X matrices  */
    C[0][0] = 0.0;
    C[0][1] = 0.0;
    C[1][0] = 0.0;
    C[1][1] = 0.0;
    X[0]    = 0.0;
    X[1]    = 0.0;

    for (i = 0; i < nPts; i++) {
      C[0][0] += v2Dot(A[i][0], A[i][0]);
      C[0][1] += v2Dot(A[i][0], A[i][1]);
      // C[1][0] += v2Dot(A[i][0], A[i][1]);
      C[1][0] = C[0][1];
      C[1][1] += v2Dot(A[i][1], A[i][1]);

      tmp = v2SubII(d[first + i],
            v2AddII(
              v2ScaleIII(d[first], B0(uPrime[i])),
            v2AddII(
                v2ScaleIII(d[first], B1(uPrime[i])),
                    v2AddII(
                          v2ScaleIII(d[last], B2(uPrime[i])),
                            v2ScaleIII(d[last], B3(uPrime[i]))))));

      X[0] += v2Dot(A[i][0], tmp);
      X[1] += v2Dot(A[i][1], tmp);
    }

    /* Compute the determinants of C and X  */
    det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
    det_C0_X  = C[0][0] * X[1]    - C[1][0] * X[0];
    det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];

    /* Finally, derive alpha values  */
    alpha_l = (det_C0_C1 === 0) ? 0.0 : det_X_C1 / det_C0_C1;
    alpha_r = (det_C0_C1 === 0) ? 0.0 : det_C0_X / det_C0_C1;

    /* If alpha negative, use the Wu/Barsky heuristic (see text) */
    /* (if alpha is 0, you get coincident control points that lead to
     * divide by zero in any subsequent newtonRaphsonRootFind() call. */
    var segLength = v2DistanceBetween2Points(d[last], d[first]);
    var epsilon = 1.0e-6 * segLength;
    if (alpha_l < epsilon || alpha_r < epsilon) {
      /* fall back on standard (probably inaccurate) formula, */
      /* and subdivide further if needed. */
      var dist = segLength / 3.0;
      bezCurve[0] = d[first];
      bezCurve[3] = d[last];
      bezCurve[1] = v2Add(bezCurve[0], v2Scale(tHat1, dist));
      bezCurve[2] = v2Add(bezCurve[3], v2Scale(tHat2, dist));
      return (bezCurve);
    }

    /*  First and last control points of the Bezier curve are */
    /*  positioned exactly at the first and last data points */
    /*  Control points 1 and 2 are positioned an alpha distance out */
    /*  on the tangent vectors, left and right, respectively */
    bezCurve[0] = d[first];
    bezCurve[3] = d[last];
    bezCurve[1] = v2Add(bezCurve[0], v2Scale(tHat1, alpha_l));
    bezCurve[2] = v2Add(bezCurve[3], v2Scale(tHat2, alpha_r));
    return (bezCurve);
  }


  /**
   * reparameterize:
   * Given set of points and their parameterization, try to find
   * a better parameterization.
   * @param {Array.<Point2>} d Array of digitized points
   * @param {Number} first Indices defining region
   * @param {Number} last Indices defining region
   * @param {Array.<Number>} u Current parameter values
   * @param {Array.<Point2>} bezCurve Current fitted curve
   * @return {Number}
   */
  function reparameterize(d, first, last, u, bezCurve) {
    //nPts = last - first + 1,
    var i,
        uPrime = [], /*  New parameter values  */
        _bezCurve = [
          new Point2(bezCurve[0].x, bezCurve[0].y),
          new Point2(bezCurve[1].x, bezCurve[1].y),
          new Point2(bezCurve[2].x, bezCurve[2].y),
          new Point2(bezCurve[3].x, bezCurve[3].y)
        ];

    for (i = first; i <= last; i++) {
      uPrime[i-first] =
        newtonRaphsonRootFind(_bezCurve, d[i], u[i-first]);
    }
    return (uPrime);
  }

  /**
   * newtonRaphsonRootFind :
   * Use Newton-Raphson iteration to find better root.
   * @param {Array.<Point2>} _Q Current fitted curve
   * @param {Point2} _P Digitized point
   * @param {Number} u Parameter value for "P"
   * @return {Number}
   */
  function newtonRaphsonRootFind(_Q, _P, u) {
    var numerator, denominator,
        Q1 = [new Point2(), new Point2(), new Point2()], /*  Q' and Q'' */
        Q2 = [new Point2(), new Point2()],
        /*u evaluated at Q, Q', & Q'' */
        Q_u = new Point2(), Q1_u = new Point2(), Q2_u = new Point2(), 
        uPrime,    /*  Improved u */
        i,
        Q = [
          new Point2(_Q[0].x, _Q[0].y),
          new Point2(_Q[1].x, _Q[1].y),
          new Point2(_Q[2].x, _Q[2].y),
          new Point2(_Q[3].x, _Q[3].y),
        ],
        P = new Point2(_P.x, _P.y);
    
    /* Compute Q(u)  */
    Q_u = bezierII(3, Q, u);
    
    /* Generate control vertices for Q'  */
    for (i = 0; i <= 2; i++) {
      Q1[i].x = (Q[i+1].x - Q[i].x) * 3.0;
      Q1[i].y = (Q[i+1].y - Q[i].y) * 3.0;
    }
    
    /* Generate control vertices for Q'' */
    for (i = 0; i <= 1; i++) {
      Q2[i].x = (Q1[i+1].x - Q1[i].x) * 2.0;
      Q2[i].y = (Q1[i+1].y - Q1[i].y) * 2.0;
    }
    
    /* Compute Q'(u) and Q''(u)  */
    Q1_u = bezierII(2, Q1, u);
    Q2_u = bezierII(1, Q2, u);
    
    /* Compute f(u)/f'(u) */
    numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
    denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
              (Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);
    if (denominator === 0.0) return u;

    /* u = u - f(u)/f'(u) */
    uPrime = u - (numerator/denominator);
    return (uPrime);
  }

  /**
   * Bezier :
   * Evaluate a Bezier curve at a particular parameter value
   * @param {Number} degree The degree of the bezier curve
   * @param {Array.<Point2>} V Array of control points
   * @param {Number} t Parametric value to find point for
   * @return {Point2}
   */
  function bezierII(degree, V, t) {
    var i, j,    
        Q,        /* Point on curve at parameter t  */
        Vtemp = [];    /* Local copy of control points    */

    /* Copy array  */
    for (i = 0; i <= degree; i++) {
      Vtemp[i] = new Point2(V[i].x, V[i].y);
    }

    /* Triangle computation  */
    for (i = 1; i <= degree; i++) {  
      for (j = 0; j <= degree-i; j++) {
        Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j+1].x;
        Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j+1].y;
      }
    }

    Q = new Point2(Vtemp[0].x, Vtemp[0].y);
    return Q;
  }

  /*
   *  B0, B1, B2, B3 :
   *  Bezier multipliers
   */
  function B0(u)
  {
    var tmp = 1.0 - u;
    return (tmp * tmp * tmp);
  }


  function B1(u)
  {
    var tmp = 1.0 - u;
    return (3 * u * (tmp * tmp));
  }

  function B2(u)
  {
    var tmp = 1.0 - u;
    return (3 * u * u * tmp);
  }

  function B3(u)
  {
    return (u * u * u);
  }

  /**
   * computeLeftTangent, computeRightTangent, computeCenterTangent :
   * Approximate unit tangents at endpoints and "center" of digitized curve
   */
  /**
   * @param {Array.<Point2>} d Digitized points.
   * @param {Number} end Index to "left" end of region.
   * @return {Vector2}
   */
  function computeLeftTangent(d, end)
  {
    var  tHat1 = new Vector2();
    tHat1 = v2SubII(d[end+1], d[end]);
    tHat1 = v2Normalize(tHat1);
    return tHat1;
  }

  /**
   * @param {Array.<Point2>} d Digitized points.
   * @param {Number} end Index to "right" end of region.
   * @return {Vector2}
   */
  function computeRightTangent(d, end)
  {
    var  tHat2 = new Vector2();
    tHat2 = v2SubII(d[end-1], d[end]);
    tHat2 = v2Normalize(tHat2);
    return tHat2;
  }

  /**
   * @param {Array.<Point2>} d Digitized points.
   * @param {Number} end Index to point inside region.
   * @return {Vector2}
   */
  function computeCenterTangent(d, center)
  {
    var V1 = new Vector2(), V2 = new Vector2(), tHatCenter = new Vector2();

    V1 = v2SubII(d[center-1], d[center]);
    V2 = v2SubII(d[center], d[center+1]);
    tHatCenter.x = (V1.x + V2.x)/2.0;
    tHatCenter.y = (V1.y + V2.y)/2.0;
    tHatCenter = v2Normalize(tHatCenter);
    return tHatCenter;
  }

  /**
   * chordLengthParameterize :
   * Assign parameter values to digitized points 
   * using relative distances between points.
   * @param {Array.<Point2>} d Array of digitized points
   * @param {Number} first Indices defining region
   * @param {Number} last Indices defining region
   * @return {Number}
   */
  function chordLengthParameterize(d, first, last)
  {
    var i, 
        u;      /*  Parameterization    */

    u = []; 

    u[0] = 0.0;
    for (i = first+1; i <= last; i++) {
      u[i-first] = u[i-first-1] +
        v2DistanceBetween2Points(d[i], d[i-1]);
    }

    for (i = first + 1; i <= last; i++) {
      u[i-first] = u[i-first] / u[last-first];
    }
    return u;
  }

  /**
   * computeMaxError :
   * Find the maximum squared distance of digitized points
   * to fitted curve.
   * @param {Array.<Point2>} d Array of digitized points
   * @param {Number} first Indices defining region
   * @param {Number} last Indices defining region
   * @param {Array.<Point2>} bezCurve Fitted Bezier curve
   * @param {Array.<Number>} u Parameterization of points
   * @param {Number} splitPoint Point of maximum error
   */
  function computeMaxError(d, first, last, bezCurve, u, splitPoint)
  {
    var i,
        maxDist,                /*  Maximum error    */
        dist,                   /*  Current error    */
        P = new Point2(),       /*  Point on curve    */
        v = new Vector2();      /*  Vector from point to curve  */

    splitPoint = (last - first + 1)/2;
    maxDist = 0.0;
    for (i = first + 1; i < last; i++) {
      P = bezierII(3, bezCurve, u[i-first]);
      v = v2SubII(P, d[i]);
      dist = v2SquaredLength(v);
      if (dist >= maxDist) {
          maxDist = dist;
          splitPoint = i;
      }
    }
    return {maxError: maxDist, splitPoint: splitPoint};
  }

  function v2AddII(a, b)
  {
    var c = new Vector2();
    c.x = a.x + b.x;  c.y = a.y + b.y;
    return (c);
  }
  function v2ScaleIII(v, s)
  {
    var result = new Vector2();
    result.x = v.x * s; result.y = v.y * s;
    return (result);
  }

  function v2SubII(a, b)
  {
    var  c = new Vector2();
    c.x = a.x - b.x; c.y = a.y - b.y;
    return c;
  }

  // include "GraphicsGems.h"          
  /* 
   * GraphicsGems.h  
   * Version 1.0 - Andrew Glassner
   * from "Graphics Gems", Academic Press, 1990
   */

  /*********************/
  /* 2d geometry types */
  /*********************/

  function Point2Struct(x, y) {  /* 2d point */
    this.x = typeof x === 'number' ? x : 0;
    this.y = typeof y === 'number' ? y : 0;
  }
  function Point2(x, y) {
    Point2Struct.call(this, x, y);  
  }
  function Vector2(x, y) {
    Point2Struct.call(this, x, y);  
  }

  /***********************/
  /* two-argument macros */
  /***********************/
  /* 
  2d and 3d Vector C Library 
  by Andrew Glassner
  from "Graphics Gems", Academic Press, 1990
  */

  /******************/
  /*   2d Library   */
  /******************/

  /* returns squared length of input vector */  
  function v2SquaredLength(a) 
  {
    return((a.x * a.x) + (a.y * a.y));
  }
    
  /* returns length of input vector */
  function v2Length(a) 
  {
    return(Math.sqrt(v2SquaredLength(a)));
  }
    
  /* negates the input vector and returns it */
  function v2Negate(v) 
  {
    var result = new Point2();
    result.x = -v.x;  result.y = -v.y;
    return(result);
  }

  /* normalizes the input vector and returns it */
  function v2Normalize(v) 
  {
    var result = new Point2(),
        len = v2Length(v);
    if (len !== 0.0) { result.x = v.x / len;  result.y = v.y / len; }
    return(result);
  }


  /* scales the input vector to the new length and returns it */
  function v2Scale(v, newlen) 
  {
    var result = new Point2(),
        len = v2Length(v);
    if (len !== 0.0) {
      result.x = v.x * newlen/len;   result.y = v.y * newlen/len;
    }
    return(result);
  }

  /* return vector sum c = a+b */
  function v2Add(a, b)
  {
    var c = new Point2();
    c.x = a.x + b.x;  c.y = a.y + b.y;
    return(c);
  }
    
  /* return the dot product of vectors a and b */
  function v2Dot(a, b) 
  {
    return((a.x * b.x) + (a.y * b.y));
  }

  /* return the distance between two points */
  function v2DistanceBetween2Points(a, b)
  {
    var dx = a.x - b.x;
    var dy = a.y - b.y;
    return(Math.sqrt((dx*dx)+(dy*dy)));
  }
})();
