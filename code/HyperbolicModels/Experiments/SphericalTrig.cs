namespace HyperbolicModels
{
	using R3.Core;
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.IO;
	using System.Linq;
	using System.Numerics;

	using Math = System.Math;

	public static class SphericalTrig
	{
		static bool IsGC( Circle3D c )
		{
			return c.Center.IsOrigin;
		}

		/// <summary>
		/// Return the normal of the great circle defined by the two vectors. 
		/// Return false if we fail to get the normal (happens if p1 = p2).
		/// </summary>
		static bool GetGCNormal( Vector3D p1, Vector3D p2, out Vector3D normal )
		{
			normal = p1.Cross( p2 );
			if( !normal.Normalize() )
			{
				normal = Vector3D.DneVector();
				return false;
			}
			return true;
		}

		/// <summary>
		/// Had to break the intersection method into separate cases (depending on when circles are great circles or not),
		/// because the spherical pythagorean theorem breaks down in GC cases.
		/// </summary>
		public static bool IntersectionSmart( Vector3D sphereCenter, Circle3D c1, Circle3D c2, out Vector3D i1, out Vector3D i2 )
		{
			i1 = i2 = Vector3D.DneVector();

			Circle3D clone1 = c1.Clone(), clone2 = c2.Clone();
			clone1.Center -= sphereCenter;
			clone2.Center -= sphereCenter;

			if( IsGC( clone1 ) && IsGC( clone2 ) )
			{
				if( !IntersectionGCGC( clone1.Normal, clone2.Normal, out i1, out i2 ) )
					return false;
			}
			else if( IsGC( clone1 ) || IsGC( clone2 ) )
			{
				bool firstIsGC = IsGC( clone1 );
				Vector3D gc = firstIsGC ? clone1.Normal : clone2.Normal;
				Circle3D c = firstIsGC ? clone2 : clone1;

				List<Vector3D> iPoints = new List<Vector3D>();
				if( !IntersectionCircleGC( c, gc, iPoints ) )
					return false;

				if( iPoints.Count != 2 )
					throw new System.NotImplementedException();

				i1 = iPoints[0];
				i2 = iPoints[1];
			}
			else
				throw new System.NotImplementedException();

			// Move us back to the sphere center.
			i1 += sphereCenter;
			i2 += sphereCenter;
			return true;
		}

		/// <summary>
		/// Returns the two intersection points of two great circles.
		/// gc1 and gc2 are the great circle normal vectors.
		/// Fails if the input circles are the same.
		/// </summary>
		public static bool IntersectionGCGC( Vector3D gc1, Vector3D gc2, out Vector3D i1, out Vector3D i2 )
		{
			i1 = i2 = Vector3D.DneVector();

			// Direction vector of the intersection point of the two great circles.
			// NOTE there are actually two antipodal intersection points, +-I
			Vector3D I = gc1.Cross( gc2 );
			if( !I.Normalize() )
				return false;

			i1 = I;
			i2 = -I;
			return true;
		}

		/// <summary>
		/// Returns intersection points between a circle and a great circle.
		/// There may be 0, 1, or 2 intersection points.
		/// Returns false if the circle is the same as gc.
		/// </summary>
		static bool IntersectionCircleGC( Circle3D c, Vector3D gc, List<Vector3D> iPoints )
		{
			double radiusCosAngle = CosAngle( c.Normal, c.PointOnCircle );
			double radiusAngle = Math.Acos( radiusCosAngle );
			double radius = radiusAngle;	// Since sphere radius is 1.

			// Find the great circle perpendicular to gc, and through the circle center.
			Vector3D gcPerp = c.Normal.Cross( gc );
			if( !gcPerp.Normalize() )
			{
				// Circles are parallel => Zero or infinity intersections.
				if( Tolerance.Equal( radius, Math.PI / 2 ) )
					return false;

				return true;
			}

			// Calculate the offset angle from the circle normal to the gc normal.
			double offsetAngle = c.Normal.AngleTo( gc );
			if( Tolerance.GreaterThan( offsetAngle, Math.PI / 2 ) )
			{
				gc *= -1;
				offsetAngle = c.Normal.AngleTo( gc );
			}
			double coAngle = Math.PI / 2 - offsetAngle;

			// No intersections.
			if( radiusAngle < coAngle )
				return true;

			// Here is the perpendicular point on the great circle.
			Vector3D pointOnGC = c.Normal;
			pointOnGC.RotateAboutAxis( gcPerp, coAngle );

			// 1 intersection.
			if( Tolerance.Equal( radiusAngle, coAngle ) )
			{
				iPoints.Add( pointOnGC );
				return true;
			}

			// 2 intersections.

			// Spherical pythagorean theorem
			// http://en.wikipedia.org/wiki/Pythagorean_theorem#Spherical_geometry
			// We know the hypotenuse and one side.  We need the third leg.
			// We do this calculation on a unit sphere, to get the result as a normalized cosine of an angle.
			double sideCosA = radiusCosAngle / Math.Cos( coAngle );
			double rot = Math.Acos( sideCosA );

			Vector3D i1 = pointOnGC, i2 = pointOnGC;
			i1.RotateAboutAxis( gc, rot );
			i2.RotateAboutAxis( gc, -rot );
			iPoints.Add( i1 );
			iPoints.Add( i2 );

			Circle3D test = new Circle3D { Normal = gc, Center = new Vector3D(), Radius = 1 };


			return true;
		}

		/// <summary>
		/// NOTE: Not general, and assumes some things we know about this problem domain, 
		/// e.g. that c1 and c2 live on the same sphere of radius 1, and have two intersection points.
		/// </summary>
		public static void IntersectionCircleCircle( Vector3D sphereCenter, Circle3D c1, Circle3D c2, out Vector3D i1, out Vector3D i2 )
		{
			// Spherical analogue of our flat circle-circle intersection.
			// Spherical pythagorean theorem for sphere where r=1: cos(hypt) = cos(A)*cos(B) 

			Circle3D clone1 = c1.Clone(), clone2 = c2.Clone();
			//clone1.Center -= sphereCenter;
			//clone2.Center -= sphereCenter;

			// Great circle (denoted by normal vector), and distance between the centers.
			Vector3D gc = clone2.Normal.Cross( clone1.Normal );
			double d = clone2.Normal.AngleTo( clone1.Normal );
			double r1 = clone1.Normal.AngleTo( clone1.PointOnCircle );
			double r2 = clone2.Normal.AngleTo( clone2.PointOnCircle );

			// Calculate distances we need.  So ugly!
			// http://www.wolframalpha.com/input/?i=cos%28r1%29%2Fcos%28r2%29+%3D+cos%28x%29%2Fcos%28d-x%29%2C+solve+for+x
			double t1 = Math.Pow( Math.Tan( d / 2 ), 2 );
			double t2 = Math.Cos( r1 ) / Math.Cos( r2 );
			double t3 = Math.Sqrt( (t1 + 1) * (t1 * t2 * t2 + 2 * t1 * t2 + t1 + t2 * t2 - 2 * t2 + 1) ) - 2 * t1 * t2;
			double x = 2 * Math.Atan( t3 / (t1 * t2 + t1 - t2 + 1) );
			double y = Math.Acos( Math.Cos( r1 ) / Math.Cos( x ) );

			i1 = clone1.Normal;
			i1.RotateAboutAxis( gc, x );
			i2 = i1;

			// Perpendicular to gc through i1.
			Vector3D gc2 = i1.Cross( gc );
			i1.RotateAboutAxis( gc2, y );
			i2.RotateAboutAxis( gc2, -y );
			i1 += sphereCenter;
			i2 += sphereCenter;

			/*
			// It would be nice to do the spherical analogue of circle-circle intersections, like here:
			// http://mathworld.wolfram.com/Circle-CircleIntersection.html
			// But I don't want to jump down that rabbit hole and am going to sacrifice some speed to use
			// my existing euclidean function.

			// Stereographic projection to the plane.  XXX - Crap, circles may become lines, and this isn't being handled well.
			Circle3D c1Plane = H3Models.BallToUHS( clone1 );
			Circle3D c2Plane = H3Models.BallToUHS( clone2 );
			if( 2 != Euclidean2D.IntersectionCircleCircle( c1Plane.ToFlatCircle(), c2Plane.ToFlatCircle(), out i1, out i2 ) )
				throw new System.Exception( "Expected two intersection points" );
			i1 = H3Models.UHSToBall( i1 ); i1 += sphereCenter;
			i2 = H3Models.UHSToBall( i2 ); i2 += sphereCenter;
			*/
		}

		/// <summary>
		/// Calculate the cosine of the angle between two vectors.
		/// </summary>
		private static double CosAngle( Vector3D p1, Vector3D p2 )
		{
			double cosA = p1.Dot( p2 ) / (p1.Abs() * p2.Abs());
			return Clamp( cosA );
		}

		/// <summary>
		/// Calculate the sine of the angle between two vectors.
		/// </summary>
		private static double SinAngle( Vector3D p1, Vector3D p2 )
		{
			double sinA = (p1.Cross( p2 )).Abs() / (p1.Abs() * p2.Abs());
			return Clamp( sinA );
		}

		/// <summary>
		/// Clamps the input to [-1,1].
		/// This is to avoid floating point issues when taking arcsin and arccos.
		/// </summary>
		private static double Clamp( double val )
		{
			if( val < -1.0 )
				val = -1.0;
			if( val > 1.0 )
				val = 1.0;
			return val;
		}
	}
}
