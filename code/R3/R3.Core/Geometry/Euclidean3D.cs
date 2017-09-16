namespace R3.Geometry
{
	using R3.Core;
	using Math = System.Math;

	public static class Euclidean3D
	{
		public static double DistancePointLine( Vector3D n1, Vector3D p1, Vector3D point )
		{
			// Check to make sure that n1 is not degenerate.
			if( Tolerance.Zero( n1.MagSquared() ) )
				return double.NaN;

			// ZZZ - Can we make this a signed distance?
			return ( ( point - p1 ).Cross( n1 ) ).Abs() / n1.Abs();
		}

		public static double DistancePointPlane( Vector3D normalVector, Vector3D planePoint, Vector3D point )
		{
			// Check to make sure that plane is not degenerate.
			if( Tolerance.Zero( normalVector.MagSquared() ) )
				return double.NaN;

			// Here is the distance (signed depending on which side of the plane we are on).
			return ( point - planePoint ).Dot( normalVector ) / normalVector.Abs();
		}

		public static Vector3D ProjectOntoLine( Vector3D nl, Vector3D pl, Vector3D point )
		{
			// http://gamedev.stackexchange.com/a/72529
			// A + dot(AP,AB) / dot(AB,AB) * AB
			Vector3D AP = point - pl;
			Vector3D AB = nl;
			return pl + AB * AP.Dot( AB ) / AB.Dot( AB );
		}

		public static Vector3D ProjectOntoPlane( Vector3D normalVector, Vector3D planePoint, Vector3D point )
		{
			if( !normalVector.Normalize() )
				throw new System.ArgumentException( "Invalid normal vector." );

			double dist = DistancePointPlane( normalVector, planePoint, point );
			normalVector *= dist;
			return point - normalVector;
		}

		/// <summary>
		/// Calculate a plane normal after a transformation function is applied
		/// to the points.
		/// </summary>
		public static Vector3D NormalFrom3Points( Vector3D p1, Vector3D p2, Vector3D p3, 
			System.Func<Vector3D, Vector3D> transform )
		{
			Vector3D p1t = transform( p1 );
			Vector3D p2t = transform( p2 );
			Vector3D p3t = transform( p3 );
			return NormalFrom3Points( p1t, p2t, p3t );
		}

		public static Vector3D NormalFrom3Points( Vector3D p1, Vector3D p2, Vector3D p3 )
		{
			Vector3D v1 = p1 - p3;
			Vector3D v2 = p2 - p3;
			Vector3D normal = v1.Cross( v2 );
			normal.Normalize();
			return normal;
		}

		public static bool Coplanar( Vector3D[] points )
		{
			throw new System.NotImplementedException();
		}
	}
}
