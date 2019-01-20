namespace R3.Geometry
{
	using R3.Math;
	using System;
	using System.Numerics;

	public enum SphericalModel
	{
		Sterographic,
		Gnomonic,
		Azimuthal_Equidistant,
		Azimuthal_EqualArea,
		Equirectangular,
		Mercator,
		Orthographic,
		Sinusoidal,
		PeirceQuincuncial,
	}

	public class SphericalModels
	{
		public static Vector3D StereoToGnomonic( Vector3D p )
		{
			Vector3D sphere = Sterographic.PlaneToSphere( p );

			// We can only represent the lower hemisphere.
			if( sphere.Z >= 0 )
			{
				sphere.Z = 0;
				sphere.Normalize();
				sphere *= Infinity.FiniteScale;
				return sphere;
			}

			double z = sphere.Z;
			sphere.Z = 0;
			return -sphere * m_gScale / z;
		}

		public static Vector3D GnomonicToStereo( Vector3D g )
		{
			g /= m_gScale;
			double dot = g.Dot( g ); // X^2 + Y^2
			double z = -1 / Math.Sqrt( dot + 1 );
			return g*z / (z - 1);
		}

		public static Vector3D StereoToEqualVolume( Vector3D p )
		{
			Vector3D result = p;
			result.Normalize();
			result *= StereoToEqualVolume( p.Abs() );
			return result;
		}

		private static double StereoToEqualVolume( double dist )
		{
			if( Infinity.IsInfinite( dist ) )
				return 1;

			double dot = dist * dist; // X^2 + Y^2 + Z^2
			double w = (dot - 1) / (dot + 1);

			w = -w;	// Because I derived formula from north pole.
			double t = Math.PI / 2 - w * Math.Sqrt( 1 - w * w ) - Math.Asin( w );
			double r = Math.Pow( t * 3 / 2, 1.0 / 3 );
			return r;
		}

		public static Vector3D StereoToEqualArea( Vector3D p )
		{
			Vector3D result = p;
			result.Normalize();
			result *= StereoToEqualArea( p.Abs() );
			return result;
		}

		/// <summary>
		/// https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
		/// </summary>
		private static double StereoToEqualArea( double dist )
		{
			if( Infinity.IsInfinite( dist ) )
				return 1;

			double dot = dist * dist; // X^2 + Y^2 + Z^2
			double w = (dot - 1) / (dot + 1);

			double r = Math.Sqrt( 2 * (1 + w) );
			return r/2;
		}

		public static Vector3D EqualAreaToStereo( Vector3D p )
		{
			Vector3D result = p;
			result.Normalize();
			result *= EqualAreaToStereo( p.Abs() );
			return result;
		}

		/// <summary>
		/// https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
		/// </summary>
		private static double EqualAreaToStereo( double dist )
		{
			if( dist > 1 )
				throw new System.ArgumentException();

			// We have dist normalized between 0 and 1, so this formula is slightly 
			// different than on Wikipedia, where dist ranges up to 2.
			Vector3D v = new Vector3D( 1, 2*Math.Acos( dist ), 0 );
			v = Sterographic.SphereToPlane( SphericalCoords.SphericalToCartesian( v ) );
			return v.Abs();
		}

		public static Vector3D StereoToEquidistant( Vector3D p )
		{
			Vector3D result = p;
			result.Normalize();
			result *= StereoToEquidistant( p.Abs() );
			return result;
		}

		private static double StereoToEquidistant( double dist )
		{
			if( Infinity.IsInfinite( dist ) )
				return 1;

			double dot = dist * dist; // X^2 + Y^2 + Z^2
			double w = (dot - 1) / (dot + 1);

			double x = Math.Sqrt( 1 - w * w );
			double r = Euclidean2D.AngleToCounterClock( new Vector3D( 0, -1 ), new Vector3D( x, w ) );
			return r / Math.PI;
		}

		public static Vector3D EquidistantToStereo( Vector3D p )
		{
			Vector3D result = p;
			result.Normalize();
			result *= EquidistantToStereo( p.Abs() );
			return result;
		}

		private static double EquidistantToStereo( double dist )
		{
			if( dist > 1 )
				throw new System.ArgumentException();

			Vector3D v = new Vector3D( 0, -1 );
			v.RotateXY( dist * Math.PI );
			v = Sterographic.SphereToPlane( new Vector3D( v.X, 0, v.Y ) );
			return v.Abs();
		}

		public static Vector3D EquirectangularToStereo( Vector3D v )
		{
			// http://mathworld.wolfram.com/EquirectangularProjection.html
			// y is the latitude
			// x is the longitude
			// Assume inputs go from -1 to 1.
			Vector3D spherical = new Vector3D( 1, Math.PI / 2 * (1 - v.Y), v.X * Math.PI );
			Vector3D onBall = SphericalCoords.SphericalToCartesian( spherical );
			return Sterographic.SphereToPlane( onBall );
		}

		public static Vector3D SinusoidalToStereo(Vector3D v)
		{
			double lat = Math.PI / 2 * ( 1 - v.Y );
			Vector3D spherical = new Vector3D( 1, lat, Math.PI * v.X / Math.Cos( lat - Math.PI / 2 ) );
			Vector3D onBall = SphericalCoords.SphericalToCartesian( spherical );
			return Sterographic.SphereToPlane( onBall );
		}

		/// <summary>
		/// 2-dimensional function.
		/// http://archive.bridgesmathart.org/2013/bridges2013-217.pdf
		/// </summary>
		public static Vector3D MercatorToStereo( Vector3D v )
		{
			v *= Math.PI;	// Input is [-1,1]
			double lat = 2 * Math.Atan( Math.Exp( v.Y ) ) - Math.PI / 2;
			double inclination = lat + Math.PI / 2;
			Vector3D spherical = new Vector3D( 1, inclination, v.X );
			Vector3D onBall = SphericalCoords.SphericalToCartesian( spherical );
			return Sterographic.SphereToPlane( onBall );
		}

		/// <summary>
		/// 2-dimensional function.
		/// ZZZ - Should make this general.
		/// </summary>
		public static Vector3D OrthographicToStereo( Vector3D v )
		{
			// We can only do the projection for half of the sphere.
			double t = v.X * v.X + v.Y * v.Y;
			if( t > 1 )
				t = 1;
			v.Z = Math.Sqrt( 1 - t );
			return Sterographic.SphereToPlane( v );
		}

		private static double m_gScale = 0.5;
	}
}