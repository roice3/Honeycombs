namespace R3.Geometry
{
	using R3.Math;
	using System;
	using System.Numerics;

	public enum SphericalModel
	{
		Sterographic,
		Gnomonic
	}

	public class SphericalModels
	{
		public static Vector3D StereoToGnomonic( Vector3D p )
		{
			Vector3D sphere = Sterographic.PlaneToSphere( p );

			// We can't only represent the lower hemisphere.
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

		private static double StereoToEqualArea( double dist )
		{
			if( Infinity.IsInfinite( dist ) )
				return 1;

			double dot = dist * dist; // X^2 + Y^2 + Z^2
			double w = (dot - 1) / (dot + 1);

			double r = Math.Sqrt( 2 * (1 + w) );
			return r/2;
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

		private static double m_gScale = 0.5;
	}
}
