namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;


	public class Honeycomb
	{
		public static Geometry GetGeometry( int p, int q, int r )
		{
			double t1 = Math.Sin( PiOverNSafe( p ) ) * Math.Sin( PiOverNSafe( r ) );
			double t2 = Math.Cos( PiOverNSafe( q ) );

			if( Tolerance.Equal( t1, t2 ) )
				return Geometry.Euclidean;

			if( Tolerance.GreaterThan( t1, t2 ) )
				return Geometry.Spherical;

			return Geometry.Hyperbolic;
		}

		/// <summary>
		/// Returns the in-radius, in the induced geometry.
		/// </summary>
		public static double InRadius( int p, int q, int r )
		{
			double pip = PiOverNSafe( p );
			double pir = PiOverNSafe( r );

			double pi_hpq = Pi_hpq( p, q );
			double inRadius = Math.Sin( pip ) * Math.Cos( pir ) / Math.Sin( pi_hpq );

			switch( GetGeometry( p, q, r ) )
			{
				case Geometry.Hyperbolic:
					return DonHatch.acosh( inRadius );
				case Geometry.Spherical:
					return Math.Acos( inRadius );
			}

			throw new System.NotImplementedException();
		}
			

		private static double Pi_hpq( int p, int q )
		{
			double pi = Math.PI;
			double pip = PiOverNSafe( p );
			double piq = PiOverNSafe( q );

			double temp = Math.Pow( Math.Cos( pip ), 2 ) + Math.Pow( Math.Cos( piq ), 2 );
			double hab = pi / Math.Acos( Math.Sqrt( temp ) );

			// Infinity safe.
			double pi_hpq = pi / hab;
			if( Infinity.IsInfinite( hab ) )
				pi_hpq = 0;

			return pi_hpq;
		}

		public static double PiOverNSafe( int n )
		{
			return n == -1 ? 0 : Math.PI / n;
		}
	}
}
