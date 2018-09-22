namespace R3.Geometry
{
	using R3.Math;
	using System;
	using System.Numerics;

	public enum HyperbolicModel
	{
		Poincare,
		Klein,
		Pseudosphere,
		Hyperboloid,
		Band,
		UpperHalfPlane,
		Orthographic,
		Square,
		InvertedPoincare,
	}

	public class HyperbolicModels
	{
		public static Vector3D PoincareToKlein( Vector3D p )
		{
			double mag = 2 / (1 + p.Dot( p ));
			return p * mag;
		}

		public static double KleinToPoincare( double magSquared )
		{
			double dot = magSquared;
			if (dot > 1)    // This avoids some NaN problems I saw.
				dot = 1;
			return (1 - Math.Sqrt( 1 - dot )) / dot;
		}

		public static Vector3D KleinToPoincare( Vector3D k )
		{
			double dot = k.Dot( k );
			return k * KleinToPoincare( dot );
		}

		public static Mobius Upper
		{
			get
			{
				Cache();
				return m_upper;
			}
		}
		public static Mobius UpperInv
		{
			get
			{
				Cache();
				return m_upperInv;
			}
		}

		/// <summary>
		/// This was needed for performance.  We don't want this Mobius transform calculated repeatedly.
		/// </summary>
		private static void Cache()
		{
			if( m_cached )
				return;

			Mobius m1 = new Mobius(), m2 = new Mobius();
			m2.Isometry( Geometry.Euclidean, 0, new Complex( 0, -1 ) );
			m1.UpperHalfPlane();
			m_upper = m2 * m1;
			m_upperInv = m_upper.Inverse();

			m_cached = true;
		}
		private static bool m_cached = false;
		private static Mobius m_upper, m_upperInv;


		public static Vector3D PoincareToUpper( Vector3D v )
		{
			v = Upper.Apply( v );
			return v;
		}

		public static Vector3D UpperToPoincare( Vector3D v )
		{
			v = UpperInv.Apply( v );
			return v;
		}

		public static Vector3D PoincareToOrtho( Vector3D v )
		{
			// This may not be correct.
			// Should probably project to hyperboloid, then remove z coord.
			return SphericalModels.StereoToGnomonic( v );
		}

		public static Vector3D OrthoToPoincare( Vector3D v )
		{
			return SphericalModels.GnomonicToStereo( v );
		}


	}
}
