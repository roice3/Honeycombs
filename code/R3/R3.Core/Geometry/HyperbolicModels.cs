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
		UpperHalfPlane
	}

	public class HyperbolicModels
	{
		public static Vector3D PoincareToKlein( Vector3D p )
		{
			double mag = 2 / (1 + p.Dot( p ));
			return p * mag;
		}

		public static Vector3D KleinToPoincare( Vector3D k )
		{
			double dot = k.Dot( k );
			if( dot > 1 )	// This avoids some NaN problems I saw.
				dot = 1;
			double mag = (1 - Math.Sqrt( 1 - dot )) / dot;
			return k * mag;
		}
	}
}
