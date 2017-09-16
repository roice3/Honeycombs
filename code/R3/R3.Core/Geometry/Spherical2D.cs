namespace R3.Geometry
{
	using Math = System.Math;
	using R3.Core;
	using System.Diagnostics;

	public static class Spherical2D
	{
		// The next two methods mimic Don's stuff for the hyperbolic plane.
		public static double
		s2eNorm( double sNorm )
		{
			//if( double.IsNaN( sNorm ) )
			//	return 1.0;
			return Math.Tan( .5 * sNorm );
		}
	}
}
