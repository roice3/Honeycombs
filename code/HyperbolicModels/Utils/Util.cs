namespace R3.Geometry
{
	using R3.Core;
	using Math = System.Math;

	public static class Util
	{
		public static double PiOverNSafe( int n )
		{
			return Honeycomb.PiOverNSafe( n );
		}

		public static Geometry GetGeometry( int p, int q, int r )
		{
			return Honeycomb.GetGeometry( p, q, r );
		}
	}
}
