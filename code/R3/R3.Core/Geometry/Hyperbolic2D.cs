namespace R3.Geometry
{
	using R3.Math;

	public static class Hyperbolic2D
	{
		/// <summary>
		/// Offsets a vector by a hyperbolic distance.
		/// </summary>
		public static Vector3D Offset( Vector3D v, double hDist )
		{
			double mag = v.Abs();
			mag = DonHatch.h2eNorm( DonHatch.e2hNorm( mag ) + hDist );
			v.Normalize();
			v *= mag;
			return v;
		}
	}
}
