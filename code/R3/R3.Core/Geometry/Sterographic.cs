namespace R3.Geometry
{
	using R3.Geometry;
	using Math = System.Math;

	public static class Sterographic
	{

		public static Vector3D R3toS3( Vector3D p )
		{
			if( Infinity.IsInfinite( p ) )
				return new Vector3D( 0, 0, 0, 1 );

			p.W = 0;
			double dot = p.Dot( p ); // X^2 + Y^2 + Z^2
			return new Vector3D(
				2 * p.X / ( dot + 1 ),
				2 * p.Y / ( dot + 1 ),
				2 * p.Z / ( dot + 1 ),
				( dot - 1 ) / ( dot + 1 ) );
		}

		public static Vector3D S3toR3( Vector3D p )
		{
			double w = p.W;
			return new Vector3D(
				p.X / ( 1 - w ),
				p.Y / ( 1 - w ),
				p.Z / ( 1 - w ) );
		}
	}
}
