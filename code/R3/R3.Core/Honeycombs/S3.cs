namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	public class S3
	{
		/// <summary>
		/// Inputs and Outputs are in R3 (stereographically projected).
		/// </summary>
		public static Vector3D[] GeodesicPoints( Vector3D v1, Vector3D v2 )
		{
			Vector3D start = Sterographic.R3toS3( v1 );
			Vector3D end = Sterographic.R3toS3( v2 );

			int div = 42;
			//int div = 56;		// 343
			//int div = 50;		// 333
			Segment seg = Segment.Line( start, end );
			Vector3D[] result = seg.Subdivide( div );
			for( int i=0; i<result.Length; i++ )
			{
				result[i].Normalize();
				result[i] = Sterographic.S3toR3( result[i] );
			}

			return result;
		}
	}
}
