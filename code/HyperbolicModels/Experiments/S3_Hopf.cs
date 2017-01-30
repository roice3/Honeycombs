namespace HyperbolicModels
{
	using R3.Core;
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	public static class S3_Hopf
	{

		public static void Experiment()
		{
			TilingConfig config = new TilingConfig( 4, 3 );
			Tiling tiling = new Tiling();
			tiling.Generate( config );

			HashSet<H3.Cell.Edge> completed = new HashSet<H3.Cell.Edge>( new H3.Cell.EdgeEqualityComparer() );

			string fileName = "hopf.pov";
			using( StreamWriter sw = File.CreateText( fileName ) )
			{
				Tile[] tiles = tiling.Tiles.ToArray();
				//foreach( Tile t in tiling.Tiles )
				foreach( Tile t in new Tile[] { tiles[0] } )
				{
					foreach( Segment seg in t.Boundary.Segments )
					{
						H3.Cell.Edge e = new H3.Cell.Edge( seg.P1, seg.P2 );
						if( completed.Contains( e ) )
							continue;

						HopfLink( sw,
							Sterographic.PlaneToSphereSafe( e.Start ),
							Sterographic.PlaneToSphereSafe( e.End ), anti: false );
						completed.Add( e );
					}
				}
			}
		}

		/// <summary>
		/// Hopf Link between two points on S^2.
		/// </summary>
		public static void HopfLink( StreamWriter sw, Vector3D s2_1, Vector3D s2_2, bool anti )
		{
			Vector3D[] circlePoints;
			string circleString;
			circlePoints = OneHopfCircleProjected( s2_1, anti );
			circleString = PovRay.EdgeSphereSweep( circlePoints, SizeFunc );
			sw.WriteLine( circleString );
			circlePoints = OneHopfCircleProjected( s2_2, anti );
			circleString = PovRay.EdgeSphereSweep( circlePoints, SizeFunc );
			sw.WriteLine( circleString );

			Mesh mesh = new Mesh();
			Vector3D[] interpolated = S3.GeodesicPoints( s2_1, s2_2 );
			for( int i = 0; i < interpolated.Length - 1; i++ )
			{
				Vector3D v1 = interpolated[i];
				Vector3D v2 = interpolated[i + 1];
				Vector3D[] p1 = OneHopfCircleProjected( v1, anti );
				Vector3D[] p2 = OneHopfCircleProjected( v2, anti );

				for( int j = 0; j < p1.Length-1; j++ )
				{
					Mesh.Triangle t1 = new Mesh.Triangle( p1[j], p1[j + 1], p2[j] );
					Mesh.Triangle t2 = new Mesh.Triangle( p2[j], p1[j + 1], p2[j + 1] );
					mesh.Triangles.Add( t1 );
					mesh.Triangles.Add( t2 );
				}
			}

			PovRay.WriteMesh( sw, mesh, append: true );
		}

		private static Vector3D[] OneHopfCircleProjected( Vector3D s2Point, bool anti )
		{
			Vector3D[] circlePoints = S3.OneHopfCircle( s2Point, anti );
			for( int i = 0; i < circlePoints.Length; i++ )
				circlePoints[i] = Sterographic.S3toR3( circlePoints[i] );
			return circlePoints;
		}

		private static Sphere SizeFunc( Vector3D v )
		{
			Vector3D c;
			double r;
			H3Models.Ball.DupinCyclideSphere( v, .04 / 2, Geometry.Spherical, out c, out r );
			return new Sphere() { Center = c, Radius = r };
		}
	}
}
