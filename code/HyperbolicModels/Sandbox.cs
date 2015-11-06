namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Drawing;
	using System.IO;
	using System.Linq;
	using System.Text.RegularExpressions;
	using R3.Core;
	using R3.Geometry;
	using R3.Math;

	using Math = System.Math;
	
	internal static class Sandbox
	{
		public static void Test()
		{
			S3.HopfOrbit();

			Mobius m = new Mobius();
			m.UpperHalfPlane();

			Vector3D test = m.Apply( new Vector3D() );
			test *= 1;
		}

		public static void ResizeImages()
		{
			Size tileSize = new Size( 500, 500 );

			IEnumerable<string> files = Directory.EnumerateFiles( @"./", "*.png" );
			foreach( string file in files )
			{
				Bitmap newImage;
				using( Bitmap image = new Bitmap( file ) )
					newImage = new Bitmap( image, tileSize );
				newImage.Save( file );
			}
		}

		internal static void TileImages()
		{
			List<HoneycombAndView> images = Program.GetImageSet().ToList();
			string[] inputImages = images.Select( i => i.FormatFilename() ).ToArray();

			ImageGrid imageGrid = new ImageGrid();
			imageGrid.Generate( new ImageGrid.Settings()
			{
				//Directory = @"C:\Users\hrn\Documents\roice\honeycombs\test pics",
				Directory = @"./",
				InputImages = inputImages
			} );
		}

		internal static void CalcSelfSimilarityScale()
		{
			double inRadius = DonHatch.h2eNorm( Honeycomb.InRadius( 4, 3, 7 ) );
			Vector3D facePoint = new Vector3D( 0, 0, -inRadius );
			Sphere s = H3Models.Ball.OrthogonalSphereInterior( facePoint );
			Vector3D facePoint2 = new Vector3D( 0, 0, inRadius );
			facePoint2 = s.ReflectPoint( facePoint2 );

			facePoint = H3Models.BallToUHS( facePoint );
			facePoint2 = H3Models.BallToUHS( facePoint2 );

			double scale = facePoint.Z / facePoint2.Z;
			scale += 0;
		}

		public static void Check_pq_Distances()
		{
			double dist = Geometry2D.GetNormalizedCircumRadius( 4, 5 );
			double offset = ( 1 - dist ) / 4;

			for( int i = 1; i <= 3; i++ )
			{
				double test = dist + i*offset;
				double min = double.MaxValue;
				int found = -1;
				for( int q = 6; q < 10000; q++ )
				{
					double compare = Geometry2D.GetNormalizedCircumRadius( 4, q );
					double diff = Math.Abs( compare - test );
					if( diff < min )
					{
						min = diff;
						found = q;
					}
				}

				System.Diagnostics.Trace.WriteLine( string.Format( "{0}:{1}:{2}", i, found, min ) );
			}
		}

		public static void Cell633()
		{
			TilingConfig config = new TilingConfig( 6, 3, maxTiles: 20000 );
			Tiling tiling = new Tiling();
			tiling.GenerateInternal( config, Polytope.Projection.VertexCentered );

			double edgeLength = Honeycomb.EdgeLength( 6, 3, 3 );

			double z = 0.25;
			double offset = H3Models.UHS.ToEHorizontal( edgeLength, z );
			double scale = offset / tiling.Tiles.First().Boundary.Segments.First().Length;
			foreach( Tile tile in tiling.Tiles )
				tile.Transform( Mobius.Scale( scale ) );

			Vector3D dummy;
			double radius;
			H3Models.UHS.Geodesic( new Vector3D( 0, 0, z ), new Vector3D( scale, 0, z ), out dummy, out radius );
			Vector3D midradius = H3Models.UHSToBall( new Vector3D( 0, 0, radius ) );
			double temp = midradius.Z;
			double temp2 = ( 1 - temp ) / 2;
			double temp3 = temp + temp2;
			double temp4 = temp3;

			Vector3D circumradius = H3Models.UHSToBall( new Vector3D( 0, 0, z ) );
			temp = circumradius.Z;
			temp2 = ( 1 - temp ) / 2;
			temp3 = temp + temp2;
			temp4 = temp3;

			// Checking
			/*
			Vector3D test = new Vector3D( offset, 0, z );
			test = H3Models.UHSToBall( test );
			double edgeLength2 = DonHatch.e2hNorm( test.Abs() );
			edgeLength2 += 0;
			*/

			HashSet<H3.Cell.Edge> edges = new HashSet<H3.Cell.Edge>();
			foreach( Tile tile in tiling.Tiles )
				foreach( Segment seg in tile.Boundary.Segments )
				{
					H3.Cell.Edge edge = new H3.Cell.Edge(
						H3Models.UHSToBall( seg.P1 + new Vector3D( 0, 0, z ) ),
						H3Models.UHSToBall( seg.P2 + new Vector3D( 0, 0, z ) ) );
					edges.Add( edge );
				}

			PovRay.WriteH3Edges( new PovRay.Parameters(), edges.ToArray(), "edges.pov" );
		}

		// https://plus.google.com/u/0/117663015413546257905/posts/BnCEkdNiTZ2
		public static void TwinDodecs()
		{
			Tiling tiling = new Tiling();
			TilingConfig config = new TilingConfig( 5, 3 );
			tiling.GenerateInternal( config, Polytope.Projection.VertexCentered );	// Vertex-centered makes infinities tricky

			Dodec dodec = new Dodec();
			foreach( Tile tile in tiling.Tiles )
			foreach( Segment seg in tile.Boundary.Segments )
			{
				Vector3D p1 = seg.P1, p2 = seg.P2;
				if( Infinity.IsInfinite( p1 ) )
					p1 = Infinity.InfinityVector;
				if( Infinity.IsInfinite( p2 ) )
					p2 = Infinity.InfinityVector;

				dodec.Verts.Add( p1 );
				dodec.Verts.Add( p2 );
				dodec.Midpoints.Add( Halfway( p1, p2 ) );
			}

			// Now recursively add more vertices.
			HashSet<Vector3D> allVerts = new HashSet<Vector3D>();
			foreach( Vector3D v in dodec.Verts )
				allVerts.Add( v );
			RecurseTwins( allVerts, dodec, 0 );

			using( StreamWriter sw = File.CreateText( "dual_dodecs_points_sphere.pov" ) )
			{
				foreach( Vector3D vert in allVerts )
				{
					Vector3D onSphere = Sterographic.PlaneToSphereSafe( vert );
					sw.WriteLine( PovRay.Sphere( new Sphere() { Center = onSphere, Radius = 0.01 } ) );

					//if( !Infinity.IsInfinite( vert ) )
					//	sw.WriteLine( PovRay.Sphere( new Sphere() { Center = vert, Radius = 0.01 } ) );
				}
			}
		}

		private static Vector3D Halfway( Vector3D v1, Vector3D v2 )
		{
			Vector3D v1_ = Sterographic.PlaneToSphereSafe( v1 );
			Vector3D v2_ = Sterographic.PlaneToSphereSafe( v2 );

			Vector3D result = ( v1_ + v2_ ) / 2;
			result.Normalize();
			return Sterographic.SphereToPlane( result );
		}

		class Dodec
		{
			public readonly HashSet<Vector3D> Verts = new HashSet<Vector3D>();
			public readonly HashSet<Vector3D> Midpoints = new HashSet<Vector3D>();
		}

		private static void RecurseTwins( HashSet<Vector3D> allVerts, Dodec dodec, int level )
		{
			level++;
			if( level > 1 )
				return;

			//foreach( Vector3D v in dodec.Midpoints )
			foreach( Vector3D v in dodec.Verts )
			{
				Dodec dual = GetDual( dodec, v );
				int count = allVerts.Count;
				foreach( Vector3D dualV in dual.Verts )
					allVerts.Add( dualV );
				if( count != allVerts.Count )
					RecurseTwins( allVerts, dual, level );
			}
		}

		// Rotates a dodec about a vertex or edge to get a dual dodec.
		private static Dodec GetDual( Dodec dodec, Vector3D rotationPoint )
		{
			//double rot = System.Math.PI / 2;	// Edge-centered
			double rot = 4 * ( - Math.Atan( ( 2 + Math.Sqrt( 5 ) - 2 * Math.Sqrt( 3 + Math.Sqrt( 5 ) ) ) / Math.Sqrt( 3 ) ) );	// Vertex-centered

			Mobius m = new Mobius();
			if( Infinity.IsInfinite( rotationPoint ) )
				m.Elliptic( Geometry.Spherical, new Vector3D(), -rot );
			else
				m.Elliptic( Geometry.Spherical, rotationPoint, rot );

			Dodec dual = new Dodec();
			foreach( Vector3D v in dodec.Verts )
			{
				Vector3D rotated = m.ApplyInfiniteSafe( v );
				dual.Verts.Add( rotated );
			}

			foreach( Vector3D v in dodec.Midpoints )
			{
				Vector3D rotated = m.ApplyInfiniteSafe( v );
				dual.Midpoints.Add( rotated );
			}

			return dual;
		}

		// http://math.stackexchange.com/questions/142112/how-to-construct-a-k-regular-graph
		// k-regular graph on n vertices.
		public static void Graph()
		{
			HashSet<Edge> edges = new HashSet<Edge>();

			int k = 24, n = 350;
			//int k = 7, n = 24;

			if( k%2 == 0 )
			{
				int increment = n / ( k + 1 );
				increment = 1;

				for( int i=0; i<n; i++ )
				for( int j=1; j<=k/2; j++ )
				{
					int t1 = Clamp( i + j*increment, n );
					int t2 = Clamp( i - j*increment, n );

					Edge e1 = new Edge( i, t1 );
					Edge e2 = new Edge( i, t2 );
					edges.Add( e1 );
					edges.Add( e2 );
				}
			}
			else
			{
				int m = k / 2;
				for( int i=0; i<n; i++ )
				{
					for( int j=1; j<=m; j++ )
					{
						int t1 = Clamp( i + j, n );
						int t2 = Clamp( i - j, n );

						Edge e1 = new Edge( i, t1 );
						Edge e2 = new Edge( i, t2 );
						edges.Add( e1 );
						edges.Add( e2 );
					}

					int t3 = Clamp( i + n/2, n );
					edges.Add( new Edge( i, t3 ) );
				}
			}

			using( StreamWriter sw = File.CreateText( "733.csv" ) )
			{
				foreach( Edge e in edges )
				{
					sw.WriteLine( string.Format( "master{0};master{1}", e.V1, e.V2 ) );
				}
			}
		}

		private static int Clamp( int i, int n )
		{
			if( i < 0 )
				i += n;
			if( i >= n )
				i -= n;
			return i;
		}

		
	}
}
