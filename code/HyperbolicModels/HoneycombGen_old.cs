﻿namespace HyperbolicModels
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;
	using R3.Math;

	/// <summary>
	/// This is an attempt to cleanup the HoneycombGen code to the latest/cleanest impl.
	/// I don't want to remove all this code because I might occasionally find it useful,
	/// but consider this outdated.
	/// </summary>
	class HoneycombGen_old
	{
		/// <summary>
		/// This generates a honeycomb by reflecting in all facets of a cell.
		/// </summary>
		public static void OneHoneycombOldCode()
		{
			H3.GenHoneycomb( EHoneycomb.H435 );
			//H3.GenHoneycomb( EHoneycomb.H337 );
			//H3.GenHoneycomb( EHoneycomb.H436 );
			//H3.GenHoneycomb( EHoneycomb.H536 );
			//H3.GenHoneycomb( EHoneycomb.H444 );
			//H3.GenHoneycomb( EHoneycomb.H363 );
			//H3.GenHoneycomb( EHoneycomb.H636 );
			//H3.GenHoneycomb( EHoneycomb.H337 );
		}

		private static H3.Cell.Edge[] Cull120Cell( H3.Cell.Edge[] edges )
		{
			Func<Vector3D, bool> passes = new Func<Vector3D, bool>( v =>
			{
				//return 
				//	Math.Pow( v.Z, 2 ) < 0.13 &&
				//	Math.Pow( v.W, 2 ) < 0.13;
				return Tolerance.Equal( v.W, 0.0 );
			} );

			H3.Cell.Edge[] result = edges.Where( e =>
			{
				Vector3D start = Sterographic.R3toS3( e.Start );
				Vector3D end = Sterographic.R3toS3( e.End );
				return passes( start ) && passes( end );
			} ).ToArray();

			// Now cull valence-2 edges.
			//result = CullValence2Edges( result );

			return result;
		}

		private static H3.Cell.Edge[] CullValence2Edges( H3.Cell.Edge[] edges )
		{
			List<H3.Cell.Edge> needRemoval = new List<H3.Cell.Edge>();

			// Info we'll need to remove dangling edges.
			Dictionary<Vector3D, int> vertexCounts = new Dictionary<Vector3D, int>();
			foreach( H3.Cell.Edge edge in edges )
			{
				CheckAndAdd( vertexCounts, edge.Start );
				CheckAndAdd( vertexCounts, edge.End );
			}

			foreach( H3.Cell.Edge e in edges )
			{
				if( vertexCounts[e.Start] == 2 ||
					vertexCounts[e.End] == 2 )
				{
					needRemoval.Add( e );
				}
			}

			return edges.Except( needRemoval ).ToArray();
		}

		private static void CheckAndAdd( Dictionary<Vector3D, int> vertexCounts, Vector3D v )
		{
			int count;
			if( vertexCounts.TryGetValue( v, out count ) )
				count++;
			else
				count = 1;

			vertexCounts[v] = count;
		}

		/// <summary>
		/// This generates a honeycomb by reflecting in 4 mirrors of the fundamental simplex.
		/// This "new" method is now old.
		/// </summary>
		public static void OneHoneycombNew( HoneycombDef imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			double thickness = 0.05;
			double thicknessSpherical = Spherical2D.s2eNorm( thickness );
			double thicknessHyperbolic = R3.Math.DonHatch.h2eNorm( thickness );
			double threshold = 1;

			H3.Cell.Edge[] edges = null;
			H3.Cell[] cellsToHighlight = null;
			Sphere[] simplex = null;
			Vector3D vertex = new Vector3D();

			Geometry g = Util.GetGeometry( p, q, r );
			if( g == Geometry.Spherical )
			{
				thickness = thicknessSpherical /*.07 for 333*/  /* 0.05for 433*/  /*.025 for 533,335*/;
				threshold = 10000;

				simplex = SimplexCalcs.MirrorsSpherical( p, q, r );
				vertex = SimplexCalcs.VertexSpherical( p, q, r );

				// Ugly special casing for 333, since it has a vertex project to infinity.
				if( p == 3 && q == 3 && r == 3 )
					SpecialCase333();
			}
			else if( g == Geometry.Euclidean )
			{
				thickness = thickness / 2;
				threshold = 1/*20*/;

				//SimplexCalcs.CalcEScale();
				simplex = SimplexCalcs.MirrorsEuclidean();
				Vector3D[] verts = SimplexCalcs.VertsEuclidean();
				vertex = verts[2];
			}
			else
			{
				thickness = thicknessHyperbolic;
				threshold = 0.01;

				simplex = SimplexCalcs.Mirrors( p, q, r );
				Vector3D[] verts = SimplexCalcs.VertsBall( p, q, r );
				vertex = verts[2];

				//Vector3D[] simplexVerts = SimplexCalcs.VertsBall( p, q, r );
				//H3.Cell.Edge edge = new H3.Cell.Edge( simplexVerts[2], simplexVerts[3] );
				//H3.Cell.Edge edge = SimplexCalcs.HoneycombEdgeBall( p, q, r );
				//H3.Cell.Edge[] startingEdges = new H3.Cell.Edge[] { edge };

				//H3.Cell.Edge[] edges = Recurse.CalcEdgesSmart2( simplex, startingEdges );

				// Vertex Centered.
				bool vertexCentered = false;
				if( vertexCentered )
				{
					Vector3D v = SimplexCalcs.VertexPointBall( p, q, r );
					v = H3Models.BallToUHS( v );
					double scale = 1.0 / v.Abs();
					edges = edges.Select( e =>
					{
						Vector3D start = H3Models.UHSToBall( H3Models.BallToUHS( e.Start ) * scale );
						Vector3D end = H3Models.UHSToBall( H3Models.BallToUHS( e.End ) * scale );
						return new H3.Cell.Edge( start, end );
					} ).ToArray();
				}

				// Code to show endpoints of 535
				/*using( StreamWriter sw = File.CreateText( "535_points.pov" ) )
				{
					HashSet<Vector3D> verts = new HashSet<Vector3D>();
					foreach( H3.Cell.Edge e in edges )
					{
						verts.Add( Sterographic.SphereToPlane( e.Start ) );
						verts.Add( Sterographic.SphereToPlane( e.End ) );
					}

					foreach( Vector3D vert in verts )
						if( !Infinity.IsInfinite( vert ) )
							sw.WriteLine( PovRay.Sphere( new Sphere() { Center = vert, Radius = 0.01 } ) );
				}*/
			}

			// Recurse
			bool dual = false;
			{
				H3.Cell.Edge[] startingEdges = null;
				if( dual )
					startingEdges = new H3.Cell.Edge[] { SimplexCalcs.DualEdgeBall( simplex ) };
				else
				{
					//startingEdges = new H3.Cell.Edge[] { SimplexCalcs.HoneycombEdgeBall( simplex, vertex ) };
					Vector3D[] verts = SimplexCalcs.VertsEuclidean();
					Vector3D v1 = verts[0] + 2*verts[2]; // adjacent cube center
					Vector3D corner = verts[3];

					startingEdges = new H3.Cell.Edge[] { new H3.Cell.Edge( v1, corner ) };
				}

				edges = Recurse.CalcEdges( simplex, startingEdges, new Recurse.Settings() { G = g, Threshold = threshold } );

				edges = edges.Where( e =>
				{
					int sum = e.Depths.Count( d => d == 0 );
					return true;
				} ).ToArray();

				//CullHalfOfEdges( ref edges );

				// No need to cull edges in spherical case.
				// This was just to generate some images for 350-cell paper.
				//edges = Cull120Cell( edges );

				Simplex tet = new Simplex();
				tet.Facets = simplex;

				if( dual )
				{
					H3.Cell.Edge[] oneDualCell = edges.Where( e => e.Depths[2] == 0 ).ToArray();
					simplex = simplex.Skip( 1 ).ToArray();
					edges = Recurse.CalcEdges( simplex, oneDualCell, new Recurse.Settings() { G = g, Threshold = threshold } );

					int[] polyMirrors = new int[] { 0, 1, 3 };
					H3.Cell startingCell = HoneycombGen.PolyhedronToHighlight( g, polyMirrors, tet, new Vector3D() );
					cellsToHighlight = Recurse.CalcCells( simplex, new H3.Cell[] { startingCell } );
					//cellsToHighlight = new H3.Cell[] { startingCell };
					//cellsToHighlight = cellsToHighlight.Skip( 7 ).ToArray();
				}
				else
				{
					int[] polyMirrors = new int[] { 1, 2, 3 };
					H3.Cell startingCell = HoneycombGen.PolyhedronToHighlight( g, polyMirrors, tet, vertex );
					//cellsToHighlight = Recurse.CalcCells( simplex, new H3.Cell[] { startingCell } );
					cellsToHighlight = new H3.Cell[] { startingCell };
				}

				// Include just one cell?
				bool includeOne = false;
				if( includeOne )
				{
					edges = edges.Where( e => e.Depths[0] == 0 ).ToArray();
					//cellsToHighlight = cellsToHighlight.Where( c => c.Depths[0] == 0 ).ToArray();
				}
			}

			// Rotate
			bool rotate = false;
			if( rotate )
			{
				CompoundOfFive24Cells( ref edges );
			}

			// Write the file
			bool pov = true;
			if( pov )
			{
				string filename = string.Format( "{0}{1}{2}.pov", p, q, r );
				PovRay.WriteEdges( new PovRay.Parameters() { AngularThickness = thickness }, g, edges,
					filename, append: false );
				//File.Delete( filename );
				//PovRay.AppendFacets( cellsToHighlight, filename );

				HashSet<Vector3D> verts = new HashSet<Vector3D>();
				foreach( H3.Cell.Edge e in edges )
				{
					verts.Add( e.Start );
					verts.Add( e.End );
				}
				/*foreach( Vector3D v in verts )
				{
					Vector3D t = v;
					t.Normalize();
					t *= 0.9;
					System.Diagnostics.Trace.WriteLine( string.Format( "light_source {{ <{0},{1},{2}> White*.2 }}", t.X, t.Y, t.Z ) );
				}*/


				/*
				// Include the standard pov stuff, so we can batch this.
				string fileName = imageData.FormatFilename( string.Empty );
				using( StreamWriter sw = File.CreateText( fileName + ".pov" ) )
				{
					sw.WriteLine( "#include \"C:\\Users\\hrn\\Documents\\roice\\povray\\paper\\H3.pov\"" );
				}

				bool dummy = true;	// Doesn't matter for Pov-Ray, just Shapeways meshes.
				H3.SaveToFile( fileName, edges, dummy, append: true );
				*/
			}
			else
			{
				if( g == Geometry.Spherical )
				{
					edges = edges.Where( e => e.Start.Valid() && e.End.Valid() && !Infinity.IsInfinite( e.Start ) && !Infinity.IsInfinite( e.End ) ).ToArray();
					S3.EdgesToStl( edges );
				}
				else
					throw new System.NotImplementedException();
			}
		}

		private static void CompoundOfFive24Cells( ref H3.Cell.Edge[] edges )
		{
			List<H3.Cell.Edge> allEdges = new List<H3.Cell.Edge>();

			Vector3D v600 = Sterographic.R3toS3( SimplexCalcs.VertexSpherical( 3, 3, 5 ) );
			Vector3D v24 = Sterographic.R3toS3( SimplexCalcs.VertexSpherical( 3, 4, 3 ) );
			Sphere[] mirrors600 = SimplexCalcs.MirrorsSpherical( 3, 3, 5 );

			double a24 = v24.AngleTo( Sterographic.R3toS3( new Vector3D() ) );
			double a600 = v600.AngleTo( Sterographic.R3toS3( new Vector3D() ) );

			Matrix4D m600 = Matrix4D.MatrixToRotateinCoordinatePlane( a600, 2, 3 );
			Matrix4D m600_ = Matrix4D.MatrixToRotateinCoordinatePlane( -a600, 2, 3 );
			Matrix4D m24 = Matrix4D.MatrixToRotateinCoordinatePlane( a24, 2, 3 );
			Matrix4D m24_ = Matrix4D.MatrixToRotateinCoordinatePlane( -a24, 2, 3 );

			double eLength = 2 * Math.PI / 10;	// 600-cell edge length
			double a_id = Math.Asin( Math.Sin( eLength / 2 ) / Math.Sin( Math.PI / 3 ) * Math.Sin( Math.PI / 2 ) );

			eLength = 1.0 / Math.Sin( 2 * Math.PI / 5 );	// icosahedron edge length
			double a_i = Math.Asin( Math.Sin( eLength / 2 ) / Math.Sin( Math.PI / 3 ) * Math.Sin( Math.PI / 5 ) );

			Func<Vector3D, Vector3D> rot600 = v =>
			{
				v = Sterographic.R3toS3( v );
				v = m600.RotateVector( v );
				v = Sterographic.S3toR3( v );
				return v;
			};

			Func<Vector3D, int, Vector3D> rotOne = ( v, idx ) =>
			{
				v = Sterographic.R3toS3( v );
				v = m24.RotateVector( v );
				v = Sterographic.S3toR3( v );
				v.RotateAboutAxis( new Vector3D( 1, 0, 0 ), -a_id );

				// Vertex to cell center.
				v = Sterographic.R3toS3( v );
				v = m600_.RotateVector( v );
				v = Sterographic.S3toR3( v );

				List<int> reflections = new List<int>();
				if( idx == 0 )
				{
					reflections.Add( 2 );
					reflections.Add( 1 );
					reflections.Add( 2 );
					reflections.Add( 0 );
					reflections.Add( 1 );
					reflections.Add( 2 );
					reflections.Add( 1 );
					reflections.Add( 2 );
				}

				if( idx != 0 )
				{
					reflections.Add( 3 );
				}
				if( idx == 2 )
				{
					reflections.Add( 1 );
					reflections.Add( 2 );
				}
				if( idx == 3 )
				{
					reflections.Add( 2 );
					reflections.Add( 1 );
				}
				if( idx == 4 )
				{
					reflections.Add( 1 );
					reflections.Add( 0 );
					reflections.Add( 1 );
					reflections.Add( 2 );
					reflections.Add( 0 );
					reflections.Add( 1 );
					reflections.Add( 0 );
					reflections.Add( 1 );
				}

				foreach( int reflection in reflections )
					v = mirrors600[reflection].ReflectPoint( v );

				v = Sterographic.R3toS3( v );
				v = m600.RotateVector( v );
				v = Sterographic.S3toR3( v );

				//v.RotateAboutAxis( new Vector3D( 0, 0, 1 ), Math.PI/3 );
				//v.RotateAboutAxis( new Vector3D( 1, 0, 0 ), -a_i*2 );
				v = Sterographic.R3toS3( v );
				//v = m24_.RotateVector( v );
				v = Sterographic.S3toR3( v );
				

				return v;
			};

			for( int i = 0; i < 5; i++ )
			{
				//if( i == 0 )
				//	continue;

				allEdges.AddRange( edges.Select( e =>
				{
					H3.Cell.Edge newEdge = new H3.Cell.Edge( rotOne( e.Start, i ), rotOne( e.End, i ) );
					//H3.Cell.Edge newEdge = new H3.Cell.Edge( rot600( e.Start ), rot600( e.End ) );
					switch( i )
					{
					case 0:
						newEdge.Color = new Vector3D( 1, 0, 0 , 1 );
						//newEdge.Color = new Vector3D( 1, 1, 1, 0 );
						break;
					case 1:
						newEdge.Color = new Vector3D( 0, 1, 0 , 2 );
						break;
					case 2:
						newEdge.Color = new Vector3D( 0, 0, 1, 3 );
						break;
					case 3:
						newEdge.Color = new Vector3D( 1, 0, 1, 4 );
						break;
					case 4:
						newEdge.Color = new Vector3D( 0, 1, 1, 5 );
						break;
					}
					return newEdge;
				} ) );
			}

			edges = allEdges.ToArray();

			//edges = edges.Where( e => Tolerance.Equal( 1, e.Start.Abs() ) || Tolerance.Equal( 1, e.End.Abs() ) ).ToArray();

			HashSet<Vector3D> uniqueVerts = new HashSet<Vector3D>();
			foreach( H3.Cell.Edge e in edges )
			{
				uniqueVerts.Add( e.Start );
				uniqueVerts.Add( e.End );
			}
			System.Diagnostics.Trace.WriteLine( "Number of verts = " + uniqueVerts.Count );

			/*edges = edges.Where( e =>
			{
				Vector3D v = Tolerance.Equal( 1, e.Start.Abs() ) ? e.End : e.Start;
				if( v.Abs() >= 0.8 || v.Abs() <= 0.7 )
					return false;

				if( Tolerance.LessThanOrEqual( v.X, 0 ) || Tolerance.GreaterThanOrEqual( v.Y, 0 ) || Tolerance.GreaterThanOrEqual( v.Z, 0 ) )
					return false;

				return true;
			} ).ToArray();

			edges = edges.OrderBy( e => Tolerance.Equal( 1, e.Start.Abs() ) ? e.End.Abs() : e.Start.Abs() ).ToArray();*/
		}

		private static void SpecialCase333()
		{
			/*
				HashSet<Vector3D> verts = new HashSet<Vector3D>();
				foreach( H3.Cell.Edge e in edges )
				{
					verts.Add( e.Start );
					verts.Add( e.End );
				}

				/// We need to be smart such that the radius does not change too much at each step,
				/// and are going to do that by adding multiple edges.
				double upperTest = 35.25; // dist where rad ~300

				System.Func<double, double> locationAtRad = rad =>
				{
					double min = 0;
					double max = upperTest;
					double current = upperTest / 2;
					double searchOffset = ( max - min ) / 4;

					Vector3D cen;
					double radTest;
					H3Models.Ball.DupinCyclideSphere( verts.First() * current, .07 / 2, g, out cen, out radTest );

					// iterate to it.
					double diff = Math.Abs( rad - radTest );
					int iterations = 1000;
					for( int i = 0; i < iterations; i++ )
					{
						double t1, t2;
						H3Models.Ball.DupinCyclideSphere( verts.First() * ( current + searchOffset ), .07 / 2, g, out cen, out t1 );
						H3Models.Ball.DupinCyclideSphere( verts.First() * ( current - searchOffset ), .07 / 2, g, out cen, out t2 );
						double d1 = Math.Abs( rad - t1 );
						double d2 = Math.Abs( rad - t2 );
						if( d1 == Math.Min( Math.Min( diff, d1 ), d2 ) )
						{
							diff = d1;
							current += searchOffset;
						}
						if( d2 == Math.Min( Math.Min( diff, d1 ), d2 ) )
						{
							diff = d2;
							current -= searchOffset;
						}

						if( Tolerance.Equal( diff, 0.0, .0001 ) )
							return current;

						searchOffset /= 2;
					}

					throw new System.Exception();
				};

				List<H3.Cell.Edge> toAdd = new List<H3.Cell.Edge>();
				foreach( Vector3D v in verts )
				{
					//toAdd.Add( new H3.Cell.Edge( v, v * 49 ) );
					//toAdd.Add( new H3.Cell.Edge( v * 45, v * 48.5 ) );	// need a lot of resolution here.
					toAdd.Add( new H3.Cell.Edge( v, v * 35.25 ) );

					double maxRad = 300;
					double changePerEdge = 5;
					double last = 1.0;
					for( double rad = changePerEdge; rad < maxRad; rad += changePerEdge )
					{
						double loc = locationAtRad( rad );
						//toAdd.Add( new H3.Cell.Edge( v * last, v * loc ) );
						last = loc;
					}
				}
				toAdd.AddRange( edges );
				edges = toAdd.ToArray();
			*/

			// Sphere sweeps just not working - We need to add in clipped tori here.
			//Vector3D s3 = Sterographic.R3toS3( new Vector3D( 0.7 / 2, 0, 0 ) );
			Vector3D s3 = Sterographic.R3toS3( new Vector3D() );
			s3 += new Vector3D();

			// r2 is major radius of torus.
			// r2 - r1 is minor radius of torus
			double r1 = .07 / 2;
			double r2;
			Vector3D cen;
			double arbitrary = 0.25;
			H3Models.Ball.DupinCyclideSphere( new Vector3D( arbitrary, 0, 0 ), r1, Geometry.Spherical, out cen, out r2 );
			double y = cen.Abs();
			double rad = (r2 * r2 - r1 * r1 - y * y) / (2 * r1 - 2 * r2);
			double x = rad + r1;
			System.Diagnostics.Trace.WriteLine( x );
		}
	}
}
