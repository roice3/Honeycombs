namespace R3.Geometry
{
	using System;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;

	public enum PQR
	{
		P,
		Q,
		R
	}

	public struct HoneycombDef
	{
		public HoneycombDef( int p, int q, int r ) 
		{ 
			P = p; 
			Q = q; 
			R = r;
			Projection = Polytope.Projection.VertexCentered;
		}

		public int P;
		public int Q;
		public int R;
		public Polytope.Projection Projection;

		public string FormatDir( PQR constant, int constantVal )
		{
			switch( constant )
			{
				case PQR.P:
					return string.Format( "{0}qr", NorI( constantVal ) );
				case PQR.Q:
					return string.Format( "p{0}r", NorI( constantVal ) );
				case PQR.R:
					return string.Format( "pq{0}", NorI( constantVal ) );
			}

			throw new System.ArgumentException();
		}

		public string FormatFilename()
		{
			return FormatFilename( "png" );
		}

		public string FormatFilename( string extension )
		{
			string projection = string.Empty;
			switch( Projection )
			{
				case Polytope.Projection.FaceCentered:
					projection = "f"; break;
				case Polytope.Projection.EdgeCentered:
					projection = "e"; break;
				case Polytope.Projection.VertexCentered:
					projection = "v"; break;
				case Polytope.Projection.CellCentered:
					projection = "c"; break;
				default:
					throw new System.ArgumentException();
			}

			string fileName = string.Format( "{0}-{1}-{2}",
				NorI( P ), NorI( Q ), NorI( R ) );
			//string fileName = string.Format( "{{{0},{1},{2}}}",
			//	NorI( P ), NorI( Q ), NorI( R ) );
			if( !string.IsNullOrEmpty( extension ) )
				fileName += "." + extension;
			return fileName;

			//return string.Format( "{0}{1}{2}_{3:F2}.png",
			//	NorI( P ), NorI( Q ), NorI( R ), scaling );
		}

		private static string NorI( int n )
		{
			return n == -1 ? "i" : n.ToString();
		}
	}

	public class HoneycombGen
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

		private static void CullHalfOfEdges( ref H3.Cell.Edge[] edges )
		{
			double thresh = -.01;
			Vector3D looking = new Vector3D( 0, 0, -1 );
			edges = edges.Where( e => e.Start.Dot( looking ) > thresh || e.End.Dot( looking ) > thresh ).ToArray();
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
		/// </summary>
		public static void OneHoneycombNew( HoneycombDef imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			double thickness = 0.1;
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
				threshold = 5/*20*/;

				SimplexCalcs.CalcEScale();
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
					startingEdges = new H3.Cell.Edge[] { SimplexCalcs.HoneycombEdgeBall( simplex, vertex ) };

				edges = Recurse.CalcEdges( simplex, startingEdges, new Recurse.Settings() { G = g, Threshold = threshold } );

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
					H3.Cell startingCell = PolyhedronToHighlight( g, polyMirrors, tet, new Vector3D() );
					cellsToHighlight = Recurse.CalcCells( simplex, new H3.Cell[] { startingCell } );
					//cellsToHighlight = new H3.Cell[] { startingCell };
					//cellsToHighlight = cellsToHighlight.Skip( 7 ).ToArray();
				}
				else
				{
					int[] polyMirrors = new int[] { 1, 2, 3 };
					H3.Cell startingCell = PolyhedronToHighlight( g, polyMirrors, tet, vertex );
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

			// Write the file
			bool pov = false;
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
				foreach( Vector3D v in verts )
				{
					Vector3D t = v;
					t.Normalize();
					t *= 0.9;
					System.Diagnostics.Trace.WriteLine( string.Format( "light_source {{ <{0},{1},{2}> White*.2 }}", t.X, t.Y, t.Z ) );
				}


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
			double rad = ( r2 * r2 - r1 * r1 - y * y ) / ( 2 * r1 - 2 * r2 );
			double x = rad + r1;
			System.Diagnostics.Trace.WriteLine( x );		
		}

		public static void GoursatSet()
		{
			int baseHue = 0;
			string baseName = "535";
			List<int[]> toRun = new List<int[]>();
			//toRun.Add( new int[] { 1, 2 } );
			//toRun.Add( new int[] { 2, 3 } );	// 1100 doesn't work
			//toRun.Add( new int[] { 0, 2 } );
			//toRun.Add( new int[] { 0, 3 } );
			//toRun.Add( new int[] { 0, 1, 2 } );
			//toRun.Add( new int[] { 0, 1, 3 } );
			//toRun.Add( new int[] { 0, 1, 2, 3 } );

			baseName = "533";
			toRun.Add( new int[] { 0 } );

			/*
			string baseName = "4333";
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 1, 2 } );
			toRun.Add( new int[] { 0, 1 } );
			toRun.Add( new int[] { 2, 3 } );
			toRun.Add( new int[] { 0, 2 } );
			toRun.Add( new int[] { 0, 1, 2 } );
			toRun.Add( new int[] { 1, 2, 3 } );
			toRun.Add( new int[] { 0, 1, 2, 3 } );
			*/

			/*
			string baseName = "5333";
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 1, 2 } );
			toRun.Add( new int[] { 0, 1 } );
			toRun.Add( new int[] { 2, 3 } );
			toRun.Add( new int[] { 0, 2 } );
			toRun.Add( new int[] { 0, 1, 2 } );
			toRun.Add( new int[] { 1, 2, 3 } );
			toRun.Add( new int[] { 0, 1, 2, 3 } );
			*/

			/*
			int baseHue = 30;
			string baseName = "4343";
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 0, 1, 2 } );
			toRun.Add( new int[] { 1, 2 } );
			toRun.Add( new int[] { 0, 1 } );
			toRun.Add( new int[] { 0, 2 } );
			toRun.Add( new int[] { 0, 1, 2, 3 } );
			*/

			/*
			int baseHue = 135;
			string baseName = "4353";
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 1, 2 } );
			toRun.Add( new int[] { 1, 3 } );
			toRun.Add( new int[] { 1, 2, 3 } );
			toRun.Add( new int[] { 0, 1, 2 } );
			toRun.Add( new int[] { 2, 3 } );
			toRun.Add( new int[] { 0, 1 } );
			toRun.Add( new int[] { 0, 1, 2, 3 } );
			*/

			/*
			int baseHue = 210;
			string baseName = "5353";
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 1, 2 } );
			toRun.Add( new int[] { 0, 1 } );
			toRun.Add( new int[] { 0, 1, 3 } );	// replacement for 0,1,2, which didn't work well
			toRun.Add( new int[] { 0, 2 } );
			toRun.Add( new int[] { 0, 1, 2, 3 } );
			*/

			/*
			int baseHue = 55;
			string baseName = "53^11";
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 2 } );
			toRun.Add( new int[] { 1, 3 } );	// was 1,2, but required alt simplex
			toRun.Add( new int[] { 0, 2 } );
			toRun.Add( new int[] { 0, 1, 2 } );
			*/
			 
			foreach( int[] active in toRun )
				OneHoneycombGoursat( active, baseName, baseHue );
		}

		public static void ParacompactSet()
		{
			List<int[]> toRun = new List<int[]>();
			toRun.Add( new int[] { 0, 1 } );
			toRun.Add( new int[] { 0, 2 } );
			toRun.Add( new int[] { 0, 3 } );
			toRun.Add( new int[] { 1, 2 } );
			toRun.Add( new int[] { 1, 3 } );
			toRun.Add( new int[] { 2, 3 } );
			toRun.Add( new int[] { 0, 1, 2 } );
			toRun.Add( new int[] { 0, 1, 3 } );
			toRun.Add( new int[] { 0, 2, 3 } );
			toRun.Add( new int[] { 1, 2, 3 } );
			toRun.Add( new int[] { 0, 1, 2, 3 } );

			//toRun.Clear();
			//toRun.Add( new int[] { 0, 1, 3 } );

			int baseHue = 0;
			HoneycombDef def;

			baseHue = 135;
			def = new HoneycombDef( 6, 3, 3 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );

			baseHue = 220;
			def = new HoneycombDef( 6, 3, 4 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );

			baseHue = 180;
			def = new HoneycombDef( 6, 3, 5 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );

			baseHue = 255;
			def = new HoneycombDef( 6, 3, 6 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );

			baseHue = 105;
			def = new HoneycombDef( 3, 6, 3 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );

			baseHue = 300;
			def = new HoneycombDef( 4, 4, 3 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );

			baseHue = 0;
			def = new HoneycombDef( 4, 4, 4 );
			foreach( int[] active in toRun )
				Paracompact( def, active, baseHue );
		}

		private static void CalcThickness( int[] active )
		{
			double thickness = 0.04;
			switch( active.Length )
			{
			case 2:
				thickness = 0.025;
				break;
			case 3:
			case 4:
				thickness = 0.02;
				break;
			}
			H3.m_settings.AngularThickness = thickness;
		}

		private static string ActiveMirrorsString( int[] active )
		{
			Func<int, string> activeToString = i => active.Contains( i ) ? "1" : "0";
			string mirrorsString = string.Format( "{0}{1}{2}{3}",
				activeToString( 0 ), activeToString( 1 ), activeToString( 2 ), activeToString( 3 ) );
			return mirrorsString;
		}

		private static void OneHoneycombGoursat( int[] active, string baseName, int baseHue )
		{
			CalcThickness( active );

			// Create the simplex.
			Simplex simplex = new Simplex();
			simplex.InitializeGoursat(); 

			// Map of labels for mirrors consistent with input scheme to Goursat function.
			// Map is from wikipedia labeling scheme to the indices our function generates.
			//
			// wiki == our index
			// 0100 == 0
			// 0001 == 1
			// 1000 == 2
			// 0010 == 3
			Func<int, int> mapMirror = i =>
			{
				switch( i )
				{
					case 0: return 2;
					case 1: return 0;
					case 2: return 3;
					case 3: return 1;
				}
				throw new System.ArgumentException();
			};

			// We need to set this up before converting the mirrors.
			string mirrorsString = ActiveMirrorsString( active );
			string suffix = "-" + mirrorsString;

			// Convert our active mirrors into the Goursat tet indices.
			int[] polyMirrors = new int[] { 1, 2, 3 };
			active = active.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();
			polyMirrors = polyMirrors.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();

			Vector3D startingPoint = IterateToStartingPoint( active, simplex );
			List<H3.Cell.Edge> startingEdges = new List<H3.Cell.Edge>();
			foreach( int a in active )
			{
				Vector3D reflected = simplex.ReflectInFacet( startingPoint, a );
				startingEdges.Add( new H3.Cell.Edge( startingPoint, reflected ) );
			}

			bool doEdges = true;
			bool doCells = false;

			// Generate the honeycomb.
			H3.Cell.Edge[] edges = null;
			if( doEdges )
				edges = Recurse.CalcEdgesSmart( simplex.Facets, startingEdges.ToArray() );

			// Highlighted cells.
			H3.Cell[] cellsToHighlight = null;
			if( doCells )
			{
				H3.Cell startingCell = PolyhedronToHighlight( Geometry.Hyperbolic, polyMirrors, simplex, startingPoint );
				cellsToHighlight = Recurse.CalcCells( simplex.Facets, new H3.Cell[] { startingCell } );
				//cellsToHighlight = new H3.Cell[] { startingCell };
			}

			// plugin Wendy's nonuniform calcs here...
			//Nonuniform.Wendy( simplex, edges );

			// Trim out half the edges (the ones we won't see in our Pov-Ray view).
			Vector3D lookFrom = new Vector3D( 1, 1, 1 ) * 0.7;
			Vector3D lookAt = new Vector3D( );	// pov-ray lookat
			double thresh = -.01;
			if( doEdges )
				edges = edges.Where( e => e.Start.Dot( lookAt ) > thresh || e.End.Dot( lookAt ) > thresh ).ToArray();
			//if( doCells )
			//	cellsToHighlight = cellsToHighlight.Where( c => c.Center.Dot( lookAt ) > thresh ).ToArray();	// I don't think this works right

			// Setup Pov-ray stuff.
			// We have 16 possible mirror states.  We'll calculate the hue by converting the binary state to decimal, and doubling.
			// So for a given family, the hue will range over 32 numbers.
			int hue = baseHue + 2 * Convert.ToInt32( mirrorsString, 2 );
			string fileName = baseName + suffix;
			using( StreamWriter sw = File.CreateText( fileName + ".pov" ) )
			{
				sw.WriteLine( string.Format( "#declare lookFrom = <{0},{1},{2}>;", lookFrom.X, lookFrom.Y, lookFrom.Z ) );
				sw.WriteLine( string.Format( "#declare lookAt = <{0},{1},{2}>;", lookAt.X, lookAt.Y, lookAt.Z ) );
				sw.WriteLine( "#include \"C:\\Users\\hrn\\Documents\\roice\\povray\\H3_uniform_faces\\H3_uniform_faces.pov\"" );
				//sw.WriteLine( string.Format( "background {{ CHSL2RGB( <{0}, 1, .3> ) }}", hue ) );
				//sw.WriteLine( string.Format( "background {{ rgb <.13,.37,.31> }}" ) ); 
				sw.WriteLine( string.Format( "background {{ rgb 1 }}" ) ); 
			}

			if( doEdges )
				H3.SaveToFile( fileName, edges, finite: true, append: true );
			if( doCells )
			{
				HashSet<H3.Cell.Edge> cellEdges = new HashSet<H3.Cell.Edge>( new H3.Cell.EdgeEqualityComparer() );
				foreach( H3.Cell cell in cellsToHighlight )
					cell.AppendAllEdges( cellEdges );
				edges = cellEdges.ToArray();
				H3.SaveToFile( fileName, edges, finite: true, append: true );
				
				H3.AppendFacets( fileName, cellsToHighlight );
			}
		}

		private static void SetupBaseHue( string fileName, string mirrorsString, int baseHue )
		{
			// Setup Pov-ray stuff.
			// We have 16 possible mirror states.  We'll calculate the hue by converting the binary state to decimal, and doubling.
			// So for a given family, the hue will range over 32 numbers.
			int hue = baseHue + 2 * Convert.ToInt32( mirrorsString, 2 );
			using( StreamWriter sw = File.CreateText( fileName + ".pov" ) )
			{
				sw.WriteLine( "#include \"C:\\Users\\hrn\\Documents\\roice\\povray\\H3_paracompact\\H3_paracompact.pov\"" );
				sw.WriteLine( string.Format( "background {{ CHSL2RGB( <{0}, 1, .1> ) }}", hue ) );
			}
		}

		private static string BaseName( HoneycombDef def )
		{
			return string.Format( "{0}{1}{2}", def.P, def.Q, def.R );
		}

		public static void Paracompact( HoneycombDef def, int[] active, int baseHue )
		{
			string baseName = BaseName( def );
			string mirrorsString = ActiveMirrorsString( active );
			string suffix = "-" + mirrorsString;
			string fileName = baseName + suffix;

			if( File.Exists( fileName + ".pov" ) )
			{
				Console.WriteLine( string.Format( "Skipping {0}", fileName ) );
				return;
			}

			Console.WriteLine( string.Format( "Building {0}", fileName ) );
			CalcThickness( active );

			// The wiki mirrors are labeled in the reverse of ours.
			Func<int, int> mapMirror = i => 3 - i;
			active = active.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();

			Simplex simplex = new Simplex();
			simplex.Facets = SimplexCalcs.Mirrors( def.P, def.Q, def.R );
			simplex.Verts = SimplexCalcs.VertsBall( def.P, def.Q, def.R );

			Vector3D startingPoint = IterateToStartingPoint( active, simplex );
			if( startingPoint.DNE )
				return;
			List<H3.Cell.Edge> startingEdges = new List<H3.Cell.Edge>();
			foreach( int a in active )
			{
				Vector3D reflected = simplex.ReflectInFacet( startingPoint, a );
				startingEdges.Add( new H3.Cell.Edge( startingPoint, reflected ) );
			}

			SetupBaseHue( fileName, mirrorsString, baseHue );
			Recurse.m_background = new Vector3D( baseHue, 1, .1 );

			H3.Cell.Edge[] edges = Recurse.CalcEdgesSmart2( simplex.Facets, startingEdges.ToArray() );
			H3.SaveToFile( fileName, edges, finite: true, append: true );
		}

		// CHEAT! (would be better to do a geometrical construction)
		// We are going to iterate to the starting point that will make all edge lengths the same.
		private static Vector3D IterateToStartingPoint( int[] activeMirrors, Simplex simplex )
		{
			if( activeMirrors.Length == 1 )
				return simplex.Verts[activeMirrors[0]];

			// We are minimizing the output of this function, 
			// because we want all edge lengths to be as close as possible.
			// Input vector should be in the Ball Model.
			Func<Vector3D, double> diffFunc = v =>
			{
				List<double> lengths = new List<double>();
				for( int i = 0; i < activeMirrors.Length; i++ )
				{
					Vector3D reflected = simplex.ReflectInFacet( v, activeMirrors[i] );
					lengths.Add( H3Models.Ball.HDist( v, reflected ) );
				}

				double result = 0;
				double average = lengths.Average();
				foreach( double length in lengths )
					result += Math.Abs( length - average );
				return result;
			};

			// So that we can leverage Euclidean barycentric coordinates, we will first convert our simplex to the Klein model.
			// We will need to take care to properly convert back to the Ball as needed.
			Vector3D[] kleinVerts = simplex.Verts.Select( v => HyperbolicModels.PoincareToKlein( v ) ).ToArray();

			// Normalizing barycentric coords amounts to making sure the 4 coords add to 1.
			Func<Vector3D, Vector3D> baryNormalize = b =>
			{
				return b / ( b.X + b.Y + b.Z + b.W );
			};

			// Bary Coords to Euclidean
			Func<Vector3D[], Vector3D, Vector3D> baryToEuclidean = ( kv, b ) =>
			{
				Vector3D result =
					kv[0] * b.X + kv[1] * b.Y + kv[2] * b.Z + kv[3] * b.W;
				return result;
			};

			// Our starting barycentric coords (halfway between all active mirrors).
			Vector3D bary = new Vector3D();
			foreach( int a in activeMirrors )
				bary[a] = 0.5;
			bary = baryNormalize( bary );

			// For each iteration, we'll shrink this search offset.
			// NOTE: The starting offset and decrease factor I'm using don't guarantee convergence, 
			// but it seems to be working pretty well (even when varying these parameters).
			//double searchOffset = 1.0 - bary[activeMirrors[0]];
			//double searchOffset = bary[activeMirrors[0]];
			double factor = 1.5;	// Adjusting this helps get some to converge, e.g. 4353-1111 
			double searchOffset = bary[activeMirrors[0]] / factor;		

			double min = double.MaxValue;
			int iterations = 1000;
			for( int i = 0; i < iterations; i++ )
			{
				min = diffFunc( HyperbolicModels.KleinToPoincare( baryToEuclidean( kleinVerts, bary ) ) );
				foreach( int a in activeMirrors )
				{
					Vector3D baryTest1 = bary, baryTest2 = bary;
					baryTest1[a] += searchOffset;
					baryTest2[a] -= searchOffset;
					baryTest1 = baryNormalize( baryTest1 );
					baryTest2 = baryNormalize( baryTest2 );

					double t1 = diffFunc( HyperbolicModels.KleinToPoincare( baryToEuclidean( kleinVerts, baryTest1 ) ) );
					double t2 = diffFunc( HyperbolicModels.KleinToPoincare( baryToEuclidean( kleinVerts, baryTest2 ) ) );
					if( t1 < min )
					{
						min = t1;
						bary = baryTest1;
					}
					if( t2 < min )
					{
						min = t2;
						bary = baryTest2;
					}
				}

				if( Tolerance.Equal( min, 0.0, 1e-14 ) )
				{
					System.Console.WriteLine( string.Format( "Converged in {0} iterations.", i ) );
					break;
				}

				searchOffset /= factor;
			}

			if( !Tolerance.Equal( min, 0.0, 1e-14 ) )
			{
				System.Console.WriteLine( "Did not converge: " + min );

				// Be a little looser before thrown an exception.
				if( !Tolerance.Equal( min, 0.0, 1e-12 ) )
				{
					System.Console.ReadKey( true );
					//throw new System.Exception( "Boo. We did not converge." );
					return Vector3D.DneVector();
				}
			}

			return HyperbolicModels.KleinToPoincare( baryToEuclidean( kleinVerts, bary ) );
		}

		private static H3.Cell PolyhedronToHighlight( Geometry g, int[] mirrors, Simplex simplex, Vector3D startingPoint )
		{
			if( mirrors.Length != 3 )
				throw new System.Exception( "We need exactly three mirrors to generate a polyhedron." );

			// In the general case, we can have 3 types of polygons, each being generated by choosing 2 of the 3 mirrors.
			// When fewer polygon types exist, GenFacet will return null, so we need to check that.
			List<H3.Cell.Facet> polyFacets = new List<H3.Cell.Facet>();
			polyFacets.Add( GenFacet( g, mirrors[0], mirrors[1], simplex, startingPoint ) );
			polyFacets.Add( GenFacet( g, mirrors[1], mirrors[2], simplex, startingPoint ) );
			polyFacets.Add( GenFacet( g, mirrors[0], mirrors[2], simplex, startingPoint ) );
			polyFacets.RemoveAll( f => f == null );

			HashSet<Vector3D> completedFacetIds = new HashSet<Vector3D>();
			foreach( H3.Cell.Facet f in polyFacets )
				completedFacetIds.Add( f.ID );

			Recurse.GenPolyhedron( mirrors.Select( m => simplex.Facets[m] ).ToArray(),
				polyFacets.ToArray(), polyFacets, completedFacetIds );

			H3.Cell cell = new H3.Cell( polyFacets.ToArray() );

			// Calc a center.  XXX - totally wrong.
			Vector3D center = new Vector3D();
			foreach( Vector3D v in cell.Verts )
				center += v;
			center /= cell.Verts.Count();
			cell.Center = center;

			// Calc the spheres.
			foreach( H3.Cell.Facet f in cell.Facets )
				f.CalcSphereFromVerts( g );

			return cell;
		}

		private static H3.Cell.Facet GenFacet( Geometry g, int mirror1, int mirror2, Simplex simplex, Vector3D startingPoint )
		{
			List<Sphere> mirrors = new List<Sphere>();
			mirrors.Add( simplex.Facets[mirror1] );
			mirrors.Add( simplex.Facets[mirror2] );

			List<H3.Cell.Edge> startingEdges = new List<H3.Cell.Edge>();
			Vector3D reflected = simplex.ReflectInFacet( startingPoint, mirror1 );
			startingEdges.Add( new H3.Cell.Edge( startingPoint, reflected ) );
			reflected = simplex.ReflectInFacet( startingPoint, mirror2 );
			startingEdges.Add( new H3.Cell.Edge( startingPoint, reflected ) );
			startingEdges.RemoveAll( e => e.Start == e.End );
			if( startingEdges.Count == 0 )
				return null;

			H3.Cell.Edge[] completedEdges = Recurse.CalcEdges( mirrors.ToArray(), startingEdges.ToArray(), new Recurse.Settings() { G = g } );
			if( completedEdges.Length == 1 )
				return null;

			List<Vector3D> facetVerts = new List<Vector3D>();
			H3.Cell.Edge edge = completedEdges.First();
			Vector3D start = edge.Start;
			Vector3D current = edge.End;
			facetVerts.Add( edge.End );
			while( current != start )
			{
				edge = completedEdges.First( e => e != edge && ( e.Start == current || e.End == current ) );
				current = edge.Start == current ? edge.End : edge.Start;
				facetVerts.Add( current );
			}

			return new H3.Cell.Facet( facetVerts.ToArray() );
		}

		private static Sphere[] SphericalCellFacetMirrors( HoneycombDef imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			double inRadius = Honeycomb.InRadius( p, q, r );
			//inRadius *= 1.4;	// Experimenting with {3,3,u}

			Tiling tiling = new Tiling();
			TilingConfig config = new TilingConfig( p, q );
			tiling.GenerateInternal( config, imageData.Projection );

			Sphere[] mirrors = H3.GenFacetSpheres( tiling, inRadius )
				.Select( f => f.Sphere ).ToArray();

			return mirrors;
		}
	}
}
