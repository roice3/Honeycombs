namespace HyperbolicModels
{
	using System;
	using System.Collections.Generic;
	using System.Drawing;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Drawing;
	using R3.Geometry;
	using R3.Math;

	public static class HoneycombPaper
	{
		public static void DoStuff( Settings settings )
		{
			HoneycombDef imageData = new HoneycombDef( settings.P, settings.Q, settings.R );

			////////////////////////////////////////////////////////////// Various things we've run over time.
			//Sandbox.CalcSelfSimilarityScale();
			//Sandbox.Check_pq_Distances();
			//HyperidealSquares();
			//S3.Hypercube();
			//R3.Geometry.Euclidean.GenEuclidean();
			//HoneycombGen.OneHoneycombOldCode();
			//AnimateCell( imageData );
			//CreateCellPovRay( imageData, "cell.pov" );
			//CreateSimplex( imageData );
			//HoneycombGen_old.OneHoneycombNew( new HoneycombDef() { P = imageData.P, Q = imageData.Q, R = imageData.R } );
			//SphericalAnimate( imageData );
			OneImage( settings );

			HoneycombDef[] scaleLarger = GetImageSet().Where( h =>
				Geometry2D.GetGeometry( h.P, h.Q ) == Geometry.Euclidean ||
				Geometry2D.GetGeometry( h.P, h.Q ) == Geometry.Spherical ).ToArray();
			int count = scaleLarger.Length;
			//foreach( HoneycombAndView h in scaleLarger )
			//	Trace.WriteLine( h.FormatFilename() );

			//BatchRun( settings );
		}

		private static int ReadArg( string arg )
		{
			if( arg == "i" || arg == "inf" )
				return -1;

			return int.Parse( arg );
		}

		public static void OneImage( Settings config, double t = 0.0 )
		{
			HoneycombDef imageData = new HoneycombDef( config.P, config.Q, config.R );
			int p = imageData.P, q = imageData.Q, r = imageData.R;

			string filename = imageData.FormatFilename();
			//if( File.Exists( filename ) )
			//	return;

			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );
			double bounds = config.UhsBoundary.Bounds;

			// Calculate the color scale.
			int size = 200;
			CoxeterImages.Settings settings = new CoxeterImages.Settings()
			{
				Honeycomb = imageData,
				Width = size,
				Height = size,
				Bounds = bounds,
				Mirrors = mirrors,
				FileName = imageData.FormatFilename(),
			};

			CoxeterImages imageCalculator = new CoxeterImages();
			imageCalculator.AutoCalcScale( settings );
			if( settings.ColorScaling < 1 )
				settings.ColorScaling = 15;

			Program.Log( "\nGenerating full image..." );
			size = 500;
			settings.Width = config.UhsBoundary.ImageWidth;
			settings.Height = config.UhsBoundary.ImageHeight;
			settings.FileName = filename;
			imageCalculator.GenImage( settings, t );
		}

		internal static CoxeterImages.Settings AutoCalcScale( HoneycombDef def, Sphere[] mirrors )
		{
			// Calculate the color scale.
			int size = 200;
			CoxeterImages.Settings settings = new CoxeterImages.Settings()
			{
				Honeycomb = def,
				Width = size,
				Height = size,
				Bounds = 1.1,
				Mirrors = mirrors,
				FileName = def.FormatFilename(),
			};

			CoxeterImages imageCalculator = new CoxeterImages();
			imageCalculator.AutoCalcScale( settings );
			if( settings.ColorScaling < 1 )
				settings.ColorScaling = 15;

			return settings;
		}

		private static void BatchRun( Settings config )
		{
			bool batchRun = true;
			if( !batchRun )
				return;

			/*HoneycombDef[] fullSet = GetFullImageSet().ToArray();
			foreach( HoneycombDef iData in fullSet )
			{
				config.Angles = new int[] { iData.P, iData.Q, iData.R };
				OneImage( config );
			}*/


			//int[] ps = new int[] { 3, 4, 5, 6, 7, 8, -1 };
			//int[] qs = new int[] { 3 };
			//int[] rs = new int[] { 3, 4, 5, 6, 7, 8, -1 };

			int[] ps = new int[] { 3, 4, 5, 6, 7, 8, -1 };
			int[] qs = new int[] { 3, 4, 5, 6, 7, 8, -1 };
			int[] rs = new int[] { 3, 4, 5, 6, 7, 8, -1 };
			//int[] rs = new int[] { 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30 };

			//int[] rs = new int[] { 3 };
			foreach( int p in ps )
			foreach( int q in qs )
			foreach( int r in rs )
			{
				//if( !( Geometry2D.GetGeometry( p, q ) == Geometry.Hyperbolic || Geometry2D.GetGeometry( q, r ) == Geometry.Hyperbolic ) )
				//	continue;
				//if( !( Geometry2D.GetGeometry( p, q ) == Geometry.Spherical || Geometry2D.GetGeometry( q, r ) == Geometry.Spherical ) )
					//continue;

				if( !((p==7&&q==3&&r==3))) continue;
				//if( !(p==3) ) continue;
				//if( !(p ==4||p==5) ) continue;
				//if( q!=7 ) continue;
				//if( r!=5 ) continue;
				//if( r!=3 ) continue;
				//if( !(r==4 || r== 5)) continue;

				//if( Geometry2D.GetGeometry( q, r ) != Geometry.Spherical )
				//	continue;

				config.Angles = new int[] { p, q, r };
				//OneImage( config );
				string fileName = config.FileName( "pov" );
				System.Console.WriteLine( fileName );
				if( File.Exists( fileName ) )
				{
					System.Console.WriteLine( "Skipping because it exists." );
					//continue;
				}

				HoneycombDef def = new HoneycombDef( p, q, r );
				CreateCellPovRay( def, fileName );
				config.Angles = new int[] { p, q, r };
				//OneImage( config );
			}
		}

		internal static IEnumerable<HoneycombDef> GetImageSet()
		{
			return 
				GetImageSetForSlice( PQR.P, 3 ).Concat(
				GetImageSetForSlice( PQR.Q, 3 ) ).Concat( 
				GetImageSetForSlice( PQR.R, 3 ) ).Concat(
				GetImageSetForSlice( PQR.P, -1 ) ).Concat(
				GetImageSetForSlice( PQR.Q, -1 ) ).Concat(
				GetImageSetForSlice( PQR.R, -1 ) ).Concat(
				GetEuclidImageSet() );
		}

		private static IEnumerable<HoneycombDef> GetFullImageSet()
		{
			int[] vals = new int[] { 3, 4, 5, 6, 7, 10, 20, -1 };
			foreach( int p in vals )
			foreach( int q in vals )
			foreach( int r in vals )
			{
				yield return new HoneycombDef()
				{
					P = p,
					Q = q,
					R = r,
				};
			}
		}

		private static IEnumerable<HoneycombDef> GetImageSetForSlice( PQR constant, int constantValue )
		{
			//int i = 8;
			//int j = 8;
			for( int i=3; i<=8; i++ )
			for( int j=3; j<=8; j++ )
			//foreach( Polytope.Projection projection in System.Enum.GetValues( typeof( Polytope.Projection ) ) )
			{
				//if( Geometry2D.GetGeometry( p, q ) != Geometry.Spherical )
				//	continue;

				//if( projection == Polytope.Projection.CellCentered )
				//	continue;

				/*if( !( Geometry2D.GetGeometry( p, q ) == Geometry.Euclidean ||
					   Geometry2D.GetGeometry( q, r ) == Geometry.Euclidean ) )
					continue;*/

				int p = constant == PQR.P ? constantValue : i;
				int q = constant == PQR.Q ? constantValue : constant == PQR.P ? i : j;
				int r = constant == PQR.R ? constantValue : j;

				// Do the last as infinity
				System.Func<int, int> iSafe = input => input == 8 ? -1 : input;
				yield return new HoneycombDef()
				{
					P = iSafe( p ),
					Q = iSafe( q ),
					R = iSafe( r ),	
					//Projection = projection
				};
			}
		}

		private static IEnumerable<HoneycombDef> GetEuclidImageSet()
		{
			for( int p=3; p<=8; p++ )
			for( int q=3; q<=8; q++ )
			for( int r=3; r<=8; r++ )
			{
				if( !( Geometry2D.GetGeometry( p, q ) == Geometry.Euclidean ||
					   Geometry2D.GetGeometry( q, r ) == Geometry.Euclidean ) )
					continue;

				// Do the last as infinity
				System.Func<int, int> iSafe = input => input == 8 ? -1 : input;
				yield return new HoneycombDef()
				{
					P = iSafe( p ),
					Q = iSafe( q ),
					R = iSafe( r )
				};
			}
		}

		private static IEnumerable<HoneycombDef> GetCoxeterSet()
		{
			yield return new HoneycombDef() { P = 5, Q = 3, R = 4 };
			yield return new HoneycombDef() { P = 4, Q = 3, R = 5 };
			yield return new HoneycombDef() { P = 5, Q = 3, R = 5 };
			yield return new HoneycombDef() { P = 3, Q = 5, R = 3 };

			for( int p=3; p<=6; p++ )
			for( int q=3; q<=6; q++ )
			for( int r=3; r<=6; r++ )
			{
				if( Geometry2D.GetGeometry( p, q ) == Geometry.Spherical &&
					Geometry2D.GetGeometry( q, r ) == Geometry.Spherical )
					continue;

				if( Geometry2D.GetGeometry( p, q ) == Geometry.Hyperbolic ||
					Geometry2D.GetGeometry( q, r ) == Geometry.Hyperbolic )
					continue;

				yield return new HoneycombDef()
				{
					P = p,
					Q = q,
					R = r
				};
			}
		}

		private static void SphericalAnimate( HoneycombDef imageData )
		{
			double colorScaling = AnimColorScaling( imageData );

			//int fps = 30;
			//int frames = 60 * fps;
			int frames = 5;
			for( int i = 0; i < frames; i++ )
			{
				double t = (double)i/frames;
				string filename = "batch/" + imageData.FormatFilename( string.Empty ) + FrameString( i ) + ".png";
				OneAnimationFrame( imageData, filename, colorScaling, t );
			}
		}

		private static string FrameString( int i )
		{
			string num = i.ToString();
			num = num.PadLeft( 3, '0' );
			return "_" + num;
		}

		private static double AnimColorScaling( HoneycombDef imageData )
		{
			return 10;
		}

		private static void OneAnimationFrame( HoneycombDef imageData, string filename, double colorScaling, double t = 0.0 )
		{
			int p = imageData.P, q = imageData.Q, r = imageData.R;
			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );

			int size = 750;
			size = 1024;
			CoxeterImages.Settings settings = new CoxeterImages.Settings()
			{
				Honeycomb = imageData,
				Width = size*2,
				Height = size,
				Bounds = 1.0,
				Mirrors = mirrors,
				FileName = filename,
			};

			CoxeterImages imageCalculator = new CoxeterImages();
			settings.ColorScaling = colorScaling;
			imageCalculator.GenImage( settings, t );
		}

		/// <summary>
		/// This is tricky to get right and has a big effect on POV-Ray output.
		/// It must match the simplex mirror construction.
		/// It would be nicer if we had a general way to calculate the interior point.
		/// </summary>
		internal static Vector3D InteriorPointBall
		{
			get
			{
				Vector3D cen = new Vector3D( 0.05, 0.01, -0.1 );	// Works for a large class of honeycombs.
				cen = new Vector3D( 0.05, 0.01, 0.95 );	// p infinity
				cen = new Vector3D( 0.05, 0.01, 0.5 );	// p infinity - work's better and for many finite p
				//cen = new Vector3D( 0.05, 0.01, -.2 );
				//cen = H3Models.UHSToBall( new Vector3D( -0.95, -0.01, 5 ) );	// q infinity

				cen.RotateXY( Math.PI / 2 );	// only if we also rotate simplex mirrors.  XXX - make a setting.
				return cen;
			}
		}

		private static void AnimateCell( HoneycombDef imageData )
		{
			int numFrames = 100;
			for( int i = 0; i < numFrames; i++ )
			{
				double t = (double)i / numFrames;
				string filename = "batch/cell" + FrameString( i ) + ".pov";
				Console.WriteLine( filename );

				System.IO.File.Delete( filename );
				using( StreamWriter sw = new StreamWriter( filename ) )
					sw.WriteLine( "#include \"C:\\Users\\hrn\\Documents\\roice\\povray\\H3\\horosphere\\633.pov\"" );

				CreateCellPovRay( imageData, filename, t );
			}
		}

		private static void CreateCellPovRay( HoneycombDef def, string filename, double t = 0 )
		{
			int p = def.P;
			int q = def.Q;
			int r = def.R;

			//Vector3D trans = new Vector3D( 1.0/3, 0 ) * (2 + 2 * Math.Sin( Math.PI / 6 )) * t;
			//double scale = 1.8;
			Vector3D trans = new Vector3D();
			double scale = 1.0;
			Vector3D[] sVerts = null;   // SimplexCalcs.VertsBall( p, q, r );

			Vector3D vUHS = H3Models.BallToUHS( SimplexCalcs.VertexPointBall( p, q, r ) );

			// Just did this for everything.  Non-general position working better and will make all heads consistent.
			scale = 2.0;

			if( Geometry2D.GetGeometry( q, r ) != Geometry.Hyperbolic ) // Vertex-centered if possible
				scale = 1.0 / vUHS.Z;
			//else if( Geometry2D.GetGeometry( p, q ) == Geometry.Hyperbolic ) // Make the biggest head somewhat smaller.
			//	scale = 2.0;

			Vector3D cen = InteriorPointBall;
			/*var kleinVerts = sVerts.Select( v => HyperbolicModels.PoincareToKlein( v ) );
			Vector3D avg = new Vector3D();
			foreach( Vector3D v in kleinVerts )
				avg += v;
			avg /= kleinVerts.Count();
			Vector3D cen = HyperbolicModels.KleinToPoincare( avg );*/
			cen = H3Models.BallToUHS( cen );
			cen += trans;
			//cen *= scale;
			cen = H3Models.UHSToBall( cen );

			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: false );

			// Apply transformations.
			simplex = simplex.Select( s =>
			{
				Sphere.TranslateSphere( s, trans );
				Sphere.ScaleSphere( s, scale );
				return H3Models.UHSToBall( s );
			} ).ToArray();

			for( int i = 0; i < 4; i++ )
				if( simplex[i].IsPointInside( cen ) )
					simplex[i].Invert = true;

			Sphere[] simplexForColorScale = SimplexCalcs.Mirrors( p, q, r, moveToBall: true );
			CoxeterImages.Settings temp = AutoCalcScale( def, simplexForColorScale );
			int maxDepth = (int)temp.ColorScaling;
			//Random rand = new Random( p+q+r );
			//int randOffset = rand.Next( maxDepth );

			bool ball = true;
			bool dual = false;
			H3.Cell[] simplicesFinal = GenCell( simplex, null, cen, ball, dual );

			using( StreamWriter sw = File.CreateText( filename ) )	// We need to reuse this StreamWriter (vs. calling AppendSimplex) for performance.
			{
				sw.WriteLine( "#include \"hyper_ball.pov\"" );

				//int[] include = new int[] { 0, 1, 2, 3 };
				int[] include = new int[] { 0 };
				if( dual )
					include = new int[] { 3 };

				// Output the facets.
				foreach( H3.Cell cell in simplicesFinal )
				{
					Sphere[] facets = cell.Facets.Select( f => f.Sphere ).ToArray();
					if( m_toKlein )
						facets = facets.Select( s => H3Models.BallToKlein( s ) ).ToArray();

					int depth = cell.Depths[0] + 1;
					Color c = Coloring.ColorAlongHexagon( maxDepth,  depth );
					if( cell.Depths.Sum() % 2 == 0 )
						c = Coloring.Inverse( c );
					PovRay.AddSimplex( sw, facets, cell.Center, include, filename, Coloring.ToVec( c ) );
				}

				/*include = new int[] { 1, 2, 3 };
				foreach( H3.Cell cell in simplicesFinal )
				{
					Sphere[] facets = cell.Facets.Select( f => f.Sphere ).ToArray();
					Color c = Color.Red;
					Vector3D cv = Coloring.ToVec( c );
					cv.W = 0.9;
					PovRay.AddSimplex( sw, facets, cell.Center, include, filename, cv );
				}*/
			}

			// Output the edges/verts.
			bool includeEdges = false;
			if( includeEdges )
			{
				sVerts = sVerts.Select( v => 
				{
					v = H3Models.BallToUHS( v );
					v += trans;
					v *= scale;
					return H3Models.UHSToBall( v );
				} ).ToArray();

				H3.Cell.Edge[] edges = Recurse.CalcEdges( simplex.Skip( 1 ).ToArray(), 
					new H3.Cell.Edge[] { new H3.Cell.Edge( sVerts[2], sVerts[3], order: false ) },
					new Recurse.Settings() { Threshold = 0.01 } );
				PovRay.WriteH3Edges( new PovRay.Parameters { AngularThickness = 0.01 }, edges, filename, append: true );

				HashSet<Vector3D> verts = new HashSet<Vector3D>();
				foreach( H3.Cell.Edge e in edges )
					verts.Add( e.End );
				PovRay.WriteVerts( new PovRay.Parameters { AngularThickness = 0.02 }, Geometry.Hyperbolic, verts.ToArray(), filename, append: true );
			}
		}

		private static bool m_toKlein = false;

		/// <summary>
		/// Used to generate a regular cell as a set of simplices and save to a Pov-ray file.
		/// This will work for non-finite cells.
		/// </summary>
		internal static H3.Cell[] GenCell( Sphere[] simplex, Mesh mesh, Vector3D cen, bool ball, bool dual )
		{
			// We don't want to include the first mirror (which reflects across cells).
			Sphere[] mirrors = simplex.Skip( 1 ).ToArray();
			if( dual )
				mirrors = new Sphere[] { simplex[0], simplex[1], simplex[2] };
			//Sphere[] allMirrors = simplex.ToArray();
			mirrors = simplex.ToArray();

			// Simplices will be the "cells" in Recurse.CalcCells.
			H3.Cell.Facet[] simplexFacets = simplex.Select( m => new H3.Cell.Facet( m ) ).ToArray();

			H3.Cell startingCell = new H3.Cell( simplexFacets );
			startingCell.Center = cen;
			startingCell.Mesh = mesh;

			//FCOrient( startingCell );

			startingCell = startingCell.Clone();	// So our mirrors don't get munged after we reflect around later.
			//H3.Cell[] simplices = Recurse.CalcCells( mirrors, new H3.Cell[] { startingCell }, new Recurse.Settings() { Ball = ball } );
			H3.Cell[] simplices = Recurse.CalcCellsSmart( mirrors, new H3.Cell[] { startingCell }, new Recurse.Settings() { Ball = ball }, (int)(2.5e6) );
			//H3.Cell[] simplices = new H3.Cell[] { startingCell };

			// Layers.
			//int layer = 0;
			//return simplices.Where( s => s.Depths[0] <= layer /*&& s.Depths[0] == 3 && s.Depths[1] == 3*/ ).ToArray();
			return simplices.ToArray();
		}

		private static void CreateSimplex( HoneycombDef imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			Vector3D cen = InteriorPointBall;
			bool ball = true;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, ref cen, moveToBall: ball );

			// Offset as we do for the boundary images.
			//Sphere s = H3Models.UHSToBall( simplex[0] );
			//s = CoxeterImages.GeodesicOffset( s, 0.02, ball: true );

			if( m_toKlein )
				simplex = simplex.Select( s => H3Models.BallToKlein( s ) ).ToArray();

			int[] include = new int[] { 0, 1, 2, 3 };	// All facets
			//int[] include = new int[] { 1 };
			File.Delete( "simplex.pov" );
			PovRay.AppendSimplex( simplex, cen, include, "simplex.pov" );

			bool includeEdges = false;
			if( includeEdges )
			{
				H3.Cell.Edge[] edges = SimplexCalcs.SimplexEdgesUHS( p, q, r );
				PovRay.WriteEdges( new PovRay.Parameters { Halfspace = true, AngularThickness = 0.03 },
					Geometry.Hyperbolic, edges, "simplex.pov", append: true );
			}
		}

		/// <summary>
		/// Face centered orientation.
		/// </summary>
		internal static void FCOrient( H3.Cell cell )
		{
			// First, need to scale so the lowest triangles are the same size,
			// then reorient so that one triangle is oriented along z axis,
			// then scale so that the triangle is flat.

			// Calculate how much we need to offset to make the cell facet flat.
			Mobius m = FCOrientMobius( cell.Facets[0].Sphere );
			foreach( H3.Cell.Facet f in cell.Facets )
				H3Models.TransformInBall2( f.Sphere, m );
			cell.Center = H3Models.TransformHelper( cell.Center, m );
		}

		/// <summary>
		/// A Mobius that will put our simplex into a face centered orientation.
		/// Meant to be used with H3Models.TransformInBall2 or H3Models.TransformHelper
		/// </summary>
		internal static Mobius FCOrientMobius( Sphere cellSphereInBall )
		{
			Mobius m = new Mobius();
			double d = cellSphereInBall.Center.Abs() - cellSphereInBall.Radius;
			m.Isometry( Geometry.Hyperbolic, 0, new Vector3D( 0, d ) );
			return m;
		}

		private static void HyperidealSquares()
		{
			Mobius rot = new Mobius();
			rot.Isometry( Geometry.Spherical, Math.PI / 4, new Vector3D() );

			List<Segment> segs = new List<Segment>();
			int[] qs = new int[] { 5, -1 };
			foreach( int q in qs )
			{
				TilingConfig config = new TilingConfig( 4, q, 1 );
				Tile t = Tiling.CreateBaseTile( config );
				List<Segment> polySegs = t.Boundary.Segments;
				polySegs = polySegs.Select( s => { s.Transform( rot ); return s; } ).ToList();
				segs.AddRange( polySegs );
			}
			
			Vector3D v1 = new Vector3D(1,0);
			v1.RotateXY( Math.PI/6 );
			Vector3D v2 = v1;
			v2.Y *= -1;
			Vector3D cen;
			double rad;
			H3Models.Ball.OrthogonalCircle( v1, v2, out cen, out rad );
			Segment seg = Segment.Arc( v1, v2, cen, false );
			rot.Isometry( Geometry.Spherical, Math.PI / 2, new Vector3D() );
			for( int i = 0; i < 4; i++ )
			{
				seg.Transform( rot );
				segs.Add( seg.Clone() );
			}

			SVG.WriteSegments( "output1.svg", segs );

			System.Func<Segment, Segment> PoincareToKlein = s =>
			{
				return Segment.Line(
					HyperbolicModels.PoincareToKlein( s.P1 ),
					HyperbolicModels.PoincareToKlein( s.P2 ) );
			};
			segs = segs.Select( s => PoincareToKlein( s ) ).ToList();

			Vector3D v0 = new Vector3D( v1.X, v1.X );
			Vector3D v3 = v0;
			v3.Y *= -1;
			Segment seg1 = Segment.Line( v0, v1 ), seg2 = Segment.Line( v2, v3 );
			Segment seg3 = Segment.Line( new Vector3D( 1, 1 ), new Vector3D( 1, -1 ) );
			for( int i = 0; i < 4; i++ )
			{
				seg1.Transform( rot ); 
				seg2.Transform( rot );
				seg3.Transform( rot );
				segs.Add( seg1.Clone() );
				segs.Add( seg2.Clone() );
				segs.Add( seg3.Clone() );
			}

			SVG.WriteSegments( "output2.svg", segs );
		}
	}
}
