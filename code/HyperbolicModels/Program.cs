namespace HyperbolicModels
{
	using System;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;
	using R3.Math;

	class Program
	{
		static void Main( string[] args )
		{
			try
			{
				HoneycombPaper();

				//HoneycombGen.GoursatSet();

				//List<HoneycombAndView> honeycombs = GetCoxeterSet().ToList();
				//foreach( HoneycombAndView h in honeycombs )
					//HoneycombGen.OneHoneycombNew( h );
				//HoneycombGen.OneHoneycombNew( new HoneycombAndView() { P = 5, Q = 3, R = 3 } );
			}
			catch( System.Exception ex )
			{
				System.Diagnostics.Trace.WriteLine( ex.Message + "\n" + ex.StackTrace );
				System.Console.WriteLine( ex.Message + "\n" + ex.StackTrace );
			}
		}

		internal static IEnumerable<HoneycombAndView> GetImageSet()
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

		private static IEnumerable<HoneycombAndView> GetFullImageSet()
		{
			int[] vals = new int[] { 3, 4, 5, 6, 7, 10, 20, -1 };
			foreach( int p in vals )
			foreach( int q in vals )
			foreach( int r in vals )
			{
				yield return new HoneycombAndView()
				{
					P = p,
					Q = q,
					R = r,
				};
			}
		}

		private static IEnumerable<HoneycombAndView> GetImageSetForSlice( PQR constant, int constantValue )
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
				yield return new HoneycombAndView()
				{
					P = iSafe( p ),
					Q = iSafe( q ),
					R = iSafe( r ),	
					//Projection = projection
				};
			}
		}

		private static IEnumerable<HoneycombAndView> GetEuclidImageSet()
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
				yield return new HoneycombAndView()
				{
					P = iSafe( p ),
					Q = iSafe( q ),
					R = iSafe( r )
				};
			}
		}

		private static IEnumerable<HoneycombAndView> GetCoxeterSet()
		{
			yield return new HoneycombAndView() { P = 5, Q = 3, R = 4 };
			yield return new HoneycombAndView() { P = 4, Q = 3, R = 5 };
			yield return new HoneycombAndView() { P = 5, Q = 3, R = 5 };
			yield return new HoneycombAndView() { P = 3, Q = 5, R = 3 };

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

				yield return new HoneycombAndView()
				{
					P = p,
					Q = q,
					R = r
				};
			}
		}

		private static void HoneycombPaper()
		{
			HoneycombAndView imageData = new HoneycombAndView()
			{
				P = 6,
				Q = 5,
				R = 4,
				Projection = Polytope.Projection.VertexCentered
			};

			////////////////////////////////////////////////////////////// Various things we've run over time.
			//Sandbox.CalcSelfSimilarityScale();
			//Sandbox.Check_pq_Distances();
			//HyperidealSquares();
			//S3.Hypercube();
			//R3.Geometry.Euclidean.GenEuclidean();
			//HoneycombGen.OneHoneycombOldCode();
			//CreateCell( imageData );
			//CreateSimplex( imageData );
			//HoneycombGen.OneHoneycombNew( new HoneycombAndView() { P = imageData.P, Q = imageData.Q, R = imageData.R } );
			OneImage( imageData );

			HoneycombAndView[] scaleLarger = GetImageSet().Where( h => 
				Geometry2D.GetGeometry( h.P, h.Q ) == Geometry.Euclidean || 
				Geometry2D.GetGeometry( h.P, h.Q ) == Geometry.Spherical ).ToArray();
			int count = scaleLarger.Length;
			//foreach( HoneycombAndView h in scaleLarger )
			//	Trace.WriteLine( h.FormatFilename() );


			bool batchRun = false;
			if( batchRun )
			{
				HoneycombAndView[] fullSet = GetFullImageSet().ToArray();
				foreach( HoneycombAndView iData in fullSet )
					OneImage( iData );

				int[] rs = new int[] { 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30 };
				foreach( int r in rs )
				{
					imageData.R = r;
					//OneImage( imageData );
				}
			}
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

		private static void CreateSimplex( HoneycombAndView imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;
			

			Vector3D cen = CellCen;
			bool ball = false;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, ref cen, moveToBall: ball );

			//for( int i=0; i<4; i++ )
			//	simplex[i] = H3Models.BallToUHS( simplex[i] );

			// Offset as we do for the boundary images.
			Sphere s = H3Models.UHSToBall( simplex[0] );
			s = CoxeterImages.GeodesicOffset( s, 0.02, ball: true );
			//simplex[0] = H3Models.BallToUHS( s );

			int[] include = new int[] { 0, 1, 2, 3 };	// All facets
			//int[] include = new int[] { 3 };
			File.Delete( "simplex.pov" );
			PovRay.AppendSimplex( simplex, cen, include, "simplex.pov" );
		}

		private static Vector3D CellCen
		{
			get
			{
				Vector3D cen = new Vector3D( 0.05, 0.01, -0.05 );		// 373, 438
				//Vector3D cen = new Vector3D( 0.05, 0.01, 100 );		// 637
				//cen.RotateXY( Math.PI / 2 );	// only if we also rotate simplex mirrors.  XXX - make a setting.
				//Vector3D cen = new Vector3D( 0.1, 0.05, -0.1 );
				return cen;
			}
		}

		private static void CreateCell( HoneycombAndView imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			// Calculate scale to make vertex centered.
			//Vector3D v = SimplexCalcs.VertexPointBall( p, q, r );
			//v = H3Models.BallToUHS( v );

			/*double t = Math.Sqrt( 1 - Math.Pow( v.Abs(), 2 ) );
			v = new Vector3D( t, 0, v.Z );
			v = H3Models.BallToUHS( v );
			Vector3D v2 = H3Models.BallToUHS( new Vector3D( 1, 0, 0 ) );
			t = v.Abs() / v2.Abs();
			v = H3Models.UHSToBall( new Vector3D( t, 0, 0 ) );*/

			double scale = Geometry2D.GetNormalizedCircumRadius( p, q );

			Vector3D cen = CellCen;
			cen = H3Models.BallToUHS( cen );

			bool ball = true;
			bool dual = false;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, ref cen, moveToBall: ball/*, scaling: 1.0/v.Abs()*/ );
			//Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball/*, scaling: 1.0/v.Abs()*/ );
			//Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball, scaling: 1.0 / scale );
			System.IO.File.Delete( "cell.pov" );
			H3.Cell[] simplicesFinal = GenCell( simplex, cen, ball, dual );

			foreach( H3.Cell cell in simplicesFinal )
			{
				//int[] include = new int[] { 0, 1, 2, 3 };
				int[] include = new int[] { 0 };
				if( dual )
					include = new int[] { 3 };
				PovRay.AppendSimplex( cell.Facets.Select( f => f.Sphere ).ToArray(), cell.Center, include, "cell.pov" );
			}
		}

		/// <summary>
		/// Used to generate a regular cell as a set of simplices and save to a Pov-ray file.
		/// This will work for non-finite cells.
		/// </summary>
		private static H3.Cell[] GenCell( Sphere[] simplex, Vector3D cen, bool ball, bool dual )
		{
			// We don't want to include the first mirror (which reflects across cells).
			Sphere[] mirrors = simplex.Skip( 1 ).ToArray();
			if( dual )
				mirrors = new Sphere[] { simplex[0], simplex[1], simplex[2] };
			Sphere[] allMirrors = simplex.ToArray();

			// Simplices will be the "cells" in Recurse.CalcCells.
			H3.Cell.Facet[] simplexFacets = simplex.Select( m => new H3.Cell.Facet( m ) ).ToArray();

			H3.Cell startingCell = new H3.Cell( simplexFacets );
			startingCell.Center = cen;

			FCOrient( startingCell );

			startingCell = startingCell.Clone();	// So our mirrors don't get munged after we reflect around later.
			H3.Cell[] simplices = Recurse.CalcCells( mirrors, new H3.Cell[] { startingCell }, new Recurse.Settings() { Ball = ball } );
			//H3.Cell[] simplices = new H3.Cell[] { startingCell };

			// Layers.
			//int layer = 0;
			//return simplices.Where( s => s.Depths[2] <= layer /*&& s.Depths[0] == 3 && s.Depths[1] == 3*/ ).ToArray();
			return simplices.ToArray();
		}

		/// <summary>
		/// Face centered orientation.
		/// </summary>
		private static void FCOrient( H3.Cell cell )
		{
			// First, need to scale so the lowest triangles are the same size,
			// then reorient so that one triangle is oriented along z axis,
			// then scale so that the triangle is flat.
 
			// Rotation - not what we want.
			if( false )
			{
				Vector3D direction = cell.Facets[0].Sphere.Center;
				Vector3D southPole = new Vector3D( 0, 0, -1 );
				Vector3D axis = direction.Cross( southPole );
				double mag = direction.AngleTo( southPole );

				foreach( Sphere s in cell.Facets.Select( f => f.Sphere ) )
					Sphere.RotateSphere( s, axis, mag );
				Vector3D newCen = cell.Center;
				newCen.RotateAboutAxis( axis, mag );
				cell.Center = newCen;
			}

			if( true )
			{
				// Calculate how much we need to offset to make the cell facet flat.
				Mobius m = new Mobius();
				Sphere s = cell.Facets[0].Sphere;
				double d = s.Center.Abs() - s.Radius;
				m.Isometry( Geometry.Hyperbolic, 0, new Vector3D( 0, d ) );
				foreach( H3.Cell.Facet f in cell.Facets )
					H3Models.TransformInBall2( f.Sphere, m );
				cell.Center = H3Models.TransformHelper( cell.Center, m );
			}
		}

		private static void OneImage( HoneycombAndView imageData )
		{
			//string filename = "batch/" + imageData.FormatFilename();
			string filename = imageData.FormatFilename();
			//if( File.Exists( filename ) )
			//	return;

			int p = imageData.P, q = imageData.Q, r = imageData.R;
			Geometry gCell = Geometry2D.GetGeometry( p, q );
			Geometry gVertex = Geometry2D.GetGeometry( q, r );
			if( !( gCell == Geometry.Hyperbolic || gVertex == Geometry.Hyperbolic ) )
			//if( !( gVertex == Geometry.Hyperbolic ) )
				return;

			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );

			//double bounds = gCell == Geometry.Spherical ? 9 : 1.1;
			double bounds =  1.1;
			//bounds = 1;
			//bounds = 5;
			//bounds = 7;
			bounds = 2.2;

			CoxeterImages imageCalculator = new CoxeterImages();

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

			imageCalculator.AutoCalcScale( settings );
			//settings.ColorScaling = 45.1;	// 437
			//settings.ColorScaling = 60; // 35.58;	// 438, 637 cathedral
			//settings.ColorScaling = 20;	// 373
			//settings.ColorScaling = 15; //464
			if( settings.ColorScaling < 1 )
				settings.ColorScaling = 15;

			size = 1000;
			//settings.Width = size * 4 / 3;
			settings.Width = size;
			settings.Height = size;
			//settings.FileName = imageData.FormatFilename( "jpg" );
			settings.FileName = filename;

			imageCalculator.GenImage( settings );
		}

		private double GetScalingOld( HoneycombAndView imageData )
		{
			// pqi, pir (3qr, p3r, pq3)
			// Q = 3 = 50, Q = 7 = 130
			double scaling = 50 + ( imageData.Q - 3 ) * 20;
			if( /*imageData.P == -1 ||*/ imageData.Q == -1 )
				scaling = 130;

			// iqr
			//scaling = 100; (blue)
			// scaling = 50; (when purple, did constant 50)

			// Reverse: What did I use this for?
			/*scaling = 130 - ( imageData.Q - 3 ) * 20;
			if( imageData.Q == -1 )
				scaling = 50;*/

			// Euclidean cells/vertex figures.
			scaling = 75;

			// Gray scale
			//scaling = 35;

			return scaling;
		}
	}
}
