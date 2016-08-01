﻿namespace HyperbolicModels
{
	using System;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;
	using R3.Math;

	public static class HoneycombPaper
	{
		public static void DoStuff( string[] args )
		{
			int p, q, r;
			if( args.Length == 3 )
			{
				p = ReadArg( args[0] );
				q = ReadArg( args[1] );
				r = ReadArg( args[2] );
			}
			else
			{
				p = 3;
				q = 3;
				r = 7;
			}

			HoneycombDef imageData = new HoneycombDef( p, q, r );

			////////////////////////////////////////////////////////////// Various things we've run over time.
			//Sandbox.CalcSelfSimilarityScale();
			//Sandbox.Check_pq_Distances();
			//HyperidealSquares();
			//S3.Hypercube();
			//R3.Geometry.Euclidean.GenEuclidean();
			//HoneycombGen.OneHoneycombOldCode();
			//CreateCellPovRay( imageData );
			//CreateCellSTL( imageData );
			//CreateSimplex( imageData );
			//HoneycombGen.OneHoneycombNew( new HoneycombDef() { P = imageData.P, Q = imageData.Q, R = imageData.R } );
			//SphericalAnimate( imageData );
			Program.OneImage( imageData );

			HoneycombDef[] scaleLarger = GetImageSet().Where( h =>
				Geometry2D.GetGeometry( h.P, h.Q ) == Geometry.Euclidean ||
				Geometry2D.GetGeometry( h.P, h.Q ) == Geometry.Spherical ).ToArray();
			int count = scaleLarger.Length;
			//foreach( HoneycombAndView h in scaleLarger )
			//	Trace.WriteLine( h.FormatFilename() );

			BatchRun();
		}

		private static int ReadArg( string arg )
		{
			if( arg == "i" || arg == "inf" )
				return -1;

			return int.Parse( arg );
		}

		private static void BatchRun()
		{
			HoneycombDef imageData;
			bool batchRun = false;
			if( batchRun )
			{
				HoneycombDef[] fullSet = GetFullImageSet().ToArray();
				foreach( HoneycombDef iData in fullSet )
					Program.OneImage( iData );

				int[] rs = new int[] { 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30 };
				foreach( int r in rs )
				{
					imageData.R = r;
					//OneImage( imageData );
				}
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

			int fps = 30;
			//int frames = 60 * fps;
			int frames = 5;
			for( int i = 0; i < frames; i++ )
			{
				string num = i.ToString();
				num = num.PadLeft( 3, '0' );

				double t = (double)i/frames;
				string filename = "batch/" + imageData.FormatFilename( string.Empty ) + "_" + num + ".png";
				OneAnimationFrame( imageData, filename, colorScaling, t );
			}
		}

		private static double AnimColorScaling( HoneycombDef imageData )
		{
			return 10;

			int p = imageData.P, q = imageData.Q, r = imageData.R;
			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );

			int size = 200;
			CoxeterImages.Settings settings = new CoxeterImages.Settings()
			{
				Honeycomb = imageData,
				Width = size*2,
				Height = size,
				Bounds = 1.0,
				Mirrors = mirrors,
				FileName = string.Empty,
			};

			CoxeterImages imageCalculator = new CoxeterImages();
			imageCalculator.AutoCalcScale( settings );
			return settings.ColorScaling;
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

		internal static Vector3D CellCenBall
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

		private static void CreateCellPovRay( HoneycombDef imageData )
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

			Vector3D cen = CellCenBall;
			cen = H3Models.BallToUHS( cen );

			bool ball = true;
			bool dual = false;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, ref cen, moveToBall: ball/*, scaling: 1.0/v.Abs()*/ );
			//Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball/*, scaling: 1.0/v.Abs()*/ );
			//Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball, scaling: 1.0 / scale );
			H3.Cell[] simplicesFinal = GenCell( simplex, null, cen, ball, dual );

			System.IO.File.Delete( "cell.pov" );
			foreach( H3.Cell cell in simplicesFinal )
			{
				//int[] include = new int[] { 0, 1, 2, 3 };
				int[] include = new int[] { 0 };
				if( dual )
					include = new int[] { 3 };

				Sphere[] facets = cell.Facets.Select( f => f.Sphere ).ToArray();
				if( m_toKlein )
					facets = facets.Select( s => H3Models.BallToKlein( s ) ).ToArray();

				PovRay.AppendSimplex( facets, cell.Center, include, "cell.pov" );
			}
		}

		private static bool m_toKlein = true;

		/// <summary>
		/// Used to generate a regular cell as a set of simplices and save to a Pov-ray file.
		/// This will work for non-finite cells.
		/// </summary>
		private static H3.Cell[] GenCell( Sphere[] simplex, Mesh mesh, Vector3D cen, bool ball, bool dual )
		{
			// We don't want to include the first mirror (which reflects across cells).
			Sphere[] mirrors = simplex.Skip( 1 ).ToArray();
			if( dual )
				mirrors = new Sphere[] { simplex[0], simplex[1], simplex[2] };
			Sphere[] allMirrors = simplex.ToArray();
			mirrors = simplex.ToArray();

			// Simplices will be the "cells" in Recurse.CalcCells.
			H3.Cell.Facet[] simplexFacets = simplex.Select( m => new H3.Cell.Facet( m ) ).ToArray();

			H3.Cell startingCell = new H3.Cell( simplexFacets );
			startingCell.Center = cen;
			startingCell.Mesh = mesh;

			//FCOrient( startingCell );

			startingCell = startingCell.Clone();	// So our mirrors don't get munged after we reflect around later.
			H3.Cell[] simplices = Recurse.CalcCells( mirrors, new H3.Cell[] { startingCell }, new Recurse.Settings() { Ball = ball } );
			//H3.Cell[] simplices = new H3.Cell[] { startingCell };

			// Layers.
			int layer = 0;
			return simplices.Where( s => s.Depths[0] <= layer /*&& s.Depths[0] == 3 && s.Depths[1] == 3*/ ).ToArray();
			//return simplices.ToArray();
		}

		/// <summary>
		/// Create an STL file for a cell.
		/// This will be a 2D surface of triangles, which will need to be thickened using an external program.
		/// Currently only works for cells with both hyperideal vertices and cells.
		/// </summary>
		private static void CreateCellSTL( HoneycombDef imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			bool ball = false;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball );
			H3.Cell.Edge[] edges;
			if( ball )
				edges = SimplexCalcs.SimplexEdgesBall( p, q, r );
			else
				edges = SimplexCalcs.SimplexEdgesUHS( p, q, r );

			// Two edges of one simplex facet.
			int div = 5;
			H3.Cell.Edge e1 = edges[2];
			H3.Cell.Edge e2 = edges[3];
			Vector3D[] points1, points2;
			if( ball )
			{
				points1 = H3Models.Ball.GeodesicPoints( e1.Start, e1.End, 2 * div );
				points2 = H3Models.Ball.GeodesicPoints( e2.Start, e2.End, 2 * div );
			}
			else
			{
				points1 = H3Models.UHS.GeodesicPoints( e1.Start, e1.End, 2 * div );
				points2 = H3Models.UHS.GeodesicPoints( e2.Start, e2.End, 2 * div );
			}

			Sphere cellSphere = simplex[0];

			// Because one vertex the facet triangle is hyperideal, it will actually look like a square.
			List<Vector3D[]> allPoints = new List<Vector3D[]>();
			for( int i = 0; i < points1.Length; i++ )
			{
				Vector3D p1 = points1[i];
				Vector3D p2 = points2[i];

				// NOTE: This arc is not generally geodesic!
				Vector3D[] arcPoints;
				if( i == points1.Length - 1 )
				{
					arcPoints = ball ?
						H3Models.Ball.GeodesicPoints( p1, p2, div ) :
						H3Models.UHS.GeodesicPoints( p1, p2, div );
				}
				else
				{
					Circle3D c = Circle3D.FromCenterAnd2Points( cellSphere.Center, p1, p2 );
					double angleTot = (p1 - c.Center).AngleTo( p2 - c.Center );
					arcPoints = Shapeways.CalcArcPoints( cellSphere.Center, cellSphere.Radius, p1, c.Normal, -angleTot, div );
				}
				//Vector3D[] arcPoints = new Vector3D[] { p1, p2 };
				allPoints.Add( arcPoints );
			}

			// Create the triangles for the patch.
			Mesh mesh = new Mesh();
			for( int i = 0; i < allPoints.Count - 1; i++ )
			{
				Vector3D[] arc1 = allPoints[i];
				Vector3D[] arc2 = allPoints[i + 1];

				for( int j = 0; j < arc1.Length - 1; j++ )
				{
					// Points of (i,j) box;
					Vector3D p1 = arc1[j];
					Vector3D p2 = arc2[j];
					Vector3D p3 = arc1[j + 1];
					Vector3D p4 = arc2[j + 1];

					mesh.Triangles.Add( new Mesh.Triangle( p1, p2, p3 ) );
					mesh.Triangles.Add( new Mesh.Triangle( p2, p4, p3 ) );
				}
			}

			// Face centered orientation.
			// XXX - We need to do this prior to mesh generation, so the mesh isn't stretched out by these transformations.
			bool faceCentered = true;

			Vector3D cen = CellCenBall;
			if( faceCentered )
				SimplexCalcs.PrepForFacetCentering( p, q, simplex, ref cen );

			Mobius mUHS = SimplexCalcs.FCOrientMobius( p, q );
			Mobius mBall = FCOrientMobius( H3Models.UHSToBall( cellSphere ) );

			simplex = simplex.Select( s =>
			{
				s = H3Models.UHSToBall( s );
				H3Models.TransformInBall2( s, mBall );
				return s;
			} ).ToArray();

			{
				for( int i = 0; i < mesh.Triangles.Count; i++ )
				{
					Mesh.Triangle tri = mesh.Triangles[i];

					if( faceCentered )
					{
						tri.a = mUHS.ApplyToQuaternion( tri.a );
						tri.b = mUHS.ApplyToQuaternion( tri.b );
						tri.c = mUHS.ApplyToQuaternion( tri.c );
					}

					tri.a = H3Models.UHSToBall( tri.a );
					tri.b = H3Models.UHSToBall( tri.b );
					tri.c = H3Models.UHSToBall( tri.c );

					if( faceCentered )
					{
						tri.a = H3Models.TransformHelper( tri.a, mBall );
						tri.b = H3Models.TransformHelper( tri.b, mBall );
						tri.c = H3Models.TransformHelper( tri.c, mBall );
					}
					mesh.Triangles[i] = tri;
				}

				if( faceCentered )
					cen = H3Models.TransformHelper( cen, mBall );
			}

			// Now we need to reflect around this fundamental patch.
			bool dual = false;
			H3.Cell[] simplicesFinal = GenCell( simplex, mesh, cen, ball, dual );

			System.IO.File.Delete( "cell.stl" );
			foreach( H3.Cell cell in simplicesFinal )
				STL.AppendMeshToSTL( cell.Mesh, "cell.stl" );
			//STL.AppendMeshToSTL( simplicesFinal[0].Mesh, "cell.stl" );
		}

		private static void CreateSimplex( HoneycombDef imageData )
		{
			int p = imageData.P;
			int q = imageData.Q;
			int r = imageData.R;

			Vector3D cen = CellCenBall;
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

			if( false )
			{
				H3.Cell.Edge[] edges = SimplexCalcs.SimplexEdgesUHS( p, q, r );
				PovRay.WriteEdges( new PovRay.Parameters { Halfspace = true, AngularThickness = 0.03 },
					Geometry.Hyperbolic, edges, "simplex.pov", append: true );
			}
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
				Mobius m = FCOrientMobius( cell.Facets[0].Sphere );
				foreach( H3.Cell.Facet f in cell.Facets )
					H3Models.TransformInBall2( f.Sphere, m );
				cell.Center = H3Models.TransformHelper( cell.Center, m );
			}
		}

		/// <summary>
		/// A Mobius that will put our simplex into a face centered orientation.
		/// Meant to be used with H3Models.TransformInBall2 or H3Models.TransformHelper
		/// </summary>
		private static Mobius FCOrientMobius( Sphere cellSphereInBall )
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

		private static double GetScalingOld( HoneycombDef imageData )
		{
			// pqi, pir (3qr, p3r, pq3)
			// Q = 3 = 50, Q = 7 = 130
			double scaling = 50 + (imageData.Q - 3) * 20;
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
