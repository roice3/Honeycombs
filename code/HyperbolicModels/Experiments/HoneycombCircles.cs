namespace HyperbolicModels
{
	using R3.Core;
	using R3.Drawing;
	using R3.Geometry;
	using System.Collections.Generic;
	using System.Drawing;
	using System.Linq;
	using Math = System.Math;

	internal class HoneycommbCircles
	{
		public static void Test()
		{
			HoneycombDef def = new HoneycombDef( 5, 3, 4 );
			Simplex simplex = new Simplex();
			simplex.Facets = SimplexCalcs.Mirrors( def.P, def.Q, def.R );

			// Simplices will be the "cells"
			H3.Cell.Facet[] simplexFacets = simplex.Facets.Select( m => new H3.Cell.Facet( m ) ).ToArray();
			H3.Cell startingCell = new H3.Cell( simplexFacets );
			startingCell.AuxPoints = SimplexCalcs.VertsBall( def.P, def.Q, def.R );
			startingCell.Center = HoneycombPaper.InteriorPointBall;

			var cells = CalcCells( simplex.Facets, new H3.Cell[] { startingCell } );

			// Get all the cell centers
			HashSet<Vector3D> centers = new HashSet<Vector3D>();
			foreach( var cell in cells )
			{
				Vector3D cellCen = cell.AuxPoints[0];
				centers.Add( cellCen );
			}

			// Colors.
			Dictionary<double, Color> colors = new Dictionary<double, Color>( new DoubleEqualityComparer() );
			System.Random rand = new System.Random( 0 );

			// Get all the in-spheres
			double inRad = startingCell.AuxPoints[1].Abs();
			//inRad *= 1.16;
			List<Sphere> inSpheres = new List<Sphere>();
			foreach( Vector3D c in centers )
			{
				Vector3D p = c;
				//SphericalModels.GnomonicToStereo( c );
				Geometry g = Geometry.Hyperbolic;

				Vector3D cen;
				double rad;
				H3Models.Ball.DupinCyclideSphere( p, inRad, g, out cen, out rad );
				Sphere i = new Sphere( cen, rad );

				Color color;
				if( !colors.TryGetValue( c.Abs(), out color ) )
				{
					Vector3D rgb = ColorUtil.CHSL2RGB( new Vector3D( rand.NextDouble()*360, .5, .5 ) );
					rgb *= 255;
					color = Color.FromArgb( 255, (int)rgb.X, (int)rgb.Y, (int)rgb.Z );
					colors[c.Abs()] = color;
				}
				i.Color = color;

				inSpheres.Add( i );
			}

			// Project sphere to unit sphere.
			List<Circle3D> circlesOnUnitSphere = new List<Circle3D>();
			foreach( Sphere i in inSpheres )
			{
				if( i.Center.IsOrigin || i.Center.DNE || Infinity.IsInfinite( i.Center ) )
					continue;

				Sphere orthogonal = new Sphere( new Vector3D(), RadiusOrthogonal( i ) );
				Circle3D c = orthogonal.Intersection( i );

				// We need to scale this based on the size of the orthogonal sphere.
				c.Center /= orthogonal.Radius;
				c.Radius /= orthogonal.Radius;
				c.Color = i.Color;
				circlesOnUnitSphere.Add( c );
			}
			Circle3D unit = new Circle3D();
			//circlesOnUnitSphere.Add( unit );

			ProjectAndSave( circlesOnUnitSphere );
		}

		/// <summary>
		/// Returns the radius of the origin-centered sphere orthogonal to s.
		/// The input sphere should not contain the origin.
		/// </summary>
		private static double RadiusOrthogonal( Sphere s )
		{
			double d = s.Center.Abs();
			double r = s.Radius;
			if( r > d )
				throw new System.ArgumentException( "Input sphere should not contain the origin." );
			return Math.Sqrt( d * d - r * r );
		}

		public static H3.Cell[] CalcCells( Sphere[] mirrors, H3.Cell[] cells )
		{
			HashSet<Vector3D> completedCellIds = new HashSet<Vector3D>( cells.Select( c => c.ID ).ToArray() );
			List<H3.Cell> completedCells = new List<H3.Cell>( cells );
			ReflectCellsRecursive( mirrors, cells, completedCells, completedCellIds );
			return completedCells.ToArray();
		}

		private static void ReflectCellsRecursive( Sphere[] simplex, H3.Cell[] cells,
			List<H3.Cell> completedCells, HashSet<Vector3D> completedCellIds )
		{
			if( 0 == cells.Length )
				return;

			List<H3.Cell> newCells = new List<H3.Cell>();

			foreach( H3.Cell cell in cells )
				for( int m = simplex.Length - 1; m >= 0; m-- )
				{
					Sphere mirror = simplex[m];

					if( completedCellIds.Count > 150000 )
						throw new System.Exception( "Maxing out cells - will result in uneven filling." );

					H3.Cell newCell = cell.Clone();
					newCell.Reflect( mirror );

					// This tracks reflections across the cell facets.
					newCell.Depths[m]++;
					newCell.LastReflection = m;

					if( !CellOk( newCell ) )
						continue;

					if( completedCellIds.Add( newCell.ID ) )
					{
						// Haven't seen this cell yet, so 
						// we'll need to recurse on it.
						newCells.Add( newCell );
						completedCells.Add( newCell );
					}
				}

			ReflectCellsRecursive( simplex, newCells.ToArray(), completedCells, completedCellIds );
		}

		internal static bool CellOk( H3.Cell cell )
		{
			foreach( H3.Cell.Facet f in cell.Facets )
			{
				if( f.Sphere.IsPlane )
					continue;

				if( f.Sphere.Radius < .05 )
					return false;
			}

			return true;
		}

			public static void Test1()
		{
			Settings settings = new Settings();
			settings.Angles = new int[] { 5, 3, 3 };
			settings.PovRay = new PovRaySettings() { Active = new int[] { 0 }, NumEdges = 10000 };

			H3.Cell.Edge[] edges = HoneycombGen.OneHoneycombOrthoscheme( settings );
			HashSet<Vector3D> centersOfDualHonycomb = new HashSet<Vector3D>();
			foreach( var e in edges )
			{
				centersOfDualHonycomb.Add( e.Start );
				centersOfDualHonycomb.Add( e.End );
			}

			List<Circle3D> outerCircles = new List<Circle3D>();
			foreach( Vector3D p in centersOfDualHonycomb )
			{
				if( p.IsOrigin )
					continue;
				Circle3D c = GetCircleForBallPoint( p );
				if( c == null )
					continue;
				outerCircles.Add( c );
			}

			ProjectAndSave( outerCircles );
		}

		private static void ProjectAndSave( List<Circle3D> circlesOnUnitSphere )
		{
			List<Circle3D> projected = new List<Circle3D>();
			foreach( Circle3D c in circlesOnUnitSphere )
			{
				Vector3D[] pp = c.RepresentativePoints.Select( p => Sterographic.SphereToPlane( p ) ).ToArray();

				Circle3D cProj = new Circle3D( pp[0], pp[1], pp[2] );
				if( Infinity.IsInfinite( cProj.Radius ) )
					continue;
				cProj.Color = c.Color;
				projected.Add( cProj );
			}

			SaveToBmp( projected );
		}

		private static void SaveToBmp( List<Circle3D> projected )
		{
			int size = 2000;
			Bitmap image = new Bitmap( size, size );
			double b = 5.0;
			ImageSpace i = new ImageSpace( size, size );
			i.XMin = -b;
			i.XMax = b;
			i.YMin = -b;
			i.YMax = b;

			float scale = 0.5f;
			using( Graphics g = Graphics.FromImage( image ) )
			{
				for( int m = 0; m < projected.Count; m++ )
				{
					using( Pen p = new Pen( projected[m].Color, scale * 3.0f ) )
					{
						Circle c = projected[m].ToFlatCircle();
						if( c.IsLine )
						{
							DrawUtils.DrawLine( -c.P2 * 25, c.P2 * 25, g, i, p );   // XXX - not general.
						}
						else
						{
							DrawUtils.DrawCircle( c, g, i, p );
						}
					}
				}
			}

			image.Save( "outerCircles.png" );
		}

		public static Circle3D GetCircleForBallPoint( Vector3D p )
		{
			Sphere ball = new Sphere();

			p = SphericalModels.GnomonicToStereo( p );

			if( Tolerance.GreaterThanOrEqual( p.Abs(), 1 ) )
				return null;

			Sphere t = H3Models.Ball.OrthogonalSphereInterior( p );
			//return H3Models.Ball.IdealCircle( t );

			// Get the corresponding point on the exterior (our inversion).
			p = HyperbolicModels.PoincareToKlein( p );
			p = HyperbolicModels.PoincareToKlein( p );
			p = ball.ReflectPoint( p );

			return GetCircle( p );
		}

		/// <summary>
		/// For a point outside the unit ball, get the ideal circle associate to the dual plane.
		/// </summary>
		public static Circle3D GetCircle( Vector3D p )
		{
			Sphere ball = new Sphere();

			// Our point defines an orthogonal cone to the unit ball.
			double d = p.Abs();
			double radius = Math.Sqrt( d * d - 1 );
			Sphere s = new Sphere( p, radius );

			return ball.Intersection( s );
		}
	}
}
