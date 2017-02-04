namespace HyperbolicModels
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
		private static void CullHalfOfEdges( ref H3.Cell.Edge[] edges )
		{
			double thresh = -.01;
			Vector3D looking = new Vector3D( 0, 0, -1 );
			edges = edges.Where( e => e.Start.Dot( looking ) > thresh || e.End.Dot( looking ) > thresh ).ToArray();
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

		/// <summary>
		/// The input time should lie in [0,1]
		/// </summary>
		/// <param name="t"></param>
		public static void ParacompactAnimationFrame( double t )
		{
			var def = new HoneycombDef( 4, 4, 4 );
			var active = new int[] { 1, 2 };
			int baseHue = -1; // Black.
							  //baseHue = 135;
			ViewPath.Time = t;
			OneHoneycombOrthoscheme( def, active, baseHue );
		}
		public static ViewPath ViewPath;

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
				OneHoneycombOrthoscheme( def, active, baseHue );

			baseHue = 220;
			def = new HoneycombDef( 6, 3, 4 );
			foreach( int[] active in toRun )
				OneHoneycombOrthoscheme( def, active, baseHue );

			baseHue = 180;
			def = new HoneycombDef( 6, 3, 5 );
			foreach( int[] active in toRun )
				OneHoneycombOrthoscheme( def, active, baseHue );

			baseHue = 255;
			def = new HoneycombDef( 6, 3, 6 );
			foreach( int[] active in toRun )
				OneHoneycombOrthoscheme( def, active, baseHue );

			baseHue = 105;
			def = new HoneycombDef( 3, 6, 3 );
			foreach( int[] active in toRun )
				OneHoneycombOrthoscheme( def, active, baseHue );

			baseHue = 300;
			def = new HoneycombDef( 4, 4, 3 );
			foreach( int[] active in toRun )
				OneHoneycombOrthoscheme( def, active, baseHue );

			baseHue = 0;
			def = new HoneycombDef( 4, 4, 4 );
			foreach( int[] active in toRun )
				OneHoneycombOrthoscheme( def, active, baseHue );
		}

		private static void CalcThickness( int[] active )
		{
			double thickness = 0.04;
			switch( active.Length )
			{
			case 2:
				thickness = 0.025;
				thickness = 0.04;
				break;
			case 3:
			case 4:
				thickness = 0.01;
				thickness = 0.02;
				//thickness = 0.04;
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

		public static void OneHoneycombGoursat( Settings settings )
		{
			OneHoneycombGoursat( null, null, 0, settings );
		}

		public static void OneHoneycombGoursat( int[] active, string baseName, int baseHue, Settings settings = null )
		{
			// Setup parameters.
			int numEdges = 250000;
			if( settings != null )
			{
				active = settings.PovRay.Active;
				numEdges = settings.PovRay.NumEdges;
				baseName = string.Join( "-", settings.Angles );
			}

			CalcThickness( active );
			if( settings != null )
				H3.m_settings.AngularThickness = settings.PovRay.EdgeWidth; // ZZZ - should really stop using that settings class.

			// Create the simplex.
			Simplex simplex = new Simplex();
			if( settings != null )
				simplex.InitializeGoursat( settings.Angles );
			else
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
			string suffix = "_" + mirrorsString;

			// Convert our active mirrors into the Goursat tet indices.
			int[] polyMirrors = new int[] { 1, 2, 3 };
			active = active.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();
			polyMirrors = polyMirrors.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();

			Vector3D startingPoint = IterateToStartingPoint( null, active, simplex );
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
				edges = Recurse.CalcEdgesSmart( simplex.Facets, startingEdges.ToArray(), numEdges );

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
			/*Vector3D lookFrom = new Vector3D( 1, 1, 1 ) * 0.7;
			Vector3D lookAt = new Vector3D();   // pov-ray lookat
			double thresh = -.01;
			if( doEdges )
				edges = edges.Where( e => e.Start.Dot( lookAt ) > thresh || e.End.Dot( lookAt ) > thresh ).ToArray();
			//if( doCells )
			//	cellsToHighlight = cellsToHighlight.Where( c => c.Center.Dot( lookAt ) > thresh ).ToArray();	// I don't think this works right
			*/

			string fileName = baseName + suffix;
			if( File.Exists( fileName + ".pov" ) )
			{
				File.Delete( fileName + ".pov" );
				//Console.WriteLine( string.Format( "Skipping {0}", fileName ) );
				//return;
			}

			SetupBaseHueGoursat( fileName, mirrorsString, baseHue );

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

		/// <summary>
		/// This was used during batching.
		/// </summary>
		private static void SetupBaseHueGoursat( string fileName, string mirrorsString, int baseHue )
		{
			return;

			// Setup Pov-ray stuff.
			// We have 16 possible mirror states.  We'll calculate the hue by converting the binary state to decimal, and doubling.
			// So for a given family, the hue will range over 32 numbers.
			int hue = baseHue + 2 * Convert.ToInt32( mirrorsString, 2 );
			using( StreamWriter sw = File.CreateText( fileName + ".pov" ) )
			{
				//sw.WriteLine( string.Format( "#declare lookFrom = <{0},{1},{2}>;", lookFrom.X, lookFrom.Y, lookFrom.Z ) );
				//sw.WriteLine( string.Format( "#declare lookAt = <{0},{1},{2}>;", lookAt.X, lookAt.Y, lookAt.Z ) );
				sw.WriteLine( "#include \"C:\\Users\\hrn\\Documents\\roice\\povray\\H3_uniform_faces\\H3_uniform_faces.pov\"" );
				//sw.WriteLine( string.Format( "background {{ CHSL2RGB( <{0}, 1, .3> ) }}", hue ) );
				//sw.WriteLine( string.Format( "background {{ rgb <.13,.37,.31> }}" ) ); 
				sw.WriteLine( string.Format( "background {{ rgb 1 }}" ) );
			}
		}

		/// <summary>
		/// This was used during batching.
		/// </summary>
		private static void SetupBaseHue( string fileName, string mirrorsString, int baseHue )
		{
			return;

			// Setup Pov-ray stuff.
			// We have 16 possible mirror states.  We'll calculate the hue by converting the binary state to decimal, and doubling.
			// So for a given family, the hue will range over 32 numbers.
			int hue = baseHue + 2 * Convert.ToInt32( mirrorsString, 2 );
			using( StreamWriter sw = File.CreateText( fileName + ".pov" ) )
			{
				if( baseHue == -1 )
					sw.WriteLine( "background { rgb 1 }" );
				else
					sw.WriteLine( string.Format( "background {{ CHSL2RGB( <{0}, 1, .1> ) }}", hue ) );
				if( ViewPath != null )
					sw.WriteLine( string.Format( "#declare lookAt = {0};", PovRay.FormatVec( ViewPath.LookAt ) ) );

				// This needs to come last, since it relies on the lookAt.
				sw.WriteLine( "#include \"D:\\roice\\povray\\H3_paracompact.pov\"" );
			}
		}

		private static string BaseName( HoneycombDef def )
		{
			return string.Format( "{0}-{1}-{2}", def.P, def.Q, def.R );
		}

		/// <summary>
		/// Ultimately, I'd like this method to work for material, paracompact, and hyperideal honeycombs.
		/// Some in the last category definitely don't work right (esp. hyperbolic cells).
		/// </summary>
		public static void OneHoneycombOrthoscheme( Settings settings )
		{
			OneHoneycombOrthoscheme( new HoneycombDef(), null, 0, settings );
		}

		public static void OneHoneycombOrthoscheme( HoneycombDef def, int[] active, int baseHue, Settings settings = null )
		{
			// Setup parameters.
			int numEdges = 250000;
			if( settings != null )
			{
				active = settings.PovRay.Active;
				def = new HoneycombDef( settings.P, settings.Q, settings.R );
				numEdges = settings.PovRay.NumEdges;
			}

			CalcThickness( active );
			if( settings != null )
				H3.m_settings.AngularThickness = settings.PovRay.EdgeWidth;	// ZZZ - should really stop using that settings class.

			string baseName = BaseName( def );
			string mirrorsString = ActiveMirrorsString( active );
			string suffix = "-" + mirrorsString;
			string fileName = baseName + suffix;
			if( ViewPath != null )
				fileName += string.Format( "_{0:D4}", ViewPath.Step );

			if( File.Exists( fileName + ".pov" ) )
			{
				File.Delete( fileName + ".pov" );
				//Console.WriteLine( string.Format( "Skipping {0}", fileName ) );
				//return;
			}

			Program.Log( string.Format( "Building {0}", fileName ) );

			// The wiki mirrors are labeled in the reverse of ours.
			Func<int, int> mapMirror = i => 3 - i;
			active = active.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();

			Simplex simplex = new Simplex();
			simplex.Facets = SimplexCalcs.Mirrors( def.P, def.Q, def.R );
			simplex.Verts = SimplexCalcs.VertsBall( def.P, def.Q, def.R );

			Vector3D startingPoint = IterateToStartingPoint( def, active, simplex );
			if( startingPoint.DNE )
				return;
			List<H3.Cell.Edge> startingEdges = new List<H3.Cell.Edge>();
			foreach( int a in active )
			{
				Vector3D reflected = simplex.ReflectInFacet( startingPoint, a );
				startingEdges.Add( new H3.Cell.Edge( startingPoint, reflected ) );
				//startingEdges.Add( new H3.Cell.Edge( simplex.Verts[0], simplex.Verts[3] ) );	// Used for Borromean Rings complement image.
			}

			// If we are doing a view path, transform our geometry.
			if( ViewPath != null )
			{
				//Vector3D p = new Vector3D( 0, 0, .5 );
				Vector3D p = new Vector3D( 0.08, 0.12, 0.07 );
				simplex.Facets = simplex.Facets.Select( f => H3Models.Transform_PointToOrigin( f, p ) ).ToArray();
				simplex.Verts = simplex.Verts.Select( v => H3Models.Transform_PointToOrigin( v, p ) ).ToArray();
				startingEdges = startingEdges.Select( e => new H3.Cell.Edge(
					H3Models.Transform_PointToOrigin( e.Start, p ),
					H3Models.Transform_PointToOrigin( e.End, p ) ) ).ToList();
			}

			SetupBaseHue( fileName, mirrorsString, baseHue );
			Recurse.m_background =  baseHue == -1 ? new Vector3D() : new Vector3D( baseHue, 1, .1 );

			H3.Cell.Edge[] edges = Recurse.CalcEdgesSmart2( simplex.Facets, startingEdges.ToArray(), numEdges );
			//edges = edges.Where( e => e.Depths[0] % 2 == 1 ).ToArray();
			H3.SaveToFile( fileName, edges, finite: true, append: true );

			bool doCells = false;
			H3.Cell[] cellsToHighlight = null;
			if( doCells )
			{
				int[] polyMirrors = new int[] { 1, 2, 3 };
				active = active.Select( i => mapMirror( i ) ).OrderBy( i => i ).ToArray();

				H3.Cell startingCell = PolyhedronToHighlight( Geometry.Hyperbolic, polyMirrors, simplex, startingPoint );
				cellsToHighlight = Recurse.CalcCells( simplex.Facets, new H3.Cell[] { startingCell } );
				H3.AppendFacets( fileName, cellsToHighlight );
			}
		}

		// CHEAT! (would be better to do a geometrical construction)
		// We are going to iterate to the starting point that will make all edge lengths the same.
		public static Vector3D IterateToStartingPoint( HoneycombDef? def, int[] activeMirrors, Simplex simplex )
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
				if( Infinity.IsInfinite( result ) )
					result = double.PositiveInfinity;
				return result;
			};

			// So that we can leverage Euclidean barycentric coordinates, we will first convert our simplex to the Klein model.
			// We will need to take care to properly convert back to the Ball as needed.
			Vector3D[] kleinVerts = simplex.Verts.Select( v => HyperbolicModels.PoincareToKlein( v ) ).ToArray();
			if( def != null )
			{
				HoneycombDef d = def.Value;
				kleinVerts[3] = SimplexCalcs.VertexPointKlein( d.P, d.Q, d.R );
			}

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

			Vector3D euclidean = baryToEuclidean( kleinVerts, bary );
			return HyperbolicModels.KleinToPoincare( euclidean );
		}

		public static H3.Cell PolyhedronToHighlight( Geometry g, int[] mirrors, Simplex simplex, Vector3D startingPoint )
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
	}
}
