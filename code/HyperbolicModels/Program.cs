namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using R3.Geometry;

	class Program
	{
		/*
		Known problems:
		* Orthoscheme code may not be working out of the box for:
			- spherical honeycombs
			- honeycomb with hyperideal cells
		*/

		static void Main( string[] args )
		{
			PovRaySettings povray = new PovRaySettings() { Active = new int[] { 0 }, NumEdges = 15000, EdgeWidth = 0.1 };
			Settings set = new Settings() { Angles = new int[] { 4, 3, 5 }, PovRay = povray };
			//HoneycombGen.OneHoneycombOrthoscheme( set );
			HoneycombPaper.DoStuff( set );

			/*for( int level = 0; level <= 8; level++ )
			{
				H3.m_maxLevel = level;
				H3.GenHoneycomb( EHoneycomb.H435 );
			}*/

			//HoneycombGen_old.OneHoneycombOldCode();
			//for( int level = 1; level <= 14; level++ )
			//	Flags.GenForTiling( 4, 6, level );

			for( int level = 18; level <= 20; level++ )
				Flags.GenForTiling( 3, 7, level );

			//Flags.GenForTiling( 3, 7, 3 );
			//Flags.GenForTilingDistanceLayers( 3, 7 );

			//for( double co = 0.0; co <= .99; co += 0.1 )
			//	StlGen.H3Helicoid( co );
			return;


			return;
			//StlGen.S3BiHelicoid();
			//HoneycombGen_old.OneHoneycombNew( new HoneycombDef() { P = 4, Q = 3, R = 4 } );
			//PointGroups pg = new PointGroups();
			//pg.Gen( 5, 3, 3 );
			R3.Math.Mobius m = new R3.Math.Mobius();
			m.UpperHalfPlane();
			R3.Math.Mobius m2 = m.Inverse();

			Vector3D v = new Vector3D( 0, 1, 0 );
			Sterographic.NormalizeToHyperboloid( ref v );

			H3Models.BallToUHS( new Vector3D( 0, 1, 0 ) );

			Sphere[] mirrorsBall = SimplexCalcs.Mirrors( 5, 3, 4, true );
			Sphere[] mirrors = mirrorsBall.Select( s => H3Models.BallToKlein( s ) ).ToArray();
			Vector3D[] verts = SimplexCalcs.VertsBall( 5, 3, 4 ).Select( p => HyperbolicModels.PoincareToKlein( p ) ).ToArray();

			Vector3D testPoint = ( verts[0] + verts[1] + verts[2] + verts[3] ) / 4;
			foreach( Sphere s in mirrors )
			{
				System.Diagnostics.Trace.WriteLine( s.IsPointInside( testPoint ) ? "inside" : "outside" );
			}

			//Vector3D t0 = new Vector3D( .45, 0, -.59 );
			Vector3D t0 = HyperbolicModels.KleinToPoincare( new Vector3D( .45, 0, -.59 ) );
			Vector3D h0 = Sterographic.PoincareBallToHyperboloid( t0 );
			Vector3D t1 = mirrorsBall[0].ReflectPoint( t0 );
			Vector3D h1 = Sterographic.PoincareBallToHyperboloid( t1 );

			//UhsBoundarySettings uhs = new UhsBoundarySettings() { Bounds = 1, ImageHeight = 1000, ImageWidth = 1000 };
			//HoneycombPaper.OneImage( new Settings() { Angles = new int[] { 5, 3, 4 }, UhsBoundary = uhs } );

			return;

			try
			{
				List<string> filenames = new List<string>();
				if( args.Length > 0 &&
					File.Exists( args[0] ) )
				{
					filenames.Add( args[0] );
				}
				else
				{
					filenames = Directory.EnumerateFiles( ".", "*.xml", SearchOption.TopDirectoryOnly ).ToList();
				}

				// Go through any settings files.
				foreach( string filename in filenames )
				{
					Settings settings = LoadSettings( filename );
					if( settings == null )
						continue;

					// Boundary images.
					if( settings.UhsBoundary != null )
					{
						Log( "\nGenerating UHS boundary image for the following honeycomb:\n" + settings.HoneycombString );
						Log( "\nSettings...\n" + settings.UhsBoundary.DisplayString );
						HoneycombPaper.OneImage( settings );
					}

					// POV-Ray definition files.
					if( settings.PovRay != null )
					{
						Log( "\nGenerating POV-Ray definition file for the following honeycomb:\n" + settings.HoneycombString );
						Log( "\nSettings...\n" + settings.PovRay.DisplayString );

						if( settings.Angles.Length == 3 )
							HoneycombGen.OneHoneycombOrthoscheme( settings );
						else if( settings.Angles.Length == 6 )
							HoneycombGen.OneHoneycombGoursat( settings );
					}
				}
			}
			catch( System.Exception ex )
			{
				Log( ex.Message + "\n" + ex.StackTrace );
			}
		}

		public static Settings LoadSettings( string filename )
		{
			//DataContractHelper.SaveToXml( Defaults, filename );
			if( !File.Exists( filename ) )
				return Defaults;

			try
			{
				return (Settings)DataContractHelper.LoadFromXml( typeof( Settings ), filename );
			}
			catch( System.Exception e )
			{
				Log( string.Format( "Failed to load settings from file '{0}', so skipping.\n{1}", e.Message ) );
				return null;
			}
		}

		public static void Log( string message )
		{
			System.Diagnostics.Trace.WriteLine( message );
			System.Console.WriteLine( message );
		}

		public static Settings Defaults
		{
			get
			{
				Settings settings = new Settings();
				settings.Angles = new int[] { 3, 3, 7 };
				settings.UhsBoundary = new UhsBoundarySettings() { Bounds = 1.0, ImageHeight = 1200, ImageWidth = 1200 };
				settings.PovRay = new PovRaySettings() { Active = new int[] { 1, 0, 0, 0 }, NumEdges = (int)5e5, EdgeWidth = 0.02 };
				return settings;
			}
		}
	}
}
