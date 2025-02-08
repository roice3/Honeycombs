namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
    using R3.Core;
    using R3.Drawing;
    using R3.Geometry;
    using R3.Math;

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
			HoneycombPaper.DoStuff( new Settings()
			{
				Angles = new[] { -1, -1, -1 },
				PovRay = new PovRaySettings()
				{
					Active = new int[] { 0 },
				}
			});
			return;

			HoneycombDef def = new HoneycombDef(7, 3, 3);
			//StlGen.HoneycombHyperidealLegs( def );
			StlGen.HoneycombFiniteVertexFig( def );
			return;

			Sections.AnimationSections( new Settings()
			{
				Angles = new int[] { 7, 3, 4 },
			} );
			return;

			HoneycombGen.OneHoneycombOrthoscheme(new Settings()
            {
                Angles = new[] { 3, 3, 6 },
                PovRay = new PovRaySettings()
                {
                    Active = new int[] { 1 },
                    NumEdges = 5000,
                    EdgeWidth = .07
                }
            });
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
