namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;

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
			bool moire = true;
			if( moire )
			{
				double bounds = 50;
				//bounds = 35;		//eisenstein
				//bounds = 450;       // used for super-wide of squared picture
				//bounds = 10;
				// linear or fractional powers need huge bounds
				//bounds = 1;	// hyperbolic
				//bounds = 20;

				//int size = 125;
				//size = 5000;
				//int size = (int)bounds * 50;

				int size = 2000;
				bounds = 25;

				var mSettings = new CoxeterImages.Settings()
				{
					Width = size,
					Height = size,
					Bounds = bounds,
					FileName = "moire.png",
					Antialias = false	// Turning this on can make all kinds of other moire effects I need to understand
				};

				Moire m = new Moire();
				m.BlockNumPixels = .25; // size / m.B needs to be constant
				m.BlockNumPixels = 8;
				//m.BlockNumPixels = 100;
				m.GenImage( mSettings );
				return;
			}

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

				List<int[]> list = new List<int[]>();
				/*list.Add( new int[] { 3, 4, 4 } );
				list.Add( new int[] { 4, 4, 3 } );*/
				//list.Add( new int[] { 4, 4, 4 } );
				//list.Add( new int[] { 3, 6, 3 } );
				//list.Add( new int[] { 6, 3, 6 } );
				/*list.Add( new int[] { 3, 3, 6 } );
				
				list.Add( new int[] { 4, 3, 6 } );
				list.Add( new int[] { 6, 3, 4 } );
				list.Add( new int[] { 5, 3, 6 } );*/
				list.Add( new int[] { 6, 3, 3 } );
				//list.Add( new int[] { 6, 3, 5 } );

				// Go through any settings files.
				foreach( int[] angles in list )
				foreach( string filename in filenames )
				{
					Settings settings = LoadSettings( filename );
					if( settings == null )
						continue;

					settings.Angles = angles;

					// Boundary images.
					if( settings.UhsBoundary != null )
					{
						Log( "\nGenerating UHS boundary image for the following honeycomb:\n" + settings.HoneycombString );
						Log( "\nSettings...\n" + settings.UhsBoundary.DisplayString );
						//HoneycombPaper.OneImage( settings );
						Sections.AnimationSections( settings );
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
