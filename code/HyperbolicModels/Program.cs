namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.IO;
	using R3.Core;
	using R3.Geometry;

	class Program
	{
		static void Main( string[] args )
		{
			try
			{
				bool wiki = false;
				if( wiki )
				{
					//HoneycombGen.GoursatSet();
					//HoneycombGen.ParacompactSet();

					ViewPath path = new ViewPath();
					Vector3D[] locations = new Vector3D[]
					{
						new Vector3D(0,.1,-.5),
						new Vector3D(-.4,0,-.1),
						new Vector3D(0,-.2,-.5),
						new Vector3D(.1,0,-.8),
						new Vector3D(0,.1,-.5),
						new Vector3D(-.4,0,-.1),
						new Vector3D(0,-.2,-.5),
						new Vector3D(.1,0,-.8),
						new Vector3D(0,.1,-.5),
					};
					List<Vector3D> lookAts = new List<Vector3D>();
					path.Initialize( locations, lookAts );
					HoneycombGen.ViewPath = path;

					for( double t = 0; t < 1; t += .0016 )
					{
						//HoneycombGen.ParacompactAnimationFrame( t );
						path.Step++;
					}

					HoneycombGen.ViewPath = null;
					HoneycombDef def = new HoneycombDef() { P = 3, Q = 4, R = 4 };
					int[] active = new int[] { 1 };
					int baseHue = 220;
					HoneycombGen.Paracompact( def, active, baseHue );
				}

				Settings settings = LoadSettings();
				if( settings.UhsBoundary != null )
				{
					Log( "Generating UHS boundary image for the following honeycomb: " + settings.HoneycombString );
					Log( "Settings...\n" + settings.UhsBoundary.DisplayString );
					HoneycombPaper.OneImage( settings );
				}
			}
			catch( System.Exception ex )
			{
				Log( ex.Message + "\n" + ex.StackTrace );
			}
		}

		public static Settings LoadSettings()
		{
			string filename = "settings.xml";
			//DataContractHelper.SaveToXml( Defaults, filename );
			if( !File.Exists( filename ) )
				return Defaults;

			try
			{
				return (Settings)DataContractHelper.LoadFromXml( typeof( Settings ), filename );
			}
			catch( System.Exception e )
			{
				Log( string.Format( "Failed to load settings. Running with defaults.\n{0}", e.Message ) );
				return Defaults;
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
				return settings;
			}
		}
	}
}
