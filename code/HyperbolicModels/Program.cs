namespace HyperbolicModels
{
	using System.IO;
	using R3.Core;
	using R3.Geometry;

	class Program
	{
		static void Main( string[] args )
		{
			try
			{
				//Sandbox.TileImages();
				//Sandbox.Polarity();

				//Hopf hopf = new Hopf();
				//hopf.GenPovRay();

				//H3Ruled ruled = new H3Ruled();
				//ruled.GenPovRay();

				//Acrohedron a = new Acrohedron();
				//a.Generate( 5, 5, 4 );

				//Sandbox.Dini();

				HoneycombPaper.DoStuff( args );
				//ThreeFifty.GenImage();

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

		public static void OneImage( HoneycombDef imageData, double t = 0.0 )
		{
			//string filename = "batch/" + imageData.FormatFilename();
			string filename = imageData.FormatFilename();
			//if( File.Exists( filename ) )
			//	return;

			int p = imageData.P, q = imageData.Q, r = imageData.R;
			Geometry gCell = Geometry2D.GetGeometry( p, q );
			Geometry gVertex = Geometry2D.GetGeometry( q, r );
			//if( !( gCell == Geometry.Hyperbolic || gVertex == Geometry.Hyperbolic ) )
			//if( !( gVertex == Geometry.Hyperbolic ) )
			//	return;

			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );

			//double bounds = gCell == Geometry.Spherical ? 9 : 1.1;
			double bounds =  1.1;
			bounds = 1;
			//bounds = 5;
			//bounds = 7;
			//bounds = 4.4;
			//bounds = 10;

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

			settings.ColorScaling = 39;
			//size = 4320;
			size = 2500;
			//settings.Width = size * 4 / 3;
			settings.Width = size * 2;
			settings.Width = size;
			settings.Height = size;
			//settings.FileName = imageData.FormatFilename( "jpg" );
			settings.FileName = filename;

			imageCalculator.GenImage( settings, t );
		}
	}
}
