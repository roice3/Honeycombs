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
						HoneycombGen.ParacompactAnimationFrame( t );
						path.Step++;
					}
				}
				else
				{ 
					HoneycombPaper.DoStuff( args );
				}
			}
			catch( System.Exception ex )
			{
				System.Diagnostics.Trace.WriteLine( ex.Message + "\n" + ex.StackTrace );
				System.Console.WriteLine( ex.Message + "\n" + ex.StackTrace );
			}
		}

		public static void OneImage( HoneycombDef imageData, double t = 0.0 )
		{
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

			double bounds =  1.1;
			bounds = 1;

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

			//size = 4320;
			size = 1000;
			//settings.Width = size * 4 / 3;
			settings.Width = size * 2;
			settings.Width = size;
			settings.Height = size;
			settings.FileName = filename;

			imageCalculator.GenImage( settings, t );
		}
	}
}
