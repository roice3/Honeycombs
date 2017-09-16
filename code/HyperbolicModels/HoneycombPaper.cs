namespace HyperbolicModels
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

		public static void OneImage( Settings config, double t = 0.0 )
		{
			HoneycombDef imageData = new HoneycombDef( config.P, config.Q, config.R );
			string filename = imageData.FormatFilename();
			//if( File.Exists( filename ) )
			//	return;

			// These boundary images don't work if the geometry of cells and vertex figures are both spherical.

			int p = imageData.P, q = imageData.Q, r = imageData.R;
			Geometry gCell = Geometry2D.GetGeometry( p, q );
			Geometry gVertex = Geometry2D.GetGeometry( q, r );
			if( gCell == Geometry.Spherical && gVertex == Geometry.Spherical )
				return;

			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );
			double bounds = config.UhsBoundary.Bounds;

			// Calculate the color scale.
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

			Program.Log( "\nGenerating full image..." );
			settings.Width = config.UhsBoundary.ImageWidth;
			settings.Height = config.UhsBoundary.ImageHeight;
			settings.FileName = filename;
			imageCalculator.GenImage( settings, t );
		}
	}
}
