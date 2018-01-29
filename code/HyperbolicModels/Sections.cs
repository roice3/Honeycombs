namespace HyperbolicModels
{
	using System;
	using System.Collections.Generic;
	using System.Drawing;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Drawing;
	using R3.Geometry;
	using R3.Math;

	public class Sections
	{
		public void AnimationSections( Settings config )
		{
			HoneycombDef imageData = new HoneycombDef( config.P, config.Q, config.R );
			int p = imageData.P, q = imageData.Q, r = imageData.R;

			string filename = imageData.FormatFilename();

			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );
			double bounds = 1.0; //config.UhsBoundary.Bounds;
			bounds = 9.0;

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
			//imageCalculator.AutoCalcScale( settings );
			if( settings.ColorScaling < 1 )
				settings.ColorScaling = 15;
			settings.ColorScaling = 11;

			Program.Log( "\nGenerating sections..." );
			size = 500;
			settings.Width = size;
			settings.Height = size;
			settings.FileName = filename;

			double max = Spherical2D.e2sNorm( 15 );
			double min = Spherical2D.e2sNorm( 1.0 / 15 );
			DonHatch.e2hNorm( max );
			int numSteps = 1800; // 1 minute
			double step = (max - min) / numSteps;
			for( int i = 0; i < 1; i++ )
			{
				Program.Log( "\nSection " + i );
				imageCalculator.m_z = 1.0 / 0.5;
				Spherical2D.s2eNorm( min + step * i );
				DonHatch.h2eNorm( step * i );
				settings.FileName = string.Format( "533_{0:D4}.png", i );
				imageCalculator.GenImage( settings );
			}
		}
	}
}
