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
		public static void AnimationSections( Settings config )
		{
			HoneycombDef imageData = new HoneycombDef( config.P, config.Q, config.R );
			int p = imageData.P, q = imageData.Q, r = imageData.R;

			string filename = imageData.FormatFilename();

			/*Simplex s = new Simplex();
			s.InitializeGoursat(new int[] { 2, 5, 3, 2, 3, 3 });
			Sphere[] mirrors = s.Facets;
			Vector3D[] verts = s.Verts;*/

			Sphere[] mirrors = SimplexCalcs.Mirrors( p, q, r );
			Vector3D[] verts = SimplexCalcs.VertsBall( p, q, r );
			double bounds = 1.0; //config.UhsBoundary.Bounds;
            //bounds = 12.0;
            //bounds = 0.25;
            //bounds = 0.75;

            // Calculate the color scale.
            int size = 200;
			CoxeterImages.Settings settings = new CoxeterImages.Settings()
			{
				Honeycomb = imageData,
				G = Honeycomb.GetGeometry( p, q, r ),
				Width = size,
				Height = size,
				Bounds = bounds,
				Mirrors = mirrors,
				Verts = verts,
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

			double inrad = Honeycomb.InRadius( p, q, r );
			double circum = Honeycomb.CircumRadius( p, q, r );

			int numSteps = 20;
			for( int i = 0; i <= numSteps; i++ )
			//int i = 20;
			//int i = 0;
			{
				Program.Log( "\nSection " + i );
				double t = Util.Smoothed( (double)i / numSteps );
				//t = Math.Pow( (double)i / numSteps, 1.0 );

				// Center and radius of cutting circle
				Vector3D cen = new Vector3D();
				double rad = 0.995;

				// 6,3,3 horosphere
				if( false )
				{
					Vector3D v = verts[3];
					double startOffset = .5;
					v = Hyperbolic2D.Offset( v, +startOffset + (5.5 + startOffset) * t );
					double vZ = v.Z;
					double horoDiameter = 1 - vZ;
					//horoDiameter += -.02;
					//horoDiameter = 1.99;
					horoDiameter = 1.75+t/4;
					cen = new Vector3D( 0, 0, 1.0 - horoDiameter / 2 );
					//cen = new Vector3D(1, 0, 0);
					rad = horoDiameter / 2;


					bool applyMob = false;
					if( applyMob )
					{
						//double z = DonHatch.h2eNorm( 2 * circum * i / numSteps );
						double z = vZ; //- DonHatch.h2eNorm( -.25 + 2.5 * t );
						Mobius m = new Mobius();
						m.Isometry( Geometry.Hyperbolic, 0, new System.Numerics.Complex( 0, z ) );
						//imageCalculator.m_z = m;
					}
				}

				// Center of ball
				//cen = new Vector3D();
				//rad = 0.95 + t * 0.49;

				// 7,3,3 cell
				if( true )
				{
					// The offset point.
					Vector3D off = verts[3];
					//off = Hyperbolic2D.Offset( off, t * 5 );
					//off = Hyperbolic2D.Offset(off, -0.5 + t);
					//off = Hyperbolic2D.Offset( off, t*3 );
					off = Hyperbolic2D.Offset( off, -0.5 + t*3 );

					// Cell hypercycle
					// This only works if verts[3] is material, i.e. magnitude < 1.
					// I should be able to get the ideal circle below though, because clearly I normalized these in the wiki images. 
					// https://en.wikipedia.org/wiki/Template:Regular_honeycomb_table Look at 73q images.
					Vector3D p1 = verts[3];
					Vector3D p2 = mirrors[3].ReflectPoint( p1 );
					Vector3D p3 = mirrors[2].ReflectPoint( p2 );
					Vector3D p4 = mirrors[1].ReflectPoint( p3 );
					Sphere hyperCycle = Sphere.From4Points( p1, p2, p3, p4 );

					// Ideal circle.
					Sphere ball = new Sphere();
					Circle3D ideal = hyperCycle.Intersection( ball );

					// New offset hypercycle
					Vector3D[] idealPoints = ideal.RepresentativePoints;
					Sphere grown = Sphere.From4Points( idealPoints[0], idealPoints[1], idealPoints[2], off );

					cen = grown.Center;
					rad = grown.Radius;
				}

				// H-plane slices
				if( false )
				{
					double startOffset = 0.001;
					Vector3D off = Hyperbolic2D.Offset( new Vector3D( 0, 0, startOffset ), startOffset + (5.5 + startOffset) * t );

					Sphere slice = H3Models.Ball.OrthogonalSphereInterior( off );
					cen = slice.Center;
					rad = slice.Radius;
				}

				//double rad = DonHatch.h2eNorm( inrad * 0.85 + 2.9 * circum * i / numSteps );
				//imageCalculator.m_r = DonHatch.h2eNorm( circum * 1.7 ); //0.85; 

				if( false )
				{
					rad = Spherical2D.s2eNorm( inrad * 0.5 + 9 * inrad * i / numSteps );
				}

				// Rob's stuff, 733
				if( false )
				{
					cen = new Vector3D( 0, 0, .5 );
					rad = 1.3;
				}
                
                imageCalculator.m_cen = cen;
				imageCalculator.m_r = rad;

				settings.FileName = string.Format( "frame_{0:D4}.png", i );
				//if( File.Exists( settings.FileName ) )
				//	continue;

				imageCalculator.GenImage( settings, 0.0 );
			}
		}
	}
}
