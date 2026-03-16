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
			Sphere[] mirrorsUHS = SimplexCalcs.Mirrors( p, q, r, moveToBall: false );
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
			size = 1000;
			int maxPossibleSize = (int)(Math.Pow( 2, 16 ) - 1); // For a 1-bit image.
			size = maxPossibleSize;
			settings.Width = size;
			settings.Height = size;
			settings.FileName = filename;

			double inrad = Honeycomb.InRadius( p, q, r );
			double circum = Honeycomb.CircumRadius( p, q, r );

			int numSteps = 1;
			for( int i = 0; i < numSteps; i++ )
			//int i = 500;
			//int i = 0;
			{
				Program.Log( "\nSection " + i );
				double t = (double)i / numSteps;

				// Smooth it
				//t = Util.Smoothed( t );

				// slow it at the end
				//t = Math.Pow( t, 0.5 );

				// ZOOM overall view.
				//settings.Bounds = 1.0 - t * .75;
				//settings.Bounds = 1.0 - t * .99;
				//settings.Bounds = .01;
				settings.Bounds = 1.0 / 2;
				//settings.Bounds = 1.0 / 3;

				// Center and radius of cutting circle
				Vector3D cen = new Vector3D(), horoCen = new Vector3D();
				double rad = 0.995;
				rad = .99;

				// 6,3,3 or 4,4,3 horosphere
				if( true )
				{
					// Our Ford circle is the in-sphere
					Vector3D cellCenter = verts[0];
					Vector3D s1 = verts[1];  // a honeycomb face center
					Vector3D s2 = s1;
					s2 = mirrors[1].ReflectPoint( s2 );
					Vector3D s3 = s2;
					s3 = mirrors[2].ReflectPoint( s3 );
					Sphere ford = Sphere.From4Points( s1, s2, s3, cellCenter );
					imageCalculator.m_ford = ford;
					Sphere ford2 = ford.Clone();
					double zOff = .05;
					Vector3D fc = ford2.Center;
					fc.Z += zOff;
					ford2.Center = fc;
					ford2.Radius -= zOff;
					imageCalculator.m_ford2 = ford2;

					// Get the Euclidean cen/rad of our horospherical cut (ball model)
					// The horosphere starting cut is the ford circle.
					Vector3D pointOnZAxis = new Vector3D( 0, 0, 1.0 - ford.Radius * 2 ); 
					double hDist = 5 * t;
					hDist = 12;
					hDist = 14;
					hDist = 2;
					hDist = 8;
					pointOnZAxis = Hyperbolic2D.Offset( pointOnZAxis, hDist );
					double vZ = pointOnZAxis.Z;

					//settings.Bounds = .01;
					System.Console.WriteLine( string.Format( "\n\ttime\t{0}\th-dist\t{1}\tbounds\t{2}", t, hDist, settings.Bounds ) );

					// Special zooming to keep a circle the same size.
					bool specialZoom = false;
					if( specialZoom )
					{
						double eDist = HyperbolicModels.PoincareToUpper( new Vector3D( 0, vZ ) ).Y;
						double uhs = eDist + 1;
						double scale = 2 * Math.Sqrt( uhs - uhs*uhs );
						settings.Bounds = scale;
						System.Console.WriteLine( string.Format( "\n\ttime\t{0}\th-dist\t{1}\tuhs\t{2}\tscale\t{3}", t, hDist, uhs, scale ) );
					}
					
					horoCen = new Vector3D( 0, 0, 1 );
					double horoDiameter = 1 - vZ;
					//horoDiameter += -.02;
					//horoDiameter = 1.99;
					//horoDiameter = 1.75+t/4;
					cen = new Vector3D( 0, 0, 1.0 - horoDiameter / 2 );
					//cen = new Vector3D( 1.0 - horoDiameter / 2, 0, 0 );
					//cen = new Vector3D(1, 0, 0);
					rad = horoDiameter / 2;

					bool transformHoneycomb = false;
					if( transformHoneycomb )
					{
						// One reflection across a cell.
						int[] word = new int[] { 0 };
						word = new int[] { 0, 2, 1, 2, 1 };
						Mobius mobius = CoxeterImages.MobiusFromReflections( mirrors, word );
						mobius = Mobius.InterpMobius( mobius, t );

						imageCalculator.m_mobiusInBall = mobius;
						settings.Bounds = 0.7;
					}

					// This approach isn't working well. Tooooo "Euclidean".
					bool transformHoro = false;
					if( transformHoro )
					{
						// Rotate to another horo center.
						Vector3D north = new Vector3D( 0, 0, 1 );
						Vector3D newHoroCen = mirrors[0].ReflectPoint( north );

						Sphere cutCopy = new Sphere( cen, rad );
						cutCopy.Reflect( mirrors[0] );
						Vector3D newHoroCenCheck = cutCopy.Center;
						newHoroCenCheck.Normalize();

						//t = 1;
						Vector3D moveTo = north + (newHoroCen - north) * t;
						moveTo.Normalize();

						horoCen = moveTo;
						cen = CoxeterImages.TransformBetweenNorthPole( moveTo, cen );

						rad = rad + (cutCopy.Radius - rad) * t;

						//Vector3D axis = new Vector3D( 1, 0, 0 );
						//double angle = t * Math.PI / 2;
						//cen.RotateAboutAxis( axis, angle );
						//horoCen.RotateAboutAxis( axis, angle );

					}

					// This approach isn't working well, which ultimately came down to the Spghere From4Points function
					// not working well for planes. It ended up working in debug and not release!
					bool transformSimplex = false;
					if( transformSimplex )
					{
						t = .1;

						// One reflection across a cell.
						int[] word = new int[] { 0 };
						Mobius mobius = CoxeterImages.MobiusFromReflections( mirrors, word );
						//mobius = Mobius.InterpMobius( mobius, t );
						mobius = Mobius.Identity();

						List<Sphere> transformedMirrors = new List<Sphere>();
						for( int m=0; m<mirrors.Length; m++ )
						{
							Sphere mirror = mirrors[m].Clone();
							if( m == 1 )
							{
								// Breaks things in Release (but not debug)!!!
								mirror.ApplyMobius( mobius );
								mirror.Radius = double.PositiveInfinity;
							}

							transformedMirrors.Add( mirror );
						}

						settings.Mirrors = transformedMirrors.ToArray();
					}

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
				if( false )
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
					double startOffset = 0.99;
					//Vector3D off = Hyperbolic2D.Offset( new Vector3D( 0, 0, startOffset ), startOffset + (5.5 + startOffset) * t );
					Vector3D off = Hyperbolic2D.Offset( new Vector3D( 0, startOffset, 0 ), startOffset + (5.5 + startOffset) * t );

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
                
                imageCalculator.m_cuttingSphereCenter = cen;
				imageCalculator.m_cuttingSphereRad = rad;

				settings.FileName = string.Format( "frame_{0:D4}.png", i );
				//if( File.Exists( settings.FileName ) )
				//	continue;

				// max size testing...
				if( false )
				{
					for( int w = 20000; w <= 50000; w += 5000 )
					{
						// got to this with iterative testing.
						//int maxPossibleSize = (int)(Math.Pow( 2, 16 ) - 1);	// For a 1-bit image.
						//w = maxPossibleSize;
						settings.Width = settings.Height = w;
						imageCalculator.TestSize( settings );
					}
				}

				imageCalculator.GenImage( settings, 0.0 );
			}
		}
	}
}
