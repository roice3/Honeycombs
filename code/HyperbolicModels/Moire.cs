namespace HyperbolicModels
{
	using R3.Core;
	using R3.Geometry;
	using System.Collections.Generic;
	using System.Drawing;
	using System.Drawing.Imaging;
	using Math = System.Math;
	using System.Numerics;
	using System.Threading.Tasks;
	using System;
	using R3.Math;

	internal class Moire
	{
		private readonly object m_lock = new object();

		private double RandInBounds( Random random, CoxeterImages.Settings settings )
		{
			double randDouble = random.NextDouble();
			double scaled = (randDouble - 0.5) * settings.Bounds * 2;
			return scaled;
		}

		private void SetupGrid( CoxeterImages.Settings settings )
		{
			m_grid = new NearTree( Metric.Euclidean );
			int id = 0;

			// Random
			double numPoints = (double)settings.Width / BlockNumPixels * settings.Height / BlockNumPixels;
			Random random = new Random();
			for( int i = 0; i < numPoints; i++ )
			{
				Complex location = new Complex( RandInBounds( random, settings ), RandInBounds( random, settings ) );
				m_grid.InsertObject( new NearTreeObject() { ID = id, Location = Vector3D.FromComplex( location ) } );
			}

			// Eisenstein
			/*
			// We don't need to make an entire tiling, which ends up slow with as many as we need.
			// We can just use Eisenstein integers.
			int range = 10;
			double scale = .2;
			Console.WriteLine( "Eisenstein Integers:" );
			for( int a = -range; a <= range; a++ )
			{
				for( int b = -range; b <= range; b++ )
				{
					id++;
					Complex omega = new Complex( -0.5, Math.Sqrt( 3 ) / 2 ); // Primitive cube root of unity
					Complex location = ( new Complex( a, 0 ) + omega * b ) * scale;
					m_grid.InsertObject( new NearTreeObject() { ID = id, Location = Vector3D.FromComplex( location ) } );
				}
			}
			*/

			/*TilingConfig config = new TilingConfig( 3, 6, maxTiles: 1000000 );	// slow and memory insane
			Tiling tiling = new Tiling();
			tiling.Generate( config );
			foreach( Tile t in tiling.Tiles )
			{
				id++;
				m_grid.InsertObject( new NearTreeObject() { ID = id, Location = t.Center } );
			}
			*/
		}

		NearTree m_grid;

		public void GenImage( CoxeterImages.Settings settings )
		{
			//SetupGrid( settings );

			int width = settings.Width;
			int height = settings.Height;
			Bitmap image = new Bitmap( width, height );

			// Cycle through all the pixels and calculate the color.
			int row = 0;
			double bounds = settings.Bounds;
			double xoff = settings.XOff;
			double yoff = settings.YOff;
			Parallel.For( 0, width, i =>
			{
				if( row++ % 20 == 0 )
					System.Console.WriteLine( string.Format( "Processing Line {0}", row ) );

				for( int j = 0; j < height; j++ )
				{
					double x = -bounds + i * xoff;
					double y = -bounds + j * yoff;

					// Doesn't seem to make any real difference here.
					// Maybe should be in round functions.
					double continuousTranslate = 0.0;
					x += continuousTranslate;
					y += continuousTranslate;

					if( settings.Antialias )
					{
						const int div = 3;
						//const int div = 1;
						List<Color> colors = new List<Color>();
						for( int k = 0; k <= div; k++ )
							for( int l = 0; l <= div; l++ )
							{
								double xa = x - xoff / 2 + k * xoff / div;
								double ya = y - yoff / 2 + l * yoff / div;
								Vector3D v = new Vector3D( xa, ya );

								Color color = CalcColor( settings, v );
								colors.Add( color );
							}

						lock( m_lock )
						{
							Color avg = CoxeterImages.AvgColor( colors );
							image.SetPixel( i, j, avg );
						}
					}
					else
					{
						lock( m_lock )
						{
							Vector3D v = new Vector3D( x, y );
							image.SetPixel( i, j, CalcColor( settings, v ) );
						}
					}
				}
			} );

			image.Save( settings.FileName, ImageFormat.Png );
		}

		private Color CalcColor( CoxeterImages.Settings settings, Vector3D v )
		{
			//v.X += v.Y * Math.Tan( System.Math.PI / 6 );

			// Round to closest grid point.
			// THIS IS WHERE WE CHOOSE SQUARE V EISENSTEIN
			//Vector3D quantized = RoundGaussian( settings, v );
			Vector3D quantized = RoundEisenstein( settings, v );

			// Neartree makes our grid. Ugh, too slow once I scale.
			/*NearTreeObject closest;
			m_grid.FindNearestNeighbor( out closest, v, 1 );
			if( closest == null )
				return Color.Black;
			Vector3D quantized = closest.Location;
			*/

			//double mag = quantized.Abs();
			double ds2 = quantized.X * quantized.X + quantized.Y * quantized.Y;    // euclidean	
			//double ds2 = quantized.X * quantized.X - quantized.Y * quantized.Y;        // lorentz metric


			double scaled = ds2;	// use with lorentz??
			//double scaled = mag * mag; // Why is squared so special? maybe because it is the square of the metric??
			//double scaled = Math.Pow( mag, .5 );
			//double scaled = Math.Pow( Math.E, mag );
			//double scaled = Math.Pow( 2.0, mag );
			//double scaled = Math.Pow( mag, .5 ) * 1000;	// fractional powers need be scaled up a ton to get moire
			//double scaled = Math.Log( mag ) * 2500;
			//double scaled = Spherical2D.e2sNorm( mag ) * 2000;
			//double scaled = DonHatch.e2hNorm( mag ) * 2000;

			double quantizedTransform = Math.PI * 0.0;

			double intensity = (1 + Math.Sin( scaled + quantizedTransform )) / 2;

			// spike it
			double spikiness = 3;	// good val is 3
			intensity = Math.Pow( intensity, spikiness );

			Color color = Color.Blue;
			return ColorUtil.AdjustL( color, intensity );
		}

		/// <summary>
		/// "block size", e.g. how may screen pixels are our virtual "pixels"?
		/// </summary>
		public double BlockNumPixels { get; set; }

		private double ScreenPixelSize( CoxeterImages.Settings settings )
		{
			// We go from -bounds to bounds
			return 2 * settings.Bounds / settings.Width;
		}

		private Vector3D RoundGaussian( CoxeterImages.Settings settings, Vector3D v )
		{
			double screenPixelSize = ScreenPixelSize( settings );

			// BlockNumPixels determines # screen pixels in in-radius of a square.
			double quantaSize = screenPixelSize * BlockNumPixels;

			Vector3D result = new Vector3D(
				Quantize( v.X, quantaSize ),
				Quantize( v.Y, quantaSize ) );
			return result;
		}

		/// <summary>
		/// Quantize an input.
		/// </summary>
		private double Quantize( double input, double quantaSize )
		{
			return Math.Round( input / quantaSize ) * quantaSize;
		}

		private Vector3D RoundGaussianDev( CoxeterImages.Settings settings, Vector3D v )
		{
			// sensitive to relative size of xoff
			// I think we want this a bit bigger than the size of a pixel.
			int awayFromZero = 1;
			//digits = (int)(-Math.Log( settings.XOff, 10 ));
			//awayFromZero = 0;

			// hmmm, maybe it is more a function of bounds?

			// Question: Does number of pixels play in?
			// NO! (at least not if we have the squares from rounding bigger than pixels)
			// I made identical 500, 1500, and 6000 pixel pictures when rounding to 0 digits.

			// Round to a fraction 1/n with this.
			// n = 10 and awaFromZero = 0 is the same as n = 1 and awayFromZero = 1.
			// n = 100 and awaFromZero = 0 is the same as n = 1 and awayFromZero = 2.
			double n = 1;   // Let's go slowly from 1 to 100;
			//double n = 30;
			//n = 5;	// good with awayFromZero = 1;

			Vector3D result = new Vector3D( 
				Math.Round( v.X * n, awayFromZero ) / n, 
				Math.Round( v.Y * n, awayFromZero ) / n );
			return result;
		}

		private Vector3D RoundEisenstein( CoxeterImages.Settings settings, Vector3D v )
		{
			double quantaSize = ScreenPixelSize( settings ) * BlockNumPixels;

			double realPart = v.X;
			double imaginaryPart = v.Y;

			// Round to nearest Eisenstein integer
			var (a, b) = RoundToNearestEisenstein( realPart, imaginaryPart, quantaSize );

			var location = ComputeEisenstein( a, b );
			return Vector3D.FromComplex( location );
		}

		// Constants for ω = -1/2 + i√3/2
		// We scale this if we want to change the grid size.
		Complex m_omega = new Complex( -0.5, Math.Sqrt( 3 ) / 2 ); // Primitive cube root of unity

		private Complex ComputeEisenstein( double a, double b )
		{	
			Complex location = (new Complex( a, 0 ) + m_omega * b);
			return location;
		}

		private (double a, double b) RoundToNearestEisenstein( double real, double imag, double quantaSize )
		{
			// Compute the Eisenstein coordinates (a, b)
			double a = real - m_omega.Real * imag / m_omega.Imaginary;
			double b = imag / m_omega.Imaginary;

			/* experimentation
			// Round to nearest... see notes above...
			int awayFromZero = 0;
			double n = 7;
			n = 5;
			//awayFromZero = 2; n = 4;  // hyperbolic (ugh, what did I mean by using this term?)  I will say, it rounds to pixels small enough when bounds are 50, that it reverts back to gaussian.
			double roundedA = Math.Round( a * n, awayFromZero ) / n;
			double roundedB = Math.Round( b * n, awayFromZero ) / n;
			*/

			double roundedA = Quantize( a, quantaSize );
			double roundedB = Quantize( b, quantaSize );

			return (roundedA, roundedB);
		}
	}
}
