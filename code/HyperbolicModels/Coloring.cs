namespace R3.Drawing
{
	using System;
	using System.Collections.Generic;
	using System.Drawing;
	using System.Linq;
	using R3.Geometry;

	internal static class Coloring
	{
		/// <summary>
		/// This will calculate the color along a hexagon on the edges of an RGB cube.
		/// incrementsUntilRepeat is the value where we return to the starting point of the hexagon (white).
		/// increments is used as the distance-along-hexagon parameter.
		/// </summary>
		public static Color ColorAlongHexagon( int incrementsUntilRepeat, double increments )
		{
			//if( 0 == increments )
			//	return Color.FromArgb( 255, 187, 23, 23 );
				//return Color.FromArgb( 0, 255, 255, 255 );

			//int temp = (increments - 2) * 125 + 80;
			/*int temp = increments * 40;
			if( temp < 0 )
				temp = 0;
			if( temp > 255 )
				temp = 255;
			Color gray = Color.FromArgb( 255, temp, temp, temp );
			Color c = increments > 2 ? Color.White : gray;
			return c;*/

			//464
			//increments = (int)( Math.Pow( (double)increments, 1.35 ) );

			// Bring to main hexagon (handle looping)
			increments = increments % incrementsUntilRepeat;

			// 0 to 6, so we can have each edge of the hexagon live in a unit interval.
			double distAlongHex = increments * 6 / incrementsUntilRepeat;

			// 464
			// Give the first 
			//double percentage = (double)increments / incrementsUntilRepeat;
			// Sigmoid!  http://en.wikipedia.org/wiki/Sigmoid_function
			//distAlongHex = ( 1 / ( 1 + Math.Exp( -percentage ) ) ) ) * 6;
			//distAlongHex = percentage * 6;

			Func<double, int> subtractive = d => (int)( 255.0 * ( 1.0 - d ) );
			Func<double, int> addative = d => (int)( 255.0 * d );

			bool blue = true;
			if( blue )
			{
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, subtractive( distAlongHex ), 255, 255 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 0, subtractive( distAlongHex ), 255 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 0, 0, subtractive( distAlongHex ) );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, addative( distAlongHex ), 0, 0 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 255, addative( distAlongHex ), 0 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 255, 255, addative( distAlongHex ) );
			}
			else
			{
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 255, 255, subtractive( distAlongHex ) );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 255, subtractive( distAlongHex ), 0 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, subtractive( distAlongHex ), 0, 0 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 0, 0, addative( distAlongHex ) );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, 0, addative( distAlongHex ), 255 );
				distAlongHex--;
				if( distAlongHex < 1 )
					return Color.FromArgb( 255, addative( distAlongHex ), 255, 255 );
			}

			throw new System.Exception( "Bad impl" );
		}
											 
		/// <summary>
		/// This method allows calculating a color based on a color scaling input.
		/// This was the way I did this for a long time and for many images, until working with Henry on coloring.
		/// Keeping this code here for reference.
		/// </summary>
		public static Color DepthColor( int depth, double colorScaling, bool invert = false )
		{
			int scaling = 50;	// 50 good. // XXX - make a setting.
			scaling = 150;
			//scaling = 130;
			//scaling = 60;
			scaling = (int)colorScaling;

			/* tuned 733
			int mag = 0;
			switch( depth )
			{
				case 1: mag = 150; break;
				case 2: mag = 350; break;
				case 3: mag = 465; break;
				case 4: mag = 400; break;
				case 5: mag = 200; break;
				case 6: mag = 100; break;
			} */

			int wraps = 0;
			int mag = depth * scaling;
			//mag = (int)( System.Math.Pow( depth, colorScaling ) * scaling );
			int c1 = mag, c2 = 0, c3 = 0;
			while( mag > 0 )	// Comment this line out for no color wrapping.
			{
				c1 = mag; c2 = c3 = 0;
				if( c1 > 255 )
				{
					c1 = 255;
					mag -= 255;
					c2 = mag;
				}
				if( c2 > 255 )
				{
					c2 = 255;
					mag -= 255;
					c3 = mag;
				}
				if( c3 > 255 )
				{
					c3 = 255;
					wraps++;
				}

				mag -= 255;
			}

			System.Func<int, int> min75 = i => 255 - i < 75 ? 75 : 255 - i;
			System.Func<int, int> max180 = i => i > 180 ? 180 : i;

			Color c = Color.White;

			//Color c = Color.FromArgb( 255, min75( c3 ), 255 - c1, 255 - c2 );	// purple red
			//Color c = Color.FromArgb( 255, min75( c3 ), 255 - c2, 255 - c1 );	// *yellow->red

			//Color c = Color.FromArgb( 255, 255 - c1, min75( c3 ), 255 - c2 );	// blue then green
			//Color c = Color.FromArgb( 255, 255 - c2, min75( c3 ), 255 - c1 ); // green->yellow

			int modulo = wraps % 2;
			switch( modulo )
			{
				case 0:
					c = Color.FromArgb( 255, 255 - c1, 255 - c2, min75( c3 ) );	// *blue
					break;
				case 1:
					c = Color.FromArgb( 255, c1, c2, max180( c3 ) );
					break;
				case 2:
					c = Color.FromArgb( 255, 255 - c2, min75( c3 ), 255 - c1 ); // green->yellow
					break;
				case 3:
					c = Color.FromArgb( 255, min75( c3 ), 255 - c2, 255 - c1 );	// *yellow->red
					break;
			}

			// c = Color.FromArgb( 255, 255 - c1, 255 - c2, min75( c3 ) );	// *blue
			//Color c = Color.FromArgb( 255, 255 - c1, 255 - c1, min75( c1 ) );	// *gray with tint of blue
			//Color c = Color.FromArgb( 255, min75( c1 ), 255 - c1, 255 - c1 );

			//Color c = Color.FromArgb( 255, 255 - c1, 255 - c2, 255 - c3 );
			//c = depth >= 4 ? Color.White : c;

			/*
			int temp = (depth-3) * 40;
			if( temp < 0 ) temp = 0;
			if( temp > 255 ) temp = 255;
			Color gray = Color.FromArgb( 255, temp, temp, temp );
			c = depth < 3 ? Color.Black : gray;
			*/

			if( invert )
				return Inverse( c );

			return c;
		}

		public static Color Inverse( Color c )
		{
			return Color.FromArgb( 255, 255 - c.R, 255 - c.G, 255 - c.B );
		}

		public static Vector3D ToVec( Color c )
		{
			return new Vector3D( (double)c.R / 255, (double)c.G / 255, (double)c.B / 255 );
		}

		public static Color AvgColor( List<Color> colors )
		{
			//if( colors.Contains( Color.White ) )
			//	return Color.White;

			int a = (int)colors.Select( c => (double)c.A ).Average();
			int r = (int)colors.Select( c => (double)c.R ).Average();
			int g = (int)colors.Select( c => (double)c.G ).Average();
			int b = (int)colors.Select( c => (double)c.B ).Average();
			return Color.FromArgb( a, r, g, b );
		}
	}
}
