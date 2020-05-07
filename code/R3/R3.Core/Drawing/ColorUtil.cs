namespace R3.Core
{
	using R3.Geometry;
	using System;
	using System.Collections.Generic;
	using System.Drawing;
	using System.Linq;
	using Math = System.Math;

	public static class ColorUtil
	{
		// Takes Hue value as input, returns RGB vector.
		// Copied from POV-Ray
		public static Vector3D CH2RGB( double H )
		{
			double R = 0, G = 0, B = 0;
			if( H >= 0 && H < 120 )
			{
				R = (120 - H) / 60;
				G = (H - 0) / 60;
				B = 0;
			}
			else if( H >= 120 && H < 240 )
			{
				R = 0;
				G = (240 - H) / 60;
				B = (H - 120) / 60;
			}
			else if( H >= 240 && H <= 360 )
			{
				R = (H - 240) / 60;
				G = 0;
				B = (360 - H) / 60;
			}

			return new Vector3D(
				Math.Min( R, 1 ),
				Math.Min( G, 1 ),
				Math.Min( B, 1 ) );
		}

		// Copied from POV-Ray
		// Putting this here for speed. It was too expensive to do this at render time in POV-Ray.
		public static Vector3D CHSL2RGB( Vector3D hsl )
		{
			Vector3D ones = new Vector3D( 1, 1, 1 );

			double H = hsl.X;
			double S = hsl.Y;
			double L = hsl.Z;
			Vector3D SatRGB = CH2RGB( H );
			Vector3D Col = 2 * S * SatRGB + (1 - S) * ones;
			Vector3D rgb;
			if( L < 0.5 )
				rgb = L * Col;
			else
				rgb = (1 - L) * Col + (2 * L - 1) * ones;

			return rgb;
		}

		public static Color AvgColor( List<Color> colors )
		{
			if( colors.Count == 0 )
				return Color.FromArgb( 0, 0, 0, 0 );

			int a = (int)colors.Select( c => (double)c.A ).Average();
			int r = (int)colors.Select( c => (double)c.R ).Average();
			int g = (int)colors.Select( c => (double)c.G ).Average();
			int b = (int)colors.Select( c => (double)c.B ).Average();
			return Color.FromArgb( a, r, g, b );
		}

		public static Color AvgColorSquare( List<Color> colors )
		{
			if( colors.Count == 0 )
				return Color.FromArgb( 0, 0, 0, 0 );

			int a = (int)Math.Sqrt( colors.Select( c => (double)c.A * c.A ).Average() );
			int r = (int)Math.Sqrt( colors.Select( c => (double)c.R * c.R ).Average() );
			int g = (int)Math.Sqrt( colors.Select( c => (double)c.G * c.G ).Average() );
			int b = (int)Math.Sqrt( colors.Select( c => (double)c.B * c.B ).Average() );
			return Color.FromArgb( a, r, g, b );
		}

		public static Color InterpColor( Color c1, Color c2, double input )
		{
			System.Func<int, int, double, int> interp = ( i1, i2, d ) =>
			{
				return (int)( (double)i1 + d * (double)( i2 - i1 ) );
			};

			int a = interp( c1.A, c2.A, input );
			int r = interp( c1.R, c2.R, input );
			int g = interp( c1.G, c2.G, input );
			int b = interp( c1.B, c2.B, input );
			return Color.FromArgb( a, r, g, b );
		}

		public static Color Inverse( Color c )
		{
			return Color.FromArgb( 255, 255 - c.R, 255 - c.G, 255 - c.B );
		}

		public static Color FromRGB( Vector3D rgb )
		{
			if( rgb.DNE )
				return Color.FromArgb( 0, 255, 255, 255 );

			rgb *= 255;
			return Color.FromArgb( 255, (int)rgb.X, (int)rgb.Y, (int)rgb.Z );
		}

		public static Color AdjustH( Color c, double h )
		{
			Vector3D hsl = new Vector3D( c.GetHue(), c.GetSaturation(), c.GetBrightness() );
			hsl.X = h;
			Vector3D rgb = CHSL2RGB( hsl );
			return FromRGB( rgb );
		}

		public static Color AdjustS( Color c, double s )
		{
			Vector3D hsl = new Vector3D( c.GetHue(), c.GetSaturation(), c.GetBrightness() );
			hsl.Y = s;
			Vector3D rgb = CHSL2RGB( hsl );
			return FromRGB( rgb );
		}

		public static Color AdjustL( Color c, double l )
		{
			if( l > 1 )
				l = 1;
			if( l < 0 )
				l = 0;

			Vector3D hsl = new Vector3D( c.GetHue(), c.GetSaturation(), c.GetBrightness() );
			hsl.Z = l;
			Vector3D rgb = CHSL2RGB( hsl );
			return FromRGB( rgb );
		}

		/// <summary>
		/// This will calculate the color along a hexagon on the edges of an RGB cube.
		/// incrementsUntilRepeat is the value where we return to the starting point of the hexagon (white).
		/// increments is used as the distance-along-hexagon parameter.
		/// </summary>
		public static Color ColorAlongHexagon( int incrementsUntilRepeat, int increments )
		{
			// Bring to main hexagon (handle looping)
			increments += (int)(.0 * incrementsUntilRepeat);    // an offset along the color hexagon
			increments = increments % incrementsUntilRepeat;

			// 0 to 6, so we can have each edge of the hexagon live in a unit interval.
			double distAlongHex = (double)increments * 6 / incrementsUntilRepeat;

			Func<double, int> subtractive = d => (int)(255.0 * (1.0 - d));
			Func<double, int> addative = d => (int)(255.0 * d);

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
	}
}
