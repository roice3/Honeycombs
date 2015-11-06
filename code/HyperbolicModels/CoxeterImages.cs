namespace R3.Geometry
{
	using System.Collections.Generic;
	using System.Drawing;
	using System.Drawing.Imaging;
	using System.IO;
	using System.Linq;
	using System.Numerics;
	using System.Threading.Tasks;
	using R3.Drawing;
	using R3.Math;
	using Math = System.Math;

	public class CoxeterImages
	{
		public class Settings
		{
			public Settings() { Antialias = true; Dual = true; }

			public HoneycombAndView Honeycomb { get; set; }
			public int Width { get; set; }
			public int Height { get; set; }
			public double Bounds { get; set; }			// y bounds (x will be scaled accordingly).
			public Sphere[] Mirrors { get; set; }		// in ball model
			public string FileName { get; set; }
			public bool Antialias { get; set; }
			public double ColorScaling { get; set; }	// Depth to cycle through one hexagon
			public bool Dual { get; set; }

			private double Aspect
			{
				get
				{
					return (double)Width / Height;
				}
			}

			public double XOff
			{
				get
				{
					return Aspect * Bounds * 2 / ( Width - 1 );
				}
			}

			public double YOff
			{
				get
				{
					return Bounds * 2 / (Height - 1);
				}
			}
		}

		Sphere[] m_whiteBoundary;

		public void GenImage( Settings settings )
		{
			int width = settings.Width;
			int height = settings.Height;
			Bitmap image = new Bitmap( width, height );

			bool drawMirrors = false;
			if( drawMirrors )
			{
				DrawMirrors( image, settings );
				image.Save( settings.FileName, ImageFormat.Png );
				return;
			}

			// Setup the boundary which will determine coloring.
			m_whiteBoundary = WhiteBoundary( settings );

			// Cycle through all the pixels and calculate the color.
			int row = 0;
			double bounds = settings.Bounds;
			double xoff = settings.XOff;
			double yoff = settings.YOff;
			Parallel.For( 0, width, i =>
			//for( int i=0; i<width; i++ )
			{
				if( row++ % 20 == 0 )
					System.Console.WriteLine( string.Format( "Processing Line {0}", row ) );

				for( int j=0; j<height; j++ )
				{
					double x = -bounds + i * xoff;
					double y = -bounds + j * yoff;

					if( settings.Antialias )
					{
						const int div = 3;
						List<Color> colors = new List<Color>();
						for( int k=0; k<=div; k++ )
						for( int l=0; l<=div; l++ )
						{
							double xa = x - xoff/2 + k * xoff/div;
							double ya = y - yoff/2 + l * yoff/div;
							Vector3D v = new Vector3D( xa, ya );

							v = ApplyTransformation( v );
							v = PlaneModelToBall( v );
							Color color = CalcColor( settings, v );
							colors.Add( color );
						}

						lock( m_lock )
						{
							image.SetPixel( i, j, AvgColor( colors ) );
						}
					}
					else
					{
						lock( m_lock )
						{
							image.SetPixel( i, j, CalcColor( settings,
								PlaneModelToBall( new Vector3D( x, y ) ) ) );
						}
					}
				}
			} );
			
			image.Save( settings.FileName, ImageFormat.Png );

			// Save as high quality jpeg.
			/*ImageCodecInfo[] codecs = ImageCodecInfo.GetImageDecoders();
			ImageCodecInfo jgpEncoder = codecs.First( c => c.FormatID == ImageFormat.Jpeg.Guid );
			EncoderParameter encoderParam = new EncoderParameter( Encoder.Quality, 100L );	// 100 is least compression
			EncoderParameters encoderParams = new EncoderParameters( 1 );
			encoderParams.Param[0] = encoderParam;
			image.Save( settings.FileName, jgpEncoder, encoderParams );*/
		}

		public struct RecursionStat
		{
			public double Percent;
			public int CellFlips;
		}

		public RecursionStat[] AutoCalcScale( Settings settings )
		{
			System.Console.WriteLine( string.Format( "Auto-calculating the color scale, yo." ) );

			int width = settings.Width;
			int height = settings.Height;
			Bitmap image = new Bitmap( width, height );

			// Setup the boundary which will determine coloring.
			m_whiteBoundary = WhiteBoundary( settings );

			List<double> flipsPerPixel = new List<double>();
			int[] binnedFlips = new int[256];

			// Cycle through all pixels.
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
					Vector3D v = new Vector3D( x, y );
					v = ApplyTransformation( v );
					v = PlaneModelToBall( v );

					lock( m_lock )
					{
						// Only track converged flips?
						// Worry that throwing out unconverged could bias.
						int cellFlips = 0, totalFlips = 0;
						if( ReflectToFundamental( settings, ref v, ref cellFlips, ref totalFlips ) )
							flipsPerPixel.Add( cellFlips );
						
						cellFlips = Math.Min( 255, cellFlips );
						image.SetPixel( i, j, Color.FromArgb( 255, cellFlips, cellFlips, cellFlips ) );
						binnedFlips[cellFlips]++;
					}
				}
			} );

			double mean = flipsPerPixel.Average();
			double median = flipsPerPixel.Median();
			//settings.ColorScaling = 115 / Math.Sqrt( mean ); // Calibrated by looking at 733, 373, 337 and iii.
			//settings.ColorScaling = 115 / Math.Sqrt( perc90 );

			System.Console.WriteLine( string.Format( "Mean: {0:G}, Median: {1:G}", mean, median ) );
			System.Console.WriteLine( string.Format( "Color Scale: {0:G}", settings.ColorScaling ) );
			//image.Save( settings.FileName, ImageFormat.Png );

			/*using( StreamWriter sw = File.AppendText( "binnedFlips.txt" ) )
			{
				string line = "H" + settings.Honeycomb.FormatFilename( string.Empty ) + ",";
				line += string.Join( ",", binnedFlips );
				//sw.WriteLine( line );
			}*/

			// Henry and I settled on this.
			// We want 1% of the pixels to go past black, halfway around the hexagon (hence the factor of 2).
			settings.ColorScaling = (int)flipsPerPixel.ElementAtPercentage( 0.99 ) * 2;
			settings.ColorScaling = 13.3 * Math.Sqrt( mean );	// equivalent to old method...  FINALLY SETTLED ON THIS!

			//////////////////// Everything below is just exploratory.
			List<RecursionStat> stats = new List<RecursionStat>();
			for( int i = 0; i < 10; i++ )
			{
				double percent = 0.9 + i * 0.01;
				int cellFlips = (int)flipsPerPixel.ElementAtPercentage( percent );
				//stats.Add( new RecursionStat() { Percent = percent, CellFlips = cellFlips } );
			}
			stats.Add( new RecursionStat() { Percent = 0.9, CellFlips = (int)flipsPerPixel.ElementAtPercentage( 0.9 ) } );
			stats.Add( new RecursionStat() { Percent = 0.95, CellFlips = (int)flipsPerPixel.ElementAtPercentage( 0.95 ) } );
			stats.Add( new RecursionStat() { Percent = 0.99, CellFlips = (int)flipsPerPixel.ElementAtPercentage( 0.99 ) } );
			stats.Add( new RecursionStat() { Percent = 0.995, CellFlips = (int)flipsPerPixel.ElementAtPercentage( 0.995 ) } );
			stats.Add( new RecursionStat() { Percent = 0.999, CellFlips = (int)flipsPerPixel.ElementAtPercentage( 0.999 ) } );
			return stats.ToArray();
		}

		/// <summary>
		/// This allows us to change the model we have on the plane.
		/// We usually want UHS, but for Pov-Ray mapping these images to sphere, we need to have it be an equirectangular projection
		/// NOTE: The bounds should be set to 1.0 for this to work! v.X and v.Y must be in-between -1 and 1. (also, don't rotate simplex mirrors)  
		/// </summary>
		private Vector3D PlaneModelToBall( Vector3D v )
		{
			bool ballMapForPovRay = false;
			if( !ballMapForPovRay )
			{
				// Normal UHS.
				return H3Models.UHSToBall( v );
			}
			else
			{
				// http://mathworld.wolfram.com/EquirectangularProjection.html
				// y is the latitude
				// x is the longitude
				// Assume inputs go from -1 to 1.
				Vector3D spherical = new Vector3D( 1, Math.PI / 2 * ( 1 - v.Y ), v.X * Math.PI );
				return SphericalCoords.SphericalToCartesian( spherical );
			}
		}

		/// <summary>
		/// Using this to move the view around in interesting ways.
		/// </summary>
		private Vector3D ApplyTransformation( Vector3D v )
		{
			bool applyNone = true;
			if( applyNone )
				return v;

			// self-similar scale for 437
			//v*= 4.259171776329806;

			Mobius m1 = new Mobius(), m2 = new Mobius(), m3 = new Mobius();

			m1.Isometry( Geometry.Hyperbolic, 0, new Vector3D( 0, -.75 ) );
			//return m1.Apply( v );

			Sphere unitSphere = new Sphere();
			v = unitSphere.ReflectPoint( v );
			//v /= 3;
			//v.RotateXY( Math.PI );
			//v.RotateXY( Math.PI/2 );
			return v;

			//v.Y = v.Y / Math.Cos( Math.PI / 6 );	// 637 repeatable
			//return v;

			// 464
			/*
			// Also, don't apply rotations during simplex generation.
			m1.UpperHalfPlane();
			m2 = Mobius.Scale( 1.3 );
			m3.Isometry( Geometry.Euclidean, 0, new Vector3D( 1.55, -1.1 ) );
			v = ( m1 * m2 * m3 ).Apply( v );
			return v;*/

			// iii
			m1.Isometry( Geometry.Hyperbolic, 0, new Complex( 0, Math.Sqrt( 2 ) - 1 ) );
			m2.Isometry( Geometry.Euclidean, -Math.PI / 4, 0 );
			m3 = Mobius.Scale( 5 );
			//v = ( m1 * m2 * m3 ).Apply( v );

			// Vertical Line
			/*v = unitSphere.ReflectPoint( v );
			m1.MapPoints( new Vector3D(-1,0), new Vector3D(), new Vector3D( 1, 0 ) );
			m2 = Mobius.Scale( .5 );
			v = (m1*m2).Apply( v ); */

			/*
			m1 = Mobius.Scale( 0.175 );
			v = unitSphere.ReflectPoint( v );
			v = m1.Apply( v );
			 * */

			// Inversion
			//v = unitSphere.ReflectPoint( v );
			//return v;

			/*Mobius m1 = new Mobius(), m2 = new Mobius(), m3 = new Mobius();
			m1.Isometry( Geometry.Spherical, 0, new Complex( 0, 1 ) );
			m2.Isometry( Geometry.Euclidean, 0, new Complex( 0, -1 ) );
			m3 = Mobius.Scale( 0.5 );
			v = (m1 * m3 * m2).Apply( v );*/

			//Mobius m = new Mobius();
			//m.Isometry( Geometry.Hyperbolic, 0, new Complex( -0.88, 0 ) );
			//m.Isometry( Geometry.Hyperbolic, 0, new Complex( 0, Math.Sqrt(2) - 1 ) );
			//m = Mobius.Scale( 0.17 );
			//m.Isometry( Geometry.Spherical, 0, new Complex( 0, 3.0 ) );
			//v = m.Apply( v );

			// 63i, 73i
			m1 = Mobius.Scale( 6.0 ); // Scale {3,i} to unit disk.
			m1 = Mobius.Scale( 1.0 / 0.14062592996431983 ); // 73i	(constant is abs of midpoint of {3,7} tiling, if we want to calc later for other tilings).
			m2.MapPoints( Infinity.InfinityVector, new Vector3D( 1, 0 ), new Vector3D() ); // swap interior/exterior
			m3.UpperHalfPlane();
			v *= 2.9;

			// iii
			/*m1.MapPoints( new Vector3D(), new Vector3D(1,0), new Vector3D( Math.Sqrt( 2 ) - 1, 0 ) );
			m2.Isometry( Geometry.Euclidean, -Math.PI / 4, 0 );
			m3 = Mobius.Scale( 0.75 );*/

			Mobius m = m3 * m2 * m1;
			v = m.Inverse().Apply( v );	// Strange that we have to do inverse here.

			return v;
		}

		private void DrawMirrors( Bitmap image, Settings settings )
		{
			double b = settings.Bounds;
			ImageSpace i = new ImageSpace( settings.Width, settings.Height );
			i.XMin = -b; i.XMax = b;
			i.YMin = -b; i.YMax = b;

			float scale = 16;

			using( Graphics g = Graphics.FromImage( image ) )
			using( Pen p = new Pen( Color.Red, scale*3.0f ) )
			//using( Pen p2 = new Pen( Color.FromArgb( 255, 255, 214, 0 ), 3.0f ) )
			using( Pen p2 = new Pen( Color.Orange, scale*3.0f ) )
			for( int m=0; m<settings.Mirrors.Length; m++ )			
			{
				Sphere s = settings.Mirrors[m];
				Circle c = H3Models.UHS.IdealCircle( s );	// XXX - not correct
				if( c.IsLine )
				{
					DrawUtils.DrawLine( -c.P2*25, c.P2*25, g, i, p );	// XXX - not general.
				}
				else
				{
					Sphere temp = H3Models.BallToUHS( s );
					DrawUtils.DrawCircle( new Circle { Center = temp.Center, Radius = temp.Radius }, g, i, m == 0 ? p2 : p );
				}

				/* // iii
				Circle c = new Circle();
				c.Radius = Math.Sqrt( 2 );
				c.Center = new Vector3D( 1, Math.Sqrt( 2 ) );
				DrawUtils.DrawCircle( c, g, i, p );
				c.Center = new Vector3D( -1, Math.Sqrt( 2 ) );
				DrawUtils.DrawCircle( c, g, i, p );
				c.Center = new Vector3D( Math.Sqrt( 2 ) - 1, 0 );
				c.Radius = 2 - Math.Sqrt( 2 );
				DrawUtils.DrawCircle( c, g, i, p );

				DrawUtils.DrawLine( new Vector3D( -2, 0 ), new Vector3D( 2, 0 ), g, i, p );
				 */
			}
		}

		private readonly object m_lock = new object();

		private Color AvgColor( List<Color> colors )
		{
			int a = (int)colors.Select( c => (double)c.A ).Average();
			int r = (int)colors.Select( c => (double)c.R ).Average();
			int g = (int)colors.Select( c => (double)c.G ).Average();
			int b = (int)colors.Select( c => (double)c.B ).Average();
			return Color.FromArgb( a, r, g, b );
		}

		private Color CalcColor( Settings settings, Vector3D v )
		{
			int cellFlips = 0, totalFlips = 0;
			if( !ReflectToFundamental( settings, ref v, ref cellFlips, ref totalFlips ) )
				return Color.White;

			//return Color.Black;
			return ColorFunc( m_whiteBoundary, v, cellFlips, settings.ColorScaling, false /*0 == totalFlips % 2*/ );
		}

		/// <summary>
		/// Somewhat based on http://commons.wikimedia.org/wiki/User:Tamfang/programs
		/// </summary>
		private bool ReflectToFundamental( Settings settings, ref Vector3D v, ref int cellFlips, ref int totalFlips )
		{
			Vector3D original = v;

			int iterationCount = 0;
			cellFlips = 0;
			totalFlips = 0;
			//int clean = 0;
			while( true && iterationCount < m_maxIterations )
			{
				/* Original way I did it (a la Anton), which ended up with potential for flips to get incremented within a cell when it shouldn't.
				for( int i=0; i<settings.Mirrors.Length; i++ )
				{
					Sphere mirror = settings.Mirrors[i];
					bool outsideFacet = mirror.IsPointInside( v );	// XXX - confusing, clean this up, especially for planes.
					if( outsideFacet )
					{
						v = mirror.ReflectPoint( v );
						//if( i == 0 ) 
							flips += 1;
						clean = 0;
					}
					else
					{
						clean++;
						if( clean >= settings.Mirrors.Length )
							return ColorFunc( m_whiteBoundary, v, flips );
					}
				}*/

				// First get things on the right side of the 3 mirrors which reflect within a cell.
				// We have to do it this way so the flips across cells are calculated appropriately.
				// (Otherwise, the cell boundary mirror can also reflect within a cell, and this avoids that).
				Sphere[] temp = settings.Mirrors.Skip( 1 ).ToArray();
				//Mirror[] temp = new Mirror[] { settings.Mirrors[0], settings.Mirrors[1], settings.Mirrors[3] };	// Experiments with using other mirrors as "cell" boundary.
				if( !ReflectAcrossMirrors( temp, ref v, ref totalFlips, ref iterationCount ) )
					continue;

				if( !ReflectAcrossMirror( settings.Mirrors[0], ref v ) )
				{
					cellFlips++;
					totalFlips++;
				}
				else
				{
					return true;
				}

				//iterationCount++;
			}

			System.Console.WriteLine( string.Format( "Did not converge at point {0}", original.ToString() ) );
			return false;
		}

		private int m_maxIterations = 4000;

		private bool ReflectAcrossMirror( Sphere mirror, ref Vector3D v )
		{
			bool outsideFacet = mirror.IsPointInside( v );
			if( outsideFacet )
			{
				v = mirror.ReflectPoint( v );
				return false;
			}

			return true;
		}

		private bool ReflectAcrossMirrors( Sphere[] mirrors, ref Vector3D v, ref int flips, ref int iterationCount )
		{
			int clean = 0;
			while( true && iterationCount < m_maxIterations )
			{
				for( int i=0; i<mirrors.Length; i++ )
				{
					Sphere mirror = mirrors[i];
					if( !ReflectAcrossMirror( mirror, ref v ) )
					{
						flips++;
						clean = 0;
					}
					else
					{
						clean++;
						if( clean >= mirrors.Length )
							return true;
					}
				}
				iterationCount++;
			}

			return false;
		}

		/// <summary>
		/// Grabs the portion of the fundamental region we want to remain white.
		/// </summary>
		private static Sphere[] WhiteBoundary( Settings settings )
		{
			// We only need to do this if we have hyperideal vertices.
			//if( Geometry2D.GetGeometry( settings.Honeycomb.Q, settings.Honeycomb.R ) != Geometry.Hyperbolic )
			//	return settings.Mirrors;

			Sphere[] fundamentalRegion = settings.Mirrors;

			double offset = 0.008;	// XXX - put in settings. // XXX - should be hyperbolic distance (doesn't scale the same for ball/UHS)
			//double offset = 0.02;	// i3i pictures

			List<Sphere> result = new List<Sphere>();
			for( int i=0; i<fundamentalRegion.Length; i++ )
			{
				Sphere s = fundamentalRegion[i];
				Sphere facetSphere;
				if( i == 0 )
				//if( i == 2 )	// Dual
				//if( true )  XXX - try with all and see how it looks!
				{
					//facetSphere = GeodesicOffset( s, offset );
					facetSphere = AlteredFacetForTrueApparent2DTilings( fundamentalRegion );
				}
				else
					facetSphere = s;

				result.Add( facetSphere );
			}

			return result.ToArray();
		}

		/// <summary>
		/// This will return an altered facet to create true apparent 2D tilings (proper bananas) on the boundary.
		/// Notes:
		///		The input simplex must be in the ball model.
		///		The first mirror of the simplex (the one that mirrors across cells) is the one we end up altering.
		/// </summary>
		public static Sphere AlteredFacetForTrueApparent2DTilings( Sphere[] simplex )
		{
			// We first need to find the size of the apparent 2D disk surrounding the leg.
			// This is also the size of the apparent cell head disk of the dual.
			// So we want to get the midsphere (insphere would also work) of the dual cell with that head, 
			// then calculate the intersection of that with the boundary.
			Sphere cellMirror = simplex[0];
			if( cellMirror.IsPlane )
				throw new System.NotImplementedException();

			// The point centered on a face is the closest point of the cell mirror to the origin.
			// This will be the point centered on an edge on the dual cell.
			Vector3D facePoint = cellMirror.ProjectToSurface( new Vector3D() );

			// Reflect it to get 3 more points on our midsphere.
			Vector3D reflected1 = simplex[1].ReflectPoint( facePoint );
			Vector3D reflected2 = simplex[2].ReflectPoint( reflected1 );
			Vector3D reflected3 = simplex[0].ReflectPoint( reflected2 );
			Sphere midSphere = Sphere.From4Points( facePoint, reflected1, reflected2, reflected3 );

			// Get the ideal circles of the cell mirror and midsphere.
			// Note: The midsphere is not geodesic, so we can't calculate it the same.
			Sphere cellMirrorUHS = H3Models.BallToUHS( cellMirror );
			Circle cellMirrorIdeal = H3Models.UHS.IdealCircle( cellMirrorUHS );

			Sphere ball = new Sphere();
			Circle3D midSphereIdealBall = ball.Intersection( midSphere ); // This should exist because we've filtered for honeycombs with hyperideal verts.
			Circle3D midSphereIdeal3D = H3Models.BallToUHS( midSphereIdealBall );
			Circle midSphereIdeal = new Circle { Center = midSphereIdeal3D.Center, Radius = midSphereIdeal3D.Radius };

			// The intersection point of our cell mirror and the disk of the apparent 2D tiling
			// gives us "ideal" points on the apparent 2D boundary. These points will be on the new cell mirror.
			Vector3D i1, i2;
			if( 2 != Euclidean2D.IntersectionCircleCircle( cellMirrorIdeal, midSphereIdeal, out i1, out i2 ) )
				throw new System.ArgumentException( "Since we have hyperideal vertices, we should have an intersection with 2 points." );

			// Transform the intersection points to a standard Poincare disk.
			// The midsphere radius is the scale of the apparent 2D tilings.
			double scale = midSphereIdeal.Radius;
			Vector3D offset = midSphereIdeal.Center;
			i1 -= offset;
			i2 -= offset;
			i1 /= scale;
			i2 /= scale;
			Circle3D banana = H3Models.Ball.OrthogonalCircle( i1, i2 );
			Vector3D i3 = H3Models.Ball.ClosestToOrigin( banana );
			i3 = Hyperbolic2D.Offset( i3, -.025 );

			// Transform back.
			i1 *= scale; i2 *= scale; i3 *= scale;
			i1 += offset; i2 += offset; i3 += offset;

			// Construct our new simplex mirror with these 3 points.
			Circle3D c = new Circle3D( i1, i2, i3 );

			Sphere result = new Sphere() { Center = c.Center, Radius = c.Radius };
			return H3Models.UHSToBall( result );
		}

		public static Sphere GeodesicOffset( Sphere s, double offset, bool ball = true )
		{
			Sphere offsetSphere;
			if( ball )
			{
				// Geodesic offset (ball).

				{	// Hyperbolic honeycomb
					double mag = s.Center.Abs() - s.Radius;
					mag = s.IsPlane ? DonHatch.h2eNorm( offset ) :
						DonHatch.h2eNorm( DonHatch.e2hNorm( mag ) - offset );

					Vector3D closestPointToOrigin = s.IsPlane ? s.Normal : s.Center;
					closestPointToOrigin.Normalize();
					closestPointToOrigin *= mag;
					offsetSphere = H3Models.Ball.OrthogonalSphereInterior( closestPointToOrigin );

					// There are multiple ultraparallel spheres.
					// This experiments with picking others.
					Mobius m = new Mobius();
					m.Isometry( Geometry.Hyperbolic, 0, new Vector3D( 0, -0.2 ) );
					//H3Models.TransformInBall2( offsetSphere, m );
				}

				{	// Spherical honeycomb
					//offset *= -1;
					double mag = -s.Center.Abs() + s.Radius;
					Spherical2D.s2eNorm( Spherical2D.e2sNorm( mag ) + offset );

					offsetSphere = s.Clone();
					offsetSphere.Radius += offset*10;
				}
			}
			else
			{
				// Geodesic offset (UHS).
				// XXX - not scaled right.
				offsetSphere = s.Clone();
				offsetSphere.Radius += offset;
			}

			return offsetSphere;
		}

		/// <summary>
		/// This is a 2D function for now.
		/// Given an input geodesic in the plane, returns an equidistant circle.
		/// Offset would be the offset at the origin.
		/// </summary>
		public static Circle3D EquidistantOffset( Segment seg, double offset )
		{
			// Get the ideal endpoints of the input.
			Vector3D i1, i2;
			H3Models.Ball.GeodesicIdealEndpoints( seg.P1, seg.P2, out i1, out i2 );

			// Get the offset amount at the location.
			Vector3D center;
			double radius;
			H3Models.Ball.OrthogonalCircle( i1, i2, out center, out radius );
			Vector3D closest = center;
			closest.Normalize();
			closest *= (center.Abs() - radius);
			H3Models.Ball.DupinCyclideSphere( closest, offset, out center, out radius );
			double mag = center.Abs() - radius;
			closest.Normalize();
			closest *= mag;

			// Resulting circle will go through i1 and i2, but will be offset slightly towards the origin.
			Circle3D result = new Circle3D( i1, closest, i2 );
			return result;
		}

		/// <summary>
		/// Same as above, but works with all geometries.
		/// </summary>
		public static Circle EquidistantOffset( Geometry g, Segment seg, double offset )
		{
			Mobius m = new Mobius();
			Vector3D direction;
			if( seg.Type == SegmentType.Line )
			{
				direction = seg.P2 - seg.P1;
				direction.RotateXY( Math.PI / 2 );
			}
			else
			{
				direction = seg.Circle.Center;
			}

			direction.Normalize();
			m.Isometry( g, 0, direction * offset );

			// Transform 3 points on segment.
			Vector3D p1 = m.Apply( seg.P1 );
			Vector3D p2 = m.Apply( seg.Midpoint );
			Vector3D p3 = m.Apply( seg.P2 );

			return new Circle( p1, p2, p3 );
		}

		private static Color ColorFunc( Sphere[] boundary, Vector3D v, int reflections, double colorScaling, bool invert = false )
		{
			foreach( Sphere facet in boundary )
			{
				bool outsideFacet = facet.IsPointInside( v );
					//&& H3Models.BallToUHS( v ).Abs() <= 1;	// To remove "black dots" (I disagree with doing this myself, but Henry wants it.)
				if( outsideFacet )
					return Color.Black;
					//return Color.White;
			}

			//return DepthColor( reflections, colorScaling, invert );
			return Coloring.ColorAlongHexagon( (int)colorScaling, reflections );
		}
	}
}
