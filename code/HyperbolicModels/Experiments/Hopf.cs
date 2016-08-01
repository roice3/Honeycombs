namespace HyperbolicModels
{
	using R3.Core;
	using R3.Drawing;
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Drawing;
	using System.IO;
	using System.Linq;
	using System.Numerics;

	using Math = System.Math;

	/// <summary>
	/// Reference for hyperbolic "fiberwise homogenous" fibrations.
	/// https://www.math.upenn.edu/grad/dissertations/NuchiThesis.pdf
	/// 
	/// Thoughts:
	///		Make meta-fibers that are collections of fibers.
	/// </summary>
	public class Hopf
	{
		public bool Ball { get { return true;  } }

		public void GenPovRay( )
		{
			//H3.Cell.Edge[] fibers = ParameterizedFibers( zParam );

			List<H3.Cell.Edge> fiberList = new List<H3.Cell.Edge>();

			bool drawAll = false;
			if( drawAll )
			{
				Vector3D[] spherePoints = SpherePoints();
				foreach( Vector3D v in spherePoints )
					fiberList.Add( CalcFiber( v ) );
				WriteFibers( fiberList, "fibers.pov" );
			}

			string baseDir = @"C:\Users\hrn\Documents\roice\povray\hopf\pulled\";

			bool drawPulledBack = true;
			if( drawPulledBack )
			{
				int frames = 1;
				for( int i = 0; i < frames; i++ )
				{
					double t = 1.0 - (double)i / frames;
					t = Smoothed( t, 1.0 );	// from 0 to 1, smoothed.

					CalcFibers( fiberList, t );

					string num = i.ToString();
					num = num.PadLeft( 3, '0' );
					string filename = Path.Combine( baseDir, string.Format( "pulled{0}.pov", num ) );
					PrependPovFile( filename );
					WriteFibers( fiberList, filename );

					AppendMesh( fiberList, filename );

					/* // Lift of a {3,7} tiling.
					var allPoints = BasePointsTiling();
					foreach( Vector3D[] array in allPoints )
					{
						fiberList.Clear();
						CalcFibers( array, fiberList, t: 1.0 );
						AppendMesh( fiberList, filename );
					}*/
				}
			}
		}

		private void AppendMesh( List<H3.Cell.Edge> fiberList, string filename )
		{
			Mesh mesh = new Mesh();
			for( int j = 0; j < fiberList.Count - 1; j++ )
				AddRow( mesh, fiberList[j], fiberList[j + 1] );
			//AddRow( mesh, fiberList.Last(), fiberList.First() );
			PovRay.WriteMesh( mesh, filename, append: true );
		}

/////////////////////////////////////////////////////// Methods involved in sequence parameterization.

		private static Complex ZParam
		{
			get
			{
				//return new Complex( 5, 4 );
				return new Complex( 2, 1 );	// 45 degree angle and symmetric about real line, seems most canonical of the twisted fibrations.
				//return new Complex( 0, 1 );
			}
		}

		void CalcFibers( List<H3.Cell.Edge> fibers, double t )
		{
			fibers.Clear();

			//t = 0.5;
			Vector3D[] basePoints = BasePointsCircle( t );
			CalcFibers( basePoints, fibers, t );
		}

		void CalcFibers( Vector3D[] basePoints, List<H3.Cell.Edge> fibers, double t )
		{
			foreach( Vector3D v in basePoints )
			{
				H3.Cell.Edge e = CalcFiberFromBase( v );
				if( e != null )
					e.Color = Color( Transform( H3Models.UHSToBall( v ), t ) );
				fibers.Add( Transform( e, t ) );
			}
		}

		H3.Cell.Edge Transform( H3.Cell.Edge edge, double t )
		{
			if( edge == null )
				return edge;

			H3.Cell.Edge result = new H3.Cell.Edge( Transform( edge.Start, t ), Transform( edge.End, t ) );
			result.Color = edge.Color;
			return result;
		}

/////////////////////////////////////////////////////// 

		void PrependPovFile( string filename )
		{
			using( StreamWriter sw = File.CreateText( filename ) )
			{
				if( Ball )
					sw.WriteLine( "#include \"" + @"C:\Users\hrn\Documents\roice\povray\hopf\hopf.pov" + "\"" );
				else
					sw.WriteLine( "#include \"" + @"C:\Users\hrn\Documents\roice\povray\hopf\hopfUHS.pov" + "\"" );
			}
		}

		void AddRow( Mesh mesh, H3.Cell.Edge edge1, H3.Cell.Edge edge2 )
		{
			int div = 60;	// Has to be the same for all edges?

			if( edge1 == null && edge2 == null )
				return;

			System.Func<Vector3D, Vector3D, int, Vector3D[]> geoPoints = null;
			if( Ball )
				geoPoints = H3Models.Ball.GeodesicPoints;
			else
				geoPoints = H3Models.UHS.GeodesicPoints;

			if( edge1 == null || edge2 == null )
			{
				H3.Cell.Edge e = edge1 == null ? edge2 : edge1;
				if( e.Start == e.End )
					return;

				Vector3D[] verts = geoPoints( e.Start, e.End, div );
				Vector3D center, dummy;
				double radius, dummyD;
				if( Ball )
					H3Models.Ball.Geodesic( e.Start, e.End, out center, out radius, out dummy, out dummyD );
				else
					H3Models.UHS.Geodesic( e.Start, e.End, out center, out radius, out dummy, out dummyD );
				for( int i = 0; i < verts.Length - 1; i++ )
				{
					mesh.Triangles.Add( new Mesh.Triangle( verts[i], verts[i + 1], center ) );	// Radius will be small, so reasonable.
				}

				return;
			}

			Vector3D[] verts1 = geoPoints( edge1.Start, edge1.End, div );
			Vector3D[] verts2 = geoPoints( edge2.Start, edge2.End, div );
			for( int i=0; i<verts1.Length-1; i++ )
			{
				mesh.Triangles.Add( new Mesh.Triangle( verts1[i], verts1[i + 1], verts2[i] ) );
				mesh.Triangles.Add( new Mesh.Triangle( verts2[i], verts1[i + 1], verts2[i+1] ) );
			}
		}

		/// <summary>
		/// XXX - Need to make this a function of clock, so it probably needs to not be static.
		/// </summary>
		/// <param name="basePointUHS"></param>
		/// <returns>A base point in the disk, possibly transformed.</returns>
		Vector3D Transform( Vector3D basePointUHS, double t )
		{
			if( !Ball )
				return basePointUHS;

			Vector3D basePointDisk = H3Models.UHSToBall( basePointUHS );

			// Apply transformations.
			// XXX - make this a mobius we pass in, so we can edit it on the outside.
			//basePointDisk.RotateAboutAxis( new Vector3D( 0, 1, 0 ), 2 * Math.PI * t );
			
			return basePointDisk;
		}

		/// <summary>
		/// Given a disk point, return an associated color with HSL components.
		/// This is a mapping from H^2 to an HSL color wheel (with L = 0.5).
		/// </summary>
		static Vector3D Color( Vector3D basePointDisk )
		{
			// Y component of input will be 0 (based on how we are choosing our base surface).
			Complex c = new Complex( basePointDisk.X, basePointDisk.Z );

			double phase = c.Phase / ( 2*Math.PI );
			if( phase < 0 )
				phase += 1;

			Vector3D color = new Vector3D(
				phase * 360,				// POV-Ray CHSL2RGB macro is messed up.  Expects hue from 0 to 360, but other entries from 0 to 1.
				c.Magnitude,
				0.5
			);
			Trace.WriteLine( color.ToStringXYZOnly() );
			return color;
		}

		static double Smoothed( double input, double max )
		{
			return (max / 2.0) * (-Math.Cos( Math.PI * input / max ) + 1);
		}

		void WriteFibers( List<H3.Cell.Edge> fiberList, string fileName )
		{
			H3.Cell.Edge[] fibers = fiberList.Where( f => f != null ).ToArray();

			if( Ball )
				PovRay.WriteH3Edges( new PovRay.Parameters { AngularThickness = 0.020 }, fibers, fileName, append: true );
			else
				PovRay.WriteH3Edges( new PovRay.Parameters { AngularThickness = 0.010, Halfspace = true }, fibers, fileName, append: true );
		}

		private static List<Vector3D[]> BasePointsTiling()
		{
			TilingConfig config = new TilingConfig( 3, 7, maxTiles: 100 );
			Tiling tiling = new Tiling();
			tiling.Generate( config );

			HashSet<H3.Cell.Edge> finished = new HashSet<H3.Cell.Edge>( new H3.Cell.EdgeEqualityComparer() );

			int numPerSeg = 25;
			List<Vector3D[]> basePoints = new List<Vector3D[]>();
			foreach( Tile t in tiling.Tiles )
			foreach( Segment s in t.Boundary.Segments )
			{
				H3.Cell.Edge e = new H3.Cell.Edge( s.P1, s.P2 );
				if( finished.Contains( e ) )
					continue;
				finished.Add( e );

				Vector3D[] points = s.Subdivide( numPerSeg ).Select( p => 
				{
					p = new Vector3D( p.X, 0, p.Y );
					return H3Models.BallToUHS( p );
				} ).ToArray();
				basePoints.Add( points );
			}

			return basePoints;
		}

		/// <summary>
		/// Calculate basePoints along a circle defined by a center on the z axis, and going through the point (0,0,1).
		/// .5 is horosphere, 1 is geodesic
		/// Result is in UHS model
		/// </summary>
		private static Vector3D[] BasePointsCircle( double centerUHS )
		{
			if( centerUHS > 1 )
				centerUHS = 1;
			if( centerUHS < 0 )
				centerUHS = 0;

			if( centerUHS == 0 )
				return new Vector3D[] { new Vector3D() };

			bool hyperbolicOffsets = true;
			if( hyperbolicOffsets )
			{
				// h-distance between each point.
				double d = 0.1;

				// We need to work in 2D first, then we'll switch to xz plane.
				Circle circle = new Circle( new Vector3D( 0, 1 ), new Vector3D( 0, 1 - 2 * centerUHS ), new Vector3D( centerUHS, 1 - centerUHS ) );

				// XXX - check to make sure total distance won't wrap around the circle.

				List<Vector3D> points = new List<Vector3D>();
				int count = 75;
				for( int i=-count; i<=count; i++ )
				{
					double currentD = d * i;

					// Angle t around a circle with center c and radius r to get an h-distance.
					// https://en.wikipedia.org/wiki/Poincar%C3%A9_half-plane_model
					// http://www.wolframalpha.com/input/?i=d+%3D+arcosh%281%2B%28%28r*sin%28t%29%29%5E2%2B%28r*cos%28t%29-r%29%5E2%29%2F%282*%28c%2Br%29*%28c%2Br*cos%28t%29%29%29%29%2C+solve+for+t
					double c = circle.Center.Y;
					double r = circle.Radius;
					double coshd = DonHatch.cosh( currentD );
					double numerator = c * c - c * (c + r) * coshd + c * r + r * r;
					double denominator = r * ((c + r) * coshd - c);
					double angle = Math.Acos( numerator / denominator );
					if( i < 0 )
						angle *= -1;
					points.Add( new Vector3D( r * Math.Sin( angle ), 0, c + r * Math.Cos( angle ) ) );

					/*
					// XXX - This was my first attempt, but this code only works for geodesics, not general arcs!
					// Equidistant lines in UHS will all be lines through the origin.
					// In the following formula, x is the angle to use to get h-spaced equidistant line a distance d away.
					// http://www.wolframalpha.com/input/?i=d+%3D+arccosh%28sec%28x%29%29%2C+solve+for+x
					double angle = Math.Acos( 1.0 / DonHatch.cosh( currentD ) );
					angle = Math.PI/2 - angle;
					if( i < 0 )
						angle *= -1;

					Vector3D p1, p2;
					Euclidean2D.IntersectionLineCircle( new Vector3D(), new Vector3D( Math.Cos(angle), Math.Sin(angle) ), circle, out p1, out p2 );
					Vector3D highest = p1.Y > p2.Y ? p1 : p2;
					points.Add( new Vector3D( highest.X, 0, highest.Y ) );
					*/
				}
				return points.ToArray();
			}
			else
			{
				// equal euclidean spacing.
				Circle3D c = new Circle3D( new Vector3D( 0, 0, 1 ), new Vector3D( 0, 0, 1 - 2 * centerUHS ), new Vector3D( centerUHS, 0, 1 - centerUHS ) );
				return c.Subdivide( 125 );
			}
		}

		private static int NumFibers
		{
			get
			{
				return 1000;
			}
		}

		/// <summary>
		/// Get a distribution of points on a sphere.
		/// The points are the vertices of a geodesic dome.
		/// </summary>
		private static Vector3D[] SpherePoints()
		{
			List<Vector3D> spherePoints = new List<Vector3D>();
			TilingConfig config = new TilingConfig( 3, 5 );
			Tiling tiling = new Tiling();
			tiling.GenerateInternal( config );

			Tile baseTile = tiling.Tiles.First();
			Vector3D[] templateTextureCoords = TextureHelper.TextureCoords( baseTile.Boundary, Geometry.Spherical, doGeodesicDome: true );
			foreach( Tile tile in tiling.Tiles )
			{
				Isometry isom = new Isometry();
				isom.CalculateFromTwoPolygons( baseTile, tile, Geometry.Spherical );
				Vector3D[] textureCoords = Isometry.TransformVertices( templateTextureCoords, isom.Inverse() );
				spherePoints.AddRange( textureCoords );
			}

			return spherePoints.Select( p => H3Models.UHSToBall( p ) ).Distinct().ToArray();
		}

		private static H3.Cell.Edge CalcFiberFromBase( Vector3D basePointUHS )
		{
			if( Tolerance.LessThan( basePointUHS.Z, 0 ) )
				return null;

			// parallel case not implemented.

			// Standard fiber.
			Vector3D z = Vector3D.FromComplex( ZParam );
			Vector3D minusi = new Vector3D( 0, -1 );

			// Its intersection with base space (plane above the real line).
			Vector3D center = (z + minusi) / 2;
			Vector3D intersection = minusi + (z - minusi) * (0 - minusi.Y) / (z.Y - minusi.Y);
			double rad = (center - minusi).Abs();
			intersection.Z = Math.Sqrt( rad*rad - Math.Pow( (intersection - center).Abs(), 2 ) );

			double dilation = basePointUHS.Z / intersection.Z;
			intersection *= dilation;
			double translation = basePointUHS.X - intersection.X;
			intersection.X += translation;
			if( intersection != basePointUHS )
				throw new System.Exception();

			z *= dilation;
			z.X += translation;
			minusi *= dilation;
			minusi.X += translation;

			return new H3.Cell.Edge( minusi, z );
		}

		private static H3.Cell.Edge CalcFiber( Vector3D boundaryPointBall )
		{
			bool parallel = false;
			if( parallel )
			{
				Vector3D v2 = new Vector3D( 0, 0, -1 );	// Really we can pick any boundary point.
				return new H3.Cell.Edge( boundaryPointBall, v2 );
			}
			else
			{
				// If the point is above the real line, it will have been connected to z.
				// If below the real line, it will have been connected to -i
				// If on the real line, it's degenerate.
				Vector3D pUHS = H3Models.BallToUHS( boundaryPointBall );
				if( Tolerance.Equal( pUHS.Y, 0 ) )
					return null;

				Vector3D z = Vector3D.FromComplex( ZParam );
				Vector3D minusi = new Vector3D( 0, -1 );

				double dilation = double.NaN;
				double translation = double.NaN;
				if( pUHS.Y > 0 )
				{
					dilation = pUHS.Y / z.Y;
					z *= dilation;
					translation = pUHS.X - z.X;
					z.X += translation;
					if( H3Models.UHSToBall( z ) != boundaryPointBall )
						throw new System.Exception();

					minusi *= dilation;
					minusi.X += translation;
					return new H3.Cell.Edge( H3Models.UHSToBall( minusi ), boundaryPointBall );
					//return null;
				}
				else
				{
					dilation = pUHS.Y / -1;
					minusi *= dilation;
					translation = pUHS.X - minusi.X;
					minusi.X += translation;
					if( H3Models.UHSToBall( minusi ) != boundaryPointBall )
						throw new System.Exception();

					z *= dilation;
					z.X += translation;
					return new H3.Cell.Edge( boundaryPointBall, H3Models.UHSToBall( z ) );
					//return null;
				}
			}
		}

		private static H3.Cell.Edge[] ParameterizedFibers( Complex z )
		{
			List<H3.Cell.Edge> fibers = new List<H3.Cell.Edge>();

			double scale = 0.1;
			Vector3D v1 = new Vector3D( 0, -1 ); // -i
			Vector3D v2 = Vector3D.FromComplex( z );

			int count = (int)Math.Round( Math.Sqrt( NumFibers ), 0 );
			for( int i=-count/2; i<count/2; i++ )
			for( int j=1; j<count; j++ )	// dilations should remain positive.
			{
				double t1 = scale * i;
				double t2 = scale * j;

				// Apply the dilation first.
				Vector3D _v1 = v1;
				Vector3D _v2 = v2;
				_v1 *= t2;
				_v2 *= t2;
				_v1.X += t1;
				_v2.X += t1;

				fibers.Add( new H3.Cell.Edge(
					H3Models.UHSToBall( _v1 ),
					H3Models.UHSToBall( _v2 ) ) );
			}

			return fibers.ToArray();
		}

		/// <summary>
		/// Gets fibers for the only fibration with parallel (vs. ultraparallel) fibers.
		/// Returns result in the ball model.
		/// </summary>
		private static H3.Cell.Edge[] ParallelFibers()
		{
			List<H3.Cell.Edge> fibers = new List<H3.Cell.Edge>();

			double scale = 0.3;

			// Just a grid of vertical fibers.
			// We could do any kind of grid we want here really (square, hexagonal, random...)
			// Each would likely produce different results.
			// It'd be nice to figure out how to space out the results near the north pole.

			int count = (int)Math.Round( Math.Sqrt( NumFibers ), 0 );
			for( int i=-count/2; i<count/2; i++ )
			for( int j=-count/2; j<count/2; j++ )
			{
				//double off1 = Math.Pow( scale*i, 3 );
				//double off2 = Math.Pow( scale*j, 3 );
				double off1 = scale*i;
				double off2 = scale*j;

				Vector3D v1 = new Vector3D( off1, off2 );
				Vector3D v2 = new Vector3D( double.PositiveInfinity, 1 );

				// Don't use Infinity.InfinityVector, because we want to distiguish this point as being on the boundary.
				fibers.Add( new H3.Cell.Edge( 
					H3Models.UHSToBall( v1 ), 
					H3Models.UHSToBall( v2 ) ) );
			}

			return fibers.ToArray();
		}
	}
}
