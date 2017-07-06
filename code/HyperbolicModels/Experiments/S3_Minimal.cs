namespace HyperbolicModels
{
	using R3.Algorithm;
    using R3.Core;
	using R3.Drawing;
    using R3.Math;
    using R3.Geometry;
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Numerics;

    using Math = System.Math;

    public class S3_Minimal
	{
		/* Questions for Peter
		 - offset coord, how does this transform when we move from s to x?  A: doesn't change.
		 - hmmm, but they go below 1. Answered in email.
		 - Notation of fluxes.  Is F_1- a shorthand?  A: Yes
		 - What is "far enough away" from the catenoid?
		*/ 

        public class Params
        {
            public int K = 3;
            public int M = 7;
            public double Scale = .03;
            //public double CatHeight = 0.06;
            //public double CatWidth = 0.02;
            public int Res = 50;
        }

        public static Params m_params = new Params();

		public class RLD_outputs
		{
			/// <summary>
			/// The profile of the RLD solution.
			/// </summary>
			public Vector3D[] profile;

			/// <summary>
			/// The height at the s_i.
			/// </summary>
			public double[] phi_i;

			/// <summary>
			/// The scale factor of the graph.
			/// </summary>
			public double scale;

			/// <summary>
			/// x_i
			/// </summary>
			public double[] x_i;

			/// <summary>
			/// The catenoidal bridge waist sizes.
			/// </summary>
			public double[] t_i;
		}

        public static void RLD_Surface()
		{
			RLD_outputs outputs;
			Mesh mesh = new Mesh();
			SurfaceInternal( out outputs );
			double scale = m_params.Scale;

			// Now add in all the catenoids.
			double mInc = Math.PI * 2 / m_params.M;
			for( int k = 1; k < outputs.x_i.Length; k++ )
			for( int m = 0; m < m_params.M; m++ )
			{
				Vector3D loc = SphericalCoords.SphericalToCartesian( new Vector3D( 1, Math.PI/2 - outputs.x_i[k], m*mInc ) );
				mesh.Append( Catenoid( scale, loc, outputs.phi_i[k], outputs.t_i[k] ) );	
			}

			PovRay.WriteMesh( mesh, "RLD.pov" );
		}

		public static void CatenoidBasedSurface()
		{
			RLD_outputs outputs;
			SurfaceInternal( out outputs );
			double scale = m_params.Scale;

			// Map a point for a given k/m from the hemihypersphere to the complex plane.
			// You can also pass in -1 for k to get a point on the equator of the hemihypersphere.
			double mInc = Math.PI * 2 / m_params.M;
			Func<RLD_outputs, int, int, Vector3D> onPlane = ( o, k, m ) =>
			{
				double theta = k == -1 ? 0 : outputs.x_i[k];
				theta += Math.PI / 2;
				return
					Sterographic.SphereToPlane(
						SphericalCoords.SphericalToCartesian(
							new Vector3D( 1, theta, m * mInc )
						)
					);
			};

			// Setup texture coords on fundamental triangle.
			// We'll use a fundamental triangle in the southern hemisphere,
			// with stereographically projected coords at (0,0), (1,0), and CCW on the unit circle depending on M.  
			Polygon p = new Polygon();
			p.Segments.Add( Segment.Line( new Vector3D(), new Vector3D( 1, 0 ) ) );
			p.Segments.Add( Segment.Arc( new Vector3D( 1, 0 ), onPlane( outputs, 1, 1 ), onPlane( outputs, -1, 1 ) ) );
			p.Segments.Add( Segment.Line( onPlane( outputs, -1, 1 ), new Vector3D() ) );
			int levels = 9;
			TextureHelper.SetLevels( levels );
			Vector3D[] coords = TextureHelper.TextureCoords( p, Geometry.Spherical, doGeodesicDome: true );
			int[] elementIndices = TextureHelper.TextureElements( 1, levels );

			// Setup a nearTree for the catenoid locations (on the plane).
			NearTree nearTree = new NearTree( Metric.Spherical );
			for( int k = 1; k < outputs.x_i.Length; k++ )
			for( int m = 0; m <= 1; m++ )
			{
				Vector3D loc = onPlane( outputs, k, m );
				nearTree.InsertObject( new NearTreeObject() { ID = k, Location = loc } );
			}

			// Given a point on the plane, find the nearest catenoid center and calculate the height of the surface based on that.
			// This also calculates the locking of the point.
			Func<Vector3D, Tuple<double,Vector3D,Vector3D>> heightAndLocking = coord =>
			{
				NearTreeObject closest;
				if( !nearTree.FindNearestNeighbor( out closest, coord, double.MaxValue ) )
					throw new System.Exception();

				Vector3D locked = new Vector3D();
				if( p.Segments[0].IsPointOn( coord ) ||
					p.Segments[2].IsPointOn( coord ) )
					locked = new Vector3D( 1, 1, 0, 0 );
				//if( p.Segments[1].IsPointOn( v ) )		// Not working right for some reason, but line below will work.
				if( Tolerance.Equal( coord.Abs(), 1 ) )
					locked = new Vector3D( 1, 1, 1, 0 );

				Vector3D vSphere = Sterographic.PlaneToSphere( coord );
				Vector3D cSphere = Sterographic.PlaneToSphere( closest.Location );
				double dist = vSphere.AngleTo( cSphere );

				int k = (int)closest.ID;
				double waist = outputs.t_i[k];
				double rld_height = outputs.phi_i[k];

				double h = waist * 3.5 * 2;						// height where catenoid will meet rld_height.
				double factor = scale * rld_height * 2 / h;		// Artifical scaling so we can see things.
				dist /= factor;

				double z = double.NaN;
				if( dist >= waist )
				{
					z = waist * DonHatch.acosh( dist / waist );
				}
				else if( dist >= 0.7 * waist )
				{
					z = 0;

					// Move the coord to the thinnest waist circle.
					Mobius m = new Mobius();
					m.Hyperbolic( Geometry.Spherical, coord.ToComplex(), waist / dist );
					coord = m.Apply( coord );
				}

				if( dist < waist * 20 )
					locked = new Vector3D( 1, 1, 1, 1 );

				return new Tuple<double, Vector3D, Vector3D>( z*factor, locked, coord );
			};

			// Calculate all the coordinates.
			Vector3D[] locks = new Vector3D[coords.Length];
			for( int i=0; i<coords.Length; i++ )
			{
				Vector3D coord = coords[i];
				var hl = heightAndLocking( coord );
				locks[i] = hl.Item2;
				coord = hl.Item3;
				coords[i] = Normal( Sterographic.PlaneToSphere( coord ), (double)hl.Item1 );
			}

			// Relax it.
			Relax( coords, elementIndices, locks );

			Mesh mesh = new Mesh();
			Sphere s = new Sphere();
			for( int i=0; i<elementIndices.Length; i+=3 )
			{
				Vector3D a = coords[elementIndices[i]];
				Vector3D b = coords[elementIndices[i + 1]];
				Vector3D c = coords[elementIndices[i + 2]];
				if( a.DNE || b.DNE || c.DNE )
					continue;

				for( int m = 0; m <= 0; m++ )
				{
					mesh.Triangles.Add( new Mesh.Triangle( a, b, c ) );
					mesh.Triangles.Add( new Mesh.Triangle(
						s.ReflectPoint( a ),
						s.ReflectPoint( b ),
						s.ReflectPoint( c ) ) );
					a.RotateXY( mInc );
					b.RotateXY( mInc );
					c.RotateXY( mInc );
				}
			}

			PovRay.WriteMesh( mesh, "RLD.pov" );
		}

		private static void Relax( Vector3D[] coordsR3, int[] elementIndices, Vector3D[] locks )
		{
			int dim = 4;
			Graph g = new Graph();
			for( int i=0; i<coordsR3.Length; i++ )
			{
				Vector3D v = coordsR3[i];
				GraphNode node = new GraphNode( new VectorND( Sterographic.R3toS3( v ) ), new VectorND( dim ) );
				node.Lock = new VectorND( locks[i] );
				g.Nodes.Add( node );
			}

			for( int i=0; i<elementIndices.Length; i+=3 )
			{
				int a = elementIndices[i];
				int b = elementIndices[i + 1];
				int c = elementIndices[i + 2];
				g.AddEdge( new GraphEdge( a, b ) );
				g.AddEdge( new GraphEdge( b, c ) );
				g.AddEdge( new GraphEdge( c, a ) );
			}

			GraphRelaxation relaxer = new GraphRelaxation();
			relaxer.Graph = g;
			relaxer.NodeRepulsion = 0;
			relaxer.EdgeAttraction = 0.5;
			relaxer.EdgeRepulsion = 0;
			relaxer.Relax( 1000 );

			for( int i = 0; i < coordsR3.Length; i++ )
			{
				GraphNode node = g.Nodes[i];
				coordsR3[i] = Sterographic.S3toR3( node.Position.ToVec3D() );
			}
		}

		/// <summary>
		/// Given an x value (0 <= x <= pi/2) and an offset normal to the hemisphere in S^3,
		/// returns the stereographically projected point in R^3, in the xz plane.
		/// </summary>
		private static Vector3D Normal( double x, double o )
        {
            Vector3D center;
            double rad;
            H3Models.Ball.DupinCyclideSphere( new Vector3D( 1, 0 ), o, Geometry.Spherical, out center, out rad );
            double mag = center.X + rad;

            // In the xz plane.
            Complex c = Complex.FromPolarCoordinates( mag, x );
            return new Vector3D( c.Real, 0, c.Imaginary );
        }

		/// <summary>
		/// Given a location and an offset normal to the hemisphere in S^3,
		/// returns the stereographically projected point in R^3, in the xz plane.
		/// </summary>
		private static Vector3D Normal( Vector3D loc, double o )
		{
			Vector3D center;
			double rad;
			H3Models.Ball.DupinCyclideSphere( new Vector3D( 1, 0 ), o, Geometry.Spherical, out center, out rad );
			double mag = center.X + rad;

			loc.Normalize();
			loc *= mag;
			return loc;
		}

		/// <summary>
		/// Generate a catenoid, then move it to a point on the hemisphere in S^3.
		/// </summary>
		private static Mesh Catenoid( double scale, Vector3D loc, double rld_height, double waist )
        {
			double h = waist * 3.5;
			//Mesh mesh = StandardCatenoid( waist, h );
			Mesh mesh = CatenoidSquared( waist, h );
			mesh.Scale( scale * rld_height * 2 / h );	// To make the catenoid meet up with the RLD solution.

            // First move to north pole, then rotate down to the location.
            Func<Vector3D, Vector3D> transform = v =>
            {
                v = H3Models.BallToUHS( v );
				Vector3D northPole = new Vector3D( 0, 0, 1 );
                Vector3D axis = loc.Cross( northPole );
                if( !axis.Normalize() ) // North or south pole?
                    return v;
				double anglXY = Euclidean2D.AngleToCounterClock( new Vector3D( 1, 0 ), new Vector3D( loc.X, loc.Y ) );
				v.RotateXY( anglXY );
				double angleDown = loc.AngleTo( northPole );
                v.RotateAboutAxis( axis, angleDown );

				if( v.DNE || Infinity.IsInfinite( v ) )
					throw new System.Exception();
				return v;
            };

            mesh.Transform( transform );
            return mesh;
        }

		/// <summary>
		/// Calculates a mesh for a standard euclidean catenoid.
		/// This will need to be transformed to the various locations later.
		/// </summary>
		private static Mesh StandardCatenoid( double waist, double height )
        {
            Mesh mesh = new Mesh();
            int res = m_params.Res*2;

            Func<double, Vector3D[]> oneCircle = z =>
            {
                double r = waist * Math.Cosh( z / waist );
                Vector3D cen = new Vector3D( 0, 0, z );
                Vector3D radius = new Vector3D( r, 0 );
                return Shapeways.Disk( cen, new Vector3D( 0, 0, 1 ), radius, res );
            };

            double inc = height/(res*2);
            for( int i=0; i<res; i++ )
            {
                double z1 = inc*i;
                double z2 = inc*(i+1);
                mesh.AddBand( oneCircle( z1 ), oneCircle( z2 ) );
                mesh.AddBand( oneCircle( -z1 ), oneCircle( -z2 ) );
            }

            return mesh;
        }

		/// <summary>
		/// Calculates a mesh for a standard euclidean catenoid.
		/// This will need to be transformed to the various locations later.
		/// 
		/// Like above, but we adjust the xy components of the mesh using one of the mappings described here:
		/// https://arxiv.org/ftp/arxiv/papers/1509/1509.06344.pdf
		/// I found that paper here:
		/// https://stackoverflow.com/questions/13211595/how-can-i-convert-coordinates-on-a-circle-to-coordinates-on-a-square
		/// This is so we can connect up to the RLD mesh later.
		/// </summary>
		private static Mesh CatenoidSquared( double waist, double height )
		{
			Mesh mesh = new Mesh();
			int res = m_params.Res * 2;
			double diskRad = waist * Math.Cosh( height/2 / waist ); ;

			// NOTE: A band is *not* a constant height slice,
			//		 so the input z value is the height at the edge midpoints of the square.
			Func<double, Vector3D[]> oneCircle = z =>
			{
				bool neg = z < 0;
				z = Math.Abs( z );

				// Radius on disk at a starting edge midpoint.
				double r = waist * Math.Cosh( z / waist );
				Vector3D start = new Vector3D( r, 0 );
				Vector3D axis = new Vector3D( 0, 0, 1 );

				List<Vector3D> points = new List<Vector3D>();
				double angleInc = 2 * Math.PI / res;
				double angle = 0;
				for( int i = 0; i < res; i++ )
				{
					Vector3D point = start;
					point.RotateAboutAxis( axis, angle );
					point = DiskToSquare( point, diskRad );

					double zi = waist * DonHatch.acosh( point.Abs() / waist );
					if( double.IsNaN( zi ) )
						zi = 0;
					if( neg )
						zi *= -1;
					Vector3D newPoint = new Vector3D( point.X, point.Y, zi );
					if( newPoint.DNE )
						throw new System.Exception();
					points.Add( newPoint );

					angle += angleInc;
				}

				return points.ToArray();
			};

			double inc = height / (res * 2);
			for( int i = 0; i < res; i++ )
			{
				double z1 = inc * i;
				double z2 = inc * (i + 1);
				mesh.AddBand( oneCircle( z1 ), oneCircle( z2 ) );
				mesh.AddBand( oneCircle( -z1 ), oneCircle( -z2 ) );
			}

			return mesh;
		}

		private static Vector3D DiskToSquare( Vector3D disk, double diskRad )
		{
			disk /= diskRad;

			double u = disk.X;
			double v = disk.Y;

			double u2mv2 = u*u - v*v;
			double ut = 2 * Math.Sqrt( 2 ) * u;
			double vt = 2 * Math.Sqrt( 2 ) * v;
			double x = .5 * Math.Sqrt( 2 + u2mv2 + ut ) - 0.5 * Math.Sqrt( 2 + u2mv2 - ut );
			double y = .5 * Math.Sqrt( 2 - u2mv2 + vt ) - 0.5 * Math.Sqrt( 2 - u2mv2 - vt );
			Vector3D result = new Vector3D( x, y );
			result *= diskRad;
			return result;
		}

        /// <summary>
        /// Calculates a mesh for an RLD surface.
        /// </summary>
        private static Mesh SurfaceInternal( out RLD_outputs outputs )
		{
			Mesh mesh = new Mesh();

			Func<double, double, double, Vector3D[]> oneCircle = ( x, o, scale ) =>
			{
                Vector3D normal = Normal( x, o * scale );
				Vector3D cen = new Vector3D( 0, 0, normal.Z );
				Vector3D radius = new Vector3D( normal.X, 0 );
				return Shapeways.Disk( cen, new Vector3D( 0, 0, 1 ), radius, m_params.Res );
			};

			outputs = RLD_Sphere( m_params.K );
			double s = m_params.Scale;
			//double s = outputs.scale;

			// Add in two bands for each segment along the profile, one for +z coords, and one for -z coords.
			Vector3D[] profile = outputs.profile;
			for( int i = 0; i < profile.Length - 1; i++ )
			{
				Vector3D p1 = profile[i];
				Vector3D p2 = profile[i + 1];
				mesh.AddBand( oneCircle( p1.X, p1.Y, s ), oneCircle( p2.X, p2.Y, s ) );
				mesh.AddBand( oneCircle( p1.X, -p1.Y, s ), oneCircle( p2.X, -p2.Y, s ) );
			}

			return mesh;
		}

        // Returns the profile for the RLD surface.
        // First component is the X value, second component is the offset.
		private static RLD_outputs RLD_Sphere( int k )
		{
			RLD_Calculator calc = new RLD_Calculator( k );
			calc.Calc();

			RLD_outputs outputs = new RLD_outputs();
			outputs.profile = calc.ProfileS3();
			outputs.phi_i = calc.s_i.Select( s => calc.Phi( s ) ).ToArray();
			outputs.x_i = calc.s_i.Select( s => CylToSphere( s ) ).ToArray();
			outputs.t_i = calc.t_i();
			outputs.scale = outputs.t_i[1] * m_params.M / (calc.Phi( calc.s_i[1] ) * calc.Flux);
			return outputs;
		}

		/// <summary>
		/// We need to calculate all the s_i values and the A_i, B_i constants.
		/// Here's the process...
		/// Start with a small s_1
		/// ...
		/// NOTE: This class works on the cylinder coordinates.
		/// </summary>
		private class RLD_Calculator
		{
			public RLD_Calculator( int k )
			{
				K = k;
				s_i = new double[k+1];
				A_i = new double[k+1];
				B_i = new double[k+1];
				s_i[0] = 0;
				A_i[0] = 1;
				B_i[0] = 0;
			}

			public int K { get; set; }
			public readonly double[] s_i;
			public readonly double[] A_i;
			public readonly double[] B_i;

			/// <summary>
			/// Returns the offsets at points along the equatorial sphere.
			/// (X is x coordinate from the paper, Y is offset).
			/// </summary>
			public Vector3D[] ProfileS3()
			{
				List<Vector3D> result = new List<Vector3D>();

				int numPoints = 250;
				double x = 0;
				double xInc = Math.PI / 2 / (numPoints - 1);
				for( int i=0; i<numPoints-1; i++ )
				{
					double offset = Phi( SphereToCyl( x ) );	// Offset is equivalent in cylindrical or spherical.
					result.Add( new Vector3D( x, offset ) );
					x += xInc;
				}

				// Add in the last one, since the calcs above result in DNEs.
				result.Add( new Vector3D( Math.PI / 2, 1 ) );

				return result.ToArray();
			}

			public double[] t_i()
			{
				double[] result = new double[K + 1];

				double C = 1;
				double m = m_params.M;
				double t_1 = C / m * Math.Exp( -m / Flux );
				result[1] = t_1;

				for( int i=2; i<=K; i++ )
					result[i] = Phi( s_i[i] ) / Phi( s_i[1] ) * t_1;

				return result;
			}

			public void Calc()
			{
				double sStart = S3_Minimal.SphereToCyl( Math.PI/2/K/10 );
				double offset = sStart / 10;

				// Using Newton's method didn't work for this (making A_i go to zero),
				// without first getting close to the roots.
				double sRunning = sStart;
				for( int i = 0; i < 1000; i++ )
				{
					sRunning += offset;
					double result = CalcOne( sRunning );
					if( result > 0 )
						continue;

					Newton n = new Newton();
					n.MaxIterations = 25;
					n.Function = s => CalcOne( s );
					double s_1 = n.FindRoot( sRunning );
					double[] x_i = s_i.Select( si => CylToSphere( si ) ).ToArray();
					return;
				}
			}

			/// <summary>
			/// Do the calculation for a starting testing s_1
			/// This returns as output the last A_i, which we want to go to zero.
			/// so that we can iterate on this until the closing condition is fulfilled.
			/// </summary>
			private double CalcOne( double sTest )
			{
				// Set s_1 so we can calculate everything.
				s_i[1] = sTest;

				for( int i=1; i<K; i++ )
				{
					CalcA_i( i );
					CalcB_i( i );
					CalcNextS( i );
				}

				// Calculate the last values.
				CalcA_i( K );
				CalcB_i( K );
				return A_i[K];
			}

			public double Phi( double s )
			{
				// Get the right index.
				int idx = 0;
				while( idx < K && s > s_i[idx+1] )
					idx++;

				return phi( s, A_i[idx], B_i[idx] );
			}

			private void CalcA_i( int i )
			{
				if( i == 0 )
				{
					A_i[0] = 1;
					return;
				}

				A_i[i] = A_i[i - 1] - 2 * phi( i ) * Flux * phiOdd( s_i[i] );
			}
			private void CalcB_i( int i )
			{
				if( i == 0 )
				{
					B_i[0] = 0;
					return;
				}

				B_i[i] = B_i[i - 1] + 2 * phi( i ) * Flux * phiEven( s_i[i] );
			}

			/// <summary>
			/// https://www.wolframalpha.com/input/?i=y+%3D+(a*(-tanh(x)+-+x*(sech(x))%5E2)+%2B+b*(sech(x))%5E2)%2F(a*(1-x*tanh(x))%2Bb*tanh(x)),+solve+for+x
			/// </summary>
			private void CalcNextS( int i )
			{
				// Find the s where the flux matches.
				Newton n = new Newton();
				n.Function = s => Flux + CalcFlux( s, A_i[i], B_i[i] );
				double nextS = n.FindRoot2( s_i[i] );
				s_i[i + 1] = nextS;
			}

			/// <summary>
			/// We can just calculate this at s_1, since by construction
			/// we make all the fluxes the same. s_1 should be set before calling this.
			/// </summary>
			public double Flux
			{
				get
				{
					double s = s_i[1];
					double A = A_i[0];
					double B = B_i[0];
					return -CalcFlux( s, A, B );
				}
			}

			private static double CalcFlux( double s, double A, double B )
			{
				return dPhi( s, A, B ) / phi( s, A, B );
			}

			/// <summary>
			/// NOTE: We calculate this from the left, so that we don't need to know A_i, B_i.
			/// </summary>
			private double phi( int i )
			{
				if( i == 0 )
					return phiEven( 0 );    // s_0 = 0

				return phi( s_i[i], A_i[i - 1], B_i[i - 1] );
			}

			private static double phi( double s, double A, double B )
			{
				return A * phiEven( s ) + B * phiOdd( s );
			}

			/// <summary>
			/// The derivative of the phi function.
			/// https://www.wolframalpha.com/input/?i=derivative+of+a*(1-x*tanh(x))%2Bb*tanh(x)
			/// </summary>
			private static double dPhi( double s, double A, double B )
			{
				double sech_squared = 1 / Math.Pow( Math.Cosh( s ), 2 );
				return A * (-Math.Tanh( s ) - s * sech_squared) + B * sech_squared;
			}

			private static double phiEven( double s )
			{
				return 1 - s * Math.Tanh( s );
			}

			private static double phiOdd( double s )
			{
				return Math.Tanh( s );
			}
		}

		public static double CylToSphere( double s )
		{
			return Math.Asin( Math.Tanh( s ) );
		}

		public static double SphereToCyl( double x )
		{
			return DonHatch.atanh( Math.Sin( x ) );
		}
	}

	/// <summary>
	/// https://en.wikipedia.org/wiki/Newton%27s_method
	/// </summary>
	internal class Newton
	{
		/// <summary>
		/// The function from x -> y
		/// </summary>
		public Func<double, double> Function { get; set; }

		/// <summary>
		/// The tolerance we want to find the root at.
		/// </summary>
		public double Tolerance = 1e-5;

		public int MaxIterations = 1000;

		public double Scale = 0.0001;

		public double FindRoot( double startX )
		{
			double root = double.NaN;
			double x = startX;

			int count = 0;
			while( count < MaxIterations )
			{
				double y = Function( x );
				if( Math.Abs( y ) < Tolerance )
					return x;

				double dOff = Scale;
				double derivative = (Function( x + dOff ) - Function( x - dOff )) / (2 * dOff);
				x -= y / derivative;
				count++;
			}

			return root;
		}

        /// <summary>
        /// Newton's method isn't working because of discontinuities.
		/// Let's go until we hit our first negative, then call the standard Newton's method.
        /// </summary>
        public double FindRoot2( double startX )
		{
			double x = startX;

			int count = 0;
			while( count < MaxIterations * 100 )
			{
				double y = Function( x );
				if( y < 0 )
					break;

				x += Scale;
				count++;
			}

			return FindRoot( x );
		}
	}
}
