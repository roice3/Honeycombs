namespace HyperbolicModels
{
	using R3.Core;
	using R3.Math;
	using R3.Geometry;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	public class GreatSphere
	{
		public GreatSphere( Vector3D pole ) { Pole = pole; }

		/// <summary>
		/// Pole, in 4D coordinates.
		/// </summary>
		public Vector3D Pole
		{
			get
			{
				return m_pole;
			}
			set
			{
				m_pole = value;

				/* Not correct (doesn't handle points on equitorial 2-sphere
				if( Tolerance.LessThan( value.W, 0 ) )
					m_pole = value * -1;
				else
					m_pole = value;*/
			}
		}
		private Vector3D m_pole;

		/// <summary>
		/// Calculates our pole from a stereographically projected sphere.
		/// Assume the input sphere is a geodesic sphere in the conformal model.
		/// </summary>
		public static GreatSphere FromSphere( Sphere s )
		{
			if( s.Center.IsOrigin )
			{
				return new GreatSphere( new Vector3D( 0, 0, 0, 1 ) );
			}

			if( s.IsPlane )
			{
				// If we are a plane, the pole is is on the unit sphere.
				Vector3D pole = s.Normal;
				pole.Normalize();
				return new GreatSphere( Sterographic.R3toS3( pole ) );
			}

			Vector3D v1 = Sterographic.R3toS3( s.Center + new Vector3D( s.Radius, 0, 0, 0 ) );
			Vector3D v2 = Sterographic.R3toS3( s.Center + new Vector3D( 0, s.Radius, 0, 0 ) );
			Vector3D v3 = Sterographic.R3toS3( s.Center + new Vector3D( 0, 0, s.Radius, 0 ) );

			Vector3D o = Orthogonal( v1, v2, v3 );
			return new GreatSphere( o );
		}

		/// <summary>
		/// Returns a stereographically sphere representing us.
		/// </summary>
		public Sphere ToSphere()
		{
			// Equatorial sphere?
			if( Pole == new Vector3D( 0, 0, 0, 1 ) )
				return new Sphere();

			// A plane?
			Vector3D poleR3 = Sterographic.S3toR3( Pole );
			if( Tolerance.Equal( poleR3.Abs(), 1 ) )
				return Sphere.Plane( poleR3 );

			// Get 4 points on the sphere.
			Vector3D e1 = new Vector3D( 1, 0, 0, 0 );
			Vector3D e2 = new Vector3D( 0, 1, 0, 0 );
			Vector3D e3 = new Vector3D( 0, 0, 1, 0 );
			Vector3D e4 = new Vector3D( 0, 0, 0, 1 );

			System.Func<Vector3D, Vector3D> one = v =>
			{
				v = Euclidean3D.ProjectOntoPlane( Pole, new Vector3D(), v );
				v.Normalize();
				return Sterographic.S3toR3( v );
			};

			return Sphere.From4Points( one( e1 ), one( e2 ), one( e3 ), one( e4 ) );	
		}

		/// <summary>
		/// Return a normalized vector orthogonal to 3 vectors.
		/// https://math.stackexchange.com/questions/904172/how-to-find-a-4d-vector-perpendicular-to-3-other-4d-vectors
		/// </summary>
		public static Vector3D Orthogonal( Vector3D a, Vector3D b, Vector3D c )
		{
			Vector3D e1 = new Vector3D( 1, 0, 0, 0 );
			Vector3D e2 = new Vector3D( 0, 1, 0, 0 );
			Vector3D e3 = new Vector3D( 0, 0, 1, 0 );
			Vector3D e4 = new Vector3D( 0, 0, 0, 1 );

			Vector3D det =
				e4 * a.Z * b.Y * c.X - e3 * a.W * b.Y * c.X - e4 * a.Y * b.Z * c.X + e2 * a.W * b.Z * c.X +
				e3 * a.Y * b.W * c.X - e2 * a.Z * b.W * c.X - e4 * a.Z * b.X * c.Y + e3 * a.W * b.X * c.Y +
				e4 * a.X * b.Z * c.Y - e1 * a.W * b.Z * c.Y - e3 * a.X * b.W * c.Y + e1 * a.Z * b.W * c.Y +
				e4 * a.Y * b.X * c.Z - e2 * a.W * b.X * c.Z - e4 * a.X * b.Y * c.Z + e1 * a.W * b.Y * c.Z +
				e2 * a.X * b.W * c.Z - e1 * a.Y * b.W * c.Z - e3 * a.Y * b.X * c.W + e2 * a.Z * b.X * c.W +
				e3 * a.X * b.Y * c.W - e1 * a.Z * b.Y * c.W - e2 * a.X * b.Z * c.W + e1 * a.Y * b.Z * c.W;

			det.Normalize();
			return det;
		}
	}

	/// <summary>
	/// A comparison that equates antipodes.
	/// </summary>
	public class AntipodeEqualityComparer : IEqualityComparer<Vector3D>
	{
		public AntipodeEqualityComparer( bool stereo )
		{
			m_stereo = stereo;
		}

		// Whether we work on stereographically projected coordinates.
		private readonly bool m_stereo;

		public bool Equals( Vector3D v1, Vector3D v2 )
		{
			if( v1 == v2 )
				return true;

			Vector3D t1 = m_stereo ? Sterographic.R3toS3( v1 ) : v1;
			Vector3D t2 = m_stereo ? Sterographic.R3toS3( v2 ) : v2;
			if( t1 == t2 * -1 )
				return true;

			return false;
		}

		public int GetHashCode( Vector3D v )
		{
			Vector3D t = m_stereo ? Sterographic.R3toS3( v ) : v;
			t = new Vector3D( Math.Abs( t.X ), Math.Abs( t.Y ), Math.Abs( t.Z ), Math.Abs( t.W ) );
			return t.GetHashCode();
		}
	}

	public class GyrationPlaneEqualityComparer : IEqualityComparer<Circle3D>
	{
		public bool Equals( Circle3D c1, Circle3D c2 )
		{
			if( c1.Center != c2.Center )
				return false;

			if( c1.Normal == c2.Normal )
				return true;
			if( c1.Normal == c2.Normal * -1 )
				return true;

			return false;
		}

		public int GetHashCode( Circle3D c )
		{
			Vector3D t = c.Normal;
			t = new Vector3D( Math.Abs( t.X ), Math.Abs( t.Y ), Math.Abs( t.Z ), Math.Abs( t.W ) );
			return c.Center.GetHashCode() ^ t.GetHashCode();
		}
	}

	public class PointGroups
	{
		public void Gen(int p, int q, int r)
		{
			Geometry g = Util.GetGeometry( p, q, r );
			if( g != Geometry.Spherical )
				throw new System.Exception( "Point group code only for spherical geometry." );

			Simplex simplex = new Simplex();
			simplex.Facets = SimplexCalcs.Mirrors( p, q, r );
			Vector3D cen = new Vector3D();
			Vector3D faceCenter = SimplexCalcs.FaceCenterSpherical( p, q, r );
			Vector3D edgeMid = SimplexCalcs.EdgeMidpointSpherical( p, q, r );
			Vector3D vertex = SimplexCalcs.VertexSpherical( p, q, r );

			List<Vector3D> startingPoles = new List<Vector3D>();
			startingPoles.Add( Sterographic.S3toR3( GreatSphere.FromSphere( simplex.Facets[0] ).Pole ) );
			GreatSphere[] spheres = CalcSpheres( simplex.Facets, startingPoles.ToArray() );

			string filename = "point_group.pov";
			using( StreamWriter sw = File.CreateText( filename ) )
			{
				//foreach( var greatSphere in spheres )
					//sw.WriteLine( PovRay.Sphere( greatSphere.ToSphere() ) );
			}

			/*
			List<H3.Cell.Edge> startingEdges = new List<H3.Cell.Edge>();
			startingEdges.Add( new H3.Cell.Edge( cen, faceCenter ) );
			H3.Cell.Edge[] edges = Recurse.CalcEdges( simplex.Facets, startingEdges.ToArray(), new Recurse.Settings() { G = Geometry.Spherical, Threshold = 0.001 } );
			//edges = edges.Where( e => !( Infinity.IsInfinite( e.Start ) || Infinity.IsInfinite( e.End ) ) ).ToArray();

			PovRay.WriteEdges( new PovRay.Parameters() { AngularThickness = 0.01 }, Geometry.Spherical, edges, filename, append: true );
			*/

			double minRad = 0;
			double thick = 0.005;
			System.Func<Vector3D,Sphere> sizeFunc = v =>
			{
				Vector3D c;
				double rad;
				H3Models.Ball.DupinCyclideSphere( v, thick / 2, g, out c, out rad );
				return new Sphere() { Center = c, Radius = Math.Max( rad, minRad ) };
			};

			// All geodesics
			List<Circle3D> startingCircles = new List<Circle3D>();
			startingCircles.Add( GeodesicFrom2Points( cen, edgeMid ) );
			Circle3D[] geodesics = CalcGeodesics( simplex.Facets, startingCircles.ToArray() );

			Vector3D color = new Vector3D( 0, 0, 1 );
			using( StreamWriter sw = File.CreateText( filename ) )
			{
				Shapeways shapeways = new Shapeways();
				foreach( Circle3D geodesic in geodesics )
				{
					Vector3D[] points;
					if( Infinity.IsInfinite( geodesic.Radius ) )
					{
						double cutoff = 15;
						Segment seg = Segment.Line( geodesic.Normal * cutoff, geodesic.Normal * -cutoff );
						points = seg.Subdivide( 42 ); 
					}
					else
					{
						List<Vector3D> tempPoints = geodesic.Subdivide( 150 ).ToList();
						tempPoints.Add( tempPoints[0] );
						tempPoints.Add( tempPoints[1] );
						points = tempPoints.ToArray();
					}

					List<Vector3D> ePoints = new List<Vector3D>();
					List<double> eRadii = new List<double>();
					foreach( Vector3D pNE in points )
					{
						Sphere sphere = sizeFunc( pNE );
						ePoints.Add( sphere.Center );
						eRadii.Add( sphere.Radius );
					}
					shapeways.AddCurve( ePoints.ToArray(), eRadii.ToArray() );

					sw.WriteLine( PovRay.EdgeSphereSweep( points, sizeFunc, color ) );
				}

				STL.SaveMeshToSTL( shapeways.Mesh, "533.stl" );
			}
		}

		public static Circle3D GeodesicFrom2Points( Vector3D a, Vector3D b )
		{
			if( a == b || a == -b )
				throw new System.Exception( "Geodesic not unique" );

			System.Func<Vector3D, Circle3D> lineFunc = p =>
			{
				Circle3D circ = new Circle3D();
				circ.Radius = double.PositiveInfinity;
				p.Normalize();
				circ.Normal = p;    // Hacky representation of a line.
				return circ;
			};

			if( a.IsOrigin || b.IsOrigin )
			{
				Vector3D p = a.IsOrigin ? b : a;
				return lineFunc( p );
			}

			double mag1 = a.Abs(), mag2 = b.Abs();
			if( Tolerance.Equal( mag1, 1 ) && Tolerance.Equal( mag2, 1 ) )
				return new Circle3D( a, b, a * -1 );

			// The antipode in S^3 of a or b will give us a 3rd point.
			Vector3D antipode = Sterographic.S3toR3( Sterographic.R3toS3( a ) * -1 );

			// If the antipode is also an an antipode in R3, the points are colinear
			// (or they are on the equatorial 2-sphere, but that is checked above).
			if( a == antipode * -1 )
				return lineFunc( a );

			return new Circle3D( a, b, antipode );
		}

		public static GreatSphere[] CalcSpheres( Sphere[] simplex, Vector3D[] poles )
		{
			HashSet<Vector3D> completedPoles = new HashSet<Vector3D>( poles, new AntipodeEqualityComparer( stereo: true ) );
			ReflectSpheresRecursive( simplex, completedPoles.ToArray(), completedPoles );
			return completedPoles.Select( p => new GreatSphere( Sterographic.R3toS3( p ) ) ).ToArray();
		}

		private static void ReflectSpheresRecursive( Sphere[] simplex, Vector3D[] poles, HashSet<Vector3D> completedPoles )
		{
			if( 0 == poles.Length )
				return;

			HashSet<Vector3D> newPoles = new HashSet<Vector3D>( new AntipodeEqualityComparer( stereo: true ) );

			foreach( Vector3D pole in poles )
			for( int m = 0; m < simplex.Length; m++ )
			{
				Sphere mirror = simplex[m];
				Vector3D reflected = mirror.ReflectPoint( pole );
				GreatSphere gs = new GreatSphere( reflected );
				if( completedPoles.Add( gs.Pole ) )
					newPoles.Add( gs.Pole );
			}

			ReflectSpheresRecursive( simplex, newPoles.ToArray(), completedPoles );
		}

		public static Circle3D[] CalcGeodesics( Sphere[] simplex, Circle3D[] geodesics )
		{
			HashSet<Circle3D> completed = new HashSet<Circle3D>( geodesics, new GyrationPlaneEqualityComparer() );
			ReflectGeodesicsRecursive( simplex, completed.ToArray(), completed );
			return completed.ToArray();
		}

		private static void ReflectGeodesicsRecursive( Sphere[] simplex, Circle3D[] geodesics, HashSet<Circle3D> completed )
		{
			if( 0 == geodesics.Length )
				return;

			HashSet<Circle3D> newGeodesics = new HashSet<Circle3D>( new GyrationPlaneEqualityComparer() );

			foreach( Circle3D geodesic in geodesics )
			for( int m = 0; m < simplex.Length; m++ )
			{
				Sphere mirror = simplex[m];

				Vector3D[] points = Infinity.IsInfinite( geodesic.Radius ) ?
					new Vector3D[] { -geodesic.Normal, new Vector3D(), geodesic.Normal } :
					geodesic.RepresentativePoints;

				Vector3D[] reflected = points.Select( p => mirror.ReflectPoint( p ) ).ToArray();
				Circle3D newCirc = new Circle3D( reflected[0], reflected[1], reflected[2] );

				if( completed.Add( newCirc ) )
					newGeodesics.Add( newCirc );
			}

			ReflectGeodesicsRecursive( simplex, newGeodesics.ToArray(), completed );
		}
	}
}
