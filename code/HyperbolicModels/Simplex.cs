namespace R3.Geometry
{
	using System.Collections.Generic;
	using System.Linq;
	using System.Numerics;
	using R3.Core;
	using R3.Geometry;
	using R3.Math;
	using Math = System.Math;

	public class Simplex
	{
		public Simplex() {}

		public Vector3D[] Verts { get; set; }
		public Sphere[] Facets { get; set; }
		
		/// <summary>
		/// Construct an edge using verts indices a -> b.
		/// </summary>
		public H3.Cell.Edge Edge( int a, int b )
		{
			return new H3.Cell.Edge( Verts[a], Verts[b] );
		}

		public void InitializeGoursat()
		{
			//Vector3D[] test = SimplexCalcs.GoursatTetrahedron( 3.5, 3.8, 3.1, 2.2, 2.01, 2.1 );	// Example values from paper.
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 4, 3, 2, 3, 3 );				// 4,3,3,3
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 5, 3, 2, 3, 3 );				// 5,3,3,3
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 4, 3, 2, 4, 3 );				// 4,3,4,3
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 4, 3, 2, 5, 3 );				// 4,3,5,3
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 5, 3, 2, 5, 3 );				// 5,3,5,3
			//Verts = SimplexCalcs.GoursatTetrahedron( 3, 5, 3, 2, 2, 2 );				// 5,3^1,1
			//Verts = SimplexCalcs.GoursatTetrahedron( 3, 2, 2, 2, 5, 3 );				// 5,3^1,1 alt (to avoid vertices at origin).

			// Paracompact doesn't work :(
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 6, 3, 2, 6, 3 );				// 6,3,6,3

			// Regular
			Verts = SimplexCalcs.GoursatTetrahedron( 2, 5, 3, 2, 4, 2 );						// 5,3,4
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 4, 3, 2, 5, 2 );						// 4,3,5
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 5, 3, 2, 5, 2 );						// 5,3,5
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 5, 2, 2, 5, 3 );						// 5,3,5 (Alt)
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 3, 5, 2, 3, 2 );						// 3,5,3

			// Spherical/Paracompact doesn't work :(
			//Verts = SimplexCalcs.GoursatTetrahedron( 2, 5, 3, 2, 3, 2 );						// 5,3,3

			Facets = SimplexCalcs.Mirrors( Verts );
		}

		public Vector3D ReflectInFacet( Vector3D v, int facet )
		{
			return Facets[facet].ReflectPoint( v );
		}
	}

	public static class SimplexCalcs
	{
		public static void ToPovRay( Sphere[] mirrors )
		{
			System.IO.File.Delete( "simplex.pov" );
			PovRay.CreateSimplex( mirrors, "simplex.pov" );
		}

		public static Vector3D VertexPointKlein( int p, int q, int r )
		{
			Geometry vertexGeometry = Geometry2D.GetGeometry( q, r );
			if( vertexGeometry != Geometry.Hyperbolic )
			{
				throw new System.NotImplementedException();
			}

			Sphere[] mirrors = Mirrors( p, q, r );
			Sphere klein = H3Models.BallToKlein( mirrors[0] );
			Vector3D off = klein.Offset;
			double h = off.Abs() / Math.Cos( off.AngleTo( new Vector3D( 0, 0, -1 ) ) );
			return new Vector3D( 0, 0, -h );
		}

		/// <summary>
		/// Calculates the point of our simplex that is at a vertex.
		/// </summary>
		public static Vector3D VertexPointBall( int p, int q, int r )
		{
			Geometry vertexGeometry = Geometry2D.GetGeometry( q, r );
			if( vertexGeometry == Geometry.Hyperbolic )
			{
				// Outside the ball, and not in a good way.  Use the Klein version.
				//throw new System.NotImplementedException();
				return new Vector3D(0,0,-1);
			}

			// Get in UHS first.
			Sphere cellFacet = Mirrors( p, q, r, moveToBall: false ).First();
			double rSquared = Math.Pow( cellFacet.Radius, 2 );
			double cSquared = Math.Pow( cellFacet.Center.Abs(), 2 );

			Vector3D uhs;
			if( Tolerance.Equal( rSquared, cSquared ) )	// e.g. 363
			{
				uhs = new Vector3D();
			}
			else
			{
				double height = Math.Sqrt( rSquared - cSquared );
				uhs = new Vector3D( 0, 0, height );
			}
			return H3Models.UHSToBall( uhs );
		}

		/// <summary>
		/// Calculates the point of our simplex that is at the middle of an edge.
		/// </summary>
		private static Vector3D MidEdgePointBall( int p, int q, int r )
		{
			// We need the mid-radius, but we have to do the calculation
			// with our Euclidean simplex mirrors (to avoid infinities that happen in the formulas).
			Circle3D edge = HoneycombEdgeUHS( p, q, r );
			if( edge.Radius == 0 )
				return edge.Center;

			Geometry cellGeometry = Geometry2D.GetGeometry( p, q );
			switch( cellGeometry )
			{
				case Geometry.Spherical:
				{
					Sphere sphereInBall = H3Models.UHSToBall( new Sphere() { Center = edge.Center, Radius = edge.Radius } );
					Vector3D mid = sphereInBall.ProjectToSurface( new Vector3D() );	// Project origin to sphere.
					return mid;
				}
				case Geometry.Euclidean:
				{
					Vector3D mid = H3Models.UHSToBall( edge.Center + new Vector3D( 0, 0, edge.Radius ) );
					return mid;
				}
				case Geometry.Hyperbolic:
				{
					throw new System.NotImplementedException();
				}
			}

			throw new System.ArgumentException();
		}

		/// <summary>
		/// Calculates the point of our simplex that is at the middle of a face.
		/// </summary>
		private static Vector3D FaceCenterBall( int p, int q, int r )
		{
			Geometry cellGeometry = Geometry2D.GetGeometry( p, q );
			switch( cellGeometry )
			{
			case Geometry.Spherical:
				{
					Vector3D cellCenter = CellCenterBall( p, q, r );
					Sphere[] mirrors = Mirrors( p, q, r, moveToBall: true );
					return mirrors[0].ProjectToSurface( cellCenter );
				}
			case Geometry.Euclidean:
				{
					Sphere[] mirrors = Mirrors( p, q, r, moveToBall: false );
					Vector3D faceCenterUHS = mirrors[0].Center;
					faceCenterUHS.Z += mirrors[0].Radius;
					return H3Models.UHSToBall( faceCenterUHS );
				}
			case Geometry.Hyperbolic:
				{
					throw new System.NotImplementedException();
				}
			}

			throw new System.ArgumentException();
		}

		/// <summary>
		/// Calculates the point of our simplex that is at the center of a cell.
		/// </summary>
		private static Vector3D CellCenterBall( int p, int q, int r )
		{
			Geometry cellGeometry = Geometry2D.GetGeometry( p, q );
			switch( cellGeometry )
			{
			case Geometry.Spherical:
				{
					return new Vector3D();
				}
			case Geometry.Euclidean:
				{
					return new Vector3D( 0, 0, 1 );
				}
			case Geometry.Hyperbolic:
				{
					//throw new System.NotImplementedException();
					return new Vector3D( 0, 0, 1 );
				}
			}

			throw new System.ArgumentException();
		}

		public static H3.Cell.Edge HoneycombEdgeBall( int p, int q, int r )
		{
			Sphere[] facets = Mirrors( p, q, r, moveToBall: true );
			Vector3D[] verts = VertsBall( p, q, r );
			return HoneycombEdgeBall( facets, verts[2] );
		}

		public static H3.Cell.Edge HoneycombEdgeBall( Sphere[] facets, Vector3D vertex )
		{
			Vector3D v1 = vertex;
			Vector3D v2 = facets[3].ReflectPoint( v1 );
			return new H3.Cell.Edge( v1, v2 );
		}

		public static H3.Cell.Edge DualEdgeBall( Sphere[] facets )
		{
			Vector3D v1 = new Vector3D();
			Vector3D v2 = facets[0].ReflectPoint( v1 );
			return new H3.Cell.Edge( v1, v2 );
		}

		private static Circle3D HoneycombEdgeUHS( int p, int q, int r )
		{
			Sphere[] simplex = Mirrors( p, q, r, moveToBall: false );
			Sphere s1 = simplex[0].Clone();
			Sphere s2 = s1.Clone();
			s2.Reflect( simplex[1] );
			Circle3D intersection = s1.Intersection( s2 );
			return intersection;
		}

		private static Segment[] BaseTileSegments( int p, int q )
		{
			Tiling tiling = new Tiling();
			TilingConfig config = new TilingConfig( p, q, 1 );
			Tile baseTile = Tiling.CreateBaseTile( config );
			//baseTile.Transform( Mobius.Scale( 2 ) );				// Only works in Euclidean case
			return baseTile.Boundary.Segments.ToArray();
		}

		/// <summary>
		/// Mirrors for Spherical geometry, in the ball model.
		/// </summary>
		public static Sphere[] MirrorsSpherical( int p, int q, int r )
		{
			// Get a {q,p} tiling on the z=0 plane.
			Segment[] baseTileSegments = BaseTileSegments( q, p );

			// This will be unit length.
			Vector3D pFaceDirection = H3Models.UHSToBall( baseTileSegments.First().P1 );

			// In-radius is in conformal model
			double inRadius = Spherical2D.s2eNorm( Honeycomb.InRadius( p, q, r ) );
			double centerOfSphereNE = ( 1 - inRadius ) / ( 1 + inRadius );
			Vector3D center;
			double radius;
			H3Models.Ball.DupinCyclideSphere( -pFaceDirection * centerOfSphereNE, 1.0 /*geodesic circle*/, Geometry.Spherical, out center, out radius );		
			Sphere cellBoundary = new Sphere() { Center = center, Radius = radius };
			//cellBoundary = H3Models.BallToUHS( cellBoundary );

			Sphere[] interior = InteriorMirrors( p, q );
			interior = interior.Select( s => H3Models.UHSToBall( s ) ).ToArray();

			Sphere[] surfaces = new Sphere[] { cellBoundary, interior[0], interior[1], interior[2] };
			surfaces[0].Invert = true;

			// Apply rotations.
			bool applyRotations = true;
			if( applyRotations )
			{
				double rotation = Math.PI / 2;
				foreach( Sphere s in surfaces )
					RotateSphere( s, rotation );
			}

			return surfaces;
		}

		public static Vector3D VertexSpherical( int p, int q, int r )
		{
			double circumRadius = Spherical2D.s2eNorm( Honeycomb.CircumRadius( p, q, r ) );
			return new Vector3D( 0, 0, -circumRadius );
		}

		public static Sphere[] MirrorsEuclidean()
		{
			int p = 4;
			int q = 3;
			//int r = 4;

			// Get a {q,p} tiling on the z=0 plane.
			Segment[] baseTileSegments = BaseTileSegments( q, p );

			// This will be unit length.
			Vector3D pFaceDirection = H3Models.UHSToBall( baseTileSegments.First().P1 );

			// NOTES: 
			//	Center is the plane normal.
			//	Cell face is only one with an offset
			//  This is constructed already in the ball.
			Sphere cellBoundary = new Sphere() { Center = -pFaceDirection, Offset = pFaceDirection * m_eScale, Radius = double.PositiveInfinity };

			Sphere[] interior = InteriorMirrors( p, q );
			interior = interior.Select( s => H3Models.UHSToBall( s ) ).ToArray();
			Sphere[] surfaces = new Sphere[] { cellBoundary, interior[0], interior[1], interior[2] };

			// Apply rotations.
			bool applyRotations = false;
			if( applyRotations )
			{
				double rotation = Math.PI / 2;
				foreach( Sphere s in surfaces )
					RotateSphere( s, rotation );
			}

			return surfaces;
		}

		public static Vector3D[] VertsEuclidean()
		{
			int p = 4;
			int q = 3;
			//int r = 4;

			// Get a {q,p} tiling on the z=0 plane.
			Segment[] baseTileSegments = BaseTileSegments( q, p );

			// These will be unit length.
			Vector3D pFaceDirection = H3Models.UHSToBall( baseTileSegments.First().P1 );
			Vector3D pMidEdgeDirection = H3Models.UHSToBall( baseTileSegments.First().Midpoint );
			Vector3D pVertexDirection  = new Vector3D( 0, 0, -1 );

			// Order is same as facets (these are the points opposite facets).
			List<Vector3D> verts = new List<Vector3D>();
			verts.Add( new Vector3D() );
			verts.Add( pFaceDirection * m_eScale );
			verts.Add( pVertexDirection * Math.Sqrt( 3 ) * m_eScale );
			verts.Add( pMidEdgeDirection * Math.Sqrt( 2 ) * m_eScale );

			// Apply rotations.
			double rotation = Math.PI / 2;
			Vector3D zAxis = new Vector3D( 0, 0, 1 );
			for( int i=0; i<4; i++ )
				verts[i].RotateAboutAxis( zAxis, rotation );

			return verts.ToArray();
		}

		public static void CalcEScale()
		{
			// Euclidean scale is arbitrary, but put it in the middle of the projections of 433 and 435.
			double r3 = Spherical2D.s2eNorm( Honeycomb.CircumRadius( 4, 3, 3 ) );
			double r5 = DonHatch.h2eNorm( Honeycomb.CircumRadius( 4, 3, 5 ) );
			m_eScale = ( r3 + r5 ) / ( 2 * Math.Sqrt(3) );
		}

		// Euclidean scale is arbitrary.
		static double m_eScale = 0.5;

		/// <summary>
		/// Return the 4 simplex vertices in the ball model.
		/// </summary>
		public static Vector3D[] VertsBall( int p, int q, int r )
		{
			// Order same as facets (these are the points opposite facets).
			Vector3D cellCenter = CellCenterBall( p, q, r );
			Vector3D faceCenter = FaceCenterBall( p, q, r );
			Vector3D edgeCenter = MidEdgePointBall( p, q, r );
			Vector3D vertexPoint = VertexPointBall( p, q, r );

			return new Vector3D[] { cellCenter, faceCenter, edgeCenter, vertexPoint };
		}

		/// <summary>
		/// Returns the 6 simplex edges in the Ball model.
		/// </summary>
		public static H3.Cell.Edge[] SimplexEdgesBall( int p, int q, int r )
		{
			H3.Cell.Edge[] edges = SimplexEdgesUHS( p, q, r );
			foreach( H3.Cell.Edge e in edges )
			{
				e.Start = H3Models.UHSToBall( e.Start );
				e.End = H3Models.UHSToBall( e.End );
			}
			return edges;
		}

		/// <summary>
		/// Returns the 6 simplex edges in the UHS model.
		/// </summary>
		public static H3.Cell.Edge[] SimplexEdgesUHS( int p, int q, int r )
		{
			// Only implemented for honeycombs with both hyperideal edges/vertices right now.
			if( !( Geometry2D.GetGeometry( p, q ) == Geometry.Hyperbolic &&
					Geometry2D.GetGeometry( q, r ) == Geometry.Hyperbolic ) )
				throw new System.NotImplementedException();

			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: false );

			Circle[] circles = simplex.Select( s => H3Models.UHS.IdealCircle( s ) ).ToArray();

			Vector3D[] defPoints = new Vector3D[6];
			Vector3D dummy;
			Euclidean2D.IntersectionLineCircle( circles[1].P1, circles[1].P2, circles[0], out defPoints[0], out dummy );
			Euclidean2D.IntersectionLineCircle( circles[2].P1, circles[2].P2, circles[0], out defPoints[1], out dummy );
			Euclidean2D.IntersectionLineCircle( circles[1].P1, circles[1].P2, circles[3], out defPoints[2], out dummy );
			Euclidean2D.IntersectionLineCircle( circles[2].P1, circles[2].P2, circles[3], out defPoints[3], out dummy );

			Circle3D c = simplex[0].Intersection( simplex[3] );

			Vector3D normal = c.Normal;
			normal.RotateXY( Math.PI / 2 );
			Vector3D intersection;
			double height, off;

			Euclidean2D.IntersectionLineLine( c.Center, c.Center + normal, circles[1].P1, circles[1].P2, out intersection );
			off = ( intersection - c.Center ).Abs();
			height = Math.Sqrt( c.Radius * c.Radius - off * off );
			intersection.Z = height;
			defPoints[4] = intersection;

			Euclidean2D.IntersectionLineLine( c.Center, c.Center + normal, circles[2].P1, circles[2].P2, out intersection );
			off = ( intersection - c.Center ).Abs();
			height = Math.Sqrt( c.Radius * c.Radius - off * off );
			intersection.Z = height;
			defPoints[5] = intersection;

			bool order = false;
			H3.Cell.Edge[] edges = new H3.Cell.Edge[]
			{
				new H3.Cell.Edge( new Vector3D(), new Vector3D( 0, 0, 10 ) ),
				new H3.Cell.Edge( defPoints[4], defPoints[5], order ),
				new H3.Cell.Edge( defPoints[0], defPoints[4], order ),
				new H3.Cell.Edge( defPoints[1], defPoints[5], order ),
				new H3.Cell.Edge( defPoints[2], defPoints[4], order ),
				new H3.Cell.Edge( defPoints[3], defPoints[5], order ),
			};

			return edges;
		}

		/// <summary>
		/// Helper to construct some points we need for calculating simplex facets for a {p,q,r} honeycomb.
		/// </summary>
		private static void TilePoints( int p, int q, out Vector3D p1, out Vector3D p2, out Vector3D p3, out Segment seg )
		{
			if( Infinite( p ) && Infinite( q ) /*&& FiniteOrInfinite( r )*/ )
			{
				p1 = new Vector3D( 1, 0, 0 );
				p2 = new Vector3D( 0, Math.Sqrt( 2 ) - 1 );
				p3 = Vector3D.DneVector();

				Circle3D arcCircle;
				H3Models.Ball.OrthogonalCircleInterior( p2, p1, out arcCircle );
				seg = Segment.Arc( p1, p2, arcCircle.Center, clockwise: true );
			}
			else
			{
				Segment[] baseTileSegments;
				if( Infinite( q ) )
					baseTileSegments = BaseTileSegments( p, q );	// Can't use dual here.
				else
					baseTileSegments = BaseTileSegments( q, p );	// Intentionally using dual.

				seg = baseTileSegments.First();

				p1 = seg.P1;
				p2 = seg.Midpoint;
				p3 = p2;
				p3.RotateXY( -Math.PI / 2 );
			}
		}

		public static Sphere[] Mirrors( int p, int q, int r, bool moveToBall = true, double scaling = -1 )
		{
			Vector3D dummy = new Vector3D();
			return Mirrors( p, q, r, ref dummy, moveToBall, scaling );
		}

		public static Sphere[] Mirrors( int p, int q, int r, ref Vector3D cellCenter, bool moveToBall = true, double scaling = -1 )
		{
			Geometry g = Util.GetGeometry( p, q, r );
			if( g == Geometry.Spherical )
				return SimplexCalcs.MirrorsSpherical( p, q, r );
			else if( g == Geometry.Euclidean )
				return SimplexCalcs.MirrorsEuclidean();

			// This is a rotation we'll apply to the mirrors at the end.
			// This is to try to make our image outputs have vertical bi-lateral symmetry and the most consistent in all cases.
			// NOTE: + is CW, not CCW. (Because the way I did things, our images have been reflected vertically, and I'm too lazy to go change this.)
			double rotation = Math.PI / 2;

			// Some construction points we need.
			Vector3D p1, p2, p3;
			Segment seg = null;
			TilePoints( p, q, out p1, out p2, out p3, out seg );

			//
			// Construct in UHS
			//

			Geometry cellGeometry = Geometry2D.GetGeometry( p, q );

			Vector3D center = new Vector3D();
			double radius = 0;
			if( cellGeometry == Geometry.Spherical )
			{
				// Finite or Infinite r

				// Spherical trig
				double halfSide = Geometry2D.GetTrianglePSide( q, p );
				double mag = Math.Sin( halfSide ) / Math.Cos( Util.PiOverNSafe( r ) );
				mag = Math.Asin( mag );

				// e.g. 43j
				//mag *= 0.95;

				// Move mag to p1.
				mag = Spherical2D.s2eNorm( mag );
				H3Models.Ball.DupinCyclideSphere( p1, mag, Geometry.Spherical, out center, out radius );
			}
			else if( cellGeometry == Geometry.Euclidean )
			{
				center = p1;
				radius = p1.Dist( p2 ) / Math.Cos( Util.PiOverNSafe( r ) );
			}
			else if( cellGeometry == Geometry.Hyperbolic )
			{
				if( Infinite( p ) && Infinite( q ) && FiniteOrInfinite( r ) )
				{
					//double iiiCellRadius = 2 - Math.Sqrt( 2 );
					//Circle3D iiiCircle = new Circle3D() { Center = new Vector3D( 1 - iiiCellRadius, 0, 0 ), Radius = iiiCellRadius };
					//radius = iiiCellRadius;	// infinite r
					//center = new Vector3D( 1 - radius, 0, 0 );

					// For finite r, it was easier to calculate cell facet in a more symmetric position, 
					// then move into position with the other mirrors via a Mobius transformation.
					double rTemp = 1 / ( Math.Cos( Util.PiOverNSafe( r ) ) + 1 );
					Mobius m = new Mobius();
					m.Isometry( Geometry.Hyperbolic, -Math.PI / 4, new Vector3D( 0, Math.Sqrt( 2 ) - 1 ) );
					Vector3D c1 = m.Apply( new Vector3D( 1 - 2 * rTemp, 0, 0 ) );
					Vector3D c2 = c1;
					c2.Y *= -1;
					Vector3D c3 = new Vector3D( 1, 0 );
					Circle3D c = new Circle3D( c1, c2, c3 );

					radius = c.Radius;
					center = c.Center;
				}
				else if( Infinite( p ) && Finite( q ) && FiniteOrInfinite( r ) )
				{
					// http://www.wolframalpha.com/input/?i=r%2Bx+%3D+1%2C+sin%28pi%2Fp%29+%3D+r%2Fx%2C+solve+for+r
					// radius = 2 * Math.Sqrt( 3 ) - 3;	// Appolonian gasket wiki page
					//radius = Math.Sin( Math.PI / q ) / ( Math.Sin( Math.PI / q ) + 1 );
					//center = new Vector3D( 1 - radius, 0, 0 );

					// For finite r, it was easier to calculate cell facet in a more symmetric position, 
					// then move into position with the other mirrors via a Mobius transformation.
					double rTemp = 1 / ( Math.Cos( Util.PiOverNSafe( r ) ) + 1 );
					Mobius m = new Mobius();
					m.Isometry( Geometry.Hyperbolic, 0, p2 );
					Vector3D findingAngle = m.Inverse().Apply( new Vector3D( 1, 0 ) );
					double angle = Math.Atan2( findingAngle.Y, findingAngle.X );

					m.Isometry( Geometry.Hyperbolic, angle, p2 );
					Vector3D c1 = m.Apply( new Vector3D( 1 - 2 * rTemp, 0, 0 ) );
					Vector3D c2 = c1;
					c2.Y *= -1;
					Vector3D c3 = new Vector3D( 1, 0 );
					Circle3D c = new Circle3D( c1, c2, c3 );

					radius = c.Radius;
					center = c.Center;
				}
				else if( Finite( p ) && Infinite( q ) && FiniteOrInfinite( r ) )
				{
					radius = p2.Abs(); // infinite r
					radius = DonHatch.asinh( Math.Sinh( DonHatch.e2hNorm( p2.Abs() ) ) / Math.Cos( Util.PiOverNSafe( r ) ) );	// hyperbolic trig

					// 4j3
					//m_jOffset = radius * 0.02;
					//radius += m_jOffset ;

					radius = DonHatch.h2eNorm( radius );
					center = new Vector3D();
					rotation *= -1;
				}
				else if( /*Finite( p ) &&*/ Finite( q ) )
				{
					// Infinite r
					//double mag = Geometry2D.GetTrianglePSide( q, p );	

					// Finite or Infinite r
					double halfSide = Geometry2D.GetTrianglePSide( q, p );
					double mag = DonHatch.asinh( Math.Sinh( halfSide ) / Math.Cos( Util.PiOverNSafe( r ) ) );	// hyperbolic trig
					H3Models.Ball.DupinCyclideSphere( p1, DonHatch.h2eNorm( mag ), out center, out radius );
				}
				else
					throw new System.NotImplementedException();
			}
			Sphere cellBoundary = new Sphere()
			{
				Center = center,
				Radius = radius
			};

			Sphere[] interior = InteriorMirrors( p, q );
			Sphere[] surfaces = new Sphere[] { cellBoundary, interior[0], interior[1], interior[2] };

			// Apply rotations.
			bool applyRotations = true;
			if( applyRotations )
			{
				foreach( Sphere s in surfaces )
					RotateSphere( s, rotation );
				p1.RotateXY( rotation );
			}

			// Apply scaling
			bool applyScaling = scaling != -1;
			if( applyScaling )
			{
				//double scale = 1.0/0.34390660467269524;
				//scale = 0.58643550768408892;
				foreach( Sphere s in surfaces )
					Sphere.ScaleSphere( s, scaling );
			}

			bool facetCentered = false;
			if( facetCentered )
			{
				PrepForFacetCentering( p, q, surfaces, ref cellCenter );
			}

			// Move to ball if needed.
			Sphere[] result = surfaces.Select( s => moveToBall ? H3Models.UHSToBall( s ) : s ).ToArray();
			cellCenter = moveToBall ? H3Models.UHSToBall( cellCenter ) : cellCenter;

			return result;
		}

		// This was used for making some images "beyond infinity" (paper appendix).
		//private static double m_jOffset;

		/// <summary>
		/// Inputs must be in UHS!
		/// </summary>
		public static void PrepForFacetCentering( int p, int q, Sphere[] spheres, ref Vector3D cellCenter )
		{
			Mobius m = FCOrientMobius( p, q );
			foreach( Sphere s in spheres )
				H3Models.TransformInUHS2( s, m );

			spheres[3].Center *= -1; // Super-hack for 437 & 737, since TransformInUHS2 function is falling short.
			//spheres[2].Center *= -1;

			cellCenter = m.ApplyToQuaternion( cellCenter );
		}

		public static Mobius FCOrientMobius( int p, int q )
		{
			Vector3D p1, p2, p3;
			Segment seg = null;
			TilePoints( p, q, out p1, out p2, out p3, out seg );
			Geometry cellGeometry = Geometry2D.GetGeometry( p, q );

			p1.RotateXY( Math.PI / 2 );	// XXX - repeated from above.  Implement a better way to sequence transformations.

			Mobius m = new Mobius();
			m.Isometry( cellGeometry, 0, -p1 );
			return m;
		}

		/// <summary>
		/// Calculates the 3 mirrors connected to the cell center.
		/// This works in all geometries and returns results in the UHS model (or the appropriate analogue).
		/// </summary>
		private static Sphere[] InteriorMirrors( int p, int q )
		{
			// Some construction points we need.
			Vector3D p1, p2, p3;
			Segment seg = null;
			TilePoints( p, q, out p1, out p2, out p3, out seg );

			Geometry cellGeometry = Geometry2D.GetGeometry( p, q );

			// XZ-plane
			Sphere s1 = new Sphere()
			{
				Center = new Vector3D( 0, 1, 0 ),
				Radius = double.PositiveInfinity
			};

			Sphere s2 = null;
			if( cellGeometry == Geometry.Euclidean )
			{
				s2 = new Sphere()
				{
					Center = -p2,
					Offset = p2,
					Radius = double.PositiveInfinity
				};
			}
			else if(
				cellGeometry == Geometry.Spherical ||
				cellGeometry == Geometry.Hyperbolic )
			{
				s2 = new Sphere()
				{
					Center = seg.Center,
					Radius = seg.Radius,
					//Invert = true
				};

				// j34
				/*double off = seg.Center.Abs() - seg.Radius;
				off *= 1.05;
				Vector3D vOff = seg.Center;
				vOff.Normalize();
				vOff *= off;
				s2 = H3Models.Ball.OrthogonalSphereInterior( vOff );*/

				// 4j3
				/*double off = seg.Center.Abs() - seg.Radius;
				off = DonHatch.h2eNorm( DonHatch.e2hNorm( off ) + m_jOffset );
				Vector3D vOff = seg.Center;
				vOff.Normalize();
				vOff *= off;
				s2 = H3Models.Ball.OrthogonalSphereInterior( vOff );*/
			}

			Sphere s3;
			if( Infinite( p ) && Infinite( q ) /*&& FiniteOrInfinite( r )*/ )
			{
				Vector3D tempCenter = seg.Center;
				tempCenter.X *= -1;
				s3 = new Sphere()
				{
					Center = tempCenter,
					Radius = seg.Radius,
					//Invert = true,
				};
			}
			else
			{
				s3 = new Sphere()
				{
					Center = p3,
					Radius = double.PositiveInfinity
				};
			}

			// The reason for the special case ordering is to make sure
			// the last mirror is always the cell-reflecting mirror of the dual honeycomb.
			// The order we want are the mirrors opposite: cell, face, edge, vertex
			// NOTE: This also puts the mirrors in the same order as in standard presentations, 
			// e.g. from Coxeter's 57-cell paper.
			Sphere[] surfaces;
			if( !Infinite( p ) && Infinite( q ) )
				surfaces = new Sphere[] { s2, s1, s3 };
			else
				surfaces = new Sphere[] { s3, s1, s2 };
			
			return surfaces;
		}

		/// <summary>
		/// Helper to rotate a sphere about the z axis.
		/// </summary>
		internal static void RotateSphere( Sphere s, double rotation )
		{
			Vector3D zAxis = new Vector3D( 0, 0, 1 );
			Sphere.RotateSphere( s, zAxis, rotation );
		}

		private static bool Infinite( int t )
		{
			return t == -1;
		}

		private static bool Finite( int t )
		{
			return !Infinite( t );
		}

		private static bool FiniteOrInfinite( int t )
		{
			return true;
		}

		/// <summary>
		/// This calculates the 4 vertices of a general (but finite) Goursat Tetrahedron.  Result is in the ball model.
		/// 
		/// The method comes from the dissertation "Hyperbolic polyhedra: volume and scissors congruence", 
		/// by Yana Zilberberg Mohanty, section 2.4, steps 1-5.
		/// 
		/// A,B,C are the three dihedral angles surrounding a vertex.
		/// A_,B_,C_ are the three oppoite dihedral angles.
		/// </summary>
		public static Vector3D[] GoursatTetrahedron( double A, double B, double C, double A_, double B_, double C_ )
		{
			// Step 1: Construct Gram matrix with reversed rows/columns.
			// NOTE: The sign of the diagonal in the paper was incorrect.
			double pi = Math.PI;
			double[,] gramMatrixData = new double[,] 
			{ 
				{ -1,				Math.Cos(pi/A_),	Math.Cos(pi/B_),	Math.Cos(pi/C) }, 
				{ Math.Cos(pi/A_),	-1,					Math.Cos(pi/C_),	Math.Cos(pi/B) }, 
				{ Math.Cos(pi/B_),	Math.Cos(pi/C_),	-1,					Math.Cos(pi/A) }, 
				{ Math.Cos(pi/C),	Math.Cos(pi/B),		Math.Cos(pi/A),		-1 }, 
			};
			Matrix4D gramMatrix = new Matrix4D( gramMatrixData );
			gramMatrix *= -1;
			Matrix4D identity = Matrix4D.Identity();

			// Step 2: Gram-Schmidt.
			Matrix4D W = GramSchmidt( identity, gramMatrix );

			// Step 3: Divide 4th row by i (already effectively done in our Gram-Schmidt routine below), and reverse order of rows.
			Matrix4D W_ = ReverseRows( W );

			// Step 4
			Matrix4D D = identity.Clone();
			D[0, 0] = -1;
			Matrix4D U = Matrix4D.Transpose( D * W_ );

			// Step 5
			Matrix4D U_ = ReverseRows( U );
			for( int i=0; i<4; i++ )
				MinkowskiNormalize( U_[i] );

			// Now move from the hyperboloid model to the ball.
			List<Vector3D> result = new List<Vector3D>();
			for( int i=0; i<4; i++ )
				result.Add( HyperboloidToBall( U_[i] ) );
			return result.ToArray();
		}

		public static Matrix4D GramSchmidt( Matrix4D input, Matrix4D innerProductValues )
		{
			Matrix4D result = input.Clone();
			for( int i = 0; i < 4; i++ )
			{
				for( int j = 0; j < i; j++ )
				{
					VectorND iVec = result[i];
					VectorND jVec = result[j];
					double inner = innerProductValues[i].Dot( jVec );
					iVec -= inner * jVec;
					result[i] = iVec;
				}

				// Normalize.  We don't use VectorND normalize because we might have timelike vectors.
				double mag2 = innerProductValues[i].Dot( result[i] );
				double abs = mag2 < 0 ? -Math.Sqrt( -mag2 ) : Math.Sqrt( mag2 );
				result[i].Divide( abs );
			}

			return result;
		}

		public static Matrix4D ReverseRows( Matrix4D m )
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
				result[i] = m[3-i];
			return result;
		}

		// Minkowski inner product.
		private static double MinkowskiInnerProduct( VectorND v1, VectorND v2 )
		{
			double inner = -v1.X[0] * v2.X[0];
			for( int i = 1; i < v1.Dimension; i++ )
				inner += v1.X[i] * v2.X[i];
			return inner;
		}

		// Minkowski normalization.
		private static VectorND MinkowskiNormalize( VectorND v )
		{
			double mag2 = MinkowskiInnerProduct( v, v );
			double abs = mag2 < 0 ? Math.Sqrt( -mag2 ) : Math.Sqrt( mag2 );
			v.Divide( abs );
			return v;
		}

		public static Vector3D HyperboloidToBall( VectorND hyperboloidPoint )
		{	
			double t = hyperboloidPoint.X[0];
			return new Vector3D(
				hyperboloidPoint.X[1] / ( 1 + t ),
				hyperboloidPoint.X[2] / ( 1 + t ),
				hyperboloidPoint.X[3] / ( 1 + t ) );
		}

		/// <summary>
		/// Given the 4 verts of a tetrahedron (must lie within ball),
		/// Calculate the faces of the tetrahedron.
		/// Input and Output in ball model.
		/// </summary>
		public static Sphere[] Mirrors( Vector3D[] verts )
		{
			// Order so we get faces opposite verts.
			Sphere s1 = H3Models.Ball.OrthogonalSphereInterior( verts[1], verts[2], verts[3] );
			Sphere s2 = H3Models.Ball.OrthogonalSphereInterior( verts[0], verts[2], verts[3] );
			Sphere s3 = H3Models.Ball.OrthogonalSphereInterior( verts[0], verts[1], verts[3] );
			Sphere s4 = H3Models.Ball.OrthogonalSphereInterior( verts[0], verts[1], verts[2] );
			return new Sphere[] { s1, s2, s3, s4 };
		}
	}
}
