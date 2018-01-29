namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	public class S3
	{
		private static void HopfFibration( Tiling tiling )
		{
			int segDivisions = 10;
			Shapeways mesh = new Shapeways();

			HashSet<Vector3D> done = new HashSet<Vector3D>();
			foreach( Tile tile in tiling.Tiles )
			foreach( Segment seg in tile.Boundary.Segments )
			{
				if( done.Contains( seg.Midpoint ) )
					continue;

				// Subdivide the segment, and project points to S2.
				Vector3D[] points = seg.Subdivide( segDivisions ).Select( v => Spherical2D.PlaneToSphere( v ) ).ToArray();
				foreach( Vector3D point in points )
				{
					Vector3D[] circlePoints = OneHopfCircle( point );
					ProjectAndAddS3Points( mesh, circlePoints, shrink: false );
				}

				done.Add( seg.Midpoint );
			}

			STL.SaveMeshToSTL( mesh.Mesh, @"D:\p4\R3\sample\out1.stl" );
		}

		public static void HopfOrbit()
		{
			List<Vector3D> s2Points = new List<Vector3D>();
			for( double theta = Math.PI * .1; theta <= Math.PI * .9; theta += Math.PI *.2 )
			for( double lon = -Math.PI; lon <= Math.PI; lon += Math.PI/10 )
			{
				s2Points.Add( SphericalCoords.SphericalToCartesian( new Vector3D( 1.0, theta, lon ) ) );
			}

			using( StreamWriter sw = File.CreateText( @".\out.pov" ) )
			{
				System.Func<Vector3D, Sphere> sizeFunc = v => new Sphere() { Center = v, Radius = 0.01 };
				foreach( Vector3D s2Point in s2Points )
				{
					Vector3D[] circlePoints = OneHopfCircle( s2Point );

					//for( int i = 0; i < circlePoints.Length; i++ )
					//	circlePoints[i] = circlePoints[i].ProjectTo3DSafe( 1.0 );

					// Note: effectively orthogonal projects here because EdgeSphereSweep doesn't write W coord.
					string circleString = PovRay.EdgeSphereSweep( circlePoints, sizeFunc );
					sw.WriteLine( circleString );
				}
			}
		}

		public static Vector3D[] OneHopfCircle( Vector3D s2Point, bool anti = false )
		{
			int circleDivisions = 125;

			// Get the hopf circle.
			// http://en.wikipedia.org/wiki/Hopf_fibration#Explicit_formulae
			double a = s2Point.X;
			double b = s2Point.Y;
			double c = s2Point.Z;
			double factor = 1 / ( Math.Sqrt( 1 + c ) );
			if( Tolerance.Equal( c, -1 ) )
				return new Vector3D[] {};

			List<Vector3D> circlePoints = new List<Vector3D>();
			double angleInc = 2 * Math.PI / circleDivisions;
			double angle = 0;
			for( int i = 0; i <= circleDivisions; i++ )
			{
				double sinTheta = Math.Sin( angle );
				double cosTheta = Math.Cos( angle );
				Vector3D point = new Vector3D(
					( 1 + c ) * cosTheta,
					anti ? - a * sinTheta - b * cosTheta : a * sinTheta - b * cosTheta,
					anti ?   a * cosTheta - b * sinTheta : a * cosTheta + b * sinTheta,
					( 1 + c ) * sinTheta );
				point.Normalize();
				circlePoints.Add( point );

				angle += angleInc;
			}

			return circlePoints.ToArray();
		}

		private static void ShapewaysPolytopes()
		{
			VEF loader = new VEF();
			loader.Load( @"C:\Users\roice\Documents\projects\vZome\VefProjector\data\24cell-cellFirst.vef" );

			int divisions = 25;

			Shapeways mesh = new Shapeways();
			//int count = 0;
			foreach( GraphEdge edge in loader.Edges )
			{
				Segment seg = Segment.Line(
					loader.Vertices[edge.V1].ConvertToReal(),
					loader.Vertices[edge.V2].ConvertToReal() );
				Vector3D[] points = seg.Subdivide( divisions );

				bool shrink = true;
				ProjectAndAddS3Points( mesh, points, shrink );

				//if( count++ > 10 )
				//	break;
			}

			STL.SaveMeshToSTL( mesh.Mesh, @"D:\p4\R3\sample\out1.stl" );
		}

		public static void EdgesToStl( H3.Cell.Edge[] edges )
		{
			Shapeways mesh = new Shapeways();
			
			int divisions = 25;
			foreach( H3.Cell.Edge edge in edges )
			{
				Segment seg = Segment.Line( 
					Sterographic.R3toS3( edge.Start ), 
					Sterographic.R3toS3( edge.End ) );
				Vector3D[] points = seg.Subdivide( divisions );

				ProjectAndAddS3Points( mesh, points );
			}

			for( int i = 0; i < mesh.Mesh.Triangles.Count; i++ )
			{
				mesh.Mesh.Triangles[i] = new Mesh.Triangle(
					SphericalModels.StereoToEqualVolume( mesh.Mesh.Triangles[i].a ),
					SphericalModels.StereoToEqualVolume( mesh.Mesh.Triangles[i].b ),
					SphericalModels.StereoToEqualVolume( mesh.Mesh.Triangles[i].c ) );
			}

			STL.SaveMeshToSTL( mesh.Mesh, @"output.stl" );
		}

		private static void ProjectAndAddS3Points( Shapeways mesh, Vector3D[] pointsS3 )
		{
			double r = 0.02;

			List<Vector3D> projected = new List<Vector3D>();
			List<double> radii = new List<double>();
			foreach( Vector3D v in pointsS3 )
			{
				v.Normalize();
				Vector3D c = v.ProjectTo3DSafe( 1.0 );

				Vector3D p;
				double d;
				H3Models.Ball.DupinCyclideSphere( c, r, Geometry.Spherical, out p, out d );
				projected.Add( p );
				radii.Add( d );
			}

			mesh.AddCurve( projected.ToArray(), radii.ToArray() );
		}

		/// <summary>
		/// Helper to project points from S3 -> S2, then add an associated curve.
		/// XXX - Not completely correct.
		/// </summary>
		private static void ProjectAndAddS3Points( Shapeways mesh, Vector3D[] pointsS3, bool shrink )
		{
			List<Vector3D> projected = new List<Vector3D>();
			foreach( Vector3D v in pointsS3 )
			{
				v.Normalize();
				Vector3D c = v.ProjectTo3DSafe( 1.0 );

				// Pull R3 into a smaller open disk.
				if( shrink )
				{
					double mag = Math.Atan( c.Abs() );
					c.Normalize();
					c *= mag;
				}

				projected.Add( c );
			}

			System.Func<Vector3D, double> sizeFunc = v =>
			{
				// Constant thickness.
				// return 0.08;

				double sphericalThickness = 0.05;

				double abs = v.Abs();
				if( shrink )
					abs = Math.Tan( abs );	// The unshrunk abs.

				// The thickness at this vector location.
				double result = Spherical2D.s2eNorm( Spherical2D.e2sNorm( abs ) + sphericalThickness ) - abs;

				if( shrink )
					result *= Math.Atan( abs ) / abs;	// shrink it back down.

				return result;
			};

			mesh.AddCurve( projected.ToArray(), sizeFunc );
		}

		public static void Hypercube()
		{
			List<Vector3D> vertices = new List<Vector3D>();
			vertices.Add( new Vector3D(  1,  1,  1,  1 ) );
			vertices.Add( new Vector3D(  1,  1,  1, -1 ) );
			vertices.Add( new Vector3D(  1,  1, -1,  1 ) );
			vertices.Add( new Vector3D(  1,  1, -1, -1 ) );
			vertices.Add( new Vector3D(  1, -1,  1,  1 ) );
			vertices.Add( new Vector3D(  1, -1,  1, -1 ) );
			vertices.Add( new Vector3D(  1, -1, -1,  1 ) );
			vertices.Add( new Vector3D(  1, -1, -1, -1 ) );
			vertices.Add( new Vector3D( -1,  1,  1,  1 ) );
			vertices.Add( new Vector3D( -1,  1,  1, -1 ) );
			vertices.Add( new Vector3D( -1,  1, -1,  1 ) );
			vertices.Add( new Vector3D( -1,  1, -1, -1 ) );
			vertices.Add( new Vector3D( -1, -1,  1,  1 ) );
			vertices.Add( new Vector3D( -1, -1,  1, -1 ) );
			vertices.Add( new Vector3D( -1, -1, -1,  1 ) );
			vertices.Add( new Vector3D( -1, -1, -1, -1 ) );

			HashSet<H3.Cell.Edge> edges = new HashSet<H3.Cell.Edge>( new H3.Cell.EdgeEqualityComparer() );
			foreach( Vector3D v1 in vertices )
			foreach( Vector3D v2 in vertices )
				if( v1.Dist( v2 ) == 2 )
					edges.Add( new H3.Cell.Edge( v1, v2 ) );

			// Radial project to S3, then stereographic to R3.
			foreach( H3.Cell.Edge edge in edges )
			{
				edge.Start.Normalize();
				edge.Start = Sterographic.S3toR3( edge.Start );
				edge.End.Normalize();
				edge.End = Sterographic.S3toR3( edge.End );
			}

			PovRay.WriteEdges( new PovRay.Parameters() { AngularThickness = .05 }, Geometry.Spherical, edges.ToArray(), "433.pov", append: false ); 
		}

		/// <summary>
		/// Inputs and Outputs are in R3 (stereographically projected).
		/// </summary>
		public static Vector3D[] GeodesicPoints( Vector3D v1, Vector3D v2 )
		{
			Vector3D start = Sterographic.R3toS3( v1 );
			Vector3D end = Sterographic.R3toS3( v2 );
			AvoidNorthPole( ref start, end );
			AvoidNorthPole( ref end, start );

			int div = 42;
			//int div = 56;		// 343
			//int div = 50;		// 333
			Segment seg = Segment.Line( start, end );
			Vector3D[] result = seg.Subdivide( div );
			for( int i=0; i<result.Length; i++ )
			{
				result[i].Normalize();
				result[i] = Sterographic.S3toR3( result[i] );
			}

			return result;
		}

		private static void AvoidNorthPole( ref Vector3D v, Vector3D direction )
		{
			if( !Tolerance.Equal( v.W, 1 ) )
				return;

			Vector3D cutEnd = v - direction;
			double abs = cutEnd.Abs();
			abs -= 0.35;
			cutEnd.Normalize();
			cutEnd *= abs;
			v = direction + cutEnd;
			v.Normalize();
		}
	}
}