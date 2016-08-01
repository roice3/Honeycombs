namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	// NOTE: Wanted to name this class R3 (parallel to H3/S3), but namespace problems happened.
	public class Euclidean
	{
		public static Vector3D[] GeodesicPoints( Vector3D v1, Vector3D v2 )
		{
			int div = 4;
			Segment seg = Segment.Line( v1, v2 );
			Vector3D[] result = seg.Subdivide( div );
			return result;
		}

		public static void GenEuclidean()
		{
			Shapeways mesh = new Shapeways();
			HashSet<H3.Cell.Edge> completed = new HashSet<H3.Cell.Edge>( new H3.Cell.EdgeEqualityComparer() );

			int count = 20;
			for( int i = -count; i < count; i++ )
			for( int j = -count; j < count; j++ )
			for( int k = -count; k < count; k++ )
			{
				// Offset
				double io = i + 0.5;
				double jo = j + 0.5;
				double ko = k + 0.5;

				// Do every edge emanating from this point.
				AddEuclideanEdge( completed, new Vector3D( io, jo, ko ), new Vector3D( io + 1, jo, ko ) );
				AddEuclideanEdge( completed, new Vector3D( io, jo, ko ), new Vector3D( io - 1, jo, ko ) );
				AddEuclideanEdge( completed, new Vector3D( io, jo, ko ), new Vector3D( io, jo + 1, ko ) );
				AddEuclideanEdge( completed, new Vector3D( io, jo, ko ), new Vector3D( io, jo - 1, ko ) );
				AddEuclideanEdge( completed, new Vector3D( io, jo, ko ), new Vector3D( io, jo, ko + 1 ) );
				AddEuclideanEdge( completed, new Vector3D( io, jo, ko ), new Vector3D( io, jo, ko - 1 ) );
			}

			//STL.SaveMeshToSTL( mesh.Mesh, "d:\\temp\\434.stl" );
			PovRay.WriteEdges( new PovRay.Parameters() { AngularThickness = .05 }, Geometry.Euclidean, completed.ToArray(), "434.pov", append: false ); 
		}

		private static void AddEuclideanEdge( HashSet<H3.Cell.Edge> completed, Vector3D start, Vector3D end )
		{
			completed.Add( new H3.Cell.Edge( start, end ) );
		}

		private static void AddEuclideanEdgeToMesh( Shapeways mesh, HashSet<H3.Cell.Edge> completed, Vector3D start, Vector3D end )
		{
			H3.Cell.Edge edge = new H3.Cell.Edge( start, end );
			if( completed.Contains( edge ) )
				return;

			Shapeways tempMesh = new Shapeways();
			Segment seg = Segment.Line( start, end );

			int div = 20 - (int)(start.Abs() * 4);
			if( div < 1 )
				div = 1;

			tempMesh.AddCurve( seg.Subdivide( div ), .05 );
			Transform( tempMesh.Mesh );

			mesh.Mesh.Triangles.AddRange( tempMesh.Mesh.Triangles );
			completed.Add( edge );
		}

		private static void Transform( Mesh mesh )
		{
			Sphere sphere = new Sphere();
			sphere.Radius = 0.1;
			for( int i = 0; i < mesh.Triangles.Count; i++ )
			{
				mesh.Triangles[i] = new Mesh.Triangle(
					sphere.ReflectPoint( mesh.Triangles[i].a ),
					sphere.ReflectPoint( mesh.Triangles[i].b ),
					sphere.ReflectPoint( mesh.Triangles[i].c ) );
			}
		}
	}
}
