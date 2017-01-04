namespace R3.Core
{
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;

	public class VEF
	{
		/// <summary>
		/// The vertices.
		/// </summary>
		public GoldenVector4D[] Vertices { get; private set; }

		/// <summary>
		/// The edges.
		/// </summary>
		public GraphEdge[] Edges { get; private set; }

		/// <summary>
		/// Load from a VEF file.
		/// </summary>
		/// <param name="fileName"></param>
		public void Load( string fileName )
		{
			IEnumerator<string> lines = File.ReadLines( fileName ).GetEnumerator();

			lines.MoveNext();
			int numVertices = int.Parse( lines.Current );

			List<GoldenVector4D> vertices = new List<GoldenVector4D>();
			for( int i = 0; i < numVertices; i++ )
			{
				lines.MoveNext();
				GoldenVector4D v = new GoldenVector4D();
				v.ReadVector( lines.Current );
				vertices.Add( v );
			}

			lines.MoveNext();
			int numEdges = int.Parse( lines.Current );
			List<GraphEdge> edges = new List<GraphEdge>();
			for( int i = 0; i < numEdges; i++ )
			{
				lines.MoveNext();
				GraphEdge e = new GraphEdge();
				e.ReadEdge( lines.Current );
				edges.Add( e );
			}

			this.Vertices = vertices.ToArray();
			this.Edges = edges.ToArray();
		}
	}
}
