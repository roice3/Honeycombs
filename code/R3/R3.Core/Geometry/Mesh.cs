namespace R3.Geometry
{
	using R3.Drawing;
	using R3.Math;
	using System.Collections.Generic;
	using System.Linq;
	using System.Diagnostics;

	public class Mesh
	{
		public Mesh()
		{
			Triangles = new List<Triangle>();
		}

		public struct Triangle
		{
			public Triangle( Vector3D _a, Vector3D _b, Vector3D _c ) { a = _a; b = _b; c = _c; color = new Vector3D(1,1,1); }
			public Vector3D a;
			public Vector3D b;
			public Vector3D c;

			// The reason we use a vector here is so the components 
			// can be interpreted in different color schemes (HLS, RGB, etc.)
			public Vector3D color;

			public Vector3D Normal
			{
				get
				{
					return (b - a).Cross( c - a );
				}
			}
		}

		public List<Triangle> Triangles { get; set; }

		public Mesh Clone()
		{
			Mesh clone = new Mesh();
			clone.Triangles = Triangles.Select( t => t ).ToList();
			return clone;
		}
	}
}
