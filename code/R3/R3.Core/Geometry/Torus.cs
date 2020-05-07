namespace R3.Geometry
{
	using Math = System.Math;
	using System.Diagnostics;
	using System.Drawing;

	/// <summary>
	/// Class to generate tori on a 3-sphere
	/// </summary>
	public class Torus
	{
		/// <summary>
		/// The things that define us.
		/// </summary>
		public class Parameters
		{
			public Parameters()
			{
				NumSegments1 = NumSegments2 = 50;
				TubeRadius1 = 0.5;
				Radius = 1.0;
			}

			/// <summary>
			/// The number of segments to generate in the first direction of the torus surface.
			/// </summary>
			public int NumSegments1 { get; set; }

			/// <summary>
			/// The number of segments to generate in the second direction of the torus surface.
			/// </summary>
			public int NumSegments2 { get; set; }

			/// <summary>
			/// The first tube radius of our torus.  
			/// NOTES: 
			///		- The second tube radius is determined by this and the 3-sphere radius.
			///		- This radius must be less than or equal the 3-sphere radius
			///		- If equal 0 or equal to the 3-sphere radius, one tube will be empty (torus will be a line).
			/// </summary>
			public double TubeRadius1 { get; set; }

			/// <summary>
			/// The radius of our 3-sphere
			/// </summary>
			public double Radius { get; set; }
		}

		public Parameters Params { get; set; }

		/// <summary>
		/// Our vertices.
		/// NOTE: Not realy a Vector3D here (need to rename my class).
		/// </summary>
		public Vector3D[][] Vertices { get; set; }

		/// <summary>
		/// Size our Vertices matrix.
		/// </summary>
		private void InitVerts()
		{
			int n1 = this.Params.NumSegments1;
			int n2 = this.Params.NumSegments2;

			Vertices = new Vector3D[n1][];
			for( int i = 0; i < n1; i++ )
				Vertices[i] = new Vector3D[n2];
		}

		/// <summary>
		/// Special case of CreateTorus for the Clifford Torus.
		/// </summary>
		public static Torus CreateClifford( Parameters parameters )
		{
			parameters.TubeRadius1 = parameters.Radius / 2;
			return CreateTorus( parameters );
		}

		/// <summary>
		/// Calculates a torus which divides the 3-sphere in two.
		/// </summary>
		public static Torus CreateTorus( Parameters parameters )
		{
			Torus t = new Torus();
			t.Params = parameters;
			t.InitVerts();

			// Shorter var names for inputs.
			int n1 = parameters.NumSegments1;
			int n2 = parameters.NumSegments2;
			double r = parameters.Radius;
			double r1 = parameters.TubeRadius1;

			// Calc r2.
			double r2 = r - r1;
			if( r2 < 0 )
				r2 = 0;

			r1 *= Math.Sqrt( 2 );
			r2 *= Math.Sqrt( 2 );
				
			double angleInc1 = 2 * Math.PI / n1;
			double angleInc2 = 2 * Math.PI / n2;

			double angle1 = 0;
			for( int i = 0; i < n1; i++ )
			{
				double angle2 = 0;
				for( int j = 0; j < n2; j++ )
				{
					t.Vertices[i][j].X = r1 * Math.Cos( angle1 );
					t.Vertices[i][j].Y = r1 * Math.Sin( angle1 );
					t.Vertices[i][j].Z = r2 * Math.Cos( angle2 );
					t.Vertices[i][j].W = r2 * Math.Sin( angle2 );
					angle2 += angleInc2;
				}
				angle1 += angleInc1;
			}

			return t;
		}
	}
}
