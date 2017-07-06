namespace R3.Math
{
	public struct GraphEdge
	{
		public GraphEdge( int v1, int v2 )
			: this()
		{
			// Keep it ordered
			if( v1 < v2 )
			{
				V1 = v1;
				V2 = v2;
			}
			else
			{
				V1 = v2;
				V2 = v1;
			}
		}

		public int V1 { get; set; }
		public int V2 { get; set; }

		/// <summary>
		/// Given a vertex index, find the vertex at the other end of the edge.
		/// </summary>
		public int Opposite( int idx )
		{
			return idx == V1 ? V2 : V1;
		}

		/// <summary>
		/// vZome VEF format.
		/// </summary>
		public void ReadEdge( string line )
		{
			//string[] split = line.Split( '\t' );
			string[] split = line.Split( new char[] { '\t', ' ' }, System.StringSplitOptions.RemoveEmptyEntries );
			V1 = int.Parse( split[0] );
			V2 = int.Parse( split[1] );
		}

		/// <summary>
		/// vZome VEF format.
		/// </summary>
		public string WriteEdge()
		{
			return V1.ToString() + "\t" + V2.ToString();
		}
	}
}
