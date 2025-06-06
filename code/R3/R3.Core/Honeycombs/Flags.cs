namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	public class Flag
	{
		public Flag( int cell, int facet, int edge, int vert )
		{
			Cell = cell;
			Facet = facet;
			Edge = edge;
			Vert = vert;
		}

		public int Cell;
		public int Facet;
		public int Edge;
		public int Vert;

		public override string ToString()
		{
			return string.Format( "{0}\t{1}\t{2}\t{3}", Cell, Facet, Edge, Vert );
		}
	}

	public static class Flags
	{
		private static int BoundaryRadiusEdges( Tiling tiling )
		{
			int p = tiling.TilingConfig.P;
			int result = tiling.Tiles.Sum( t => p - t.EdgeIncidences.Count );
			return result;
		}

		public static Vector3D[] BoundaryVerts( Tiling tiling )
		{
			// Count all vertices without q surrounding tiles.
			Dictionary<Vector3D, int> vCounts = new Dictionary<Vector3D, int>();
			foreach( Tile t in tiling.Tiles )
			{
				foreach( Vector3D v in t.Boundary.Vertices )
				{
					int count;
					if( !vCounts.TryGetValue( v, out count ) )
						count = 0;
					vCounts[v] = ++count;
				}
			}

			int q = tiling.TilingConfig.Q;
			//int result = vCounts.Where( kvp => kvp.Value < q ).Count();

			return vCounts.Where( kvp => kvp.Value < q ).Select( kvp => kvp.Key ).ToArray();
		}

		public static void GenForTiling( int p, int q, int levels )
		{
			s_zeroIndexed = true;
			s_baseFn = string.Format( "{0}-{1}_layer_{2}_{3}-indexed", p, q, levels, "0" );

			TilingConfig config = new TilingConfig( p, q, maxTiles: 1 );
			config.Levels = levels;
			Tiling tiling = new Tiling();
			tiling.Generate( config );

			int bCount1 = BoundaryRadiusEdges( tiling );
			BoundaryVerts( tiling );

			OutputTiling( tiling, levels );
		}

		private static bool s_zeroIndexed;
		private static string s_baseFn;

		public static void GenForTilingDistanceLayers( int p, int q )
		{
			TilingConfig config1 = new TilingConfig( p, q, maxTiles: (int)0.1e6 );
			Tiling tiling1 = new Tiling();
			tiling1.Generate( config1 );
			var map1 = SortIntoDistanceLayers( tiling1 );

			/*
			TilingConfig config2 = new TilingConfig( p, q, maxTiles: (int)1.6e4 );
			Tiling tiling2 = new Tiling();
			tiling2.Generate( config2 );
			var map2 = SortIntoDistanceLayers( tiling2 );
			*/

			int layers = 10;
			//for( int i = 1; i <= layers; i++ )
			int i = 6;
			{
				double linkLength = 2 * Geometry2D.GetTrianglePSide( p, q );
				double cutoff = linkLength * i;
				cutoff = linkLength * 5.909;//4.29;
				double cutoffE = DonHatch.h2eNorm( cutoff );

				IEnumerable<Tile> tiles = tiling1.Tiles.Where( t => t.Center.Abs() <= cutoffE );

				tiling1.SetTiles( tiles );
				tiling1.FillOutIncidences();

				BoundaryRadiusEdges( tiling1 );

				// Check that we recursed deep enough.
				foreach( Tile tile in tiles )
				{
					int incidenceCount = ( q - 2 ) * p;
					if( tile.VertexIndicences.Count != incidenceCount )
						throw new System.Exception( "didn't recurse deep enough!" );
				}

				//OutputTiling( tiles, p, q, i );
			}
		}

		private static Dictionary<double, List<Tile>> SortIntoDistanceLayers( Tiling tiling )
		{
			Dictionary<double, List<Tile>> result = new Dictionary<double, List<Tile>>( new DoubleEqualityComparer( Tolerance.ThresholdStrict ) );

			foreach( Tile tile in tiling.Tiles )
			{
				List<Tile> list;
				double dist = tile.Center.Abs();
				if( !result.TryGetValue( dist, out list ) )
				{
					list = new List<Tile>();
					result[dist] = list;
				}

				list.Add( tile );
			}

			return result;
		}

		private static void OutputTiling( Tiling tiling, int levels )
		{
			IEnumerable<Tile> tiles = tiling.Tiles;
			int p = tiling.TilingConfig.P;
			int q = tiling.TilingConfig.Q;

			// Convert to Honeycomb structs.
			List<H3.Cell.Facet> facets = new List<H3.Cell.Facet>();
			foreach( Tile tile in tiles )
			{
				H3.Cell.Facet facet = new H3.Cell.Facet( tile.Boundary.Vertices );
				facets.Add( facet );
			}
			H3.Cell dummy = new H3.Cell( facets.ToArray() );

			AppendCounts( "counts.txt", "Level " + levels );
			string filename = string.Format( "{0}_{1}_level_{2}", p, q, levels );
			Gen( new H3.Cell[] { dummy }, filename, tiling );
		}

		private static void AppendCounts( string filename, string line )
		{
			using( StreamWriter sw = File.AppendText( filename ) )
				sw.WriteLine( line );
		}

		public static void Gen( H3.Cell[] completedCells, string filename = "flags.txt", Tiling tiling = null )
		{
			Dictionary<H3.Cell, int> allCells = new Dictionary<H3.Cell, int>();
			Dictionary<H3.Cell.Facet, int> allFacets = new Dictionary<H3.Cell.Facet, int>( new H3.Cell.ElementEqualityComparer() );
			Dictionary<H3.Cell.Edge, int> allEdges = new Dictionary<H3.Cell.Edge, int>( new H3.Cell.ElementEqualityComparer() );
			Dictionary<Vector3D, int> allVerts = new Dictionary<Vector3D, int>();
			HashSet<Flag> flags = new HashSet<Flag>();
			int cNum = 1, fNum = 1, eNum = 1, vNum = 1;
			foreach( H3.Cell c in completedCells )
			{
				if( c.ID.DNE )
					throw new System.Exception( "Cell tolerance exceeded." );

				c.ElementId = cNum;
				allCells[c] = cNum++;
				foreach( H3.Cell.Facet f in c.Facets )
				{
					if( f.ID.DNE )
						throw new System.Exception( "Facet tolerance exceeded." );

					int eId;
					if( !allFacets.TryGetValue( f, out eId ) )
					{
						f.ElementId = fNum;
						allFacets[f] = fNum++;
					}
					else
					{
						f.ElementId = eId;
					}

					for( int i = 0; i < f.Verts.Length; i++ )
					{
						Vector3D s = f.Verts[i];
						Vector3D e = i == f.Verts.Length - 1 ? f.Verts[0] : f.Verts[i + 1];
						H3.Cell.Edge e_ = new H3.Cell.Edge( s, e );
						if( !allEdges.TryGetValue( e_, out eId ) )
						{
							e_.ElementId = eNum;
							allEdges[e_] = eNum++;
						}
						else
						{
							e_.ElementId = eId;
						}

						if( e_.ID.DNE )
							throw new System.Exception( "Edge tolerance exceeded." );

						Vector3D[] eVerts = new Vector3D[] { s, e };
						foreach( Vector3D eVert in eVerts )
						{
							if( eVert.DNE )
								throw new System.Exception( "Vertex tolerance exceeded." );

							int vElementId;
							if( !allVerts.TryGetValue( eVert, out vElementId ) )
							{
								vElementId = vNum;
								allVerts[eVert] = vNum++;
							}

							Flag flag = new Flag( c.ElementId, f.ElementId, e_.ElementId, vElementId );
							flags.Add( flag );
						}
					}
				}
			}

			Dictionary<int, HashSet<int>> vertAdj = new Dictionary<int, HashSet<int>>();
			foreach( H3.Cell.Edge e in allEdges.Keys )
			{
				int sId = allVerts[e.Start];
				int eId = allVerts[e.End];

				{
					HashSet<int> adj;
					if( !vertAdj.TryGetValue( sId, out adj ) )
						adj = new HashSet<int>();
					adj.Add( eId );
					vertAdj[sId] = adj;
				}

				{
					HashSet<int> adj;
					if( !vertAdj.TryGetValue( eId, out adj ) )
						adj = new HashSet<int>();
					adj.Add( sId );
					vertAdj[eId] = adj;
				}
			}

			var test = vertAdj.Where( kvp => kvp.Value.Count() < 8 ).Select( kvp => kvp.Key ).ToArray();	// Does not give all boundary verts
			Vector3D[] bVertsArray = BoundaryVerts( tiling );
			foreach( Vector3D v in bVertsArray )
			{
				int vId = allVerts[v];
				if( !test.Contains( vId ) )
				{
					//throw new System.Exception();
				}
			}

			System.Func<int, int> applyZeroIndexing = i => s_zeroIndexed ? i - 1 : i;


			using( StreamWriter sw = new StreamWriter( s_baseFn + "_vert_incidence.txt" ) )
			{
				foreach( var kvp in vertAdj.OrderBy( kvp => kvp.Key ) )
				{
					int[] vAdj = kvp.Value.OrderBy( i => i ).Select( i => applyZeroIndexing( i ) ).ToArray();
					string line = string.Join( "\t", vAdj );
					sw.WriteLine( line );
				}
			}

			using( StreamWriter sw = new StreamWriter( s_baseFn + "_vert_locations.txt" ) )
			{
				foreach( var kvp in allVerts.OrderBy( kvp => kvp.Value ) )
				{
					string line = string.Format( "{0}\t{1}\t{2}", applyZeroIndexing( kvp.Value ), kvp.Key.X, kvp.Key.Y );
					sw.WriteLine( line );
				}
			}

			using( StreamWriter sw = new StreamWriter( s_baseFn + "_boundary_vert_IDs.txt" ) )
			{
				foreach( Vector3D v in bVertsArray )
				{
					int id = allVerts[v];
					string line = string.Format( "{0}", applyZeroIndexing( id ) );
					sw.WriteLine( line );
				}
			}

			using( StreamWriter sw = new StreamWriter( s_baseFn + "_boundary_verts.txt" ) )
			{
				foreach( Vector3D v in bVertsArray )
				{
					string line = string.Format( "{0}\t{1}", v.X, v.Y );
					sw.WriteLine( line );
				}
			}

			return;

			//CalcH3BoundaryValues( allCells, allFacets, allEdges, allVerts );

			// Interactive testing
			//var t1 = flags.Where( f => f.Vert == 20 ).GroupBy( f => f.Cell ).ToArray();

			filename = string.Format( "L{0}", m_level );
			using( StreamWriter sw = new StreamWriter( filename + "_flags.txt" ) )
			{
				foreach( Flag flag in flags )
					sw.WriteLine( flag.ToString() );
			}

			return;

			string fn = "counts.txt";
			AppendCounts( fn, "\tNum tiles: " + allFacets.Count );
			AppendCounts( fn, "\tNum edges: " + allEdges.Count );
			AppendCounts( fn, "\tNum verts: " + allVerts.Count );

			// Geometry tiles.
			using( StreamWriter sw = new StreamWriter( filename + "_tile_geometry.txt" ) )
			{
				foreach( H3.Cell.Facet f in allFacets.Keys )
				{
					string line = "";
					foreach( Vector3D v in f.Verts )
						line += string.Format( "{0}\t{1}\t", v.X, v.Y );
					sw.WriteLine( line );
				}
			}

			// Geometry edges.
			using( StreamWriter sw = new StreamWriter( filename + "_edge_geometry.txt" ) )
			{
				foreach( H3.Cell.Edge e in allEdges.Keys )
				{
					string line = string.Format( "{0}\t{1}\t{2}\t{3}\t", e.Start.X, e.Start.Y, e.End.X, e.End.Y );
					sw.WriteLine( line );
				}
			}
		}

		private static void CalcH3BoundaryValues(
			Dictionary<H3.Cell, int> allCells,
			Dictionary<H3.Cell.Facet, int> allFacets,
			Dictionary<H3.Cell.Edge, int> allEdges,
			Dictionary<Vector3D, int> allVerts )
		{
			Dictionary<H3.Cell.Facet, int> facetCellCount = new Dictionary<H3.Cell.Facet, int>( new H3.Cell.ElementEqualityComparer() );
			foreach( H3.Cell cell in allCells.Keys )
			foreach( H3.Cell.Facet facet in cell.Facets )
			{
				int count;
				if( !facetCellCount.TryGetValue( facet, out count ) )
					count = 0;
				facetCellCount[facet] = ++count;
			}

			// All boundary facets have a cell count of 1.
			IEnumerable<H3.Cell.Facet> boundaryFacets = facetCellCount.Where( kvp => kvp.Value == 1 ).Select( kvp => kvp.Key );

			// Boundary facets have boundary elements (note: this is not true for cells).
			HashSet<H3.Cell.Edge> boundaryEdges = new HashSet<H3.Cell.Edge>( new H3.Cell.ElementEqualityComparer() );
			HashSet<Vector3D> boundaryVerts = new HashSet<Vector3D>();
			foreach( H3.Cell.Facet f in boundaryFacets )
			{
				for( int i = 0; i < f.Verts.Length; i++ )
				{
					Vector3D s = f.Verts[i];
					Vector3D e = i == f.Verts.Length - 1 ? f.Verts[0] : f.Verts[i + 1];
					H3.Cell.Edge e_ = new H3.Cell.Edge( s, e );
					boundaryEdges.Add( e_ );
					boundaryVerts.Add( s );
					boundaryVerts.Add( e );
				}
			}

			string fn = "boundaryCounts.txt";
			AppendCounts( fn, "Level " + m_level );
			AppendCounts( fn, "\tNum boundary facets: " + boundaryFacets.Count() );
			AppendCounts( fn, "\tNum boundary edges: " + boundaryEdges.Count );
			AppendCounts( fn, "\tNum boundary verts: " + boundaryVerts.Count );
		}

		public static int m_level;

		private static void Sandbox()
		{
			Dictionary<H3.Cell.Edge, int> completedEdges = new Dictionary<H3.Cell.Edge, int>();

			Dictionary<Vector3D, int> vertices = new Dictionary<Vector3D, int>();
			int count = 1;
			foreach( var e in completedEdges )
			{
				AddVert( vertices, e.Key.Start, ref count );
				AddVert( vertices, e.Key.End, ref count );
			}

			Dictionary<int, HashSet<int>> vertexIndicences = new Dictionary<int, HashSet<int>>();
			foreach( var e in completedEdges )
			{
				AddVertIncidence( vertices, vertexIndicences, e.Key.Start, e.Key.End );
				AddVertIncidence( vertices, vertexIndicences, e.Key.End, e.Key.Start );
			}

			using( StreamWriter sw = new StreamWriter( "vert_incidences.txt" ) )
			{
				foreach( var v in vertexIndicences.OrderBy( vi => vi.Key ) )
				{
					int[] incident = v.Value.OrderBy( i => i ).ToArray();
					sw.WriteLine( string.Format( "{0}:{1}", v.Key, string.Join( ",", incident ) ) );
				}
			}
		}

		private static void AddVertIncidence( Dictionary<Vector3D, int> vertices, Dictionary<int, HashSet<int>> vertexIndicences,
			Vector3D v1, Vector3D v2 )
		{
			int v1_ = vertices[v1];
			int v2_ = vertices[v2];
			HashSet<int> incident;
			if( !vertexIndicences.TryGetValue( v1_, out incident ) )
				vertexIndicences[v1_] = incident = new HashSet<int>();
			incident.Add( v2_ );
		}

		private static void AddVert( Dictionary<Vector3D, int> vertices, Vector3D v, ref int count )
		{
			int dummy;
			if( vertices.TryGetValue( v, out dummy ) )
				return;

			vertices[v] = count;
			count++;
		}
	}
}
