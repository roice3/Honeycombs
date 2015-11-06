namespace HyperbolicModels
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;

	internal class Nonuniform
	{
		/// <summary>
		/// Wendy's 77
		/// </summary>
		public static void Wendy( Simplex simplex, H3.Cell.Edge[] edges )
		{
			H3.Cell startingCell = null;

			Vector3D start = startingCell.Verts.First();

			Func<Vector3D, Vector3D> findAntipode = input =>
			{
				Vector3D antipode = new Vector3D();
				double max = double.MinValue;
				foreach( Vector3D v in startingCell.Verts )
				{
					double d = H3Models.Ball.HDist( v, input );
					if( d > max )
					{
						max = d;
						antipode = v;
					}
				}
				return antipode;
			};

			H3.Cell.Edge[] diagonals = new H3.Cell.Edge[] { new H3.Cell.Edge( start, findAntipode( start ) ) };
			diagonals = Recurse.CalcEdges( simplex.Facets, diagonals, new Recurse.Settings() { Threshold = 0.9983 } );

			// diagonals includes too much at this point (it includes all icosahedra diagonals, but we only want one diagonal from each cell).
			// We need to begin at 4 start points, and branch out from each to find the ones we want.

			var vertsToDiagonals = FindConnectedEdges( diagonals );
			var connectedEdges = FindConnectedEdges( edges );

			// Get all edges (not diagonals) connected to start.
			List<H3.Cell.Edge> connectedToStart = connectedEdges[start];
			Vector3D startOpp = connectedToStart[0].Opp( start );
			List<H3.Cell.Edge> connectedToStartOpp = connectedEdges[startOpp];

			// We need to pick 4 of these edges, arranged in a tetrahedron for our starting points.
			List<Vector3D> startingPoints = new List<Vector3D>();

			List<double> distances = new List<double>();

			// View1
			//startingPoints.Add( start );
			//foreach( Vector3D v in connectedToStartOpp.Select( e => e.Opp( startOpp ) ) )
			//	distances.Add( H3Models.Ball.HDist( startingPoints.First(), v ) );
			//startingPoints.Add( connectedToStartOpp[10].Opp( startOpp ) );
			//startingPoints.Add( connectedToStartOpp[13].Opp( startOpp ) );
			//startingPoints.Add( connectedToStartOpp[14].Opp( startOpp ) );

			// View2
			startingPoints.Add( startOpp );
			foreach( Vector3D v in connectedToStart.Select( e => e.Opp( start ) ) )
				distances.Add( H3Models.Ball.HDist( startingPoints.First(), v ) );
			startingPoints.Add( connectedToStart[10].Opp( start ) );
			startingPoints.Add( connectedToStart[13].Opp( start ) );
			startingPoints.Add( connectedToStart[14].Opp( start ) );

			distances.Clear();
			distances.Add( H3Models.Ball.HDist( startingPoints[1], startingPoints[2] ) );
			distances.Add( H3Models.Ball.HDist( startingPoints[1], startingPoints[3] ) );
			distances.Add( H3Models.Ball.HDist( startingPoints[2], startingPoints[3] ) );
			distances.Add( H3Models.Ball.HDist( startingPoints[0], startingPoints[1] ) );
			distances.Add( H3Models.Ball.HDist( startingPoints[0], startingPoints[2] ) );
			distances.Add( H3Models.Ball.HDist( startingPoints[0], startingPoints[3] ) );
			double dist = 3.097167;

			Func<Vector3D[], H3.Cell.Edge[]> RemoveVerts = starting =>
			{
				List<H3.Cell.Edge> keepers = new List<H3.Cell.Edge>();
				HashSet<Vector3D> removedVerts = new HashSet<Vector3D>();
				Recurse.BranchAlongVerts( starting.ToArray(), vertsToDiagonals, removedVerts );
				foreach( H3.Cell.Edge e in edges )
				{
					if( removedVerts.Contains( e.Start ) || removedVerts.Contains( e.End ) )
						continue;
					keepers.Add( e );
				}
				return keepers.ToArray();
			};

			edges = RemoveVerts( startingPoints.ToArray() );

			bool done = false;
			while( !done )
			{
				done = true;
				var newConnectedEdges = FindConnectedEdges( edges );
				foreach( Vector3D v in newConnectedEdges.Keys )
				{
					List<H3.Cell.Edge> oldEdgeList = connectedEdges[v];
					List<H3.Cell.Edge> newEdgeList = newConnectedEdges[v];

					// Only work edges that were full originally.
					if( oldEdgeList.Count != 20 )
						continue;

					// We need at least two to find the rest.
					int newCount = newEdgeList.Count;
					if( newCount > 16 && newCount < 19 )
					{
						List<H3.Cell.Edge> removed = oldEdgeList.Except( newEdgeList, new H3.Cell.EdgeEqualityComparer() ).ToList();

						H3.Cell.Edge[] toTrim = newEdgeList.FindAll( e =>
						{
							foreach( H3.Cell.Edge alreadyRemoved in removed )
							{
								double d = H3Models.Ball.HDist( alreadyRemoved.Opp( v ), e.Opp( v ) );
								if( !Tolerance.Equal( dist, d, 0.00001 ) )
									return false;
							}

							return true;
						} ).ToArray();

						edges = RemoveVerts( toTrim.Select( e => e.Opp( v ) ).ToArray() );
						done = false;
					}

					if( newCount == 20 )
						done = false;
				}
			}
		}

		private static Dictionary<Vector3D, List<H3.Cell.Edge>> FindConnectedEdges( H3.Cell.Edge[] edges )
		{
			var result = new Dictionary<Vector3D, List<H3.Cell.Edge>>();
			System.Action<Vector3D, H3.Cell.Edge> addOne = ( v, e ) =>
			{
				List<H3.Cell.Edge> edgeList;
				if( !result.TryGetValue( v, out edgeList ) )
					edgeList = new List<H3.Cell.Edge>();

				edgeList.Add( e );
				result[v] = edgeList;
			};

			foreach( H3.Cell.Edge edge in edges )
			{
				addOne( edge.Start, edge );
				addOne( edge.End, edge );
			}

			return result;
		}
	}
}
