namespace R3.Algorithm
{
	using Math = System.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;
	using R3.Math;

	public class GraphNode
	{
		public GraphNode( VectorND pos, VectorND vel )
		{
			Position = pos;
			Velocity = vel;
			Lock = new VectorND( pos.Dimension );	// No locking.
		}

		public VectorND Position { get; set; }
		public VectorND Velocity { get; set; }

		/// <summary>
		/// Used to control how the node will move during relaxation.
		/// To lock a coordinate, set that coordinate to a positive number.
		/// </summary>
		public VectorND Lock { get; set; }

		public void FullLock()
		{
			for( int i = 0; i < Lock.Dimension; i++ )
				Lock.X[i] = 1;
		}
	}

	public class Graph
	{
		public Graph()
		{
			Nodes = new List<GraphNode>();
			Edges = new List<GraphEdge>();
			m_connections = new Dictionary<int,List<int>>();
		}

		public List<GraphNode> Nodes { get; set; }
		public List<GraphEdge> Edges { get; set; }
		
		public static Graph FromMesh( Mesh mesh )
		{
			throw new System.NotImplementedException();
		}

		public Mesh ToMesh()
		{
			throw new System.NotImplementedException();
		}

		/// <summary>
		/// Setup a complete graph with n nodes.
		/// Positions are randomly distributed.
		/// </summary>
		public void SetupCompleteGraph( int n )
		{
			Nodes.Clear();
			Edges.Clear();

			GraphNode southPole = new GraphNode( new VectorND( new double[] { 0, 0, 0, -1 } ), new VectorND( 4 ) );
			southPole.FullLock();
			Nodes.Add( southPole );

			System.Random rand = new System.Random( 0 );
			for( int i = 1; i < n; i++ )
			{ 
				VectorND randPosition = new VectorND( new double[] {
					rand.NextDouble() - 0.5,
					rand.NextDouble() - 0.5,
					rand.NextDouble() - 0.5,
					rand.NextDouble() - 0.5 } );
				randPosition.Normalize();
				Nodes.Add( new GraphNode( randPosition, new VectorND( 4 ) ) );
			}

			for( int i = 0; i < n; i++ )
			for( int j = i + 1; j < n; j++ )
			{ 
				AddEdge( new GraphEdge( i, j ) );
			}
		}

		public void AddEdge( GraphEdge e )
		{
			List<int> vals;

			if( m_connections.TryGetValue( e.V1, out vals ) )
			{
				if( !vals.Contains( e.V2 ) )
					vals.Add( e.V2 );
			}
			else
				m_connections[e.V1] = new List<int>( new int[] { e.V2 } );

			/*if( m_connections.TryGetValue( e.V2, out vals ) )
			{
				if( !vals.Contains( e.V1 ) )
					vals.Add( e.V1 );
			}
			else
				m_connections[e.V2] = new List<int>( new int[] { e.V1 } );*/
			
			Edges.Add( e );
		}

		public Dictionary<int,List<int>> m_connections;

		/// <summary>
		/// ZZZ - slow impl.  Should calc this once and cache.
		/// </summary>>
		public bool Connected( int n1, int n2 )
		{
			//return Edges.Any( e => ( e.V1 == n1 && e.V2 == n2 ) || ( e.V2 == n1 && e.V1 == n2 ) );
			List<int> vals;
			if( m_connections.TryGetValue( n1, out vals ) && vals.Contains( n2 ) )
				return true;
			if( m_connections.TryGetValue( n2, out vals ) && vals.Contains( n1 ) )
				return true;
			return false;
		}

		public void Normalize()
		{
			var ordered = Nodes.OrderByDescending( n => n.Position.Abs );
			GraphNode largest = ordered.First();

			double mag = largest.Position.Abs;
			foreach( GraphNode node in Nodes )
			{
				node.Position /= mag;
				node.Position *= (2 + 1 * Golden.tau);
			}
		}
	}

	public class GraphRelaxation
	{
		/// <summary>
		/// The graph we're operating on.
		/// </summary>
		public Graph Graph { get; set; }

		/// <summary>
		/// An attractive force between nodes connected by an edge.  
		/// This is like a rubber band along the axis of the edge.
		/// </summary>
		public double EdgeAttraction { get; set; }

		/// <summary>
		/// A repulsion force between all nodes.
		/// Nodes act like similarly charged ions.
		/// </summary>
		public double NodeRepulsion { get; set; }

		/// <summary>
		/// A repulsion force between all edges.
		/// </summary>
		public double EdgeRepulsion { get; set; }

		public void Relax( int steps )
		{
			// ZZZ - add convergence criterion, rather than hard coded number of steps.
			System.Console.WriteLine( "Relaxing" );
			for( int i=0; i<steps; i++ )
			{
				if( i%20 == 0 )
					System.Console.WriteLine( i );

				VectorND[] accelerations = CalcAccelerations();
				for( int j = 0; j < Graph.Nodes.Count; j++ )
					UpdatePositionAndVelocity( Graph.Nodes[j], accelerations[j] );
			}
		}

		private int m_dim = 4;

		private VectorND[] CalcAccelerations()
		{
			int count = Graph.Nodes.Count;
			VectorND[] accelerations = new VectorND[count];
			for( int i = 0; i < count; i++ )
				accelerations[i] = new VectorND( m_dim );

			bool nodeRepulse = !Tolerance.Zero( NodeRepulsion );
			bool edgeRepulse = !Tolerance.Zero( EdgeRepulsion );

			/*
			for( int i = 0; i < count; i++ )
			for( int j = i + 1; j < count; j++ )
			{
				if( nodeRepulse )
				{
					VectorND nodeForce = CalculateForce( Graph.Nodes[i], Graph.Nodes[j], NodeRepulsion, square: true );
					accelerations[i] -= nodeForce;	// Repulsive.
					accelerations[j] += nodeForce;
				}

				if( Graph.Connected( i, j ) )
				{
					VectorND edgeForce = CalculateForce( Graph.Nodes[i], Graph.Nodes[j], EdgeAttraction, square: false );
					accelerations[i] += edgeForce;	// Attractive.
					accelerations[j] -= edgeForce;
				}
			}*/

			foreach( var kvp in Graph.m_connections )
			{
				int i = kvp.Key;
				foreach( int j in kvp.Value )
				{
					VectorND edgeForce = CalculateForce( Graph.Nodes[i], Graph.Nodes[j], EdgeAttraction, square: false );
					accelerations[i] += edgeForce;  // Attractive.
					accelerations[j] -= edgeForce;
				}
			}

			if( edgeRepulse )
			{ 
				count = Graph.Edges.Count;
				for( int i = 0; i < count; i++ )
				for( int j = i + 1; j < count; j++ )
				{
					// Rather than mess with torques and doing this "right" (like it was two charged rod segments),
					// We'll calculate the effect on the two COMs, and give half the force to each node.

					int n1 = Graph.Edges[i].V1;
					int n2 = Graph.Edges[i].V2;
					int n3 = Graph.Edges[j].V1;
					int n4 = Graph.Edges[j].V2;
					GraphNode center1 = new GraphNode( ( Graph.Nodes[n1].Position + Graph.Nodes[n2].Position ) / 2, new VectorND( m_dim ) );
					GraphNode center2 = new GraphNode( ( Graph.Nodes[n3].Position + Graph.Nodes[n4].Position ) / 2, new VectorND( m_dim ) );

					VectorND force = CalculateForce( center1, center2, EdgeRepulsion, square: true ) / 2;

					accelerations[n1] -= force;
					accelerations[n2] -= force;
					accelerations[n3] += force;
					accelerations[n3] += force;
				}
			}

			return accelerations;
		}

		private VectorND CalculateForce( GraphNode node1, GraphNode node2, double strength, bool square )
		{
			double distance = node1.Position.Dist( node2.Position );

			// Here is the direction vector of the force.
			VectorND force = node2.Position - node1.Position;
			if( !force.Normalize() )
			{
				Debug.Assert( false );
				return new VectorND( node1.Position.Dimension );
			}

			// Calculate the magnitude.
			double mag = 0;
			if( square )
				mag = strength / Math.Pow( distance, 2 );
			else
			{
				mag = strength * distance;	// http://en.wikipedia.org/wiki/Hooke's_law
				
				// Try to make all edges a specific length.
				
				/*double length = 0.1;
				double diff = distance - length;
				mag = strength * diff;*/

				//if( diff < 0 )
					//mag = -strength / Math.Pow( diff, 2 );
				//else
					//mag = strength / Math.Pow( diff, 2 );
			}

			if( mag > 100 )
				mag = 100;

			return force * mag;
		}

		private void UpdatePositionAndVelocity( GraphNode node, VectorND acceleration )
		{
			//if( allLocked )
			//	return;

			VectorND position = node.Position;
			VectorND velocity = node.Velocity;
			//if( position.IsOrigin )
			//	return;

			// Leapfrog method.
			double timeStep = 1;
			velocity += acceleration * timeStep;
			velocity *= .5; // Damping.
			ZeroComponents( node.Lock, ref velocity );
			position += velocity * timeStep;
			//position.Normalize(); position *= 5;

			//if( position.MagSquared > 1 )
			//{
			//	position.Normalize();
			//	velocity = new VectorND( 3 );
			//}

			position.Normalize();

			node.Position = position;
			node.Velocity = velocity;
			//node.Acceleration = acceleration;  Any reason to store this?
		}

		private void ZeroComponents( VectorND components, ref VectorND input )
		{
			for( int i = 0; i < 4; i++ )
				if( components.X[i] > 0 )
					input.X[i] = 0;
		}
	}
}
