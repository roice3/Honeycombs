namespace HyperbolicModels
{
	// TODO
	// EdgeEqualityComparer requirement everywhere is buggy.  Perhaps just make the edge class comparable itself, so this always happens.

	using R3.Core;
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.IO;
	using System.Linq;
	using System.Numerics;

	using Math = System.Math;
	using Edge = R3.Geometry.H3.Cell.Edge;
	using VertexMap = System.Collections.Generic.Dictionary<R3.Geometry.Vector3D, System.Collections.Generic.List<R3.Geometry.Polygon>>;
	using EdgeMap = System.Collections.Generic.Dictionary<R3.Geometry.H3.Cell.Edge, System.Collections.Generic.List<R3.Geometry.Polygon>>;
	using VertexToEdgeMap = System.Collections.Generic.Dictionary<R3.Geometry.Vector3D, System.Collections.Generic.HashSet<R3.Geometry.H3.Cell.Edge>>;

	public class Acrohedron
	{
		// Variables to track various things about our polyhedron.
		private List<Polygon> m_polygons = new List<Polygon>();
		private readonly VertexMap m_vertexMap = new VertexMap();
		private readonly EdgeMap m_edgeMap = new EdgeMap( new H3.Cell.EdgeEqualityComparer() );
		private readonly VertexToEdgeMap m_vertexToEdgeMap = new VertexToEdgeMap();

		public void Generate( int X, int Y, int Z )
		{
			// The starting 3 polygons, connected up at a vertex.
			Polygon[] starting = GenerateStarting( X, Y, Z );
			starting.ToList().ForEach( p => AddPoly( p ) );

			for( int i = 0; i < 20; i++ )
				OneChain();

			//m_polygons = m_polygons.Skip( 55 ).Take(3).ToList();

			// Write it out so we can see it.
			File.Delete( "test.pov" );
			PovRay.AppendEuclideanPolygons( m_polygons.ToArray(), "test.pov" );
		}

		private static Polygon[] GenerateStarting( int X, int Y, int Z )
		{
			Geometry g = Geometry.Euclidean;

			// First generate the 3 starting polygons.
			Polygon poly1 = Polygon.CreateEuclidean( X );
			Polygon poly2 = Polygon.CreateEuclidean( Y );
			Polygon poly3 = Polygon.CreateEuclidean( Z );
			//double a1 = Polygon.InteriorAngle( X );
			//double a2 = Polygon.InteriorAngle( Y );
			//double a3 = Polygon.InteriorAngle( Z );

			// Normalize polygons to have edge length = 1.
			poly1.Scale( 1.0 / poly1.Segments.First().Length );
			poly2.Scale( 1.0 / poly2.Segments.First().Length );
			poly3.Scale( 1.0 / poly3.Segments.First().Length );

			// Transform them so their first vertex is at the origin,
			// And poly1/poly2 share an edge, as well as poly1/poly3
			// poly2 and poly3 may or may not share an edge.
			Complex origin = new Complex( 0, 0 );
			poly1.Transform( Mobius.CreateFromIsometry( g, 0, -poly1.Vertices.First() ) );
			poly2.Transform( Mobius.CreateFromIsometry( g, 0, -poly2.Vertices.First() ) );
			poly3.Transform( Mobius.CreateFromIsometry( g, 0, -poly3.Vertices.First() ) );
			poly2.Transform( Mobius.CreateFromIsometry( g, poly1.Vertices.Last().AngleTo( poly2.Vertices.Skip( 1 ).First() ), origin ) );
			poly3.Transform( Mobius.CreateFromIsometry( g, -poly1.Vertices.Last().AngleTo( poly3.Vertices.Skip( 1 ).First() ), origin ) );

			// Now connect poly2/poly3.
			ConnectPolys( poly1, poly2, poly3 );

			return new Polygon[] { poly1, poly2, poly3 };
		}

		private void OneChain()
		{
			// Here are all unattached edges.
			var edgeKvps = m_edgeMap.Where( kvp => kvp.Value.Count == 1 ).ToArray();
			Edge[] edges = edgeKvps.Select( kvp => kvp.Key ).ToArray();

			Trace.WriteLine( string.Format( "Chain length {0}", edges.Length ) );

			// Find all the vertices connected to two polygons.
			// This will necessarily lie on the edges above.
			// These are vertices that are candidates for filling in with a polygon.
			Vector3D[] vertices = m_vertexMap.Where( kvp => kvp.Value.Count > 1 ).Select( kvp => kvp.Key ).ToArray();
			HashSet<Vector3D> edgeVerts = new HashSet<Vector3D>();
			foreach( Edge e in edges )
			{
				edgeVerts.Add( e.Start );
				edgeVerts.Add( e.End );
			}
			//Vector3D[] vertices = edgeVerts.ToArray();

			// XXX - Things that are hard.
			//	- How can we tell when a vertex is done?
			//  - How to avoid duplicating edges and vertices?

			// Attempt to fill those verts.
			foreach( Vector3D v in vertices )
			{
				Edge[] connectedEdges = m_vertexToEdgeMap[v].Intersect( edges, new H3.Cell.EdgeEqualityComparer() ).ToArray();

				// Is it on the chain?
				if( connectedEdges.Length != 2 )
					continue;

				Edge e1 = connectedEdges[0];
				Edge e2 = connectedEdges[1];
				Vector3D v1 = e1.Opp( v ) - v;
				Vector3D v2 = e2.Opp( v ) - v;
				double angle = v1.AngleTo( v2 );
				int nSides= PolyWithAngle( angle );
				if( nSides != 0 )
				{
					Vector3D normal = v1.Cross( v2 );
					Polygon poly = Polygon.CreateEuclidean( nSides, v, e1.Opp( v ), normal );

					// Things may have changed, so we must double check.
					if( m_edgeMap[e1].Count == 1 &&
						m_edgeMap[e2].Count == 1 )
					{
						AddPoly( poly );
						Trace.WriteLine( string.Format( "Added {0}-gon", nSides ) );
					}
				}
			}
		}

		private int m_max = 10;

		private int PolyWithAngle( double angle )
		{
			for( int i=3; i<m_max; i++ )
			{
				if( Tolerance.Equal( angle, Polygon.InteriorAngle( i ) ) )
					return i;
			}

			return 0;
		}

		/// <summary>
		/// Connects two polygons by rotating them about edges of a static polygon
		/// The three polygons will meet at a single vertex afterward (this vertex will be convex).
		/// </summary>
		private static void ConnectPolys( Polygon staticPoly, Polygon poly1, Polygon poly2 )
		{
			// All three polygons should share one vertex.
			// The edges connected to this vertex (that are in the static polyhedron)
			// will be the hinges we rotate around to connect poly1 and poly2.
			// These are also the edges the static poly shares with the other two.

			Vector3D[] iVerts = staticPoly.Vertices.Intersect( poly1.Vertices ).Intersect( poly2.Vertices ).ToArray();
			if( iVerts.Length != 1 )
				throw new System.ArgumentException();
			Vector3D vertex = iVerts.First();

			// Get our hinges and the two edges we want to merge.
			Edge hinge1 = Hinge( staticPoly, poly1 );
			Edge hinge2 = Hinge( staticPoly, poly2 );
			Edge e1 = AdjacentEdge( poly1, hinge1, vertex );
			Edge e2 = AdjacentEdge( poly2, hinge2, vertex );

			// Construct circles through the edges they will rotate through when hinged.
			Circle3D c1 = ConstructRotationCircle( hinge1, vertex, e1.Opp( vertex ) );
			Circle3D c2 = ConstructRotationCircle( hinge2, vertex, e2.Opp( vertex ) );

			// These circles both lie on the same sphere (with center = vertex and radius 1), 
			// and will intersect in two points. Those two points will determine the two 
			// locations we can hinge our polys to.
			Vector3D i1, i2;
			SphericalTrig.IntersectionSmart( vertex, c1, c2, out i1, out i2 );

			// XXX - Need to think more about how to choose among options.
			double a1 = (e1.Opp( vertex ) - c1.Center).AngleTo( i1 - c1.Center );
			double a2 = (e2.Opp( vertex ) - c2.Center).AngleTo( i1 - c2.Center );
			RotatePoly( poly1, EdgeAxis( hinge1 ), -a1 );
			RotatePoly( poly2, EdgeAxis( hinge2 ), a2 );
		}

		private static void RotatePoly( Polygon poly, Vector3D axis, double rot )
		{
			foreach( Segment seg in poly.Segments )
			{
				Vector3D p1 = seg.P1, p2 = seg.P2;
				p1.RotateAboutAxis( axis, rot );
				p2.RotateAboutAxis( axis, rot );
				seg.P1 = p1;
				seg.P2 = p2;
			}
		}

		private static Circle3D ConstructRotationCircle( Edge hinge, Vector3D hingePoint, Vector3D polyPoint )
		{
			Vector3D hingeAxis = EdgeAxis( hinge );
			Circle3D c = new Circle3D
			{
				Center = Euclidean3D.ProjectOntoLine( hingeAxis, hingePoint, polyPoint ),
				Normal = hingeAxis,
				Radius = Euclidean3D.DistancePointLine( hingeAxis, hingePoint, polyPoint )
			};

			if( c.Normal.Dot( c.Center ) < 0 )
				c.Normal *= -1;

			return c;
		}

		private static Vector3D EdgeMidpoint( Edge e )
		{
			return (e.Start + e.End) / 2;
		}

		private static Vector3D EdgeAxis( Edge e )
		{
			return e.End - e.Start;
		}

		/// <summary>
		/// Returns the edge adjacent to e and connected to v.
		/// </summary>
		private static Edge AdjacentEdge( Polygon poly, Edge e, Vector3D v )
		{
			Edge[] edges = poly.Segments.Select( s => SegToEdge( s ) ).ToArray();
			return AdjacentEdge( edges, e, v );
		}

		private static Edge AdjacentEdge( Edge[] edges, Edge e, Vector3D v )
		{
			H3.Cell.EdgeEqualityComparer comparer = new H3.Cell.EdgeEqualityComparer();
			foreach( Edge edge in edges )
			{
				if( comparer.Equals( e, edge ) )
					continue;
				if( edge.Start == v || edge.End == v )
					return edge;
			}

			return null;
		}

		private static Edge Hinge( Polygon p1, Polygon p2 )
		{
			Segment[] intersection = p1.Segments.Intersect( p2.Segments, new SegmentEqualityComparer() ).ToArray();
			if( intersection.Length != 1 )
				throw new System.ArgumentException();
			return SegToEdge( intersection.First() );
		}

		private static Edge SegToEdge( Segment s )
		{
			return new Edge( s.P1, s.P2 );
		}

		internal class SegmentEqualityComparer : IEqualityComparer<Segment>
		{
			public bool Equals( Segment s1, Segment s2 )
			{
				if( s1.Type == SegmentType.Arc || s2.Type == SegmentType.Arc )
					throw new System.NotImplementedException();

				Edge e1 = new Edge( s1.P1, s1.P2 );
				Edge e2 = new Edge( s2.P1, s2.P2 );
				return m_comparer.Equals( e1, e2 );
			}

			public int GetHashCode( Segment s )
			{
				Edge e = SegToEdge( s );
				return m_comparer.GetHashCode( e );
			}

			private readonly H3.Cell.EdgeEqualityComparer m_comparer = new H3.Cell.EdgeEqualityComparer();
		}

		/// <summary>
		/// Adds a polygon to our polyhedron.
		/// This will update our data structures, and so should only be called once the polygon
		/// is in its final position!
		/// NOTE: Don't call multiple times with the same polygon.  This isn't smart about that.
		/// </summary>
		private void AddPoly( Polygon poly )
		{
			m_polygons.Add( poly );

			// Update all the maps.
			UpdateVertexMap( m_vertexMap, poly );
			UpdateEdgeMap( m_edgeMap, poly );
			UpdateVertexToEdgeMap( m_vertexToEdgeMap, poly );
		}

		private static void UpdateVertexMap( VertexMap vMap, Polygon poly )
		{
			foreach( Vector3D v in poly.Vertices )
			{
				List<Polygon> polys;
				if( !vMap.TryGetValue( v, out polys ) )
				{
					polys = new List<Polygon>();
					vMap[v] = polys;
				}

				polys.Add( poly );
			}
		}

		private static void UpdateEdgeMap( EdgeMap eMap, Polygon poly )
		{
			foreach( Segment s in poly.Segments )
			{
				Edge e = SegToEdge( s );
				List<Polygon> polys;
				if( !eMap.TryGetValue( e, out polys ) )
				{
					polys = new List<Polygon>();
					eMap[e] = polys;
				}

				polys.Add( poly );
			}
		}

		private static void UpdateVertexToEdgeMap( VertexToEdgeMap map, Polygon poly )
		{
			foreach( Segment s in poly.Segments )
			{
				Edge e = SegToEdge( s );
				Vector3D[] verts = new Vector3D[] { e.Start, e.End };
				foreach( Vector3D v in verts )
				{
					HashSet<Edge> edges;
					if( !map.TryGetValue( v, out edges ) )
					{
						edges = new HashSet<Edge>( new H3.Cell.EdgeEqualityComparer() );
						map[v] = edges;
					}

					edges.Add( e );
				}
			}
		}
	}
}
