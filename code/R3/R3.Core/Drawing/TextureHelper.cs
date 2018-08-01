namespace R3.Drawing
{
	using R3.Geometry;
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using Math = System.Math;

	public class TextureHelper
	{
		public TextureHelper()
		{
			SetLevels( 3 );
		}

		public static int SetLevels( int levels )
		{
			m_maxSubdivisions = (int)Math.Pow( 2, levels );
			return m_maxSubdivisions;
		}

		/// <summary>
		/// Stores the triangle element indices for different levels of detail.
		/// There are 4 entries in the list, and the first entry will have the least detail.
		/// The arrays specify indices into the texture coords, and represent triangle elements.
		/// </summary>
		public List<int[]> ElementIndices { get; set; }

		/// <summary>
		/// Sets up our list of element indices
		/// </summary>
		public void SetupElementIndices( Polygon poly )
		{
			//int numBaseTriangles = poly.Segments.Count == 3 ? 1 : poly.Segments.Count;	// For geodesic saddles.
			int numBaseTriangles = poly.Segments.Count;

			ElementIndices = new List<int[]>();
			for( int i=0; Math.Pow( 2, i ) <= m_maxSubdivisions; i++ )
				ElementIndices.Add( TextureElements( numBaseTriangles, i ) );
		}

		private static int m_maxSubdivisions = 8;	// Must be a power of 2.


		///////////////////////////////////////////////////////////////// PLAYING AROUND WITH GEODESIC SADDLES

		private static Vector3D[] CalcPointsUsingTwoSegments( Segment seg1, Segment seg2, int divisions, Geometry g )
		{
			List<Vector3D> points = new List<Vector3D>();
			Vector3D[] s1 = SubdivideSegmentInGeometry( seg1.P1, seg1.P2, divisions, g );
			Vector3D[] s2 = SubdivideSegmentInGeometry( seg2.P2, seg2.P1, divisions, g );
			for( int i = 0; i < divisions; i++ )
				points.AddRange( SubdivideSegmentInGeometry( s1[i], s2[i], divisions - i, g ) );

			points.Add( seg1.P2 );
			return points.ToArray();
		}

		public static Vector3D[] CalcViaProjections( Vector3D p1, Vector3D p2, Vector3D p3, int divisions, Geometry g )
		{
			if( g == Geometry.Euclidean )
				throw new System.NotImplementedException();

			Vector3D h1 = new Vector3D(), h2 = new Vector3D(), h3 = new Vector3D();
			if( g == Geometry.Hyperbolic )
			{
				h1 = Sterographic.PlaneToHyperboloid( p1 );
				h2 = Sterographic.PlaneToHyperboloid( p2 );
				h3 = Sterographic.PlaneToHyperboloid( p3 );
			}
			else if( g == Geometry.Spherical )
			{
				h1 = Sterographic.PlaneToSphereSafe( p1 );
				h2 = Sterographic.PlaneToSphereSafe( p2 );
				h3 = Sterographic.PlaneToSphereSafe( p3 );
			}

			List<Vector3D> temp = new List<Vector3D>();
			Segment seg1 = Segment.Line( h1, h2 );
			Segment seg2 = Segment.Line( h3, h2 );
			Vector3D[] s1 = seg1.Subdivide( divisions );
			Vector3D[] s2 = seg2.Subdivide( divisions );
			for( int i = 0; i < divisions; i++ )
			{
				Segment seg = Segment.Line( s1[i], s2[i] );
				temp.AddRange( seg.Subdivide( divisions - i ) );
			}
			temp.Add( h2 );

			List<Vector3D> result = new List<Vector3D>();
			foreach( Vector3D v in temp )
			{
				Vector3D copy = v;
				if( g == Geometry.Hyperbolic )
				{
					Sterographic.NormalizeToHyperboloid( ref copy );
					result.Add( Sterographic.HyperboloidToPlane( copy ) );
				}
				else if( g == Geometry.Spherical )
				{
					copy.Normalize();
					result.Add( Sterographic.SphereToPlane( copy ) );
				}
			}
			return result.ToArray();
		}

		private static Vector3D FindClosestPoint( Vector3D v, Vector3D[] list )
		{
			Vector3D result = new Vector3D();

			double dist = double.MaxValue;
			foreach( Vector3D t in list )
			{
				double abs = ( v - t ).Abs();
				if( abs < dist )
				{
					dist = abs;
					result = t;
				}
			}

			return result;
		}

		/////////////////////////////////////////////////////////////////

		

		/// <summary>
		/// Helper to generate a set of texture coordinates.
		/// </summary>
		public static Vector3D[] TextureCoords( Polygon poly, Geometry g, bool doGeodesicDome = false )
		{
			int divisions = m_maxSubdivisions;

			List<Vector3D> points = new List<Vector3D>();
			if( 0 == poly.Segments.Count )
				return points.ToArray();

			// ZZZ - Should we do this different handling of triangles?
			// I think no, this was just for investigating "geodesic saddles".
			bool doGeodesicSaddles = doGeodesicDome;
			if( 3 == poly.Segments.Count && doGeodesicSaddles )
			{
				return CalcViaProjections( poly.Segments[0].P1, poly.Segments[1].P1, poly.Segments[2].P1, divisions, g );
			}
			else
			{
				// We make a triangle lattice for each segment.
				// Think of the segment and the poly center making one big triangle,
				// which is subdivided into smaller triangles.
				foreach( Segment s in poly.Segments )
				{
					Vector3D[] s1 = SubdivideSegmentInGeometry( s.P1, poly.Center, divisions, g );
					Vector3D[] s2 = SubdivideSegmentInGeometry( s.P2, poly.Center, divisions, g );
					for( int i = 0; i < divisions; i++ )
						points.AddRange( SubdivideSegmentInGeometry( s1[i], s2[i], divisions - i, g ) );

					points.Add( poly.Center );
				}
			}

			return points.ToArray();
		}

		/// <summary>
		/// Subdivides a segment from p1->p2 with the two endpoints not on the origin, in the respective geometry.
		/// </summary>
		public static Vector3D[] SubdivideSegmentInGeometry( Vector3D p1, Vector3D p2, int divisions, Geometry g )
		{
			// Handle this specially, so we can keep things 3D if needed.
			if( g == Geometry.Euclidean )
			{
				Segment seg = Segment.Line( p1, p2 );
				return seg.Subdivide( divisions );
			}

			Mobius p1ToOrigin = new Mobius();
			p1ToOrigin.Isometry( g, 0, -p1 );
			Mobius inverse = p1ToOrigin.Inverse();

			Vector3D newP2 = p1ToOrigin.Apply( p2 );
			Segment radial = Segment.Line( new Vector3D(), newP2 );
			Vector3D[] temp = SubdivideRadialInGeometry( radial, divisions, g );

			List<Vector3D> result = new List<Vector3D>();
			foreach( Vector3D v in temp )
				result.Add( inverse.Apply( v ) );

			return result.ToArray();
		}

		/// <summary>
		/// Equally subdivides a segment with a startpoint at the origin, in the respective geometry.
		/// </summary>
		private static Vector3D[] SubdivideRadialInGeometry( Segment radial, int divisions, Geometry g )
		{
			List<Vector3D> result = new List<Vector3D>();
			if( radial.Type != SegmentType.Line )
			{
				Debug.Assert( false );
				return result.ToArray();
			}

			switch( g )
			{
				case Geometry.Spherical:
				{
					double eLength = radial.Length;
					double sLength = Spherical2D.e2sNorm( eLength );
					double divLength = sLength / divisions;

					for( int i = 0; i <= divisions; i++ )
					{
						double temp = Spherical2D.s2eNorm( divLength * i );
						result.Add( radial.P2 * temp / eLength );
					}

					break;
				}
				case Geometry.Euclidean:
					return radial.Subdivide( divisions );

				case Geometry.Hyperbolic:
				{
					double eLength = radial.Length;
					double hLength = DonHatch.e2hNorm( eLength );
					double divLength = hLength / divisions;

					for( int i = 0; i <= divisions; i++ )
					{
						double temp = DonHatch.h2eNorm( divLength * i );
						result.Add( radial.P2 * temp / eLength );
					}

					break;
				}
			}

			return result.ToArray();
		}

		/// <summary>
		/// Returns the sum of all the integers up to and including n.
		/// </summary>
		private static int TriangularNumber( int n )
		{
			return n * (n + 1) / 2;
		}

		/// <summary>
		/// Grabs an array of indices into the coordinate array for TextureCoords.
		/// The array represents individual triangles (each set of 3 is one triangle).
		/// </summary>
		public static int[] TextureElements( int numBaseTriangles, int LOD )
		{
			int divisions = m_maxSubdivisions;
			int stride = divisions / (int)Math.Pow( 2, LOD );

			// 9 + 8 + 7 + 6 + 5 + 4 + 3 + 2 + 1
			int numVertsPerSegment = TriangularNumber( divisions + 1 );

			List<int> result = new List<int>();
			int offset = 0;
			for( int count = 0; count < numBaseTriangles; count++ )
			{
				// Make the triangles.
				int start1 = offset, start2 = offset;
				for( int i = 0; i < divisions; i += stride )
				{
					start1 = start2;

					int temp = divisions - i + 1;
					for( int j = 0; j < stride; j++ )
					{
						start2 += temp;
						temp--;
					}

					for( int j = 0; j < divisions - i; j += stride )
					{
						int idx1 = start1 + j;
						int idx2 = start1 + j + stride;
						int idx3 = start2 + j;

						result.Add( idx1 );
						result.Add( idx2 );
						result.Add( idx3 );
					}

					for( int j = 0; j < divisions - i - stride; j += stride )
					{
						int idx1 = start2 + j;
						int idx2 = start1 + j + stride;
						int idx3 = start2 + j + stride;

						result.Add( idx1 );
						result.Add( idx2 );
						result.Add( idx3 );
					}
				}

				offset += numVertsPerSegment;
			}

			return result.ToArray();
		}

		public static Dictionary<int, int[]> ElementGraph( int numBaseTriangles, int LOD )
		{
			Dictionary<int, int[]> result = new Dictionary<int, int[]>();

			// Brute force.
			int[] textureElements = TextureElements( numBaseTriangles, LOD );
			Dictionary<GraphEdge, List<int>> edgeToTriangles = new Dictionary<GraphEdge, List<int>>();
			for( int i = 0; i < textureElements.Length / 3; i++ )
			{
				int idx1 = i * 3;
				int idx2 = i * 3 + 1;
				int idx3 = i * 3 + 2;

				System.Action<GraphEdge, int> addEdge = (e,idx) =>
				{
					List<int> tris;
					if( !edgeToTriangles.TryGetValue( e, out tris ) )
						tris = new List<int>();
					tris.Add( idx );
					edgeToTriangles[e] = tris;
				};

				addEdge( new GraphEdge( textureElements[idx1], textureElements[idx2] ), i );
				addEdge( new GraphEdge( textureElements[idx2], textureElements[idx3] ), i );
				addEdge( new GraphEdge( textureElements[idx3], textureElements[idx1] ), i );
			}

			Dictionary<int, List<int>> temp = new Dictionary<int, List<int>>();
			System.Action<int, int> addIncident = ( idx1, idx2 ) =>
			{
				List<int> incident;
				if( !temp.TryGetValue( idx1, out incident ) )
					incident = new List<int>();
				incident.Add( idx2 );
				temp[idx1] = incident;
			};

			foreach( var tris in edgeToTriangles.Values )
			{
				if( tris.Count == 1 )
					continue;
				else if( tris.Count == 2 )
				{
					addIncident( tris[0], tris[1] );
					addIncident( tris[1], tris[0] );
				}
				else
					throw new System.Exception();
			}

			int divisions = m_maxSubdivisions;
			for( int i = 0; i < divisions; i++ )
				addIncident( i, -1 );

			foreach( var kvp in temp )
				result[kvp.Key] = kvp.Value.ToArray();
			
			return result;
		}
	}
}
