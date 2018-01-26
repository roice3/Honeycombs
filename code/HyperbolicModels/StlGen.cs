namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Drawing;
	using R3.Geometry;
	using R3.Math;
	using Math = System.Math;

	public class StlGen
	{
		static int m_div = 7;
		static double m_thresh = 0.1;//0.004;

		public static void Helicoid()
		{
			Mesh thinMesh = new Mesh();

			H3Ruled ruled = new H3Ruled();
			H3.Cell.Edge[] edgesBall = ruled.Helicoid();
			List<Vector3D[]> boundaryPoints = new List<Vector3D[]>();
			List<Vector3D> starts = new List<Vector3D>();
			List<Vector3D> ends = new List<Vector3D>();
			for( int i = 0; i < edgesBall.Length - 1; i++ )
			{
				int idx1 = i;
				int idx2 = i + 1;
				H3.Cell.Edge e1 = edgesBall[idx1];
				H3.Cell.Edge e2 = edgesBall[idx2];

				int div = 50;
				Vector3D[] points1 = H3Models.Ball.GeodesicPoints( e1.Start, e1.End, div );
				Vector3D[] points2 = H3Models.Ball.GeodesicPoints( e2.Start, e2.End, div );

				thinMesh.AddBand( points1, points2, close: false );

				starts.Add( e1.Start );
				ends.Add( e1.End );
				if( idx1 == 0 )
					boundaryPoints.Add( points1 );
				if( idx2 == edgesBall.Length - 1 )
				{
					boundaryPoints.Add( points2 );
					starts.Add( e2.Start );
					ends.Add( e2.End );
				}
			}
			starts.Reverse();
			boundaryPoints.Add( starts.ToArray() );
			boundaryPoints.Add( ends.ToArray() );

			// Build a normal map.
			Dictionary<Vector3D, Vector3D> normalMap = new Dictionary<Vector3D, Vector3D>();
			{
				Dictionary<Vector3D, List<Mesh.Triangle>> nodesToTriangles = new Dictionary<Vector3D, List<Mesh.Triangle>>();
				System.Action<Vector3D,Mesh.Triangle> doOne = ( v, t ) =>
				{
					List<Mesh.Triangle> connected;
					if( !nodesToTriangles.TryGetValue( v, out connected ) )
						connected = new List<Mesh.Triangle>();
					connected.Add( t );
					nodesToTriangles[v] = connected;
				};

				foreach( Mesh.Triangle tri in thinMesh.Triangles )
				{
					doOne( tri.a, tri );
					doOne( tri.b, tri );
					doOne( tri.c, tri );
				}

				foreach( Vector3D v in nodesToTriangles.Keys )
				{
					var tris = nodesToTriangles[v];
					Vector3D n = new Vector3D();
					foreach( Mesh.Triangle t in tris )
						n += t.Normal;
					n.Normalize();
					normalMap[v] = -n;
				}
			}

			Mesh fullMesh = new Mesh();
			foreach( Mesh.Triangle tri in thinMesh.Triangles )
			{
				Mesh.Triangle[] thickened = ThickenInBallSimple( tri, normalMap );
				fullMesh.Triangles.AddRange( thickened );
			}
			//fullMesh.Append( thinMesh );

			System.Func<Vector3D, System.Tuple<Vector3D, Vector3D>> thickenFn = v => ThickenInBallSimple( v, normalMap );
			fullMesh.Append( ThickenBoundary( boundaryPoints[0], thickenFn ) );
			fullMesh.Append( ThickenBoundary( boundaryPoints[2], thickenFn ) );
			fullMesh.Append( ThickenBoundary( boundaryPoints[3], thickenFn ) );
			var temp = ThickenBoundary( boundaryPoints[1], thickenFn );
			ReverseTris( temp );
			fullMesh.Append( temp );

			string filename = "helicoid.stl";
			System.IO.File.Delete( filename );
			using( StreamWriter sw = File.AppendText( filename ) )
			{
				STL.AppendMeshToSTL( fullMesh, sw );

				Vector3D aStart = H3Ruled.Transform( new Vector3D( 0, 0, -1 ) );
				Vector3D aEnd = H3Ruled.Transform( new Vector3D( 0, 0, 1 ) );
				Mesh m3 = new Mesh();
				AddEdge( m3, aStart, aEnd );
				STL.AppendMeshToSTL( m3, sw );
			}
		}

		public static void HoneycombFiniteVertexFig( HoneycombDef def )
		{
			// This will be used to remove duplicates.
			// Our faces will be doubled-up.  We'll make a hash from one of the interior meshPoints.
			Dictionary<Vector3D, H3.Cell> complete = new Dictionary<Vector3D, H3.Cell>();

			m_thresh = 0.07;
			HoneycombFiniteVertexFig( def, 3, complete );

			m_thresh = 0.02;
			HoneycombFiniteVertexFig( def, 2, complete );

			m_thresh = 0.007;
			HoneycombFiniteVertexFig( def, 1, complete );

			//m_thresh = 0.005;
			//CreateHoneycombSTL( def, 0, complete );

			string filename = "cell.stl";
			System.IO.File.Delete( filename );
			using( StreamWriter sw = File.AppendText( filename ) )
			{
				foreach( H3.Cell cell in complete.Values )
				{
					//STL.AppendMeshToSTL( cell.Mesh, sw );

					bool reverse = cell.Depths.Sum() % 2 == 1;

					Mesh m = new Mesh();
					Sphere normal = cell.Facets[0].Sphere;
					foreach( Mesh.Triangle tri in cell.Mesh.Triangles )
					{
						Mesh.Triangle[] thickened = ThickenInBallSimple( tri, normal );
						m.Triangles.AddRange( thickened );
					}

					if( reverse )
						ReverseTris( m );

					STL.AppendMeshToSTL( m, sw );

					System.Func<Vector3D, System.Tuple<Vector3D, Vector3D>> thickenFn = v => ThickenInBallSimple( v, normal );
					int stride = (int)Math.Sqrt( cell.Mesh.Triangles.Count ) + 1;
					Vector3D[] e1 = cell.AuxPoints.Skip( 1 + 0 * stride ).Take( stride ).ToArray();
					Vector3D[] e2 = cell.AuxPoints.Skip( 1 + 1 * stride ).Take( stride ).ToArray();
					Vector3D[] e3 = cell.AuxPoints.Skip( 1 + 2 * stride ).Take( stride ).ToArray();
					Mesh m1 = ThickenBoundary( e1, thickenFn ), m2 = ThickenBoundary( e2, thickenFn ), m3 = ThickenBoundary( e3, thickenFn );
					if( reverse )
					{
						ReverseTris( m1 );
						ReverseTris( m2 );
						ReverseTris( m3 );
					}

					STL.AppendMeshToSTL( m1, sw );
					STL.AppendMeshToSTL( m2, sw );
					STL.AppendMeshToSTL( m3, sw );
				}
			}
		}

		private static void HoneycombFiniteVertexFig( HoneycombDef def, int lod, Dictionary<Vector3D, H3.Cell> complete )
		{
			int p = def.P;
			int q = def.Q;
			int r = def.R;

			double scale = 1.0;
			Vector3D vUHS = H3Models.BallToUHS( SimplexCalcs.VertexPointBall( p, q, r ) );
			if( Geometry2D.GetGeometry( q, r ) != Geometry.Hyperbolic ) // Vertex-centered if possible
				scale = 1.0 / vUHS.Z;
			System.Func<Vector3D, Vector3D> trans = v =>
			{
				v = H3Models.BallToUHS( v );
				v *= scale;
				v = H3Models.UHSToBall( v );
				return v;
			};

			bool ball = true;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball );
			simplex = simplex.Select( s => 
			{
				s = H3Models.BallToUHS( s );
				Sphere.ScaleSphere( s, scale );
				s = H3Models.UHSToBall( s );
				return s;
			} ).ToArray();
			H3.Cell.Edge[] edges = SimplexCalcs.SimplexEdgesBall( p, q, r );

			// Two edges of the simplex facet.
			// NOTE: This contruction only works for material triangles, and matches the construction in the TextureHelper.
			m_div = TextureHelper.SetLevels( lod );
			int[] elementIndices = TextureHelper.TextureElements( 1, lod );
			List<Vector3D> points = new List<Vector3D>();
			H3.Cell.Edge e1 = edges[2];
			H3.Cell.Edge e2 = edges[3];
			Vector3D p1 = trans( e1.Start ), p2 = trans( e1.End ), p3 = trans( e2.End );
			Vector3D[] points1 = H3Models.Ball.GeodesicPoints( p2, p1, m_div );
			Vector3D[] points2 = H3Models.Ball.GeodesicPoints( p3, p1, m_div );
			for( int i = 0; i < m_div; i++ )
				points.AddRange( H3Models.Ball.GeodesicPoints( points1[i], points2[i], m_div - i ) );
			points.Add( p1 );

			Mesh mesh = new Mesh();
			for( int i = 0; i < elementIndices.Length / 3; i++ )
			{
				int idx1 = i * 3;
				int idx2 = i * 3 + 1;
				int idx3 = i * 3 + 2;
				Vector3D v1 = points[elementIndices[idx1]];
				Vector3D v2 = points[elementIndices[idx2]];
				Vector3D v3 = points[elementIndices[idx3]];
				mesh.Triangles.Add( new Mesh.Triangle( v1, v2, v3 ) );
			}

			// AuxPoints will be used for multiple things.
			// - The first is a definition point for a face, so we can check for duplicates.
			// - We'll also store the points for the 3 edges of our fundamental triangle.
			List<Vector3D> auxPoints = new List<Vector3D>();
			{
				auxPoints.Add( (p1 + p2 + p3) / 3 );
				auxPoints.AddRange( points1 );
				auxPoints.AddRange( points2.Reverse() );
				auxPoints.AddRange( H3Models.Ball.GeodesicPoints( points2[0], points1[0], m_div ) );
			}

			Vector3D cen = HoneycombPaper.InteriorPointBall;
			H3.Cell[] simplices = GenCell( simplex, mesh, cen, auxPoints.ToArray(), ball );

			// Existing cells take precedence.
			foreach( H3.Cell c in simplices )
			{
				Vector3D t = c.AuxPoints[0];
				H3.Cell dummy;
				if( !complete.TryGetValue( t, out dummy ) )
					complete[t] = c;
			}
		}

		public static void HoneycombHyperidealLegs( HoneycombDef def )
		{
			// This will be used to avoid duplicates.
			// The key is the cell center.
			Dictionary<Vector3D, H3.Cell> complete = new Dictionary<Vector3D, H3.Cell>();

			m_thresh = 0.05;
			//m_thresh = 0.07;
			HoneycombHyperidealLegs( def, 3, complete );

			m_thresh = 0.01;
			//m_thresh = 0.02;
			HoneycombHyperidealLegs( def, 2, complete );

			m_thresh = 0.004;
			//m_thresh = 0.007;
			HoneycombHyperidealLegs( def, 1, complete );

			string filename = "cell.stl";
			System.IO.File.Delete( filename );
			using( StreamWriter sw = File.AppendText( filename ) )
			{
				HashSet<H3.Cell.Edge> edgesToMesh = new HashSet<H3.Cell.Edge>( new H3.Cell.EdgeEqualityComparer() );
				foreach( H3.Cell cell in complete.Values )
				{
					int depth = cell.Depths.Sum();

					Mesh m = new Mesh();
					Sphere normal = cell.Facets[0].Sphere;
					foreach( Mesh.Triangle tri in cell.Mesh.Triangles )
					{
						Mesh.Triangle[] thickened = Thicken( tri, normal );
						m.Triangles.AddRange( thickened.Select( t => Transform( t ) ) );
					}

					List<object> boundary = new List<object>();
					int skip = 2;
					int stride = (int)Math.Sqrt( cell.Mesh.Triangles.Count )/2 + 1;
					int num = stride;
					boundary.Add( cell.AuxPoints.Skip( skip ).Take( num ) );
					skip += num;
					boundary.Add( cell.AuxPoints.Skip( skip ).Take( num ) );
					skip += num;
					num = 2 * m_div + 1;
					boundary.Add( cell.AuxPoints.Skip( skip ).Take( num ) );
					skip += num;
					boundary.Add( cell.AuxPoints.Skip( skip ).Take( num ) );
					skip += num;
					foreach( object e in boundary )
					{
						var enumerable = (IEnumerable<Vector3D>)e;
						//if( depth % 2 == 0 )
						//	enumerable = enumerable.Reverse();

						Mesh m2 = ThickenBoundary( enumerable.ToArray(), normal );

						if( depth % 2 == 0 )
							ReverseTris( m2 );

						m.Triangles.AddRange( m2.Triangles.Select( t => Transform( t ) ) );
					}

					STL.AppendMeshToSTL( m, sw );
					edgesToMesh.Add( new H3.Cell.Edge( cell.AuxPoints[0], cell.AuxPoints[1] ) );
				}

				/*foreach( H3.Cell.Edge e in edgesToMesh )
				{
					Mesh m3 = new Mesh();
					AddEdge( m3, Transform( e.Start ), Transform( e.End ) );
					STL.AppendMeshToSTL( m3, sw );
				}*/
			}
		}

		/// <summary>
		/// Create an STL file for a cell.
		/// Currently only works for cells with both hyperideal vertices and cells.
		/// </summary>
		public static void HoneycombHyperidealLegs( HoneycombDef def, int lod, Dictionary<Vector3D, H3.Cell> complete )
		{
			int p = def.P;
			int q = def.Q;
			int r = def.R;

			m_div = TextureHelper.SetLevels( lod );

			bool ball = false;
			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: ball );
			H3.Cell.Edge[] edges;
			if( ball )
				edges = SimplexCalcs.SimplexEdgesBall( p, q, r );
			else
				edges = SimplexCalcs.SimplexEdgesUHS( p, q, r );

			// Two edges of one simplex facet.
			H3.Cell.Edge e1 = edges[2];
			H3.Cell.Edge e2 = edges[3];
			Vector3D[] points1, points2;
			if( ball )
			{
				points1 = H3Models.Ball.GeodesicPoints( e1.Start, e1.End, 2 * m_div );
				points2 = H3Models.Ball.GeodesicPoints( e2.Start, e2.End, 2 * m_div );
			}
			else
			{
				points1 = H3Models.UHS.GeodesicPoints( e1.Start, e1.End, 2 * m_div );
				points2 = H3Models.UHS.GeodesicPoints( e2.Start, e2.End, 2 * m_div );
			}

			Sphere cellSphere = simplex[0];
			Sphere vertexSphere = simplex[3];

			// Because one vertex the facet triangle is hyperideal, it will actually look like a square.
			List<Vector3D[]> allPoints = new List<Vector3D[]>();
			for( int i = 0; i < points1.Length; i++ )
			{
				Vector3D p1 = points1[i];
				Vector3D p2 = points2[i];

				Vector3D[] arcPoints;
				if( i == points1.Length - 1 )
				//if( false )
				{
					// NOTE: This arc is not generally geodesic!
					// Or is it?
					arcPoints = ball ?
						H3Models.Ball.GeodesicPoints( p1, p2, m_div ) :
						H3Models.UHS.GeodesicPoints( p1, p2, m_div );

					/*Circle3D arc = cellSphere.Intersection( vertexSphere );
					double angleTot = (p1 - arc.Center).AngleTo( p2 - arc.Center );
					arcPoints = Shapeways.CalcArcPoints( arc.Center, arc.Radius, p1, arc.Normal, -angleTot, div );*/
				}
				else
				{
					Circle3D c = Circle3D.FromCenterAnd2Points( cellSphere.Center, p1, p2 );
					double angleTot = (p1 - c.Center).AngleTo( p2 - c.Center );
					arcPoints = Shapeways.CalcArcPoints( cellSphere.Center, cellSphere.Radius, p1, c.Normal, -angleTot, m_div );
				}
				//Vector3D[] arcPoints = new Vector3D[] { p1, p2 };
				allPoints.Add( arcPoints );
			}

			// Create the triangles for the patch.
			Mesh mesh = new Mesh();
			for( int i = 0; i < allPoints.Count - 1; i++ )
			{
				Vector3D[] arc1 = allPoints[i];
				Vector3D[] arc2 = allPoints[i + 1];

				for( int j = 0; j < arc1.Length - 1; j++ )
				{
					// Points of (i,j) box;
					Vector3D p1 = arc1[j];
					Vector3D p2 = arc2[j];
					Vector3D p3 = arc1[j + 1];
					Vector3D p4 = arc2[j + 1];

					Mesh.Triangle tri1 = new Mesh.Triangle( p1, p2, p3 );
					Mesh.Triangle tri2 = new Mesh.Triangle( p2, p4, p3 );

					// We need to thicken after reflecting around, otherwise we can't apply a min thickness.
					/*Sphere normal = cellSphere;
					Mesh.Triangle[] thickened1 = Thicken( tri1, normal );
					Mesh.Triangle[] thickened2 = Thicken( tri2, normal );
					mesh.Triangles.AddRange( thickened1 );
					mesh.Triangles.AddRange( thickened2 );*/

					mesh.Triangles.Add( tri1 );
					mesh.Triangles.Add( tri2 );
				}
			}

			// AuxPoints will be used for multiple things.
			// - The first two points are for an an that will fill the gap where there is a missing face.
			// - We'll also store the points for the 4 edges of our fundamental triangle.
			List<Vector3D> auxPoints = new List<Vector3D>();
			{
				var edge1 = allPoints.First();
				var edge2 = allPoints.Last();
				List<Vector3D> edge3 = new List<Vector3D>(), edge4 = new List<Vector3D>();
				for( int i = 0; i < allPoints.Count; i++ )
				{
					edge3.Add( allPoints[i][0] );
					edge4.Add( allPoints[i][allPoints[i].Length - 1] );
				}
				edge4.Reverse();

				auxPoints.Add( e1.Start );
				auxPoints.Add( e1.End );
				auxPoints.AddRange( edge1.Reverse() );
				auxPoints.AddRange( edge2 );
				auxPoints.AddRange( edge3 );
				auxPoints.AddRange( edge4 );
			}

			Vector3D cen = HoneycombPaper.InteriorPointBall;

			/* Reorientation code.  Move this elsewhere.
			 
			// Face centered orientation.
			bool faceCentered = false;
			if( faceCentered )
				SimplexCalcs.PrepForFacetCentering( p, q, simplex, ref cen );

			Mobius mUHS = SimplexCalcs.FCOrientMobius( p, q );
			Mobius mBall = HoneycombPaper.FCOrientMobius( H3Models.UHSToBall( cellSphere ) );

			simplex = simplex.Select( s =>
			{
				s = H3Models.UHSToBall( s );
				//H3Models.TransformInBall2( s, mBall );
				return s;
			} ).ToArray();

			
			{
				for( int i = 0; i < mesh.Triangles.Count; i++ )
				{
					Mesh.Triangle tri = mesh.Triangles[i];

					if( faceCentered )
					{
						tri.a = mUHS.ApplyToQuaternion( tri.a );
						tri.b = mUHS.ApplyToQuaternion( tri.b );
						tri.c = mUHS.ApplyToQuaternion( tri.c );
					}

					tri.a = H3Models.UHSToBall( tri.a );
					tri.b = H3Models.UHSToBall( tri.b );
					tri.c = H3Models.UHSToBall( tri.c );

					if( faceCentered )
					{
						tri.a = H3Models.TransformHelper( tri.a, mBall );
						tri.b = H3Models.TransformHelper( tri.b, mBall );
						tri.c = H3Models.TransformHelper( tri.c, mBall );
					}
					mesh.Triangles[i] = tri;
				}

				if( faceCentered )
					cen = H3Models.TransformHelper( cen, mBall );
			}
			*/

			// Now we need to reflect around this fundamental patch.
			H3.Cell[] simplices = GenCell( simplex, mesh, cen, auxPoints.ToArray(), ball );

			// Existing cells take precedence.
			foreach( H3.Cell c in simplices )
			{
				Vector3D t = c.Center;
				H3.Cell dummy;
				if( !complete.TryGetValue( t, out dummy ) )
					complete[t] = c;
			}
		}

		private static void ReverseTris( Mesh m )
		{
			for( int i = 0; i < m.Triangles.Count; i++ )
			{
				Mesh.Triangle tri = m.Triangles[i];    // Ugh, damn you structs!
				tri.ChangeOrientation();
				m.Triangles[i] = tri;
			}
		}

		private static Vector3D Transform( Vector3D v )
		{
			//v *= 1.35;
			//v *= 2;
			v = H3Models.UHSToBall( v );
			return v;
		}

		private static Mesh.Triangle Transform( Mesh.Triangle tri )
		{
			tri.a = Transform( tri.a );
			tri.b = Transform( tri.b );
			tri.c = Transform( tri.c );
			return tri;
		}

		private static Mesh ThickenBoundary( Vector3D[] edge, System.Func<Vector3D, System.Tuple<Vector3D,Vector3D>> thickenFn )
		{
			List<Vector3D> side1 = new List<Vector3D>();
			List<Vector3D> side2 = new List<Vector3D>();
			foreach( Vector3D v in edge )
			{
				var thickened = thickenFn( v );
				side1.Add( thickened.Item1 );
				side2.Add( thickened.Item2 );
			}

			Mesh m = new Mesh();
			m.AddBand( side1.ToArray(), side2.ToArray(), close: false );
			return m;
		}

		//
		// Used to connect up the four thickened edges of the region.
		//
		private static Mesh ThickenBoundary( Vector3D[] edge, Sphere normal )
		{
			System.Func<Vector3D, System.Tuple<Vector3D, Vector3D>> fn = v => Thicken( v, normal );
			return ThickenBoundary( edge, fn );
		}

		/// <summary>
		/// Ball
		/// </summary>
		private static void AddEdge( Mesh mesh, Vector3D v1, Vector3D v2 )
		{
			double quality = v1.Dist( v2 ) * 8;	// roughly 0 to 1.
			int div = (int)(quality * 16);
			if( div > 16 )
				div = 16;
			if( div < 6 )
				div = 6;
			div = 50;

			Vector3D[] points = H3Models.Ball.GeodesicPoints( v1, v2, div );
			List<Vector3D> ePoints = new List<Vector3D>();
			List<double> eRadii = new List<double>();
			foreach( Vector3D pNE in points )
			{
				Sphere sphere = SphereFuncBall( pNE, false );
				ePoints.Add( sphere.Center );
				eRadii.Add( sphere.Radius );
			}

			Shapeways shapeways = new Shapeways();
			shapeways.AddCurve( ePoints.ToArray(), eRadii.ToArray() );
			mesh.Append( shapeways.Mesh );
		}

		internal static Sphere SphereFuncBall( Vector3D v, bool simple = true )
		{
			//bool simple = true;
			if( simple )
			{
				double size = 125;
				double thickCen = 2.0 / size;
				//double thickEdge = .7 / size;
				double thickEdge = 1.0 / size;
				double thick = thickCen - v.Abs() * (thickCen - thickEdge);
				return new Sphere() { Center = v, Radius = thick };
			}
			else
			{
				double thick = .1;
				//double minRad = 0.8 / 100;
				double minRad = 1.0 / 125;
				Vector3D c;
				double r;
				H3Models.Ball.DupinCyclideSphere( v, thick / 2, Geometry.Hyperbolic, out c, out r );
				return new Sphere() { Center = c, Radius = Math.Max( r, minRad ) };
			}
		}

		internal static Sphere SphereFuncUHS( Vector3D v )
		{
			return H3Models.BallToUHS( SphereFuncBall( H3Models.UHSToBall( v ) ) );
		}

		/// <summary>
		/// A simple thickening. This is not accurate, but meant to save cost by minimizing wall thickness everywhere.
		/// The input vector should live on normal.
		/// </summary>
		private static System.Tuple<Vector3D, Vector3D> ThickenInBallSimple( Vector3D v, Sphere normal )
		{
			double size = 125;
			double thickCen = 2.0 / size;
			//double thickEdge = .7 / size;
			double thickEdge = 1.0 / size;
			double thick = thickCen - v.Abs() * (thickCen - thickEdge);
			Vector3D direction = v - normal.Center;
			return ThickenInBallSimple( v, direction, thick );
		}

		private static System.Tuple<Vector3D,Vector3D> ThickenInBallSimple( Vector3D v, Vector3D direction, double thickness )
		{
			direction.Normalize();
			direction *= thickness;
			return new System.Tuple<Vector3D, Vector3D>( v + direction, v - direction );
		}

		private static System.Tuple<Vector3D, Vector3D> ThickenInBallSimple( Vector3D v, Dictionary<Vector3D, Vector3D> normalMap )
		{
			Sphere s = SphereFuncBall( v );
			return ThickenInBallSimple( v, normalMap[v], s.Radius );
		}

		private static Mesh.Triangle[] ThickenInBallSimple( Mesh.Triangle tri, Sphere normal )
		{
			List<Mesh.Triangle> result = new List<Mesh.Triangle>();

			System.Tuple<Vector3D, Vector3D> a = ThickenInBallSimple( tri.a, normal );
			System.Tuple<Vector3D, Vector3D> b = ThickenInBallSimple( tri.b, normal );
			System.Tuple<Vector3D, Vector3D> c = ThickenInBallSimple( tri.c, normal );

			result.Add( new Mesh.Triangle( a.Item1, c.Item1, b.Item1 ) );   // For consistent orientations
			result.Add( new Mesh.Triangle( a.Item2, b.Item2, c.Item2 ) );

			return result.ToArray();
		}

		private static Mesh.Triangle[] ThickenInBallSimple( Mesh.Triangle tri, Dictionary<Vector3D,Vector3D> normalMap )
		{
			List<Mesh.Triangle> result = new List<Mesh.Triangle>();

			System.Tuple<Vector3D, Vector3D> a = ThickenInBallSimple( tri.a, normalMap );
			System.Tuple<Vector3D, Vector3D> b = ThickenInBallSimple( tri.b, normalMap );
			System.Tuple<Vector3D, Vector3D> c = ThickenInBallSimple( tri.c, normalMap );

			result.Add( new Mesh.Triangle( a.Item1, c.Item1, b.Item1 ) );   // For consistent orientations
			result.Add( new Mesh.Triangle( a.Item2, b.Item2, c.Item2 ) );

			return result.ToArray();
		}

		/// <summary>
		/// Assumes we are in the UHS model.
		/// Normal to the input sphere means we are along an arc orthogonal to that sphere.
		/// </summary>
		private static System.Tuple<Vector3D, Vector3D> Thicken( Vector3D v, Sphere normal )
		{
			if( Tolerance.Zero( v.Z ) )
			{
				Sphere s = SphereFuncUHS( v );
				Vector3D direction = v - normal.Center;
				direction.Normalize();
				direction *= s.Radius;
				return new System.Tuple<Vector3D, Vector3D>( v + direction, v - direction );
			}
			else
			{
				// We need to get the circle orthogonal to "normal" and z=0 and passing through v.
				// v is on normal.
				Vector3D z = new Vector3D( 0, 0, -1 );
				Vector3D hypot = normal.Center - v;
				double angle = hypot.AngleTo( z );
				double d2 = v.Z / Math.Tan( angle );

				Vector3D cen = v;
				cen.Z = 0;
				Vector3D direction = cen - normal.Center;
				direction.Normalize();
				direction *= d2;
				cen += direction;
				double rad = v.Dist( cen );
				Vector3D n = direction;
				n.Normalize();
				n.RotateXY( Math.PI / 2 );

				// Sphere with non-euclidean center v;
				Sphere s = SphereFuncUHS( v );

				// Two points where this sphere intersects circle.
				var result = s.Intersection( new Circle3D() { Center = cen, Radius = rad, Normal = n } );

				// Our artificial thickening can cause this when s is too big relative to the circle.
				// In this case, we are very close to the boundary, so we can do a different calc.
				if( result == null )
				{
					direction = v - normal.Center;
					direction.Normalize();
					direction *= s.Radius;
					result = new System.Tuple<Vector3D, Vector3D>( v + direction, v - direction );
				}

				return result;
			}
		}

		/// <summary>
		/// Assumes we are in the UHS model.
		/// </summary>
		private static Mesh.Triangle[] Thicken( Mesh.Triangle tri, Sphere normal )
		{
			List<Mesh.Triangle> result = new List<Mesh.Triangle>();

			System.Tuple<Vector3D, Vector3D> a = Thicken( tri.a, normal );
			System.Tuple<Vector3D, Vector3D> b = Thicken( tri.b, normal );
			System.Tuple<Vector3D, Vector3D> c = Thicken( tri.c, normal );

			if( a == null || b == null || c == null ||
				a.Item1.DNE || a.Item2.DNE ||
				b.Item1.DNE || b.Item2.DNE ||
				c.Item1.DNE || c.Item2.DNE )
			{
				System.Console.WriteLine( "bah" );
			}

			result.Add( new Mesh.Triangle( a.Item1, c.Item1, b.Item1 ) );   // For consistent orientations
			result.Add( new Mesh.Triangle( a.Item2, b.Item2, c.Item2 ) );

			return result.ToArray();
		}

		private static H3.Cell[] GenCell( Sphere[] simplex, Mesh mesh, Vector3D cen, Vector3D[] auxPoints, bool ball )
		{
			//Sphere[] mirrors = simplex.Skip(1).ToArray();
			Sphere[] mirrors = simplex.ToArray();
			H3.Cell.Facet[] simplexFacets = simplex.Select( m => new H3.Cell.Facet( m ) ).ToArray();
			H3.Cell startingCell = new H3.Cell( simplexFacets );
			startingCell.Center = cen;
			startingCell.Mesh = mesh;
			startingCell.AuxPoints = auxPoints;
			startingCell = startingCell.Clone();    // So our mirrors don't get munged after we reflect around later.
			H3.Cell[] simplices = Recurse.CalcCells( mirrors, new H3.Cell[] { startingCell }, new Recurse.Settings() { Ball = ball, Threshold = m_thresh } );
			//return simplices.Where( s => s.Depths[0] <= layer /*&& s.Depths[0] == 3 && s.Depths[1] == 3*/ ).ToArray();
			return simplices.ToArray();
		}
	}
}
