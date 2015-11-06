namespace HyperbolicModels
{
	using R3.Geometry;
	using System.Collections.Generic;
	using System.Linq;
	using Math = System.Math;
	using R3.Core;

	public class Lawson
	{
		public Lawson( int m, int k )
		{
			m_m = m;
			m_k = k;
		}

		private int m_m;
		private int m_k;

		public void Gen()
		{
			Vector3D[] verts = Verts( m_m, m_k );
			Sphere[] tet = Tetrahedron( verts );
			Quad quad = new Quad() { Verts = verts };

			// We need to avoid infinities.
			//tet = tet.Select( s => H3Models.UHSToBall( s ) ).ToArray();
			//for( int i=0; i<quad.Verts.Length; i++ )
			//	quad.Verts[i] = H3Models.UHSToBall( quad.Verts[i] );

			// Reflect it around.
			//Quad[] quads = new Quad[] { quad };
			Quad[] quads = CalcQuads( tet, quad );

			List<H3.Cell.Edge> edges = new List<H3.Cell.Edge>();
			foreach( Quad q in quads )
			{
				q.R3toS3();
				edges.AddRange( q.GenEdges() );
			}

			string filename = string.Format( "lawson_{0}_{1}.pov", m_m, m_k );
			PovRay.WriteEdges( new PovRay.Parameters() { AngularThickness = 0.01 }, 
				Geometry.Spherical, edges.ToArray(), filename, append: false );
		}

		private class Quad
		{
			public Vector3D[] Verts;

			public Vector3D ID
			{
				get
				{
					if( Verts == null || Verts.Length != 4 )
						throw new System.Exception( "Quad not initialized." );

					Vector3D ID = new Vector3D();
					for( int i = 0; i < Verts.Length; i++ )
					{
						Vector3D vert = Verts[i];
						if( Infinity.IsInfinite( vert ) )
							ID += new Vector3D( 100, 100, 100 );
						else
							ID += vert;
					}
					return ID;
				}
			}

			public Quad Clone()
			{
				Quad newQuad = new Quad();
				newQuad.Verts = (Vector3D[])Verts.Clone();
				return newQuad;
			}

			public void R3toS3()
			{
				for( int i = 0; i < Verts.Length; i++ )
					Verts[i] = Sterographic.R3toS3( Verts[i] );
			}

			public H3.Cell.Edge[] GenEdges()
			{
				List<H3.Cell.Edge> result = new List<H3.Cell.Edge>();
				result.AddRange( GenEdgesInternal( Verts[0], Verts[1], Verts[3], Verts[2] ) );	// tricky to figure out indices
				result.AddRange( GenEdgesInternal( Verts[0], Verts[3], Verts[1], Verts[2] ) );
				return result.ToArray();
			}

			private H3.Cell.Edge[] GenEdgesInternal( Vector3D v1, Vector3D v2, Vector3D v3, Vector3D v4 )
			{
				int div = 5;
				List<H3.Cell.Edge> result = new List<H3.Cell.Edge>();
				for( int i = 0; i <= div; i++ )
				{
					Vector3D start = v1 + ( v2 - v1 ) * i / div;
					Vector3D end = v3 + ( v4 - v3 ) * i / div;
					start.Normalize();
					end.Normalize();

					start = Sterographic.S3toR3( start );
					end = Sterographic.S3toR3( end );

					if( Infinity.IsInfinite( start ) )
						start = end * 10;
					if( Infinity.IsInfinite( end ) )
						end = start * 10;

					result.Add( new H3.Cell.Edge( start, end ) );
				}
				return result.ToArray();
			}
		}

		/// <summary>
		/// Calculate the 4 points defining the fundamental geodesic quadrilateral.
		/// </summary>
		private static Vector3D[] Verts( int m, int k )
		{
			double dist1 = Math.PI / ( m + 1 );
			double dist2 = Math.PI / ( k + 1 );

			Vector3D p4 = new Vector3D( 1, 0 );
			p4.RotateXY( dist2 );

			return new Vector3D[]
			{
				new Vector3D(),
				new Vector3D( 1, 0 ),
				new Vector3D( 0, 0, Spherical2D.s2eNorm( dist1 ) ),
				p4
			};
		}

		/// <summary>
		/// Calculate the surfaces of the tetrahedron for the quadrilateral
		/// </summary>
		private static Sphere[] Tetrahedron( Vector3D[] quad )
		{
			// NOTE: For planes, the convention is that the normal vector points "outward"

			// The only non-planar sphere is incident with the last two points of the quad,
			// as well as the antipode of the last point. Fig 1 in Lawson's paper.
			Circle3D c = new Circle3D( quad[3], quad[2], -quad[3] );

			Vector3D plane3 = quad[3];
			plane3.RotateXY( Math.PI / 2 );

			return new Sphere[]
			{
				Sphere.Plane( new Vector3D( 0, -1 ) ),
				Sphere.Plane( new Vector3D( 0, 0, -1 ) ),
				Sphere.Plane( plane3 ),
				new Sphere( c.Center, c.Radius )
			};
		}

		private static Quad[] CalcQuads( Sphere[] mirrors, Quad start )
		{
			List<Quad> allQuads = new List<Quad>();
			allQuads.Add( start );
			HashSet<Vector3D> completed = new HashSet<Vector3D>( new Vector3D[] { start.ID } );
			ReflectEdgesRecursive( mirrors, new Quad[] { start }, allQuads, completed );
			return allQuads.ToArray();
		}

		private static void ReflectEdgesRecursive( Sphere[] mirrors, Quad[] quads, 
			List<Quad> allQuads, HashSet<Vector3D> completed )
		{
			if( 0 == quads.Length )
				return;

			List<Quad> newQuads = new List<Quad>();

			foreach( Quad quad in quads )
			//foreach( Sphere mirror in mirrors )
			{
				Sphere mirror = mirrors[3];
				Quad newQuad = quad.Clone();

				for( int i = 0; i < newQuad.Verts.Length; i++ )
					newQuad.Verts[i] = mirror.ReflectPoint( newQuad.Verts[i] );

				if( completed.Add( newQuad.ID ) )
				{
					// Haven't seen this yet, so 
					// we'll need to recurse on it.
					allQuads.Add( newQuad );
					newQuads.Add( newQuad );
				}
			}

			//ReflectEdgesRecursive( mirrors, newQuads.ToArray(), allQuads, completed );
		}
	}
}
