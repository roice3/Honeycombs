namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;

	/// <summary>
	/// A class to play around with H3 honecombs
	/// </summary>
	public class H3
	{
		/// <summary>
		/// We can track a cell by its ideal vertices.
		/// We'll work with these in the plane.
		/// </summary>
		public class Cell
		{
			public Cell( Facet[] facets ) : this( -1, facets )
			{
			}

			public Cell( int p, Facet[] facets )
			{
				P = p;
				Facets = facets;
				Depths = new int[4];
			}
			
			public int P; // Number of edges in polygon
			public Facet[] Facets;
			public Vector3D Center;

			// Not necessary.
			public Mesh Mesh;

			/// <summary>
			/// Used to track recursing depth of reflections across various mirrors.
			/// </summary>
			public int[] Depths;

			public int LastReflection = -1;

			public bool IdealVerts
			{
				get
				{
					if( this.Verts.Count() == 0 )
						return false;
					return Tolerance.Equal( this.Verts.First().MagSquared(), 1 );
				}
			}

			public void AppendAllEdges( HashSet<Edge> edges )
			{
				foreach( Facet f in Facets )
					f.AppendAllEdges( edges );
			}

			/// <summary>
			/// In Ball model.
			/// </summary>
			public void CalcCenterFromFacets()
			{
				Vector3D center = new Vector3D();
				foreach( Sphere s in this.Facets.Select( f => f.Sphere ) )
				{
					if( s.IsPlane )
						continue;

					Vector3D sCenter = s.Center;
					double abs = s.Center.Abs();
					sCenter.Normalize();
					sCenter *= abs - s.Radius;
					center += sCenter;
				}

				center /= this.Facets.Length;
				this.Center = center;
			}

			public class Facet
			{
				public Facet( Vector3D[] verts ) { Verts = verts; }
				public Facet( Sphere sphere ) { Sphere = sphere; }

				/// <summary>
				/// The facet vertices.
				/// May live on the plane, in the ball model, etc. as needed.
				/// </summary>
				public Vector3D[] Verts;

				/// <summary>
				/// This is an alternate way to track facets, and *required* 
				/// for Lorentzian honeycombs, because all the vertices are hyperideal.
				/// It is expected that these are always in the Ball model.
				/// </summary>
				public Sphere Sphere { get; set; }
				private Vector3D CenterInBall 
				{ 
					get
					{
						if( Infinity.IsInfinite( Sphere.Radius ) )
							return new Vector3D();

						// Calcs based on orthogonal circles.
						// http://mathworld.wolfram.com/OrthogonalCircles.html
						double d = Math.Sqrt( 1 + Sphere.Radius * Sphere.Radius );
						Vector3D center = Sphere.Center;
						center.Normalize();
						center *= ( d - Sphere.Radius );
						return center;
					}
				}

				public void CalcSphereFromVerts( Geometry g )
				{
					switch( g )
					{
						case Geometry.Spherical:

							Sphere = new Sphere();
							if( Verts.Where( v => v == new Vector3D() ).Count() > 0 )	// XXX - not general, I'm so lazy.
							{
								Vector3D[] nonZero = Verts.Where( v => v != new Vector3D() ).ToArray();
								Sphere.Radius = double.PositiveInfinity;
								Sphere.Center = nonZero[0].Cross( nonZero[1] );
							}
							else
							{
								// The sphere intersects the unit-sphere at a unit-circle (orthogonal to the facet center direction).
								Vector3D direction = new Vector3D();
								foreach( Vector3D v in Verts )
									direction += v;
								direction /= Verts.Length;
								direction.Normalize();

								Vector3D p1 = Euclidean3D.ProjectOntoPlane( direction, new Vector3D(), Verts[0] );
								p1.Normalize();

								Circle3D c = new Circle3D( p1, Verts[0], -p1 );
								Sphere.Center = c.Center;
								Sphere.Radius = c.Radius;
							}

							break;

						case Geometry.Euclidean:

							Sphere = new Sphere();
							Sphere.Radius = double.PositiveInfinity;
							Vector3D v1 = Verts[0], v2 = Verts[1], v3 = Verts[2];
							Sphere.Center = ( v2 - v1 ).Cross( v3 - v1 );
							Sphere.Offset = Euclidean3D.ProjectOntoPlane( Sphere.Center, v1, new Vector3D() );
							break;

						case Geometry.Hyperbolic:
							Sphere = H3Models.Ball.OrthogonalSphereInterior( Verts[0], Verts[1], Verts[2] );
							break;
					}
				}
				
				public Facet Clone()
				{
					Facet newFacet = new Facet( Verts == null ? null : (Vector3D[])Verts.Clone() );
					newFacet.Sphere = Sphere == null ? null : Sphere.Clone();
					return newFacet;
				}

				public void Reflect( Sphere sphere )
				{
					if( Verts != null )
					{
						for( int i=0; i<Verts.Length; i++ )
							Verts[i] = sphere.ReflectPoint( Verts[i] );
					}

					if( Sphere != null )
						Sphere.Reflect( sphere );
				}

				public Vector3D ID
				{
					get
					{
						Vector3D result = new Vector3D();
						if( Verts != null )
						{
							foreach( Vector3D v in Verts )
								result += v;
							return result;
						}

						if( Sphere != null )
							return Sphere.ID;

						throw new System.ArgumentException();
					}
				}

				public void AppendAllEdges( HashSet<Edge> edges )
				{
					// We can only do this if we have vertices.
					if( Verts == null )
						return;

					for( int i=0; i<Verts.Length; i++ )
					{
						int idx1 = i;
						int idx2 = i == Verts.Length - 1 ? 0 : i + 1;
						edges.Add( new Edge( Verts[idx1], Verts[idx2] ) );
					}
				}
			}

			public class Edge
			{
				public Edge( Vector3D v1, Vector3D v2, bool order = true )
				{
					// Keep things "ordered", so we can easily compare edges.
					if( order )
					{
						Vector3D[] orderedVerts = new Vector3D[] { v1, v2 };
						orderedVerts = orderedVerts.OrderBy( v => v, new Vector3DComparer() ).ToArray();
						Start = orderedVerts[0];
						End = orderedVerts[1];
					}
					else
					{
						Start = v1;
						End = v2;
					}

					Depths = new int[4];
				}

				public Vector3D Start;
				public Vector3D End;

				// The reason we use a vector here is so the components 
				// can be interpreted in different color schemes (HLS, RGB, etc.)
				public Vector3D Color = new Vector3D( 1, 1, 1 );

				/// <summary>
				/// Used to track recursing depth of reflections across various mirrors.
				/// </summary>
				public int[] Depths;

				public Edge Clone()
				{
					return (Edge)MemberwiseClone();
				}

				public Vector3D ID
				{
					get
					{
						return Start + End;
					}
				}

				public void CopyDepthsFrom( Edge e )
				{
					Depths = (int[])e.Depths.Clone();
				}
			}

			public class EdgeEqualityComparer : IEqualityComparer<Edge>
			{
				public bool Equals( Edge e1, Edge e2 )
				{
					return
						e1.Start.Compare( e2.Start, m_tolerance ) &&
						e1.End.Compare( e2.End, m_tolerance );
				}

				public int GetHashCode( Edge e )
				{
					return e.Start.GetHashCode() ^ e.End.GetHashCode();
				}

				private double m_tolerance = 0.0001;
			}

			public bool HasVerts
			{
				get
				{
					foreach( Facet f in Facets )
						if( f.Verts == null )
							return false;
					return true;
				}
			}

			public IEnumerable<Vector3D> Verts
			{
				get
				{
					foreach( Facet facet in Facets )
					{
						if( facet.Verts != null )
							foreach( Vector3D v in facet.Verts )
								yield return v;
					}
				}
			}

			/// <summary>
			/// Additional points (could be whatever) that we want reflected around with this cell.
			/// </summary>
			public Vector3D[] AuxPoints { get; set; }

			public Vector3D ID
			{
				get 
				{
					Vector3D result = new Vector3D();
					//foreach( Vector3D v in Verts )
						//result += Sterographic.PlaneToSphereSafe( v );	// XXX - what about when not working in plane.

					if( HasVerts )
					{
						foreach( Vector3D v in Verts )
							result += v;
					}
					else
					{
						// Intentionally just use the center.
					}
					result += Center;
					return result;
				}
			}

			public Cell Clone()
			{
				List<Facet> newFacets = new List<Facet>();
				foreach( Facet facet in Facets )
					newFacets.Add( facet.Clone() );

				Cell clone = new Cell( P, newFacets.ToArray() );
				clone.Center = Center;
				clone.Depths = (int[])Depths.Clone();
				clone.LastReflection = LastReflection;

				if( AuxPoints != null )
					clone.AuxPoints = (Vector3D[])AuxPoints.Clone();

				if( Mesh != null )
					clone.Mesh = Mesh.Clone();

				return clone;
			}

			/// <summary>
			/// Scales cell to circumsphere.
			/// NOTE: We should already be on a sphere, not the plane.
			/// </summary>
			public void ScaleToCircumSphere( double r )
			{
				foreach( Facet facet in Facets )
					for( int i=0; i<P; i++ )
					{
						Vector3D axis = facet.Verts[i];
						axis.Normalize();
						facet.Verts[i] = axis * r;
					}
			}

			/// <summary>
			/// Apply a Mobius transformation to us (meaning of Mobius transform is on boundary of UHS model).
			/// </summary>
			public void ApplyMobius( Mobius m )
			{
				foreach( Facet facet in Facets )
					for( int i = 0; i < P; i++ )
						facet.Verts[i] = H3Models.Ball.ApplyMobius( m, facet.Verts[i] );
				Center = H3Models.Ball.ApplyMobius( m, Center );
			}

			public void Reflect( Sphere sphere )
			{
				foreach( Facet facet in Facets )
					facet.Reflect( sphere );
				Center = sphere.ReflectPoint( Center );

				if( AuxPoints != null )
				{
					for( int i = 0; i < AuxPoints.Length; i++ )
						AuxPoints[i] = sphere.ReflectPoint( AuxPoints[i] );
				}
				
				if( this.Mesh != null )
				{
					for( int i=0; i<Mesh.Triangles.Count; i++ )
					{
						Mesh.Triangle tri = Mesh.Triangles[i];
						tri.a = sphere.ReflectPoint( tri.a );
						tri.b = sphere.ReflectPoint( tri.b );
						tri.c = sphere.ReflectPoint( tri.c );
						Mesh.Triangles[i] = tri;
					}
				}
			}
		}

		public enum Output
		{
			POVRay
		}

		public class Settings
		{
			public double Scale = 50;	// 10cm ball by default.

			public bool Halfspace = false;
			public int MaxCells = 150000;

			// Ball Model
			public double Ball_MaxLength = 3;
			//public double Ball_MinLength = 0.075;
			public double Ball_MinLength = 0.15;
			//public double Ball_MinLength = 0.05;
			//private static double Ball_MinLength = 0.45;	// lamp
			public double Ball_Cutoff = 0.95;

			// UHS
			//public double UHS_MinEdgeLength = .09;
			//public double UHS_MaxBounds = 6.5;
			public double UHS_MinEdgeLength = 0.03;
			public double UHS_MaxBounds = 2;
			public double UHS_Horocycle = 0.25;

			// Bananas
			public bool ThinEdges = false;
			public double AngularThickness = 0.06;	// an angle (the slope of the banana)
			//public double AngularThickness = 0.04;
			//public double AngularThickness = 0.25;

			// Position and associated Mobius to apply
			public Polytope.Projection Position = Polytope.Projection.CellCentered;
			public Mobius Mobius = Mobius.Identity();

			public Output Output = Output.POVRay;
		}

		public static Settings m_settings = new Settings();

		public static string m_baseDir = @"./";

		private static int[] m_levelsToKeep = new int[] { 1, 2, 3, 4, 5 };
		private static double[] m_rangeToKeep = new double[] { 0.9, 0.96 };
		private static bool CheckRange( Vector3D v )
		{
			double abs = v.Abs();
			return
				abs > m_rangeToKeep[0] &&
				abs < m_rangeToKeep[1];
		}

		public static void SaveToFile( string honeycombString, Cell.Edge[] edges, bool finite, bool append = false )
		{
			SaveToFile( honeycombString, edges.ToDictionary( e => e, e => 1 ), finite, append );
		}

		public static void SaveToFile( string honeycombString, Dictionary<Cell.Edge, int> edges, bool finite, bool append = false )
		{
		PovRay.WriteH3Edges( new PovRay.Parameters() 
					{ 
						AngularThickness = m_settings.AngularThickness,
						Halfspace = m_settings.Halfspace,
						//Halfspace = true,
						ThinEdges = m_settings.ThinEdges,
					},
					edges.Keys.ToArray(), m_baseDir + honeycombString + ".pov", append ); 
		}

		public static void AppendFacets( string honeycombString, H3.Cell[] cells )
		{
			PovRay.AppendFacets( cells, m_baseDir + honeycombString + ".pov" );
		}

		private static void CheckAndAdd( Dictionary<Vector3D, int> vertexCounts, Vector3D v )
		{
			int count;
			if( vertexCounts.TryGetValue( v, out count ) )
				count++;
			else
				count = 1;

			vertexCounts[v] = count;
		}


		public static double SizeFuncConst( Vector3D v )
		{
			return H3Models.SizeFuncConst( v, m_settings.Scale );
		}
	}
}
