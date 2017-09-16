namespace R3.Geometry
{
	using Math = System.Math;
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Linq;

	public enum SegmentType
	{
		Line,
		Arc
	}

	public class Segment : ITransformable
	{
		public SegmentType Type { get; set; }
		public Vector3D P1 { get; set; }
		public Vector3D P2 { get; set; }

		// These only apply to arc segments.
		public Vector3D Center { get; set; }
		public bool Clockwise { get; set; }

		public Segment Clone()
		{
			Segment newSeg = new Segment();
			newSeg.Type = Type;
			newSeg.P1 = P1;
			newSeg.P2 = P2;
			newSeg.Center = Center;
			newSeg.Clockwise = Clockwise;
			return newSeg;
		}

		public static Segment Line( Vector3D start, Vector3D end )
		{
			Segment newSeg = new Segment();
			newSeg.Type = SegmentType.Line;
			newSeg.P1 = start;
			newSeg.P2 = end;
			return newSeg;
		}

		public static Segment Arc( Vector3D start, Vector3D end, Vector3D center, bool clockwise )
		{
			Segment newSeg = new Segment();
			newSeg.Type = SegmentType.Arc;
			newSeg.P1 = start;
			newSeg.P2 = end;
			newSeg.Center = center;
			newSeg.Clockwise = clockwise;
			return newSeg;
		}

		public static Segment Arc( Vector3D start, Vector3D mid, Vector3D end )
		{
			Segment newSeg = new Segment();
			newSeg.Type = SegmentType.Arc;
			newSeg.P1 = start;
			newSeg.P2 = end;

			Circle c = new Circle();
			c.From3Points( start, mid, end );
			newSeg.Center = c.Center;

			// Obtain vectors from center point of circle (as if at the origin)
			Vector3D startOrigin = start - c.Center;
			Vector3D midOrigin = mid - c.Center;
			Vector3D endOrigin = end - c.Center;

			// Calculate the normal vector and angle to traverse.
			// ZZZ - worry about failure of cross product here.
			Vector3D normalVector = startOrigin.Cross( endOrigin );
			newSeg.Clockwise = normalVector.Z < 0;
			double angleToTraverse = startOrigin.AngleTo( endOrigin );

			// The normal vector might need to be reversed and the angleToTraverse adjusted.
			// This happens depending on the location of the midpoint relative to the start and end points.
			double compareAngle = startOrigin.AngleTo( midOrigin ) + midOrigin.AngleTo( endOrigin );
			bool reverse = !Tolerance.Equal( angleToTraverse, compareAngle );
			if( reverse )
				newSeg.Clockwise = !newSeg.Clockwise;

			return newSeg;
		}

		public double Radius
		{
			get
			{
				Debug.Assert( SegmentType.Arc == Type );
				return ( ( P1 - Center ).Abs() );
			}
		}

		public double Angle
		{
			get
			{
				if( SegmentType.Arc != Type )
				{
					Debug.Assert( false );
					return 0;
				}

				Vector3D v1 = P1 - Center;
				Vector3D v2 = P2 - Center;
				return Clockwise ?
					Euclidean2D.AngleToClock( v1, v2 ) :
					Euclidean2D.AngleToCounterClock( v1, v2 );
			}
		}

		public Circle Circle
		{
			get
			{
				Debug.Assert( SegmentType.Arc == Type );

				// Avoiding allocations of new circles,
				// (Memory profiling showed this was responsible
				// for many allocations.)
				if( m_circle != null )
				{
					if( m_circle.Center == this.Center &&
						m_circle.Radius == this.Radius )
						return m_circle;
				}

				m_circle = new Circle();
				m_circle.Center = Center;
				m_circle.Radius = Radius;
				return m_circle;
			}
		}
		private Circle m_circle;

		public double Length
		{
			get
			{
				if( SegmentType.Arc == Type )
				{
					return Radius * Angle;
				}
				else
				{
					return ( P2 - P1 ).Abs();
				}
			}
		}

		public Vector3D Midpoint
		{
			get
			{
				if( SegmentType.Arc == Type )
				{
					double a = Angle / 2;
					Vector3D ret = P1 - Center;
					ret.RotateXY( Clockwise ? -a : a );
					ret += Center;
					return ret;
				}
				else
				{
					return ( P1 + P2 ) / 2;
				}
			}
		}

		public void Reverse() 
		{
			SwapPoints();
			if( SegmentType.Arc == Type )
				Clockwise = !Clockwise;
		}

		/// <summary>
		/// Return the vertices from subdividing ourselves.
		/// </summary>
		public Vector3D[] Subdivide ( int numSegments )
		{
			List<Vector3D> ret = new List<Vector3D>();
			if( numSegments < 1 )
			{
				Debug.Assert( false );
				return ret.ToArray();
			}

			if( Type == SegmentType.Arc )
			{
				Vector3D v = P1 - Center;
				double angle = this.Angle / numSegments;
				for( int i=0; i<numSegments; i++ )
				{
					ret.Add( Center + v );
					v.RotateXY( Clockwise ? -angle : angle );
				}
			}
			else
			{
				Vector3D v = P2 - P1;
				v.Normalize();
				for( int i=0; i<numSegments; i++ )
					ret.Add( P1 + v * i * Length / numSegments );
			}

			// Add in the last point and return.
			ret.Add( P2 );
			return ret.ToArray();
		}

		public void SwapPoints()
		{
			Vector3D t = P1;
			P1 = P2;
			P2 = t;
		}

		public void Reflect( Segment s )
		{
			// NOTES:
			// Arcs can go to lines, and lines to arcs.
			// Rotations may reverse arc directions as well.
			// Arc centers can't be transformed directly.

			// NOTE: We must calc this before altering the endpoints.
			Vector3D mid = Midpoint;
			if( Infinity.IsInfinite( mid ) )
				mid = Infinity.IsInfinite( s.P1 ) ? s.P2 * Infinity.FiniteScale : s.P1 * Infinity.FiniteScale;

			P1 = s.ReflectPoint( P1 );
			P2 = s.ReflectPoint( P2 );
			mid = s.ReflectPoint( mid );

			// Can we make a circle out of the reflected points?
			Circle temp = new Circle();
			if( !Infinity.IsInfinite( P1 ) && !Infinity.IsInfinite( P2 ) && !Infinity.IsInfinite( mid ) &&
				temp.From3Points( P1, mid, P2 ) )
			{
				Type = SegmentType.Arc;
				Center = temp.Center;

				// Work out the orientation of the arc.
				Vector3D t1 = P1 - Center;
				Vector3D t2 = mid - Center;
				Vector3D t3 = P2 - Center;
				double a1 = Euclidean2D.AngleToCounterClock( t2, t1 );
				double a2 = Euclidean2D.AngleToCounterClock( t3, t1 );
				Clockwise = a2 > a1;
			}
			else
			{
				// The circle construction fails if the points
				// are colinear (if the arc has been transformed into a line).
				Type = SegmentType.Line;

				// XXX - need to do something about this.
				// Turn into 2 segments?
				//if( isInfinite( mid ) )
				// Actually the check should just be whether mid is between p1 and p2.
			}
		}

		public void Transform( Mobius m )
		{
			TransformInternal( m );
		}

		public void Transform( Isometry i )
		{
			TransformInternal( i );
		}

		/// <summary>
		/// Apply a transform to us.
		/// </summary>
		private void TransformInternal<T>( T transform ) where T : ITransform
		{
			// NOTES:
			// Arcs can go to lines, and lines to arcs.
			// Rotations may reverse arc directions as well.
			// Arc centers can't be transformed directly.

			// NOTE: We must calc this before altering the endpoints.
			Vector3D mid = Midpoint;
			if( Infinity.IsInfinite( mid ) )
				mid = Infinity.IsInfinite( P1 ) ? P2 * Infinity.FiniteScale : P1 * Infinity.FiniteScale;

			P1 = transform.Apply( P1 );
			P2 = transform.Apply( P2 );
			mid = transform.Apply( mid );

			// Can we make a circle out of the transformed points?
			Circle temp = new Circle();
			if( !Infinity.IsInfinite( P1 ) && !Infinity.IsInfinite( P2 ) && !Infinity.IsInfinite( mid ) &&
				temp.From3Points( P1, mid, P2 ) )
			{
				Type = SegmentType.Arc;
				Center = temp.Center;

				// Work out the orientation of the arc.
				Vector3D t1 = P1 - Center;
				Vector3D t2 = mid - Center;
				Vector3D t3 = P2 - Center;
				double a1 = Euclidean2D.AngleToCounterClock( t2, t1 );
				double a2 = Euclidean2D.AngleToCounterClock( t3, t1 );
				Clockwise = a2 > a1;
			}
			else
			{
				// The circle construction fails if the points
				// are colinear (if the arc has been transformed into a line).
				Type = SegmentType.Line;

				// XXX - need to do something about this.
				// Turn into 2 segments?
				//if( isInfinite( mid ) )
				// Actually the check should just be whether mid is between p1 and p2.
			}
		}


		public Vector3D ReflectPoint( Vector3D input ) 
		{
			if( SegmentType.Arc == Type )
			{
				Circle c = this.Circle;
				return c.ReflectPoint( input );
			}
			else
			{
				return Euclidean2D.ReflectPointInLine( input, P1, P2 );
			}
		}
	}

	public class Polygon : ITransformable
	{
		public Polygon()
		{
			Segments = new List<Segment>();
			Center = new Vector3D();
		}

		public Vector3D Center { get; set; }
		public List<Segment> Segments { get; set; }

		public void Clear()
		{
			Segments.Clear();
		}

		public Polygon Clone()
		{
			Polygon newPoly = new Polygon();
			//newPoly.Segments = new List<Segment>( Segments );
			foreach( Segment s in Segments )
				newPoly.Segments.Add( s.Clone() );
			newPoly.Center = Center;
			return newPoly;
		}

		public void CreateRegular( int numSides, int q )
		{
			int p = numSides;

			Segments.Clear();
			List<Vector3D> points = new List<Vector3D>();

			Geometry g = Geometry2D.GetGeometry( p, q );
			double circumRadius = Geometry2D.GetNormalizedCircumRadius( p, q );

			double angle = 0;
			for( int i=0; i<p; i++ )
			{
				Vector3D point = new Vector3D();
				point.X = ( circumRadius * Math.Cos( angle ) );
				point.Y = ( circumRadius * Math.Sin( angle ) );
				points.Add( point );
				angle += Utils.DegreesToRadians( 360.0 / p );
			}

			// Turn this into segments.
			for( int i=0; i<points.Count; i++ )
			{
				int idx1 = i;
				int idx2 = i == points.Count - 1 ? 0 : i+1;
				Segment newSegment = new Segment();
				newSegment.P1 = points[idx1];
				newSegment.P2 = points[idx2];

				if( g != Geometry.Euclidean )
				{
					newSegment.Type = SegmentType.Arc;

					if( 2 == p )
					{
						// Our magic formula below breaks down for digons.
						double factor = Math.Tan( Math.PI / 6 );
						newSegment.Center = newSegment.P1.X > 0 ?
							new Vector3D( 0,-circumRadius,0 ) * factor :
							new Vector3D( 0, circumRadius,0 ) * factor;
					}
					else
					{
						// Our segments are arcs in Non-Euclidean geometries.
						// Magically, the same formula turned out to work for both.
						// (Maybe this is because the Poincare Disc model of the
						// hyperbolic plane is stereographic projection as well).

						double piq = q == -1 ? 0 : Math.PI / q;	// Handle q infinite.
						double t1 = Math.PI / p;
						double t2 = Math.PI / 2 - piq - t1;
						double factor = ( Math.Tan( t1 ) / Math.Tan( t2 ) + 1 ) / 2;
						newSegment.Center = ( newSegment.P1 + newSegment.P2 ) * factor;
					}
				
					newSegment.Clockwise = Geometry.Spherical == g ? false : true;
				}

				// XXX - Make this configurable?
				// This is the color of cell boundary lines.
				//newSegment.m_color = CColor( 1, 1, 0, 1 );
				Segments.Add( newSegment );
			}
		}

		public int NumSides
		{
			get { return Segments.Count; }
		}

		public Vector3D[] EdgePoints
		{
			get
			{
				List<Vector3D> points = new List<Vector3D>();
				double arcResolution = Utils.DegreesToRadians( 4.5 );

				for( int i = 0; i < NumSides; i++ )
				{
					Segment s = Segments[i];

					// First point.
					// ZZZ - getting lazy
					//Debug.Assert( ! (isInfinite( s.m_p1 ) && isInfinite( s.m_p2 )) );
					Vector3D p1 = Infinity.IsInfinite( s.P1 ) ? s.P2 * Infinity.FiniteScale : s.P1;
					points.Add( p1 );

					// For arcs, add in a bunch of extra points.
					if( SegmentType.Arc == s.Type )
					{
						double maxAngle = s.Angle;
						Vector3D vs = s.P1 - s.Center;
						int numSegments = (int)(maxAngle / (arcResolution));
						if( numSegments < 10 )	// ZZZ - arbitrary.
							numSegments = 10;
						double angle = maxAngle / numSegments;
						for( int j = 1; j < numSegments; j++ )
						{
							vs.RotateXY( s.Clockwise ? -angle : angle );
							points.Add( vs + s.Center );
						}
					}

					// Last point.
					Vector3D p2 = Infinity.IsInfinite( s.P2 ) ? s.P1 * Infinity.FiniteScale : s.P2;
					points.Add( p2 );
				}

				return points.ToArray();
			}
		}

		public CircleNE CircumCircle
		{
			get
			{
				CircleNE result = new CircleNE();
				if( Segments.Count > 2 )
					result.From3Points( Segments[0].P1, Segments[1].P1, Segments[2].P1 );
				result.CenterNE = this.Center;
				return result;
			}
		}

		public void Reverse()
		{
			// Reverse all our segments and swap the order of them.
			foreach( Segment s in Segments )
				s.Reverse();

			Segments.Reverse();
		}

		/// <summary>
		/// Apply a Mobius transform to us.
		/// </summary>
		public void Transform( Mobius m )
		{
			foreach( Segment s in this.Segments )
				s.Transform( m );
			Center = m.Apply( Center );
		}

		/// <summary>
		/// Apply an isometry to us.
		/// </summary>
		public void Transform( Isometry isometry )
		{
			foreach( Segment s in this.Segments )
				s.Transform( isometry );
			Center = isometry.Apply( Center );
		}
	}
}