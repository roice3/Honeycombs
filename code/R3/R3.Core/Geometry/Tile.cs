namespace R3.Geometry
{
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Numerics;

	public class Tile
	{
		public Tile()
		{ 
		}

		public Tile( Polygon boundary, Polygon drawn, Geometry geometry )
			: this()
		{
			Boundary = boundary;
			Drawn = drawn;
			Geometry = geometry;

			// Make the vertex circle.
			VertexCircle = boundary.CircumCircle;

			// ZZZ - we shouldn't do this here (I did it for the slicing study page).
			//VertexCircle.Radius = 1.0;

			//
			// Below are experimentations with different vertex circle sizes.
			//

			//VertexCircle.Radius *= (1+1.0/9);

			// cuts adjacent cells at midpoint
			// Math.Sqrt(63)/6 for {3,6}
			// (1 + 1.0/5) for {3,7} 
			// (1 + 1.0/9) for {3,8}
			// (1 + 1.0/20) for {3,9}

			// cuts at 1/3rd
			// 2/Math.Sqrt(3) for {3,6}
		}

		public Polygon Boundary { get; set; }
		public Polygon Drawn { get; set; }
		public CircleNE VertexCircle { get; set; }
		public Geometry Geometry { get; set; }

		/// <summary>
		/// This will trim back the tile using an equidistant curve.
		/// It assumes the tile is at the origin.
		/// </summary>
		internal static void ShrinkTile( ref Tile tile, double shrinkFactor )
		{
			// This code is not correct in non-Euclidean cases!
			// But it works reasonable well for small shrink factors.
			// For example, you can easily use this function to grow a hyperbolic tile beyond the disk.
			Mobius m = new Mobius();
			m.Hyperbolic( tile.Geometry, new Vector3D(), shrinkFactor );
			tile.Drawn.Transform( m );
			return;

			/*
			// ZZZ
			// Wow, all the work I did below was subsumed by 4 code lines above!
			// I can't bring myself to delete it yet.

			switch( tile.Geometry )
			{
				case Geometry.Spherical:
				{
					List<Tile> clipped = new List<Tile>();
					clipped.Add( tile );

					Polygon original = tile.Drawn.Clone();
					foreach( Segment seg in original.Segments )
					{
						Debug.Assert( seg.Type == SegmentType.Arc );

						if( true )
						{
							// Unproject to sphere.
							Vector3D p1 = Spherical2D.PlaneToSphere( seg.P1 );
							Vector3D p2 = Spherical2D.PlaneToSphere( seg.P2 );

							// Get the poles of the GC, and project them to the plane.
							Vector3D pole1, pole2;
							Spherical2D.GreatCirclePole( p1, p2, out pole1, out pole2 );
							pole1 = Spherical2D.SphereToPlane( pole1 );
							pole2 = Spherical2D.SphereToPlane( pole2 );

							// Go hyperbolic, dude.
							double scale = 1.065;	// ZZZ - needs to be configurable.
							Complex fixedPlus = pole1;
							Mobius hyperbolic = new Mobius();
							hyperbolic.Hyperbolic( tile.Geometry, fixedPlus, scale );
							Vector3D newP1 = hyperbolic.Apply( seg.P1 );
							Vector3D newMid = hyperbolic.Apply( seg.Midpoint );
							Vector3D newP2 = hyperbolic.Apply( seg.P2 );

							Circle trimmingCircle = new Circle();
							trimmingCircle.From3Points( newP1, newMid, newP2 );

							Slicer.Clip( ref clipped, trimmingCircle, true );
						}
						else
						{
							// I think this block has logic flaws, but strangely it seems to work,
							// so I'm leaving it in commented out for posterity.

							Vector3D p1 = seg.P1;
							Vector3D mid = seg.Midpoint;
							Vector3D p2 = seg.P2;

							//double offset = .1;
							double factor = .9;
							double f1 = Spherical2D.s2eNorm( (Spherical2D.e2sNorm( p1.Abs() ) * factor) );
							double f2 = Spherical2D.s2eNorm( (Spherical2D.e2sNorm( mid.Abs() ) * factor) );
							double f3 = Spherical2D.s2eNorm( (Spherical2D.e2sNorm( p2.Abs() ) * factor) );
							p1.Normalize();
							mid.Normalize();
							p2.Normalize();
							p1 *= f1;
							mid *= f2;
							p2 *= f3;

							Circle trimmingCircle = new Circle();
							trimmingCircle.From3Points( p1, mid, p2 );

							Slicer.Clip( ref clipped, trimmingCircle, true );
						}
					}

					Debug.Assert( clipped.Count == 1 );
					tile = clipped[0];
					return;
				}
				case Geometry.Euclidean:
				{
					double scale = .95;

					Mobius hyperbolic = new Mobius();
					hyperbolic.Hyperbolic( tile.Geometry, new Vector3D(), scale );

					tile.Drawn.Transform( hyperbolic );

					return;
				}
				case Geometry.Hyperbolic:
				{
					List<Tile> clipped = new List<Tile>();
					clipped.Add( tile );

					Circle infinity = new Circle();
					infinity.Radius = 1.0;

					Polygon original = tile.Drawn.Clone();
					foreach( Segment seg in original.Segments )
					{
						Debug.Assert( seg.Type == SegmentType.Arc );
						Circle segCircle = seg.GetCircle();

						// Get the intersection points with the disk at infinity.
						Vector3D p1, p2;
						int count = Euclidean2D.IntersectionCircleCircle( infinity, segCircle, out p1, out p2 );
						Debug.Assert( count == 2 );

						Vector3D mid = seg.Midpoint;
						//mid *= 0.75;	// ZZZ - needs to be configurable.
							
						double offset = .03;
						double f1 = DonHatch.h2eNorm( DonHatch.e2hNorm( mid.Abs() ) - offset );
						mid.Normalize();
						mid *= f1;
							
						Circle trimmingCircle = new Circle();
						trimmingCircle.From3Points( p1, mid, p2 );

						Slicer.Clip( ref clipped, trimmingCircle, false );
					}

					Debug.Assert( clipped.Count == 1 );
					tile = clipped[0];
					return;
				}
			}
			*/
		}
	}
}
