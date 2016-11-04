namespace R3.Geometry
{
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;
	using Math = System.Math;
	using Edge = H3.Cell.Edge;

	public static class Recurse
	{
		public class Settings
		{
			public Geometry G = Geometry.Hyperbolic;
			public int MaxEdges = 2500000;
			public double Threshold = DefaultThresh();
			public bool Ball = true;
		}

		/// <summary>
		/// Attempts to calculate approx 1.3M edges when the threshold is a distance from origin in the ball model.
		/// This works for honeycombs with finite cells.
		/// </summary>
		public static Edge[] CalcEdgesSmart( Sphere[] simplex, Edge[] edges )
		{
			Settings s = new Settings();

			// The number of cells increase exponentially with hyperbolic distance,
			// so linear on a log scale.
			// We'll do to test runs to get the line, then run at the extrapolated value.
			double hDist = 5;
			s.Threshold = DonHatch.h2eNorm( hDist );
			Edge[] result = CalcEdges( simplex, edges, s );
			int count1 = result.Length;

			hDist = 5.5;
			s.Threshold = DonHatch.h2eNorm( hDist );
			result = CalcEdges( simplex, edges, s );
			int count2 = result.Length;

			double slope = ( Math.Log( count2 ) - Math.Log( count1 ) ) / 0.5;

			// Why 1.3M?  We'll get 650k after we half this.
			double desiredCount = 1.0e4;		// for testing
			//double desiredCount = 1.3e6;
			//double desiredCount = 0.4e6;	// Mid-range
			double logDesiredCount = Math.Log( desiredCount );
			hDist = 5.5 + ( logDesiredCount - Math.Log( count2 ) ) / slope;

			s.Threshold = DonHatch.h2eNorm( hDist );
			return CalcEdges( simplex, edges, s );
		}

		/// <summary>
		/// Attempts to calculate approx 1.3M edges when the threshold is a minimum edge length.
		/// This is required for honeycombs with ideal or ultra-ideal cells
		/// </summary>
		public static Edge[] CalcEdgesSmart2( Sphere[] simplex, Edge[] edges )
		{
			Settings s = new Settings();

			// I found that log(1/thresh)/log(count) was relatively constant,
			// so we'll extrapolate that to get close to the right number of edges.
			double OneOverThresh = 60;
			s.Threshold = 1 / OneOverThresh;
			Edge[] result = CalcEdges( simplex, edges, s );
			int count1 = result.Length;

			OneOverThresh = 80;
			s.Threshold = 1 / OneOverThresh;
			result = CalcEdges( simplex, edges, s );
			int count2 = result.Length;

			double slope = ( Math.Log( count2 ) - Math.Log( count1 ) ) / ( Math.Log( 80 ) - Math.Log( 60 ) );

			// Why 1.3M?  We'll get 650k after we half this.
			double desiredCount = 2e6;
			//double desiredCount = 3e4;	// For testing
			double logDesiredCount = Math.Log( desiredCount );
			double temp = Math.Log( 80 ) + ( logDesiredCount - Math.Log( count2 ) ) / slope;

			s.Threshold = 1 / Math.Exp( temp );
			return CalcEdges( simplex, edges, s );
		}

		public static Edge[] CalcEdges( Sphere[] simplex, Edge[] edges )
		{
			Settings settings = new Settings();
			return CalcEdges( simplex, edges, settings );
		}

		public static Edge[] CalcEdges( Sphere[] simplex, Edge[] edges, Settings settings )
		{
			HashSet<Edge> completedEdges = new HashSet<Edge>( edges, new H3.Cell.EdgeEqualityComparer() );
			ReflectEdgesRecursive( simplex, completedEdges.ToArray(), settings, completedEdges );
			return completedEdges.ToArray();
		}

		private static void ReflectEdgesRecursive( Sphere[] simplex, Edge[] edges, Settings settings,
			HashSet<Edge> completedEdges )
		{
			if( 0 == edges.Length )
				return;

			HashSet<Edge> newEdges = new HashSet<Edge>( new H3.Cell.EdgeEqualityComparer() );

			foreach( Edge edge in edges )
			//foreach( Sphere mirror in simplex )
			for( int m=0; m<simplex.Length; m++ )
			{
				Sphere mirror = simplex[m];

				if( completedEdges.Count > settings.MaxEdges )
					throw new System.Exception( "Maxing out edges - will result in uneven filling." );

				Vector3D r1 = mirror.ReflectPoint( edge.Start );
				Vector3D r2 = mirror.ReflectPoint( edge.End );

				Edge newEdge = new Edge( r1, r2 );
				newEdge.CopyDepthsFrom( edge );
				if( !EdgeOk( newEdge, settings ) )
					continue;

				// This tracks reflections across the cell facets.
				newEdge.Depths[m]++;

				// Edge color.
				// Make the threshold length black, or the background color.
				double percentWhite = ( r1.Dist( r2 ) - settings.Threshold ) / 0.015;
				if( percentWhite < 0 )
					percentWhite = 0;
				if( percentWhite > 1 )
					percentWhite = 1;
				//newEdge.Color = new Vector3D( percentWhite, percentWhite, percentWhite );
				newEdge.Color = m_background;
				newEdge.Color.Z = 0.1 + 0.9 * percentWhite;

				if( completedEdges.Add( newEdge ) )
				{
					// Haven't seen this edge yet, so 
					// we'll need to recurse on it.
					newEdges.Add( newEdge );
				}
			}

			ReflectEdgesRecursive( simplex, newEdges.ToArray(), settings, completedEdges );
		}

		public static Vector3D m_background;

		public static H3.Cell[] CalcCells( Sphere[] mirrors, H3.Cell[] cells )
		{
			Settings settings = new Settings();
			return CalcCells( mirrors, cells, settings );
		}

		public static H3.Cell[] CalcCells( Sphere[] mirrors, H3.Cell[] cells, Settings settings )
		{
			HashSet<Vector3D> completedCellIds = new HashSet<Vector3D>( cells.Select( c => c.ID ).ToArray() );
			List<H3.Cell> completedCells = new List<H3.Cell>( cells );
			ReflectCellsRecursive( mirrors, cells, settings, completedCells, completedCellIds );
			return completedCells.ToArray();
		}

		private static void ReflectCellsRecursive( Sphere[] simplex, H3.Cell[] cells, Settings settings,
			List<H3.Cell> completedCells, HashSet<Vector3D> completedCellIds )
		{
			if( 0 == cells.Length )
				return;

			List<H3.Cell> newCells = new List<H3.Cell>();

			foreach( H3.Cell cell in cells )
			//foreach( Sphere mirror in simplex )
			for( int m=0; m<simplex.Length; m++ )
			{
				Sphere mirror = simplex[m];
				//if( m == 2 )
				//	continue;

				//if( completedCellIds.Count > 1000 )
				//	return;

				//if( completedCellIds.Count > settings.MaxEdges/5 )
				if( completedCellIds.Count > settings.MaxEdges / 20 )
					throw new System.Exception( "Maxing out cells - will result in uneven filling." );

				H3.Cell newCell = cell.Clone();	
				newCell.Reflect( mirror );
				if( !CellOk( newCell, settings ) )
					continue;

				// This tracks reflections across the cell facets.
				newCell.Depths[m]++;

				if( completedCellIds.Add( newCell.ID ) )
				{
					// Haven't seen this cell yet, so 
					// we'll need to recurse on it.
					newCells.Add( newCell );
					completedCells.Add( newCell );
				}
			}

			ReflectCellsRecursive( simplex, newCells.ToArray(), settings, completedCells, completedCellIds );
		}

/////////////////////////////////////////////////////////
// Hacking around.  I should remove this or make CalcCells() configurable enough to deal with more cases.
		public static H3.Cell[] CalcCells2( Sphere[] mirrors, H3.Cell[] cells )
		{
			Settings settings = new Settings();
			return CalcCells2( mirrors, cells, settings );
		}

		public static H3.Cell[] CalcCells2( Sphere[] mirrors, H3.Cell[] cells, Settings settings )
		{
			HashSet<Vector3D> completedCellIds = new HashSet<Vector3D>( cells.Select( c => c.ID ).ToArray() );
			List<H3.Cell> completedCells = new List<H3.Cell>( cells );
			ReflectCellsRecursive2( mirrors, cells, settings, completedCells, completedCellIds );
			return completedCells.ToArray();
		}

		private static void ReflectCellsRecursive2( Sphere[] simplex, H3.Cell[] cells, Settings settings,
			List<H3.Cell> completedCells, HashSet<Vector3D> completedCellIds )
		{
			if( 0 == cells.Length )
				return;

			List<H3.Cell> newCells = new List<H3.Cell>();

			foreach( H3.Cell cell in cells )
			//foreach( Sphere mirror in simplex )
			for( int m = 0; m < simplex.Length; m++ )
			{
				Sphere mirror = simplex[m];
				if( completedCellIds.Count > 250000 )
					return;

				H3.Cell newCell = cell.Clone();
				newCell.Reflect( mirror );
				//if( !CellOk( newCell, settings ) )
				bool cellOk = true;
				foreach( H3.Cell.Facet f in cell.Facets )
					if( f.Sphere.Radius < 0.002 )
						cellOk = false;
				if( !cellOk )
					continue;

				// This tracks reflections across the cell facets.
				newCell.Depths[m]++;

				if( completedCellIds.Add( newCell.ID ) )
				{
					// Haven't seen this cell yet, so 
					// we'll need to recurse on it.
					newCells.Add( newCell );
					completedCells.Add( newCell );
				}
			}

			ReflectCellsRecursive2( simplex, newCells.ToArray(), settings, completedCells, completedCellIds );
		}
/////////////////////////////////////////////////////////

		private static double DefaultThresh()
		{
			//double thresh = 0.995;	// 435,353
			//double thresh = 0.994;	// 534
			//double thresh = 0.996;	// 336
			//double thresh = 0.997;	// 535, 344
			//double thresh = 0.9975;	// 436
			//double thresh = 0.9985;	// 536
			//double thresh = 0.993;	// 633
			//double thresh = 0.995;	// 443
			//double thresh = 0.9977;	// 635
			//double thresh = 0.9965;	// 634
			//double thresh = 0.9977;	// 444
			//double thresh = 0.9982;	// 636
			//double thresh = 0.9973;	// 363

			//double thresh = 0.999;	// 555
			//double thresh = 0.9982;	// 637, 445
			//thresh = 0.9988;	// 637 scale 2

			//double thresh = 0.995;	// 4333
			//double thresh = 0.997;	// 5333
			//double thresh = 0.998;	// 4343
			//double thresh = 0.998;	// 4343
			//double thresh = 0.999;	// 4353, 5353

			//double thresh = 0.991;	// 534 uniform family .991-4 range
			//double thresh = 0.996;	// 535 uniform family
			//double thresh = 0.994;	// 353 uniform family

			double thresh = 0.997; // Wendy's [77]
			//double thresh = 0.99;	// low res

			//thresh = 0.9992;	// catacombs

			//thresh = 0.97;
			//thresh = 0.6;
			thresh = 0.03;
			return thresh;
		}

		private static bool EdgeOk( Edge edge, Settings s )
		{
			if( s.G == Geometry.Spherical )
				return true;

			double thresh = s.Threshold;

			bool useEdgeLength = s.G == Geometry.Hyperbolic;
			if( useEdgeLength )
			{
				// This will also work for ideal edges.
				return edge.Start.Dist( edge.End ) > thresh;
			}
			else
			{
				return
					edge.Start.Abs() < thresh &&
					edge.End.Abs() < thresh;
			}
		}

		internal static bool CellOk( H3.Cell cell, Settings s )
		{
			if( s.G == Geometry.Spherical )
				return true;

			double thresh = s.Threshold;
			if( !cell.HasVerts )
			{
				foreach( H3.Cell.Facet f in cell.Facets )
				{
					//bool ball = s.Ball;
					bool radiusCutoff = true;
					if( radiusCutoff )
					{
						if( f.Sphere.IsPlane )
							continue;

						//double closestToOrigin = f.Sphere.Center.Abs() - f.Sphere.Radius;
						//if( closestToOrigin > thresh )
						//	return false;

						//if( f.Sphere.Radius < 0.001 )
						//if( f.Sphere.Radius < 0.01 )
						//if( f.Sphere.Radius < 0.02 )
						//if( f.Sphere.Radius < 0.05 )
						if( f.Sphere.Radius < 0.1 )
						//if( f.Sphere.Radius < 0.25 )
							return false;
					}
					else
					{
						double max = 20;
						if( f.Sphere.IsPlane )
						{
							if( f.Sphere.Offset.Abs() > max )
								return false;
						}
						else
						{
							if( f.Sphere.Center.Abs() > max )
								return false;
						}
					}
				}

				return true;
			}

			// Any vertex < threshold makes us ok.
			foreach( Vector3D v in cell.Verts )
				if( v.Abs() < thresh )
					return true;	

			return false;
		}

		/// <summary>
		/// This generates a polyhedron using recursion.  It needs to be finite 
		/// (There are not any other breakouts of the recursion, other than all facets having been generated.)
		/// </summary>
		public static void GenPolyhedron( Sphere[] mirrors, H3.Cell.Facet[] facets,
			List<H3.Cell.Facet> completedFacets, HashSet<Vector3D> completedFacetIds )
		{
			if( 0 == facets.Length )
				return;

			List<H3.Cell.Facet> newFacets = new List<H3.Cell.Facet>();

			foreach( H3.Cell.Facet facet in facets )
				foreach( Sphere mirror in mirrors )
			{
				H3.Cell.Facet newFacet = facet.Clone();
				newFacet.Reflect( mirror );
				if( completedFacetIds.Add( newFacet.ID ) )
				{
					// Haven't seen this facet yet, so 
					// we'll need to recurse on it.
					newFacets.Add( newFacet );
					completedFacets.Add( newFacet );
				}
			}

			GenPolyhedron( mirrors, newFacets.ToArray(), completedFacets, completedFacetIds );
		}

		public static void BranchAlongVerts( Vector3D[] starting, Dictionary<Vector3D, List<H3.Cell.Edge>> vertsToEdges,
			HashSet<Vector3D> foundVertices )
		{
			if( 0 == starting.Length )
				return;

			List<Vector3D> newVerts = new List<Vector3D>();
			foreach( Vector3D v in starting )
			{
				List<H3.Cell.Edge> connectedEdges = vertsToEdges[v];
				foreach( H3.Cell.Edge e in connectedEdges )
				{
					Vector3D opp = e.Start == v ? e.End : e.Start;
					if( foundVertices.Add( opp ) )
						newVerts.Add( opp );
				}
			}

			BranchAlongVerts( newVerts.ToArray(), vertsToEdges, foundVertices );
		}
	}
}
