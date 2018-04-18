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
		public enum EdgeThreshType
		{
			Radial,
			Length
		}

		public class Settings
		{
			public Geometry G = Geometry.Hyperbolic;
			public int MaxEdges = 2500000;
			public double Threshold = DefaultThresh();
			public EdgeThreshType ThreshType = EdgeThreshType.Length;
			public bool Ball = true;
		}

		/// <summary>
		/// Attempts to calculate approx 1.3M edges when the threshold is a distance from origin in the ball model.
		/// This works for honeycombs with finite cells.
		/// </summary>
		public static Edge[] CalcEdgesSmart( Sphere[] simplex, Edge[] edges, int desiredCount )
		{
			Settings s = new Settings();
			s.ThreshType = EdgeThreshType.Radial;

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
			double logDesiredCount = Math.Log( desiredCount );
			hDist = 5.5 + ( logDesiredCount - Math.Log( count2 ) ) / slope;

			s.Threshold = DonHatch.h2eNorm( hDist );
			return CalcEdges( simplex, edges, s );
		}

		/// <summary>
		/// Attempts to calculate approx 1.3M edges when the threshold is a minimum edge length.
		/// This is required for honeycombs with ideal or ultra-ideal cells
		/// </summary>
		public static Edge[] CalcEdgesSmart2( Sphere[] simplex, Edge[] edges, int desiredCount )
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
			double logDesiredCount = Math.Log( (double)desiredCount );
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
				// This also controls resolution of the edges, and can have a big effect on file size.
				// Make the threshold cutoff black, or the background color.
				double percentWhite = 1;
				if( settings.ThreshType == EdgeThreshType.Length )
					percentWhite = (r1.Dist( r2 ) - settings.Threshold) / 0.015;
				else
				{
					double closestToOrigin = Math.Min( r1.Abs(), r2.Abs() );	// Mainly ranges from 0 to 1
					if( closestToOrigin < 0.9 )
						percentWhite = 1.0;
					else
						percentWhite = 1.0 - Math.Pow( closestToOrigin - 0.9, 1.3 ) / 0.1;
				}
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

		public static H3.Cell[] CalcCellsSmart( Sphere[] mirrors, H3.Cell[] cells, Settings settings, int desiredCount )
		{
			double t1 = 80;
			double t2 = 130;

			// I found that log(1/thresh)/log(count) was relatively constant,
			// so we'll extrapolate that to get close to the right number of edges.
			double OneOverThresh = t1;
			settings.Threshold = 1 / OneOverThresh;
			H3.Cell[] result = CalcCells( mirrors, cells, settings );
			int count1 = result.Length;
			System.Console.WriteLine( string.Format( "count1: {0}", count1 ) );

			OneOverThresh = t2;
			settings.Threshold = 1 / OneOverThresh;
			result = CalcCells( mirrors, cells, settings );
			int count2 = result.Length;
			System.Console.WriteLine( string.Format( "count2: {0}", count2 ) );

			double slope = (Math.Log( count2 ) - Math.Log( count1 )) / (Math.Log( t2 ) - Math.Log( t1 ));
			double logDesiredCount = Math.Log( (double)desiredCount );
			double temp = Math.Log( t2 ) + (logDesiredCount - Math.Log( count2 )) / slope;

			settings.Threshold = 1 / Math.Exp( temp );
			System.Console.WriteLine( string.Format( "Setting threshold to: {0}", settings.Threshold ) );
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
			//for( int m=0; m<simplex.Length; m++ )
			for( int m=simplex.Length-1; m>=0; m-- )
			{
				Sphere mirror = simplex[m];
				//if( m == 2 )
				//	continue;

				//if( completedCellIds.Count > 1000 )
				//	return;

				//if( completedCellIds.Count > settings.MaxEdges/5 )
				if( completedCellIds.Count > settings.MaxEdges * 3 )
					throw new System.Exception( "Maxing out cells - will result in uneven filling." );

				H3.Cell newCell = cell.Clone();	
				newCell.Reflect( mirror );
				
				// This tracks reflections across the cell facets.
				newCell.Depths[m]++;
				newCell.LastReflection = m;

				//if( newCell.Depths[0] > 2 )
				//	continue;

				if( !CellOk( newCell, settings ) )
					continue;

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

		private static double DefaultThresh()
		{
			return 0.03;
		}

		private static bool EdgeOk( Edge edge, Settings s )
		{
			if( s.G == Geometry.Spherical )
				return true;

			double thresh = s.Threshold;

			switch( s.ThreshType )
			{
			case EdgeThreshType.Length:

				// This will also work for ideal edges.
				return edge.Start.Dist( edge.End ) > thresh;

			case EdgeThreshType.Radial:

				return
					edge.Start.Abs() < thresh &&
					edge.End.Abs() < thresh;
			}

			return false;
		}

		internal static bool CellOk( H3.Cell cell, Settings s )
		{
			if( s.G == Geometry.Spherical )
				return true;

			double thresh = s.Threshold;
			if( !cell.HasVerts )
			{
				//return cell.Depths.Sum() <= 20;

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

						if( f.Sphere.Radius < thresh )
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
