namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.Linq;
	using R3.Core;
	using R3.Geometry;

	public class Catacombs
	{
		/// <summary>
		/// This is like the GenCell method, but super hacked up for the Catacombs image with Henry.
		/// </summary>
		internal static void GenCellCatacombs( Sphere[] simplex, bool ball )
		{
			// We don't want to include the first mirror (which reflects across cells).
			Sphere[] mirrors = simplex.Skip( 1 ).ToArray();
			Sphere[] allMirrors = simplex.ToArray();

			// Simplices will be the "cells" in Recurse.CalcCells.
			H3.Cell.Facet[] simplexFacets = simplex.Select( m => new H3.Cell.Facet( m ) ).ToArray();

			// Offset cell boundary ever so slightly, to avoid artifacts of adjacent cells.
			Sphere toReflectLater = simplexFacets[0].Sphere.Clone();
			//simplexFacets[0].Sphere = CoxeterImages.GeodesicOffset( simplexFacets[0].Sphere, ball ? -1e-6 : 1e-7, ball );

			H3.Cell startingCell = new H3.Cell( simplexFacets );
			startingCell = startingCell.Clone();	// So our mirrors don't get munged after we reflect around later.
			Vector3D cen = new Vector3D( 0.05, 0.01, -0.05 );		//	373, 438
			//Vector3D cen = new Vector3D( 0.05, 0.01, 100 );		//	637
			//cen.RotateXY( Math.PI / 2 );	// only if we also rotate simplex mirrors.  XXX - make a setting.
			startingCell.Center = cen;
			H3.Cell[] simplices = Recurse.CalcCells( mirrors, new H3.Cell[] { startingCell }, new Recurse.Settings() { Ball = ball } );

			List<H3.Cell> simplicesFinal = new List<H3.Cell>();
			List<int[]> reflectionSets = new List<int[]>();
			// 1 reflects in x-axis
			// 3, 1 rotates right
			// 1, 3 rotates left

			reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 0, 1, 3 } );
			reflectionSets.Add( new int[] { 0, 3, 1, 2, 3, 1, 0, 3, 1, 2, 0, 3 } );

			// 2
			reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 3, 1, 0, 3, 1, 3, 1, 3, 1, 2, 3, 1, 3, 2, 3, 1 } );
			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 3, 1, 0,		3, 1, 3, 1, 0,	3, 1, 3, 1, 3, 1, 2,	3, 1, 3, 2, 3, 1 } );

			// 3
			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 3, 1, 0, 3, 1, 3, 1, 0, 3, 1, 2 } );
			reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 3, 0, 3, 1, 3, 1, 0, 2 } );
			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 3, 1, 0, 3, 1, 2 } );

			// 5
			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 3, 1, 0, 1, 2, 3, 1 } );

			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 0, 2, 3, 1 } );
			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 0, 3, 1, 3, 1, 3, 1, 2, 3, 1, 3, 1, 0, 1, 3 } );	// baby
			//reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 0, 3, 1, 3, 1, 3, 1, 2 } );	// maybe
			//reflectionSets.Add( new int[] { 0, 3, 1, 2, 3, 1, 0, 2, 1, 3 } );
			//reflectionSets.Add( new int[] { 0, 3, 1, 2, 3, 1, 0, 3, 1, 2 } );		// not great orientation
			// reflectionSets.Add( new int[] { 0, 3, 1, 3, 1, 2 } );	// big

			bool ceiling = true;
			if( ceiling )
				simplicesFinal = simplices.ToList();
			else
			{
				foreach( int[] set in reflectionSets )
				{
					List<H3.Cell> copy = simplices.Select( s => s.Clone() ).ToList();
					foreach( int r in set )
					{
						foreach( H3.Cell cell in copy )
							cell.Reflect( allMirrors[r] );
					}
					simplicesFinal.AddRange( copy );
				}
			}

			/*
			// A second cell.
			//toReflectLater = simplices[45].Facets[0].Sphere.Clone();
			//toReflectLater = simplices.First( s => s.Depth == 2 ).Facets[0].Sphere.Clone();
			foreach( H3.Cell cell in simplices )
				cell.Reflect( toReflectLater );

			// A third cell.
			toReflectLater = simplices[40].Facets[0].Sphere.Clone();
			//toReflectLater = simplices.First( s => s.Depth == 4 ).Facets[0].Sphere.Clone();
			foreach( H3.Cell cell in simplices )
				cell.Reflect( toReflectLater );
			
			foreach( H3.Cell cell in simplices )
				cell.Depths = new int[4];
			List<H3.Cell> simplicesFinal = Recurse.CalcCells2( mirrors, simplices ).ToList();
			simplicesFinal = simplicesFinal.Where( s => s.Depths[0] % 3 == 1 && s.Depths[1] % 2 == 0 && s.Depths[2] % 2 == 1 ).ToList();
			*/

			/*
			List<H3.Cell> simplicesFinal = new List<H3.Cell>();
			//for( int d = 0; d < 1; d+=2 )
			int d = 0;
			{
				//Sphere toReflect = simplices.First( s => s.Depth == d ).Facets[0].Sphere.Clone();
				//Sphere toReflect = simplices.Where( s => s.Depth == d ).Skip(1).Take(1).First().Facets[0].Sphere.Clone();
				List<H3.Cell> reflectionCells = simplices.Where( s => s.Depths[1] == d && s.Depths[0] % 2 == 0 ).Skip(0).Take(1).ToList();
				foreach( Sphere toReflect in reflectionCells.Select( c => c.Facets[0].Sphere.Clone() ) )
				{
					List<H3.Cell> thisCell = new List<H3.Cell>();
					foreach( H3.Cell cell in simplices )
					{
						H3.Cell clone = cell.Clone();
						clone.Reflect( toReflect );
						thisCell.Add( clone );
					}

					//Sphere toReflect2 = thisCell.First( s => s.Depth1 == d + 3 && s.Depth0 % 2 == 0 ).Facets[0].Sphere.Clone();
					//List<H3.Cell> reflectionCellsTemp = simplices.Where( s => Math.Abs( s.Depths[1] - d ) == 2 && s.Depths[0] % 2 == 0 ).ToList();
					List<H3.Cell> reflectionCellsTemp = simplices.Where( s => s.Depths[1] == 2 && s.Depths[1] == s.Depths[0] + s.Depths[2] ).ToList();
					List<H3.Cell> reflectionCells2 = reflectionCellsTemp;//.Where( ( x, i ) => i % 3 == 0 ).ToList(); // .Skip( 5 ).Take( 5 ).ToList();
					foreach( Sphere toReflect2 in reflectionCells2.Select( c => c.Facets[0].Sphere.Clone() ) )
					//Sphere toReflect2 = toReflectLater;
					{
						foreach( H3.Cell cell in thisCell )
						{
							H3.Cell clone = cell.Clone();
							clone.Reflect( toReflect2 );
							simplicesFinal.Add( clone );
						}
					}
				}
			}*/

			int count = 0;
			foreach( H3.Cell cell in simplicesFinal )
			{
				count++;

				//if( count % 2 == 0 )
				//	continue;
				/*if( count < 1 )
					continue;
				if( count > 30 )
					return;
				 */
				//int[] include = new int[] { 0, 1, 2, 3 };
				int[] include = new int[] { 0 };
				PovRay.AppendSimplex( cell.Facets.Select( f => f.Sphere ).ToArray(), cell.Center, include, "cell.pov" );
			}
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
	}
}
