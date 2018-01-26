namespace HyperbolicModels
{
	using System.Drawing;
	using System.IO;
	using System.Linq;
	using R3.Core;
	using R3.Drawing;
	using R3.Geometry;

	internal class Tree
	{
		public static void Create( HoneycombDef def, string filename)
		{
			int p = def.P;
			int q = def.Q;
			int r = def.R;

			double scale = 5.0;
			Vector3D cen = HoneycombPaper.InteriorPointBall;

			Sphere[] simplex = SimplexCalcs.Mirrors( p, q, r, moveToBall: false );

			// Apply transformations.
			simplex = simplex.Select( s =>
			{
				Sphere.ScaleSphere( s, scale );
				return H3Models.UHSToBall( s );
			} ).ToArray();

			for( int i = 0; i < 4; i++ )
				if( simplex[i].IsPointInside( cen ) )
					simplex[i].Invert = true;

			Sphere[] simplexForColorScale = SimplexCalcs.Mirrors( p, q, r, moveToBall: true );
			CoxeterImages.Settings temp = HoneycombPaper.AutoCalcScale( def, simplexForColorScale );
			int maxDepth = (int)temp.ColorScaling;

			bool ball = true;
			bool dual = false;
			H3.Cell[] simplicesFinal = HoneycombPaper.GenCell( simplex, null, cen, ball, dual );

			simplicesFinal = simplicesFinal.Where( s => s.Depths[0] < 1 ).ToArray();
			//simplicesFinal = simplicesFinal.Where( s => s.)

			// Output the facets.
			using( StreamWriter sw = File.CreateText( filename ) )  // We need to reuse this StreamWriter (vs. calling AppendSimplex) for performance.
			{
				sw.WriteLine( "#include \"hyper_ball.pov\"" );
				int[] include = new int[] { 0 };
				foreach( H3.Cell cell in simplicesFinal )
				{
					Sphere[] facets = cell.Facets.Select( f => f.Sphere ).ToArray();
					int depth = cell.Depths[0] + 1;
					Color c = Coloring.ColorAlongHexagon( maxDepth, depth );
					PovRay.AddSimplex( sw, facets, cell.Center, include, filename, Coloring.ToVec( c ) );
				}
			}
		}
	}
}
