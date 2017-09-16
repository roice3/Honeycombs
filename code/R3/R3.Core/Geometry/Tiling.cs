namespace R3.Geometry
{
	using Math = System.Math;
	using R3.Core;
	using R3.Math;
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Linq;

	/// <summary>
	/// Information we need for a tiling.
	/// </summary>
	public struct TilingConfig
	{
		public TilingConfig( int p, int q, int maxTiles ) : this() 
		{
			SetupConfig( p, q, maxTiles );
		}

		public TilingConfig( int p, int q ) : this()
		{
			if( Geometry2D.GetGeometry( p, q ) != Geometry.Spherical )
				throw new System.ArgumentException();

			SetupConfig( p, q, PlatonicSolids.NumFacets( p, q ) );
		}

		private void SetupConfig( int p, int q, int maxTiles )
		{
			P = p;
			Q = q;
			m.Unity();
			MaxTiles = maxTiles;
			Shrink = 1.0;
		}

		public int P { get; set; }
		public int Q { get; set; }

		/// <summary>
		/// The induced geometry.
		/// </summary>
		public Geometry Geometry { get { return Geometry2D.GetGeometry( this.P, this.Q ); } }

		/// <summary>
		/// A Mobius transformation to apply while creating the tiling.
		/// </summary>
		public Mobius M { get { return m; } set { m = value; } }
		private Mobius m;

		/// <summary>
		/// The max number of tiles to include in the tiling.
		/// </summary>
		public int MaxTiles { get; set; }

		/// <summary>
		/// A shrinkage to apply to the drawn portion of a tile.
		/// Default is 1.0 (no shrinkage).
		/// </summary>
		public double Shrink { get; set; }
	}

	public class Tiling
	{
		public Tiling()
		{
			m_tiles = new List<Tile>();
			this.TilePositions = new Dictionary<Vector3D, Tile>();
		}

		/// <summary>
		/// The tiling configuration.
		/// </summary>
		public TilingConfig TilingConfig { get; set; }

		/// <summary>
		/// Our tiles.
		/// </summary>
		private List<Tile> m_tiles;

		/// <summary>
		/// A dictionary from tile centers to tiles.
		/// </summary>
		public Dictionary<Vector3D, Tile> TilePositions { get; set; }

		public static Tile CreateBaseTile( TilingConfig config )
		{
			Polygon boundary = new Polygon(), drawn = new Polygon();
			boundary.CreateRegular( config.P, config.Q );
			drawn = boundary.Clone();

			//boundary.CreateRegular( 3, 10 );
			//drawn.CreateRegular( 3, 8 );
			//boundary.CreateRegular( 3, 7 );
			//drawn = Heart();

			//for( int i=0; i<drawn.NumSides; i++ )
			//	drawn.Segments[i].Center *= 0.1;

			// Good combos:
			// ( 5, 5 ), ( 10, 10 )
			// ( 3, 10 ), ( 3, 9 )
			// ( 6, 4 ), ( 6, 8 )
			// ( 7, 3 ), ( 7, 9 )

			Tile tile = new Tile( boundary, drawn, config.Geometry );
			Tile.ShrinkTile( ref tile, config.Shrink );
			return tile;
		}

		/// <summary>
		/// The number of tiles.
		/// </summary>
		public int Count
		{
			get { return m_tiles.Count; }
		}

		/// <summary>
		/// Access to all the tiles.
		/// </summary>
		public IEnumerable<Tile> Tiles
		{
			get { return m_tiles; }
		}

		/// <summary>
		/// Retrieve all the polygons in this tiling that we want to draw.
		/// </summary>
		public IEnumerable<Polygon> Polygons 
		{
			get
			{
				return m_tiles.Select( t => t.Drawn );
			}
		}

		/// <summary>
		/// Retreive all the (non-Euclidean) vertex circles in this tiling.
		/// </summary>
		public IEnumerable<CircleNE> Circles
		{
			get
			{
				return m_tiles.Select( t => t.VertexCircle );
			}
		}
	}
}
