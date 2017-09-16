namespace R3.Geometry
{
	using System.Collections.Generic;
	using System.Diagnostics;
	using System.Linq;
	using Math = System.Math;
	using R3.Core;
	using R3.Math;

	public static class Polytope
	{
		// The various projections we can do on a polytope.
		public enum Projection
		{
			CellCentered,
			FaceCentered,
			EdgeCentered,
			VertexCentered
		}
	}
}