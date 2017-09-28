namespace R3.Geometry
{
	using System.Numerics;

	/// <summary>
	/// Class with some hackish methods for dealing with points projected to infinite.
	/// </summary>
	public static class Infinity
	{
		public static Vector3D InfinityVector = new Vector3D( double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity );
		public static Vector3D LargeFiniteVector = new Vector3D( FiniteScale, FiniteScale, FiniteScale );

		public const double FiniteScale = 10000;
		public const double InfiniteScale = 500000;

		public static bool IsFinite( double input )
		{
			return ( -InfiniteScale <= input ) && ( input <= InfiniteScale );
		}

		public static bool IsInfinite( Vector3D input )
		{
			return
				!(IsFinite (input.X) &&
				  IsFinite (input.Y) &&
			      IsFinite (input.Z) &&
			      IsFinite (input.W));
		}

		public static bool IsInfinite( Complex input )
		{
			return
				IsInfinite( input.Real ) ||
				IsInfinite( input.Imaginary );
		}

		public static bool IsInfinite( double input )
		{
			return !IsFinite (input);
		}

		public static Vector3D InfinitySafe( Vector3D input )
		{
			if( Infinity.IsInfinite( input ) )
				return Infinity.LargeFiniteVector;
			return input;
		}
	}
}
