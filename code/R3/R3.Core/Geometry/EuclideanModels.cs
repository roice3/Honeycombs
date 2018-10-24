namespace R3.Geometry
{
	using R3.Geometry;
	using System.Numerics;
	using Math = System.Math;

	public enum EuclideanModel
	{
		Isometric,
		Conformal,
		Disk,
		UpperHalfPlane,
		Spiral,
		Loxodromic,
	}

	public class EuclideanModels
	{
		public static Vector3D DiskToIsometric( Vector3D v )
		{
			// ZZZ - Check that this is correct (it's quite close if not!)
			return SphericalModels.StereoToGnomonic( v );
		}

		public static Vector3D UpperHalfPlaneToIsometric( Vector3D v )
		{
			v = HyperbolicModels.UpperToPoincare( v );
			v = SphericalModels.StereoToGnomonic( v );
			return v;	
		}

		public static Vector3D SpiralToIsometric( Vector3D v )
		{
			//Mobius mob = Mobius.CreateFromIsometry( Geometry.Spherical, 0, new System.Numerics.Complex( 1, 0 ) );
			//v = mob.Apply( v );

			Complex vc = v.ToComplex();


			v = new Vector3D( Math.Log( vc.Magnitude ), vc.Phase );
			//vc = Complex.Exp( vc );
			//v = Vector3D.FromComplex( vc );

			int m_ = 7;
			int n_ = 3;
			double scale = Math.Sqrt( m_ * m_ + n_ * n_ );
			double a = Euclidean2D.AngleToClock( new Vector3D( 0, 1 ), new Vector3D( m_, n_ ) );

			v.RotateXY( a ); // Rotate
			v *= scale; // Scale

			//v *= 4;	// Make the grid more dense.

			v *= Math.Sqrt( 2 ) * Geometry2D.EuclideanHypotenuse / ( 2 * Math.PI );
			v.RotateXY( Math.PI / 4 );
			return v;
		}

		public static Vector3D LoxodromicToIsometric( Vector3D v )
		{
			throw new System.NotImplementedException();
		}
	}
}
