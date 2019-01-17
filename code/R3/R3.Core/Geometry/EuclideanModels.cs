namespace R3.Geometry
{
	using R3.Geometry;
	using R3.Math;
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

		public static Vector3D SpiralToIsometric( Vector3D v, int p, int m, int n )
		{
			Complex vc = v.ToComplex();
			v = new Vector3D( Math.Log( vc.Magnitude ), vc.Phase );

			Vector3D e1 = new Vector3D( 0, 1 );
			Vector3D e2;
			switch( p )
			{
				case 3:
					e2 = new Vector3D(); break;
				case 4:
					e2 = new Vector3D(); break;
				case 6:
					e2 = new Vector3D(); break;
				default:
					throw new System.ArgumentException();
			}

			double scale = Math.Sqrt( m * m + n * n );
			double a = Euclidean2D.AngleToClock( new Vector3D( 0, 1 ), new Vector3D( m, n ) );

			v.RotateXY( a ); // Rotate
			v *= scale; // Scale

			v *= Math.Sqrt( 2 ) * Geometry2D.EuclideanHypotenuse / ( 2 * Math.PI );
			v.RotateXY( Math.PI / 4 );
			return v;
		}

		public static Vector3D LoxodromicToIsometric( Vector3D v, int p, int m, int n )
		{
			Mobius mob = Mobius.CreateFromIsometry( Geometry.Spherical, 0, new System.Numerics.Complex( 1, 0 ) );
			v = mob.Apply( v );
			return SpiralToIsometric( v, p, m, n );
		}
	}
}
